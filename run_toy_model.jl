using JuMP, Gurobi, DataFrames, CSV, Plots, Printf

# 1. Load your empirical Ashburn, VA single data center profile
df_dc = CSV.read("single_dc_profile.csv", DataFrame)

# Pick 2 representative periods from Day 1:
# Period 1 = Hour 4 (3:00 AM - Cool night, 0% solar, low grid load)
# Period 2 = Hour 15 (2:00 PM - Hot afternoon, high cooling PUE, high solar, peak grid congestion)
idx_t = [4, 15]
L_dc = df_dc.DC_Load_MW[idx_t]      # Empirical Data Center Load (MW)
CF_sol = df_dc.Solar_CF[idx_t]      # Empirical Solar Capacity Factor [0, 1]

println("Selected Toy Model Periods from single_dc_profile.csv:")
for (k, h) in enumerate(idx_t)
    @printf("  Period %d (%s): Temp = %.1f°C, PUE = %.3f, DC Load = %.2f MW, Solar CF = %.3f\n",
            k, df_dc.Timestamp[h], df_dc.Temp_C[h], df_dc.Cooling_PUE[h], L_dc[k], CF_sol[k])
end

# 2. Solve 3-Bus Toy Model Across 3 Regimes:
#    Regime A: Islanded BTM Bridge (No Grid Connection: K_imp = 0, K_exp = 0)
#    Regime B: Grid-Connected Zero-Export Baseline (K_imp > 0, K_exp = 0)
#    Regime C: Grid-Connected Export-Enabled / On-Demand Contract (K_imp > 0, K_exp >= 0)
function solve_dc_regime(regime::String)
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    N = 1:3
    T = 1:2
    lines = [(1,2), (1,3), (2,3)]

    # Cost Parameters (Anchored in NREL ATB)
    c_u1 = 25.0          # Bus 1 Utility Baseload ($/MWh)
    c_u3 = 120.0         # Bus 3 Utility Peaker ($/MWh)
    c_gas = 45.0         # Bus 2 Data Center BTM Gas/LNG marginal cost ($/MWh)
    cap_gas = 15.0       # Amortized BTM Gas CAPEX ($/MW over horizon)
    cap_sol = 12.0       # Amortized BTM Solar CAPEX ($/MW over horizon)
    kappa_exp = 10.0     # Amortized Export Interconnection Upgrade cost ($/MW)

    # Grid Parameters
    D_grid = Dict((1,1)=>0.0, (1,2)=>0.0, (2,1)=>0.0, (2,2)=>0.0, (3,1)=>100.0, (3,2)=>220.0)
    B = Dict((1,2)=>1000.0, (1,3)=>1000.0, (2,3)=>1000.0)
    F_max = Dict((1,2)=>200.0, (1,3)=>120.0, (2,3)=>200.0) # Line 1-3 binds at 120 MW in Period 2

    # Decision Variables
    @variable(model, 0 <= p_u1[t in T] <= 300.0) # Bus 1 Utility Gen (MW)
    @variable(model, 0 <= p_u3[t in T] <= 200.0) # Bus 3 Peaker Gen (MW)
    
    # Data Center Sizing & Dispatch Variables at Bus 2
    @variable(model, 0 <= K_gas <= 150.0)        # Installed BTM Gas Capacity (MW)
    @variable(model, 0 <= K_sol <= 80.0)         # Installed BTM Solar Capacity (MW)
    @variable(model, g_gas[t in T] >= 0)         # BTM Gas generation (MW)
    @variable(model, g_sol[t in T] >= 0)         # BTM Solar generation (MW)
    @variable(model, p_imp[t in T] >= 0)         # Grid Import to Data Center (MW)
    @variable(model, p_exp[t in T] >= 0)         # Grid Export from Data Center (MW)
    @variable(model, K_exp >= 0)                 # Export Interconnection Capacity (MW)

    @variable(model, theta[i in N, t in T])
    @variable(model, f[l in lines, t in T])

    # Regime-specific Interconnection Rules
    if regime == "Islanded_BTM"
        for t in T
            fix(p_imp[t], 0.0; force=true)
            fix(p_exp[t], 0.0; force=true)
        end
        fix(K_exp, 0.0; force=true)
    elseif regime == "Grid_Zero_Export"
        fix(K_gas, 150.0; force=true) # Retains BTM fleet built during bridge phase
        fix(K_sol, 50.0; force=true)
        for t in T
            fix(p_exp[t], 0.0; force=true)
        end
        fix(K_exp, 0.0; force=true)
    elseif regime == "Grid_Export_Enabled"
        fix(K_gas, 150.0; force=true) # Retains BTM fleet & unlocks grid export
        fix(K_sol, 50.0; force=true)
    end

    # Objective: Minimize Total System Cost (Utility + Data Center CAPEX/OPEX + Interconnection)
    @objective(model, Min,
        sum(c_u1 * p_u1[t] + c_u3 * p_u3[t] + c_gas * g_gas[t] for t in T) +
        cap_gas * K_gas + cap_sol * K_sol + kappa_exp * K_exp
    )

    # Constraints
    @constraint(model, [t in T], g_gas[t] <= K_gas)
    @constraint(model, [t in T], g_sol[t] <= CF_sol[t] * K_sol)
    @constraint(model, dc_bal[t in T], g_gas[t] + g_sol[t] + p_imp[t] - p_exp[t] == L_dc[t])
    @constraint(model, [t in T], p_exp[t] <= K_exp)

    for (i, j) in lines, t in T
        @constraint(model, f[(i,j), t] == B[(i,j)] * (theta[i, t] - theta[j, t]))
        @constraint(model, -F_max[(i,j)] <= f[(i,j), t] <= F_max[(i,j)])
    end
    @constraint(model, [t in T], theta[1, t] == 0.0)

    @constraint(model, bal1[t in T], D_grid[(1,t)] + f[(1,2), t] + f[(1,3), t] - p_u1[t] == 0.0)
    @constraint(model, bal2[t in T], D_grid[(2,t)] + p_imp[t] - p_exp[t] - f[(1,2), t] + f[(2,3), t] == 0.0)
    @constraint(model, bal3[t in T], D_grid[(3,t)] - f[(1,3), t] - f[(2,3), t] - p_u3[t] == 0.0)

    optimize!(model)

    return DataFrame(
        Regime = fill(regime, 2),
        Period = [1, 2],
        DC_Load_MW = L_dc,
        K_Gas_MW = fill(round(value(K_gas), digits=1), 2),
        K_Sol_MW = fill(round(value(K_sol), digits=1), 2),
        K_Exp_MW = fill(round(value(K_exp), digits=1), 2),
        DC_Gas_Gen = [round(value(g_gas[t]), digits=1) for t in T],
        DC_Sol_Gen = [round(value(g_sol[t]), digits=1) for t in T],
        DC_Import = [round(value(p_imp[t]), digits=1) for t in T],
        DC_Export = [round(value(p_exp[t]), digits=1) for t in T],
        Bus3_Peaker = [round(value(p_u3[t]), digits=1) for t in T],
        LMP_Bus2 = [round(abs(dual(bal2[t])), digits=2) for t in T],
        LMP_Bus3 = [round(abs(dual(bal3[t])), digits=2) for t in T],
        Total_Cost = fill(round(objective_value(model), digits=2), 2)
    )
end

df_all = vcat(
    solve_dc_regime("Islanded_BTM"),
    solve_dc_regime("Grid_Zero_Export"),
    solve_dc_regime("Grid_Export_Enabled")
)

CSV.write("toy_example_results.csv", df_all)
println("\n=== THREE-REGIME TOY MODEL COMPARISON ===")
println(df_all)
# =============================================================================
#  toy_example.jl — proposal toy example: 3 buses x 2 hours (Section 5)
#
#  Same formulation as the proposal, Section 3, equations (1)-(12); metrics (13)-(15).
#  Network: triangle, every line x = 0.1 p.u. (B = 1,000 MW/rad), limit 100 MW.
#    Bus 1  G1: 300 MW at $20/MWh (reference bus)
#    Bus 2  G2: 70 MW peaker at $80/MWh + town load 40 MW (hour 1), 120 MW (hour 2)
#    Bus 3  town load 20 / 50 MW + the data center: flat 100 MW, 50 MW solar
#  Hour 1 = night (no sun), hour 2 = evening peak (solar CF 0.4).
#  Scenarios: S0 (sigma = 0, beta = 0) plus sigma in {0,25,50,100}% x beta in {25,50,100}%.
#  Last section: value of lost compute V swept over V_SWEEP (S0, S2, S7, S12).
#  Needs: JuMP, HiGHS, DataFrames, CSV.   Run:  julia toy_example.jl
# =============================================================================
using JuMP, HiGHS, DataFrames, CSV, Printf

# ---- sets
N, T = 1:3, 1:2
lines = [(1, 2), (1, 3), (2, 3)]                 # (i, j)
gen_bus = [1, 2]                                 # G1 at bus 1, G2 at bus 2
# ---- parameters
c_g   = [20.0, 80.0]                             # $/MWh
P_bar = [300.0, 70.0]                            # MW
B     = 1000.0                                   # MW/rad  (100 MVA / x = 0.1)
F_max = 100.0                                    # MW, every line
D     = [0.0 0.0; 40.0 120.0; 20.0 50.0]         # non-DC demand D[i,t] (MW)
dc    = 3                                        # data-center bus
L     = [100.0, 100.0]                           # DC load (MW), flat
L_hat = maximum(L)
CF    = [0.0, 0.4]; K_sol = 50.0                 # solar capacity factor, solar MW
K_imp = L_hat                                    # substation import limit
c_gas = 9.717 * 3.50 + 6.94                      # $40.95/MWh
V     = 3000.0                                   # value of lost compute ($/MWh): GPU rental revenue per MWh (H100/B200 neocloud)
V_SWEEP = [600.0, 1000.0, 3000.0, 6600.0]        # V sensitivity: spot H100 ... AWS B200 rental prices
R_g   = [100.0, 80.0]                            # ramp limits (MW/h)
# BTM capital: NREL ATB 2024, F-frame CT, Moderate, 2025: CAPEX $1,319.40/kW, FOM $25.70/kW-yr,
# annualized with the ATB real WACC 5.36% over 30 years, then spread over the hours of a year
r, n  = 0.0535848, 30
CRF   = r / (1 - (1 + r)^-n)
C_btm_hour = (CRF * 1319.40 + 25.70) * 1000 / 8760   # $/MW per hour of the horizon

function solve_toy(sigma, beta; V = V)
    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, 0 <= p[g = 1:2, t = T] <= P_bar[g])            # (11)
    @variable(m, 0 <= g_gas[t = T] <= beta * L_hat)             # (6)
    @variable(m, 0 <= g_sol[t = T] <= CF[t] * K_sol)            # (7)
    @variable(m, 0 <= p_imp[t = T] <= K_imp)                    # import <= load (implied by (3))
    @variable(m, u[t = T] >= 0)                                 # u >= 0
    @variable(m, -F_max <= f[l = 1:3, t = T] <= F_max)          # (9)
    @variable(m, theta[N, T])
    @expression(m, c[t = T], L[t] - g_sol[t] - p_imp[t])        # (4) curtailment
    # (1) objective
    @objective(m, Min, sum(c_g[g] * p[g, t] for g in 1:2, t in T) + sum(c_gas * g_gas[t] + V * u[t] for t in T))
    # (2) nodal balance: generation + flows in - flows out - DC import = demand  (dual = LMP)
    @constraint(m, bal[i = N, t = T],
        sum(p[g, t] for g in 1:2 if gen_bus[g] == i) +
        sum(f[l, t] for (l, (a, b)) in enumerate(lines) if b == i) -
        sum(f[l, t] for (l, (a, b)) in enumerate(lines) if a == i) -
        (i == dc ? p_imp[t] : 0.0) == D[i, t])
    @constraint(m, [t = T], g_gas[t] + g_sol[t] + p_imp[t] + u[t] == L[t])        # (3)
    @constraint(m, [t = T], 0 <= c[t] - u[t] <= sigma * L[t])                     # (5)
    @constraint(m, [(l, (a, b)) = enumerate(lines), t = T], f[l, t] == B * (theta[a, t] - theta[b, t]))  # (8)
    @constraint(m, [(l, (a, b)) = enumerate(lines), t = T], -pi/6 <= theta[a, t] - theta[b, t] <= pi/6)  # (10)
    @constraint(m, [t = T], theta[1, t] == 0)                                     # (10) reference bus
    @constraint(m, [g = 1:2, t = 2:length(T)], -R_g[g] <= p[g, t] - p[g, t-1] <= R_g[g])  # (12) ramp limits
    optimize!(m)
    termination_status(m) == MOI.OPTIMAL || error("not optimal")
    return (Z = objective_value(m), p = value.(p), gas = value.(g_gas), sol = value.(g_sol),
            imp = value.(p_imp), u = value.(u), c = value.(c), f = value.(f), lmp = dual.(bal))
end

cases = [("S0", 0.0, 0.0)]
k = 0
for b in (0.25, 0.5, 1.0), s in (0.0, 0.25, 0.5, 1.0)
    global k += 1; push!(cases, ("S$k", s, b))
end
res = Dict(n => solve_toy(s, b) for (n, s, b) in cases)
Z0 = res["S0"].Z

rows = DataFrame()
for (n, s, b) in cases, t in T
    r = res[n]
    push!(rows, (scenario = n, sigma = s, beta = b, btm_gas_MW = b * L_hat, hour = t,
                 G1_MW = r.p[1, t], G2_peaker_MW = r.p[2, t], dc_import_MW = r.imp[t], solar_MW = r.sol[t],
                 btm_gas_MW_used = r.gas[t], curtailed_MW = r.c[t], lost_compute_MW = r.u[t],
                 flow_1_3_MW = r.f[2, t], lmp_bus1 = r.lmp[1, t], lmp_bus2 = r.lmp[2, t], lmp_bus3 = r.lmp[3, t],
                 total_cost = r.Z, VoF = Z0 - r.Z,
                 NSB = Z0 - r.Z - C_btm_hour * length(T) * b * L_hat); promote = true)
end
rows = mapcols(x -> x isa AbstractVector{<:Real} ? round.(x; digits = 2) .+ 0.0 : x, rows)
mkpath(joinpath(@__DIR__, "results"))
CSV.write(joinpath(@__DIR__, "results", "toy_example.csv"), rows)
println("Hourly results, hour 2 (evening peak):")
show(rows[rows.hour .== 2, :]; allcols = true, allrows = true); println()
println("\nTotal cost and value of flexibility (2 hours):")
for (n, s, b) in cases
    @printf("%-4s sigma %4.0f%%  beta %4.0f%%  Z = \$%9.2f  VoF = \$%8.2f  NSB = \$%8.2f\n", n, 100s, 100b, res[n].Z,
            Z0 - res[n].Z, Z0 - res[n].Z - C_btm_hour * length(T) * b * L_hat)
end

# ---- V sensitivity: value of lost compute
vrows = DataFrame()
for Vs in V_SWEEP
    z0 = solve_toy(0.0, 0.0; V = Vs).Z
    for (n, sg, b) in (("S0", 0.0, 0.0), ("S2", 0.25, 0.25), ("S7", 0.5, 0.5), ("S12", 1.0, 1.0))
        r = solve_toy(sg, b; V = Vs)
        push!(vrows, (V = Vs, scenario = n, total_cost = r.Z, VoF = z0 - r.Z,
                      NSB = z0 - r.Z - C_btm_hour * length(T) * b * L_hat, lost_compute_MWh = sum(r.u),
                      lmp_bus3_hour2 = r.lmp[3, 2]))
    end
end
vrows = mapcols(x -> x isa AbstractVector{<:Real} ? round.(x; digits = 2) .+ 0.0 : x, vrows)
CSV.write(joinpath(@__DIR__, "results", "toy_V_sensitivity.csv"), vrows)
println("\nV sensitivity (value of lost compute, \$/MWh):")
show(vrows; allrows = true); println()

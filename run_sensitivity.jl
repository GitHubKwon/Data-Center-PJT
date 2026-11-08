# =============================================================================
#  run_sensitivity.jl — Regime C scenarios over the 7-day horizon, one location
#
#  1. Scenario grid: baseline S0 (sigma = beta = 0) plus sigma x beta (S1-S12),
#     with VoF = Z(S0) - Z(sigma, beta) and NSB = VoF - BTM capital cost
#  Two sensitivities (all other inputs at base values):
#  2. Congestion: line ratings at 75% (base), 80%, 90% of the rating, S0 vs S12
#  3. Value of lost compute V: S0, S2, S7 and S12 at each V in V_SWEEP
#  Run:  julia run_sensitivity.jl
# =============================================================================

include(joinpath(@__DIR__, "setup.jl"))
solve_case(par; ls = LINE_SCALE) = summarize(solve_dispatch(net, D_week, L_week, CF_week, par;
                                                            line_scale = ls, solver = SOLVER))
btm_cost_k(b) = BTM_COST_PER_MW_HOUR * length(L_week) * b * DC_PEAK / 1e3     # BTM capital over the horizon ($k)

# ---- 1. sigma x beta grid
grid = scenario_grid(net, D_week, L_week, CF_week, BASE; sigmas = SIGMAS, betas = BETAS,
                     line_scale = LINE_SCALE, solver = SOLVER)
CSV.write(joinpath(RESULTS, "sensitivity_sigma_beta.csv"), grid)
@printf("\nScenario grid (7 days, V = \$%.0f/MWh):\n", BASE.v_dc)
show(select(grid, :scenario, :sigma, :beta, :btm_gas_MW, :total_cost_k, :VoF_k, :dVoF_k, :NSB_k, :v_per_MWh,
            :curtail_hours, :curtailed_MWh, :gas_MWh, :dc_lost_MWh, :co2_t, :lmp_dc_mean, :lmp_dc_max);
     allrows = true, allcols = true)
println()

# ---- 2. congestion: line ratings (S0 and S12 at each setting)
settings = [("line rating 75% (base)", BASE, LINE_SCALE),
            ("line rating 80%", BASE, 0.80),
            ("line rating 90%", BASE, 0.90)]
oat = DataFrame()
for (label, par, ls) in settings
    r0  = solve_case(update(par; sigma = 0.0, beta = 0.0); ls = ls)
    r12 = solve_case(par; ls = ls)
    append!(oat, DataFrame([merge((setting = label, line_scale = ls, scenario = "S0"), r0),
                            merge((setting = label, line_scale = ls, scenario = "S12"), r12)]))
end
oat.VoF_k = [s == "S12" ? oat.total_cost_k[i-1] - oat.total_cost_k[i] : 0.0 for (i, s) in enumerate(oat.scenario)]
oat.NSB_k = [s == "S12" ? oat.VoF_k[i] - btm_cost_k(1.0) : 0.0 for (i, s) in enumerate(oat.scenario)]
CSV.write(joinpath(RESULTS, "sensitivity_congestion.csv"), oat)
println("\nCongestion sensitivity (S0 = inflexible, S12 = sigma 100%, beta 100%):")
show(select(oat, :setting, :scenario, :total_cost_k, :VoF_k, :NSB_k, :curtail_hours, :curtailed_MWh,
            :gas_MWh, :dc_lost_MWh, :co2_t, :lmp_dc_mean, :lmp_dc_max); allrows = true, allcols = true)
println()
println("Note: line ratings below 75% are infeasible without load shedding (no VOLL).")

# ---- 3. value of lost compute V (S0, S2, S7, S12 at each V)
vcases = [("S0", 0.0, 0.0), ("S2", 0.25, 0.25), ("S7", 0.5, 0.5), ("S12", 1.0, 1.0)]
vsweep = DataFrame()
for V in V_SWEEP
    rows = [merge((V = V, scenario = n), solve_case(update(BASE; v_dc = V, sigma = s, beta = b))) for (n, s, b) in vcases]
    d = DataFrame(rows)
    d.VoF_k = d.total_cost_k[1] .- d.total_cost_k
    d.NSB_k = d.VoF_k .- btm_cost_k.(d.beta)
    d.avoided_lost_k = V .* (d.dc_lost_MWh[1] .- d.dc_lost_MWh) ./ 1e3   # part of VoF from avoided lost compute
    append!(vsweep, d)
end
CSV.write(joinpath(RESULTS, "sensitivity_V.csv"), vsweep)
println("\nValue of lost compute V (\$/MWh):")
show(select(vsweep, :V, :scenario, :total_cost_k, :VoF_k, :NSB_k, :avoided_lost_k, :dc_lost_MWh, :lmp_dc_mean, :lmp_dc_max);
     allrows = true, allcols = true)
println()
share = combine(groupby(vsweep, :V), [:scenario, :VoF_k] => ((s, x) -> x[s .== "S2"][1] / x[s .== "S12"][1]) => :first_25pct_share)
println("\nShare of the full-flexibility saving captured by the first 25% (S2 / S12):"); show(share; allrows = true); println()

# ---- figure: VoF, NSB, lost compute and curtailment vs sigma (one line per beta), plus VoF vs V
plots = []
for (col, ylab) in ((:VoF_k, "Value of flexibility (\$k / week)"), (:NSB_k, "Net social benefit (\$k / week)"),
                    (:dc_lost_MWh, "Lost compute (MWh)"), (:curtailed_MWh, "Curtailed DC load (MWh)"))
    p = plot(; xlabel = "Curtailment share sigma", ylabel = ylab, legend = col == :VoF_k ? :bottomright : false,
             xticks = (SIGMAS, string.(Int.(100 .* SIGMAS), "%")))
    for b in BETAS
        d = grid[grid.beta .== b, :]
        plot!(p, d.sigma, d[!, col]; marker = :circle, label = "beta = $(Int(100b))% ($(round(Int, b * DC_PEAK)) MW)")
    end
    push!(plots, p)
end
pv = plot(; xlabel = "Value of lost compute V (\$/MWh)", ylabel = "VoF (\$k / week)", legend = :topleft)
for (n, _, _) in vcases[2:end]
    d = vsweep[vsweep.scenario .== n, :]
    plot!(pv, d.V, d.VoF_k; marker = :circle, label = n)
end
push!(plots, pv)
fig = plot(plots...; layout = (3, 2), size = (1100, 1100), left_margin = 5Plots.mm)
savefig(fig, joinpath(RESULTS, "sensitivity_sigma_beta.png"))
isinteractive() && display(fig)

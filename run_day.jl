# =============================================================================
#  run_day.jl — one day, hourly: the ISO schedules the data center (Regime C)
#
#  Highest-load day of the week. Inflexible baseline S0 (sigma = 0, beta = 0)
#  vs the fully flexible case S12 (sigma = 100%, beta = 100% of the 100 MW load).
#  Run:  julia run_day.jl
# =============================================================================

include(joinpath(@__DIR__, "setup.jl"))

hrs = day_hours(HOT_DAY)
D, L, CF, tm = D_week[:, hrs], L_week[hrs], CF_week[hrs], time_week[hrs]
@printf("Day: %s\n", Date(tm[1]))

cases = [("S0", INFLEX), ("S12", BASE)]
res = Dict(name => solve_dispatch(net, D, L, CF, par; line_scale = LINE_SCALE, solver = SOLVER) for (name, par) in cases)

day_summary = DataFrame([merge((scenario = name,), summarize(res[name])) for (name, _) in cases])
show(select(day_summary, :scenario, :sigma, :beta, :total_cost_k, :grid_gen_cost_k, :curtail_hours, :curtailed_MWh,
            :curtail_pct, :gas_MWh, :dc_lost_MWh, :co2_t, :lmp_dc_mean, :lmp_dc_max); allcols = true)
println()
CSV.write(joinpath(RESULTS, "day_summary.csv"), day_summary)

println("\nHourly schedule, S12:")
hourly = hourly_table(res["S12"], tm)
show(hourly; allrows = true)
println()
CSV.write(joinpath(RESULTS, "day_hourly_S0.csv"), hourly_table(res["S0"], tm))
CSV.write(joinpath(RESULTS, "day_hourly_S12.csv"), hourly)

# ---- figure: DC supply by hour and LMP at the DC bus
h = 0:23
plots = []
for (name, title) in (("S0", "S0 inflexible (sigma = 0, beta = 0)"), ("S12", "S12 flexible (sigma = 100%, beta = 100%)"))
    r = res[name]
    p = areaplot(h, hcat(r.imp, r.gas, r.solar, r.lost);
                 label = ["Grid import" "BTM gas" "BTM solar" "Lost compute"],
                 fillalpha = 0.75, linewidth = 0, ylabel = "MW", legend = :bottomleft, xlims = (-0.5, 23.5),
                 title = "Data-center supply, $(title)")
    plot!(p, h, r.L; color = :black, linewidth = 1.5, label = "DC load")
    push!(plots, p)
end
p_lmp = plot(h, max.(res["S0"].lmp_dc, 1.0); label = "S0", xlabel = "Hour", ylabel = "\$/MWh",
             title = "LMP at the data-center bus (log scale)", marker = :circle, yscale = :log10)
plot!(p_lmp, h, max.(res["S12"].lmp_dc, 1.0); label = "S12", marker = :square)
hline!(p_lmp, [gas_cost(BASE)]; linestyle = :dash, color = :gray, label = "BTM gas cost")
fig = plot(plots..., p_lmp; layout = (3, 1), size = (900, 950), left_margin = 5Plots.mm)
savefig(fig, joinpath(RESULTS, "day_dispatch.png"))
isinteractive() && display(fig)

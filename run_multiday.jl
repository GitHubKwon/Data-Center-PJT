# =============================================================================
#  run_multiday.jl — several days in one optimization (default: all 7 days)
#
#  Same two cases as run_day.jl (S0 inflexible, S12 fully flexible), solved over
#  a multi-day horizon so ramp limits link the days. Change NDAYS to use fewer days.
#  Run:  julia run_multiday.jl
# =============================================================================

include(joinpath(@__DIR__, "setup.jl"))

NDAYS = 7
hrs = 1:24 * NDAYS
D, L, CF, tm = D_week[:, hrs], L_week[hrs], CF_week[hrs], time_week[hrs]

cases = [("S0", INFLEX), ("S12", BASE)]
res = Dict(name => solve_dispatch(net, D, L, CF, par; line_scale = LINE_SCALE, solver = SOLVER) for (name, par) in cases)
multi = DataFrame([merge((scenario = name,), summarize(res[name])) for (name, _) in cases])
multi.VoF_k = multi.total_cost_k[1] .- multi.total_cost_k
println("\n$(NDAYS)-day results:")
show(select(multi, :scenario, :sigma, :beta, :total_cost_k, :grid_gen_cost_k, :btm_fuel_cost_k, :lost_compute_cost_k,
            :VoF_k, :curtail_hours, :curtailed_MWh, :curtail_pct, :gas_MWh, :dc_lost_MWh, :co2_t,
            :lmp_dc_mean, :lmp_dc_max, :solve_s); allcols = true)
println()
CSV.write(joinpath(RESULTS, "multiday_summary.csv"), multi)
CSV.write(joinpath(RESULTS, "multiday_hourly_S0.csv"), hourly_table(res["S0"], tm))
CSV.write(joinpath(RESULTS, "multiday_hourly_S12.csv"), hourly_table(res["S12"], tm))

# curtailment by day (S12) and lost compute by day (S0)
byday = DataFrame(day = Date.(tm[1:24:end]), load_factor = day_factor[1:NDAYS],
                  curtailed_MWh_S12 = [sum(res["S12"].curt[day_hours(d)]) for d in 1:NDAYS],
                  lost_MWh_S0 = [sum(res["S0"].lost[day_hours(d)]) for d in 1:NDAYS],
                  lmp_dc_max_S12 = [maximum(res["S12"].lmp_dc[day_hours(d)]) for d in 1:NDAYS])
println("\nBy day:")
show(byday; allrows = true)
println()

hh = 1:length(hrs)
r = res["S12"]
p1 = areaplot(hh, hcat(r.imp, r.gas, r.solar, r.lost);
              label = ["Grid import" "BTM gas" "BTM solar" "Lost compute"], fillalpha = 0.75,
              linewidth = 0, ylabel = "MW", legend = :bottomleft, title = "Data-center supply, S12 (fully flexible)")
plot!(p1, hh, L; color = :black, linewidth = 1, label = "DC load")
p2 = plot(hh, max.(res["S0"].lmp_dc, 1.0); label = "LMP at DC bus, S0", ylabel = "\$/MWh", yscale = :log10,
          xlabel = "Hour from $(Date(tm[1]))", title = "Price at the data-center bus (log scale)")
plot!(p2, hh, max.(r.lmp_dc, 1.0); label = "LMP at DC bus, S12")
hline!(p2, [gas_cost(BASE)]; linestyle = :dash, color = :gray, label = "BTM gas cost")
fig = plot(p1, p2; layout = (2, 1), size = (1100, 750), left_margin = 5Plots.mm)
savefig(fig, joinpath(RESULTS, "multiday_dispatch.png"))
isinteractive() && display(fig)

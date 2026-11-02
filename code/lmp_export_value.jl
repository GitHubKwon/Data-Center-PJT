using CSV, DataFrames, Dates, Printf

# Project root = parent of the code/ folder, so the script runs from any directory
const ROOT = normpath(joinpath(@__DIR__, ".."))

# Real-price check: value of BTM export capacity at the ASHBURN 35 kV node
# (PJM Data Miner 2, rt_hrl_lmps, pnode 34886139, Dominion zone).
# Price-taker: the data center is paid the observed LMP and does not move it,
# so these values are upper bounds on what export capacity is worth.
#
# Run from the PJT folder:  julia code/lmp_export_value.jl

const FMT   = dateformat"m/d/yyyy I:M:S p"
const K_GAS = 150.0                 # fixed BTM gas capacity (MW)
const K_SOL = 50.0                  # fixed BTM solar capacity (MW)
const C_GAS = [35.0, 45.0, 75.0]    # BTM gas marginal cost scenarios ($/MWh)

function load_lmp(file)
    df = CSV.read(file, DataFrame;
                  types = Dict(:datetime_beginning_utc => String,
                               :datetime_beginning_ept => String,
                               :row_is_current => String))
    df = df[lowercase.(df.row_is_current) .== "true", :]
    df.utc = DateTime.(df.datetime_beginning_utc, FMT)
    df.ept = DateTime.(df.datetime_beginning_ept, FMT)
    df = unique(sort(df, :utc), :utc)   # key on UTC; drops the repeated DST hour
    return df
end

# ---------------------------------------------------------------------------
# 1. Annual break-even export cost (Nov 1, 2025 - Nov 1, 2026)
# ---------------------------------------------------------------------------
yr = load_lmp(joinpath(ROOT, "data", "pjm_lmp", "LMP_Ashburn_Year Long Data.csv"))
@printf("Year file: %d hours, %s to %s (EPT)\n", nrow(yr), yr.ept[1], yr.ept[end])

summary = DataFrame(c_gas = Float64[], hours_profitable = Int[],
                    margin_per_MW_yr = Float64[], margin_no_congestion = Float64[],
                    congestion_share = Float64[])
for c in C_GAS
    m  = max.(yr.total_lmp_rt .- c, 0.0)                               # with congestion
    m0 = max.(yr.total_lmp_rt .- yr.congestion_price_rt .- c, 0.0)     # congestion removed
    push!(summary, (c_gas = c, hours_profitable = count(>(0), m),
                    margin_per_MW_yr = round(sum(m), digits = 0),
                    margin_no_congestion = round(sum(m0), digits = 0),
                    congestion_share = round(1 - sum(m0) / sum(m), digits = 3)))
end
println("\n=== Annual export margin per MW of export capacity (Ashburn) ===")
println("(sum over hours of max(LMP - c_gas, 0); compare with annualized kappa_exp)")
println(summary)
CSV.write(joinpath(ROOT, "results", "lmp_breakeven_summary.csv"), summary)

# ---------------------------------------------------------------------------
# 2. July 15-19, 2025: export value curve using the actual load profile
# ---------------------------------------------------------------------------
jl = load_lmp(joinpath(ROOT, "data", "pjm_lmp", "LMP_ashburn_July15-21-2025.csv"))
dc = CSV.read(joinpath(ROOT, "data", "single_dc_profile.csv"), DataFrame; types = Dict(:Timestamp => String))
dc.ept = DateTime.(dc.Timestamp, dateformat"yyyy-mm-ddTHH:MM")
wk = innerjoin(jl, dc, on = :ept)
wk.surplus = K_GAS .+ K_SOL .* wk.Solar_CF .- wk.DC_Load_MW     # MW available to export
@printf("\nJuly 2025 week: %d hours matched (%s to %s), surplus %.1f-%.1f MW\n",
        nrow(wk), minimum(wk.ept), maximum(wk.ept), minimum(wk.surplus), maximum(wk.surplus))

c = 45.0
spread = max.(wk.total_lmp_rt .- c, 0.0)
curve = DataFrame(K_exp = Float64[], marginal_value_per_MW = Float64[], total_margin = Float64[])
for K in 0.0:5.0:80.0
    mv  = sum(spread .* (wk.surplus .>= K))                   # value of the K-th MW of export capacity
    tot = sum(spread .* min.(K, max.(wk.surplus, 0.0)))       # total margin with K MW of export capacity
    push!(curve, (K_exp = K, marginal_value_per_MW = round(mv, digits = 0),
                  total_margin = round(tot, digits = 0)))
end
println("\n=== July 15-19, 2025 export value curve (c_gas = \$45/MWh) ===")
println("Optimal K_exp is where marginal_value_per_MW falls to kappa_exp for this period.")
println(curve)
CSV.write(joinpath(ROOT, "results", "lmp_july2025_export_value.csv"), curve)

println("\nSaved: results/lmp_breakeven_summary.csv, results/lmp_july2025_export_value.csv")

# =============================================================================
#  setup.jl — data and base-case settings shared by every run_*.jl script
#
#  Data source is chosen automatically:
#   * HW2 mode  (if data/case_118_uc.xlsx exists): HW2 generators (costs, ramps)
#     and HW2 24-hour bus loads; extra days = HW2 day x daily factor.
#   * HW1 mode  (otherwise): MATPOWER network + HW1 gen_cost_curves.csv; bus loads
#     = MATPOWER Pd x hourly PJM Dominion load shape x KAPPA (kept at 1.0; not a sensitivity).
#  Network: data/case_118_network.m (your homework file) if present, else data/case118.m.
# =============================================================================

include(joinpath(@__DIR__, "src", "dc_dispatch.jl"))
using Plots

const DATA    = joinpath(@__DIR__, "data")
const RESULTS = joinpath(@__DIR__, "results")
mkpath(RESULTS)

# Solver: HiGHS (no licence needed). For Gurobi: using Gurobi; const SOLVER = Gurobi.Optimizer
const SOLVER = HiGHS.Optimizer

# ---- data files
const NETWORK_FILE = isfile(joinpath(DATA, "case_118_network.m")) ?
                     joinpath(DATA, "case_118_network.m") : joinpath(DATA, "case118.m")
const UC_FILE    = joinpath(DATA, "case_118_uc.xlsx")       # HW2 unit-commitment data
const COST_FILE  = joinpath(DATA, "gen_cost_curves.csv")    # HW1 offer curves
const USE_HW2    = isfile(UC_FILE)

# ---- scenario knobs
const LINE_SCALE = 0.75                  # F^max = 75% of rating (HW2 convention: s = 75 * s_pu)
KAPPA            = USE_HW2 ? 1.0 : 1.2   # non-DC load multiplier (1.0 for HW2; not varied)
const DC_PEAK    = 100.0                 # L̂: data-center load (MW)
const FLAT_LOAD  = true                  # true: flat DC_PEAK every hour; false: hourly shape from single_dc_profile.csv
const V_LOST     = 3000.0                # V: value of lost compute ($/MWh), GPU rental revenue per MWh
const V_SWEEP    = [600.0, 1000.0, 3000.0, 6600.0]   # V sensitivity: spot H100 ... AWS B200 rental prices
const DC_BUS     = USE_HW2 ? 15 : 86     # one location: bus 15, congested west (highest no-DC LMP at peak)
const SIGMAS     = [0.0, 0.25, 0.5, 1.0] # σ: max curtailable share of hourly DC load
const BETAS      = [0.25, 0.5, 1.0]      # β: BTM gas capacity / peak DC load (S1-S12; S0 is sigma = beta = 0)
BASE   = DCParams(dc_bus = DC_BUS, sigma = 1.0, beta = 1.0, dc_peak = DC_PEAK, btm_pv = 50.0,
                  v_dc = V_LOST)                                                                 # S12: fully flexible
INFLEX = update(BASE; sigma = 0.0, beta = 0.0)                                                   # S0: inflexible baseline

# ---- hourly data: 7 days, 15-21 July 2025
time_week, L_raw, CF_week, sigma_week = load_week(DATA)
L_week = FLAT_LOAD ? fill(DC_PEAK, length(L_raw)) : L_raw .* (DC_PEAK / maximum(L_raw))
day_factor = daily_factors(sigma_week)           # daily peak / weekly peak
net0 = load_matpower(NETWORK_FILE)
if USE_HW2
    net, D24 = load_hw2_uc(net0, UC_FILE)
    D_week = KAPPA .* multiday(D24, day_factor)  # N x 168
else
    net = hw1_gen_costs(net0, COST_FILE)
    D_week = KAPPA .* (net.Pd * sigma_week')     # N x 168
end
const HOT_DAY = argmax(day_factor)               # highest-load day of the week
day_hours(d) = (d - 1) * 24 + 1:d * 24

@printf("Data: %s | network %s | %d buses, %d generators, %d lines | load %.0f-%.0f MW (kappa %.2f)\n",
        USE_HW2 ? "HW2 UC data" : "HW1 data", basename(NETWORK_FILE), length(net.Pd), length(net.gen_bus),
        length(net.fr), minimum(sum(D_week, dims = 1)), maximum(sum(D_week, dims = 1)), KAPPA)
@printf("Data center: bus %d, %s %.0f MW, BTM gas %.0f MW (beta %.2f), sigma %.2f, solar %.0f MW, gas cost %.2f \$/MWh, V %.0f \$/MWh\n",
        BASE.dc_bus, FLAT_LOAD ? "flat" : "peak", DC_PEAK, btm_capacity(BASE), BASE.beta, BASE.sigma, BASE.btm_pv,
        gas_cost(BASE), BASE.v_dc)

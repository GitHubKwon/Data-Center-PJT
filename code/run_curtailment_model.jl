# =============================================================================
# run_curtailment_model.jl
#
# Flexible (curtailable) vs. firm grid connection for a data center with a fixed
# behind-the-meter (BTM) plant. Multi-period DC-OPF LP on the IEEE 14-bus case,
# 168 hours (July 15-21, 2025), following Section 3 of the report:
#
#   min  sum_ij kappa_up_ij * dF_ij + sum_t ( sum_g c_util_gt p_gt + c_gas g_gas_t )
#
#   Regime A  islanded bridge   : p_imp = 0
#   Regime B  firm connection   : p_imp = L_DC, BTM plant off (backup only)
#   Regime C  flexible connection: c_t >= L_DC_t - p_imp_t,  sum c_t <= alpha * sum L_DC_t
#
# Run from the PJT folder:   julia code/run_curtailment_model.jl
# Packages:  ] add JuMP Gurobi CSV DataFrames      (or HiGHS instead of Gurobi)
# To use HiGHS instead of Gurobi, set SOLVER below to "HiGHS" (or set the
# environment variable PJT_SOLVER=HiGHS before starting Julia).
# =============================================================================

const SOLVER = get(ENV, "PJT_SOLVER", "Gurobi")      # "Gurobi" or "HiGHS"

using JuMP, CSV, DataFrames, Dates, Printf, Statistics
if SOLVER == "Gurobi"
    using Gurobi
else
    using HiGHS
end

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RES  = joinpath(ROOT, "results")
mkpath(RES)

# Shared Gurobi environment so the license banner prints once, not once per solve
const GRB_ENV = SOLVER == "Gurobi" ? Gurobi.Env() : nothing
new_model() = SOLVER == "Gurobi" ? Model(() -> Gurobi.Optimizer(GRB_ENV)) : Model(HiGHS.Optimizer)

# =============================================================================
# 0. ASSUMPTIONS AND SOURCED CONSTANTS  (change these, not the code below)
# =============================================================================
const HOURS_PER_YEAR = 8760.0
const T_START = DateTime(2025, 7, 15, 0)        # first modeled hour (EPT)
const NT      = 168                             # one week, hourly
const WEEK    = NT / HOURS_PER_YEAR             # share of a year in the horizon (168/8760)

# --- Fuel and generator costs ($/MWh) ---------------------------------------
const GAS_PRICE = 3.94      # $/MMBtu, EIA STEO Oct 2026, electric power sector, 2026
const HR_CC_FLEET = 7.535   # MMBtu/MWh, EIA Electric Power Annual Table 8.2 (2025), NG combined cycle
const HR_GT_FLEET = 10.975  # MMBtu/MWh, EIA Electric Power Annual Table 8.2 (2025), NG gas turbine
const HR_CT_NEW   = 9.717   # MMBtu/MWh, NREL ATB 2024, F-frame CT (the BTM turbines)
const VOM_CC      = 2.17    # $/MWh, NREL ATB 2024, 2x1 F-frame CC
const VOM_CT      = 6.94    # $/MWh, NREL ATB 2024, F-frame CT
const C_CC_UTIL   = HR_CC_FLEET * GAS_PRICE + VOM_CC     # ~ $31.9/MWh, existing utility CC (bus 2)
const C_PEAKER    = HR_GT_FLEET * GAS_PRICE + VOM_CT     # ~ $50.2/MWh, existing utility peakers
const C_GAS_BASE  = 45.0                                 # BTM gas ($/MWh); 9.717*3.94+6.94 = 45.2

# Existing gas peaker placed at bus 3 (a synchronous condenser in the original case; assumed size, MW).
# No generation sits in the 69 kV area (buses 6-14), so it is a load pocket fed through the
# 4-7, 4-9 and 5-6 transformers, like a data center cluster at the edge of the grid.
const PEAKERS = Dict(3 => 100.0)

# --- Data center and BTM plant ------------------------------------------------
const DC_BUS  = 14          # bus hosting the data center (end of the 69 kV area)
const K_GAS   = 150.0       # MW, fixed BTM gas
const K_SOL   = 50.0        # MW, fixed BTM solar
const K_BESS  = 100.0       # MWh, fixed BTM battery  ([Y] in the report; assumed 4-hour, 25 MW)
const GAMMA   = 0.25        # battery power-to-energy ratio (1/h)
const ETA_CH  = 0.92        # charge efficiency  (0.92^2 ~ 85% round trip, NREL ATB)
const ETA_DIS = 0.92        # discharge efficiency
const K_IMP   = 150.0       # MW, substation import limit (>= ~118 MW peak load)
const E_GAS   = HR_CT_NEW * 119.0 / 2204.62   # tCO2/MWh = 9.717 MMBtu/MWh * 119 lb/MMBtu (~0.524)

# --- Curtailment budget (Norris et al., 2025: 0.25%, 0.5%, 1%) ---------------
const ALPHA_BASE = 0.01
# How BTM solar is treated (it costs nothing to run, so the choice matters a lot):
# true (default) = solar is a load reduction in every regime. Regime B: p_imp = L_DC - g_sol
#                  (gas and battery stay off as backup); Regime C: c >= L_DC - g_sol - p_imp,
#                  so the curtailment budget applies only to what solar cannot cover.
#                  B and C then differ only in flexibility, not in who gets the free solar.
# false          = Section 3.5 as currently written: solar is off in Regime B, and in Regime C
#                  any solar use counts as curtailment and uses up most of the budget.
#                  Run as a sensitivity below.
const SOLAR_AS_LOAD_REDUCTION = true

# --- Transmission upgrade cost (MISO MTEP26 Cost Estimation Workbook, 2026 $) --
# Table 4.1-3: reconductor cost per circuit-mile ($M/mile); Table 3.1-5: rating after (MVA)
const RECOND_M_PER_MILE = Dict(69 => 0.4229628745588701, 138 => 0.498761783170008)
const RATING_MVA        = Dict(69 => 140.0, 138 => 394.0)
# Table 2.3-7: power transformer cost, 138 kV / 69 kV ($/MVA)
const XFMR_PER_MVA = 5999.696183250002
# MISO Guide Table 5.2 (Kentucky, nearest MISO state to Virginia):
# depreciation 2.5% + return 8.07% + expense 3.08% = annual cost as a share of capital cost
const FCR = 0.025 + 0.0807 + 0.0308
# Typical series reactance used to back out line length from the case's per-unit reactance
const X_OHM_PER_MILE = Dict(69 => 0.75, 138 => 0.80)
# Multiplier on kappa_up (1 = MISO exploratory cost; Northern Virginia is likely higher)
const KAPPA_MULT_BASE = 1.0

# --- Firm transmission service (PJM OATT, DOM zone NITS rate, CAPS handbook) ----
const NITS_PER_MW_YR = 75876.81
const KAPPA_FIRM     = NITS_PER_MW_YR * WEEK          # ~ $1,455/MW for the 168 h week
const INCLUDE_NITS_BASE = false                         # if true: + KAPPA_FIRM * peak import

# --- Network congestion knob (scales every rateA) ---------------------------
# At the PGLib ratings (1.0) the 14-bus case can serve the data center without any upgrade.
# The baseline uses 70% of rateA so the data center's arrival creates congestion, standing in
# for a capacity-constrained area; 1.0, 0.85 and 0.6 are run as sensitivities.
const FMAX_SCALE_BASE = 0.7

# =============================================================================
# 1. INPUT DATA
# =============================================================================
const PJM_FMT = dateformat"m/d/yyyy I:M:S p"

function read_pjm(file)
    df = CSV.read(file, DataFrame; types = Dict(:datetime_beginning_utc => String,
                                                :datetime_beginning_ept => String))
    if hasproperty(df, :row_is_current)
        df = df[lowercase.(string.(df.row_is_current)) .== "true", :]
    end
    df.utc = DateTime.(df.datetime_beginning_utc, PJM_FMT)
    df.ept = DateTime.(df.datetime_beginning_ept, PJM_FMT)
    df = unique(sort(df, :utc), :utc)
    hrs = T_START .+ Hour.(0:NT-1)
    df = df[in.(df.ept, Ref(Set(hrs))), :]
    nrow(df) == NT || error("$(basename(file)): expected $NT hours from $T_START, found $(nrow(df))")
    return sort(df, :ept)
end

# 1a. PJM real-time LMP at ASHBURN 35 kV (system energy price = PJM-wide marginal energy cost)
lmp = read_pjm(joinpath(ROOT, "data", "pjm_lmp", "LMP_ashburn_July14-23-2025.csv"))
# 1b. PJM metered hourly load, Dominion zone -> hourly shape for regional demand D_i,t
dom = read_pjm(joinpath(ROOT, "data", "pjm_load", "DOM_hrl_load_metered_July14-22-2025.csv"))
load_shape = dom.mw ./ maximum(dom.mw)
# 1c. Data center load and solar capacity factor (Section 4.2)
dc = CSV.read(joinpath(ROOT, "data", "single_dc_profile.csv"), DataFrame; types = Dict(:Timestamp => String))
dc.ept = DateTime.(dc.Timestamp, dateformat"yyyy-mm-ddTHH:MM")
dc = sort(dc[in.(dc.ept, Ref(Set(lmp.ept))), :], :ept)
nrow(dc) == NT || error("single_dc_profile.csv does not cover the 168 modeled hours")
@assert dc.ept == lmp.ept == dom.ept
L_DC   = Float64.(dc.DC_Load_MW)
CF_SOL = Float64.(dc.Solar_CF)
P_MKT  = Float64.(lmp.system_energy_price_rt)     # $/MWh, price of power bought from PJM at bus 1

# 1d. IEEE 14-bus network (PGLib-OPF)
function matpower_block(txt, name)
    m = match(Regex("mpc\\.$name\\s*=\\s*\\[(.*?)\\];", "s"), txt)
    rows = Vector{Vector{Float64}}()
    for line in split(m.captures[1], '\n')
        line = strip(replace(split(line, '%')[1], ';' => ""))
        isempty(line) && continue
        push!(rows, parse.(Float64, split(line)))
    end
    return rows
end
case     = read(joinpath(ROOT, "data", "pglib_opf_case14_ieee.m"), String)
BASE_MVA = parse(Float64, match(r"mpc\.baseMVA\s*=\s*([\d.]+)", case).captures[1])
busrows  = matpower_block(case, "bus")
genrows  = matpower_block(case, "gen")
brrows   = matpower_block(case, "branch")

BUSES   = Int.(getindex.(busrows, 1))
REF_BUS = Int(busrows[findfirst(r -> r[2] == 3, busrows)][1])
D_PEAK  = Dict(Int(r[1]) => r[3] for r in busrows)       # MW, treated as the weekly peak
D = Dict((i, t) => D_PEAK[i] * load_shape[t] for i in BUSES, t in 1:NT)

# Branches: buses 1-5 = 138 kV, buses 6-14 = 69 kV; zero-resistance branches are transformers
kv(i) = i <= 5 ? 138 : 69
struct Branch
    f::Int; t::Int; x::Float64; rate::Float64; kind::String; kv::Int; miles::Float64
    capex_per_MW::Float64; kappa::Float64
end
BR = Branch[]
for r in brrows
    f, t, x, rate, res = Int(r[1]), Int(r[2]), r[4], r[6], r[3]
    if res == 0.0
        kind, kvb, miles = "transformer", 138, 0.0
        capex = XFMR_PER_MVA
    else
        kind, kvb = "line", kv(f)
        zbase = kvb^2 / BASE_MVA
        miles = x * zbase / X_OHM_PER_MILE[kvb]
        capex = RECOND_M_PER_MILE[kvb] * 1e6 / RATING_MVA[kvb] * miles
    end
    push!(BR, Branch(f, t, x, rate, kind, kvb, miles, capex, capex * FCR * WEEK))
end
LINES = [(b.f, b.t) for b in BR]
BRD   = Dict((b.f, b.t) => b for b in BR)

# Utility generators: bus 1 = purchases from PJM at the hourly system energy price,
# bus 2 = existing CC, synchronous-condenser buses = existing gas peakers
struct Gen; bus::Int; pmax::Float64; name::String; end
GENS = Gen[]
for r in genrows
    b = Int(r[1])
    if b == REF_BUS
        push!(GENS, Gen(b, r[9], "PJM market (bus 1)"))
    elseif r[9] > 0
        push!(GENS, Gen(b, r[9], "Utility CC (bus $b)"))
    elseif haskey(PEAKERS, b)
        push!(GENS, Gen(b, PEAKERS[b], "Peaker (bus $b)"))
    end
end
gen_cost(g::Gen, t) = g.bus == REF_BUS ? P_MKT[t] : (g.name[1:7] == "Utility" ? C_CC_UTIL : C_PEAKER)

# Save the inputs so every number in the report can be traced
CSV.write(joinpath(RES, "curtailment_kappa_up.csv"),
    DataFrame(from = [b.f for b in BR], to = [b.t for b in BR], type = [b.kind for b in BR],
              kV = [b.kv for b in BR], rateA_MW = [b.rate for b in BR],
              est_length_miles = round.([b.miles for b in BR], digits = 1),
              capital_cost_per_MW = round.([b.capex_per_MW for b in BR], digits = 0),
              kappa_up_per_MW_week = round.([b.kappa for b in BR], digits = 2)))
CSV.write(joinpath(RES, "curtailment_inputs_hourly.csv"),
    DataFrame(hour = 1:NT, ept = lmp.ept, DOM_load_MW = dom.mw, load_shape = round.(load_shape, digits = 4),
              regional_load_MW = round.([sum(D[i, t] for i in BUSES) for t in 1:NT], digits = 2),
              DC_load_MW = L_DC, solar_CF = CF_SOL, PJM_energy_price = P_MKT,
              Ashburn_total_LMP = lmp.total_lmp_rt))

println("=== Inputs ===")
@printf("Hours: %s to %s (EPT)\n", lmp.ept[1], lmp.ept[end])
@printf("Regional load: %.1f-%.1f MW (DOM shape x 259 MW peak); DC load %.1f-%.1f MW at bus %d\n",
        minimum(sum(D[i, t] for i in BUSES) for t in 1:NT), maximum(sum(D[i, t] for i in BUSES) for t in 1:NT),
        minimum(L_DC), maximum(L_DC), DC_BUS)
@printf("PJM energy price (bus 1): mean \$%.2f, max \$%.2f/MWh; CC \$%.2f; peakers \$%.2f; BTM gas \$%.2f\n",
        mean(P_MKT), maximum(P_MKT), C_CC_UTIL, C_PEAKER, C_GAS_BASE)
@printf("kappa_up: \$%.2f-\$%.2f per MW per week (x%.1f); NITS \$%.2f per MW per week\n",
        minimum(b.kappa for b in BR), maximum(b.kappa for b in BR), KAPPA_MULT_BASE, KAPPA_FIRM)
@printf("Network: rateA x %.2f; data center at bus %d; solar as load reduction = %s; alpha = %.2f%%\n",
        FMAX_SCALE_BASE, DC_BUS, SOLAR_AS_LOAD_REDUCTION, 100 * ALPHA_BASE)

# =============================================================================
# 2. MODEL  (Section 3 of the report)
# =============================================================================
Base.@kwdef struct Scenario
    name::String       = "baseline"
    regime::Char       = 'C'
    c_gas::Float64     = C_GAS_BASE
    alpha::Float64     = ALPHA_BASE
    kappa_mult::Float64 = KAPPA_MULT_BASE
    fmax_scale::Float64 = FMAX_SCALE_BASE
    co2_cap::Float64   = Inf             # tCO2 over the horizon (Inf = no cap)
    include_nits::Bool = INCLUDE_NITS_BASE
    solar_lr::Bool     = SOLAR_AS_LOAD_REDUCTION
end

function solve(s::Scenario)
    m = new_model(); set_silent(m)
    T = 1:NT
    G = 1:length(GENS)

    @variable(m, dF[l in LINES] >= 0)                                 # MW of added thermal capacity
    @variable(m, 0 <= p[g in G, t in T] <= GENS[g].pmax)              # utility generation
    @variable(m, 0 <= g_gas[t in T] <= K_GAS)                         # BTM gas
    @variable(m, 0 <= g_sol[t in T] <= CF_SOL[t] * K_SOL)             # BTM solar
    @variable(m, 0 <= p_ch[t in T] <= GAMMA * K_BESS)                 # battery charge
    @variable(m, 0 <= p_dis[t in T] <= GAMMA * K_BESS)                # battery discharge
    @variable(m, 0 <= soc[t in T] <= K_BESS)                          # state of charge
    @variable(m, 0 <= p_imp[t in T] <= K_IMP)                         # grid import
    @variable(m, c[t in T] >= 0)                                      # curtailed import
    @variable(m, f[l in LINES, t in T])                               # branch flow
    @variable(m, theta[i in BUSES, t in T])                           # voltage angle
    @variable(m, K_firm >= 0)                                         # peak import (for NITS)

    # Objective: upgrades + utility energy + BTM fuel (+ optional NITS on peak import)
    @objective(m, Min,
        sum(s.kappa_mult * BRD[l].kappa * dF[l] for l in LINES) +
        sum(gen_cost(GENS[g], t) * p[g, t] for g in G, t in T) +
        sum(s.c_gas * g_gas[t] for t in T) +
        (s.include_nits ? KAPPA_FIRM * K_firm : 0.0))

    # Nodal balance (dual = LMP):  generation + inflow - outflow = load + DC import
    @constraint(m, bal[i in BUSES, t in T],
        sum(p[g, t] for g in G if GENS[g].bus == i) +
        sum(f[l, t] for l in LINES if l[2] == i) - sum(f[l, t] for l in LINES if l[1] == i) ==
        D[i, t] + (i == DC_BUS ? p_imp[t] : 0.0))

    # Data center microgrid balance
    @constraint(m, [t in T], g_gas[t] + g_sol[t] + p_dis[t] - p_ch[t] + p_imp[t] == L_DC[t])

    # DC power flow, thermal limits with upgrades, angle limits
    @constraint(m, [l in LINES, t in T], f[l, t] == BASE_MVA / BRD[l].x * (theta[l[1], t] - theta[l[2], t]))
    @constraint(m, [l in LINES, t in T],  f[l, t] <= s.fmax_scale * BRD[l].rate + dF[l])
    @constraint(m, [l in LINES, t in T], -f[l, t] <= s.fmax_scale * BRD[l].rate + dF[l])
    @constraint(m, [l in LINES, t in T], -pi / 6 <= theta[l[1], t] - theta[l[2], t] <= pi / 6)
    @constraint(m, [t in T], theta[REF_BUS, t] == 0)

    # Battery: starts and ends half full
    @constraint(m, [t in T], soc[t] == (t == 1 ? 0.5 * K_BESS : soc[t-1]) + ETA_CH * p_ch[t] - p_dis[t] / ETA_DIS)
    @constraint(m, soc[NT] == 0.5 * K_BESS)

    # Optional carbon cap on BTM gas
    isfinite(s.co2_cap) && @constraint(m, sum(E_GAS * g_gas[t] for t in T) <= s.co2_cap)

    # Peak import (only priced if include_nits)
    @constraint(m, [t in T], K_firm >= p_imp[t])

    # Regimes
    if s.regime == 'A'
        @constraint(m, [t in T], p_imp[t] == 0)
        @constraint(m, [t in T], c[t] == 0)
    elseif s.regime == 'B'
        if s.solar_lr
            @constraint(m, [t in T], p_imp[t] == L_DC[t] - g_sol[t])   # solar used, rest firm import
        else
            @constraint(m, [t in T], p_imp[t] == L_DC[t])
            @constraint(m, [t in T], g_sol[t] == 0)
        end
        @constraint(m, [t in T], g_gas[t] == 0)
        @constraint(m, [t in T], p_ch[t] == 0)
        @constraint(m, [t in T], p_dis[t] == 0)
        @constraint(m, [t in T], c[t] == 0)
    else   # 'C'
        if s.solar_lr
            @constraint(m, [t in T], c[t] >= L_DC[t] - g_sol[t] - p_imp[t])
        else
            @constraint(m, [t in T], c[t] >= L_DC[t] - p_imp[t])
        end
        @constraint(m, sum(c[t] for t in T) <= s.alpha * sum(L_DC))
    end

    optimize!(m)
    if !is_solved_and_feasible(m)
        return (scenario = s, status = string(termination_status(m)), model = m)
    end

    upg  = sum(s.kappa_mult * BRD[l].kappa * value(dF[l]) for l in LINES)
    util = sum(gen_cost(GENS[g], t) * value(p[g, t]) for g in G, t in T)
    fuel = sum(s.c_gas * value(g_gas[t]) for t in T)
    lmp_dc = [dual(bal[DC_BUS, t]) for t in T]
    peak_imp = maximum(Vector(value.(p_imp)))
    return (scenario = s, status = "OPTIMAL", model = m,
            total = objective_value(m), upgrade = upg, utility = util, fuel = fuel,
            nits = s.include_nits ? KAPPA_FIRM * value(K_firm) : 0.0,
            dF = Dict(l => value(dF[l]) for l in LINES),
            dF_total = sum(value(dF[l]) for l in LINES),
            imp = Vector(value.(p_imp)), gas = Vector(value.(g_gas)), sol = Vector(value.(g_sol)),
            ch = Vector(value.(p_ch)), dis = Vector(value.(p_dis)), soc = Vector(value.(soc)),
            curt = Vector(value.(c)), lmp_dc = lmp_dc,
            pk = [value(p[g, t]) for g in G, t in T],
            peak_imp = peak_imp,
            gas_MWh = sum(value.(g_gas)), co2 = E_GAS * sum(value.(g_gas)),
            curt_MWh = sum(value.(c)))
end

# =============================================================================
# 3. EXPERIMENTS  (Section 4.3)
# =============================================================================
rows = DataFrame()
function record!(r, sweep, value_label)
    s = r.scenario
    push!(rows, (sweep = sweep, setting = value_label, regime = string(s.regime), status = r.status,
        c_gas = s.c_gas, alpha = s.alpha, kappa_mult = s.kappa_mult, fmax_scale = s.fmax_scale,
        co2_cap_t = isfinite(s.co2_cap) ? round(s.co2_cap, digits = 1) : missing, include_nits = s.include_nits,
        solar_as_load_reduction = s.solar_lr,
        total_cost = r.status == "OPTIMAL" ? round(r.total, digits = 2) : missing,
        upgrade_cost = r.status == "OPTIMAL" ? round(r.upgrade, digits = 2) : missing,
        utility_cost = r.status == "OPTIMAL" ? round(r.utility, digits = 2) : missing,
        btm_fuel_cost = r.status == "OPTIMAL" ? round(r.fuel, digits = 2) : missing,
        nits_cost = r.status == "OPTIMAL" ? round(r.nits, digits = 2) : missing,
        upgrades_MW = r.status == "OPTIMAL" ? round(r.dF_total, digits = 2) : missing,
        peak_import_MW = r.status == "OPTIMAL" ? round(r.peak_imp, digits = 2) : missing,
        curtailed_MWh = r.status == "OPTIMAL" ? round(r.curt_MWh, digits = 2) : missing,
        btm_gas_MWh = r.status == "OPTIMAL" ? round(r.gas_MWh, digits = 2) : missing,
        btm_co2_t = r.status == "OPTIMAL" ? round(r.co2, digits = 2) : missing,
        avg_LMP_dc_bus = r.status == "OPTIMAL" ? round(mean(r.lmp_dc), digits = 2) : missing); promote = true)
    return r
end
run!(sweep, label; kw...) = record!(solve(Scenario(; kw...)), sweep, label)

println("\n=== Solving ===")
# Baseline: all three regimes
base = Dict(R => run!("baseline", "baseline"; regime = R) for R in ('A', 'B', 'C'))

# One-at-a-time sweeps around the baseline (Regimes B and C; A where it changes)
for k in (0.5, 2.0, 5.0, 10.0, 20.0, 50.0), R in ('B', 'C')
    run!("kappa_up multiplier", "x$(k)"; regime = R, kappa_mult = k)
end
for cg in (35.0, 75.0), R in ('A', 'C')
    run!("BTM gas cost", "\$$(Int(cg))/MWh"; regime = R, c_gas = cg)
end
for fs in (1.0, 0.85, 0.6), R in ('B', 'C')
    run!("congestion (rateA scale)", "$(fs)"; regime = R, fmax_scale = fs)
end
for a in (0.0025, 0.005)
    run!("curtailment budget alpha", "$(100a)%"; regime = 'C', alpha = a)
end
for R in ('A', 'C')                     # carbon cap: half of the uncapped BTM emissions, then zero gas
    r0 = base[R]
    r0.status == "OPTIMAL" || continue
    run!("carbon cap", "50% of uncapped"; regime = R, co2_cap = 0.5 * r0.co2)
    run!("carbon cap", "zero gas"; regime = R, co2_cap = 0.0)
end
for R in ('B', 'C')
    run!("firm service charge (NITS)", "on"; regime = R, include_nits = true)
end
for R in ('B', 'C')
    run!("solar treatment", "Section 3.5 as written"; regime = R, solar_lr = false)
end

CSV.write(joinpath(RES, "curtailment_scenarios.csv"), rows)

# Baseline hourly results and upgrades
hourly = DataFrame(hour = 1:NT, ept = lmp.ept, DC_load_MW = L_DC)
for R in ('A', 'B', 'C')
    r = base[R]; r.status == "OPTIMAL" || continue
    hourly[!, "$(R)_import"]    = round.(r.imp, digits = 3)
    hourly[!, "$(R)_btm_gas"]   = round.(r.gas, digits = 3)
    hourly[!, "$(R)_btm_solar"] = round.(r.sol, digits = 3)
    hourly[!, "$(R)_batt_net"]  = round.(r.dis .- r.ch, digits = 3)
    hourly[!, "$(R)_soc"]       = round.(r.soc, digits = 3)
    hourly[!, "$(R)_curtailed"] = round.(r.curt, digits = 3)
    hourly[!, "$(R)_LMP_dc_bus"] = round.(r.lmp_dc, digits = 2)
end
CSV.write(joinpath(RES, "curtailment_hourly_baseline.csv"), hourly)

upg = DataFrame(from = [l[1] for l in LINES], to = [l[2] for l in LINES],
                rateA_MW = [BRD[l].rate for l in LINES],
                kappa_up_per_MW_week = round.([BRD[l].kappa for l in LINES], digits = 2))
for R in ('A', 'B', 'C')
    base[R].status == "OPTIMAL" && (upg[!, "$(R)_dF_MW"] = round.([base[R].dF[l] for l in LINES], digits = 3))
end
CSV.write(joinpath(RES, "curtailment_upgrades_baseline.csv"), upg)

# =============================================================================
# 4. SUMMARY
# =============================================================================
println("\n=== Baseline (one week, July 15-21, 2025) ===")
@printf("%-8s %14s %12s %14s %12s %10s %10s %10s\n",
        "Regime", "Total (\$)", "Upgrade", "Utility", "BTM fuel", "dF (MW)", "Curt MWh", "Gas MWh")
for R in ('A', 'B', 'C')
    r = base[R]
    if r.status != "OPTIMAL"
        println("$R  $(r.status)"); continue
    end
    @printf("%-8s %14.2f %12.2f %14.2f %12.2f %10.2f %10.2f %10.2f\n",
            R, r.total, r.upgrade, r.utility, r.fuel, r.dF_total, r.curt_MWh, r.gas_MWh)
end
if base['B'].status == "OPTIMAL" && base['C'].status == "OPTIMAL"
    @printf("\nFlexible connection saves \$%.2f vs. firm (B - C) this week.\n", base['B'].total - base['C'].total)
    println("Upgraded branches (MW):")
    for l in LINES
        b, cc = base['B'].dF[l], base['C'].dF[l]
        (b > 1e-6 || cc > 1e-6) && @printf("  %2d-%-2d  B %7.2f   C %7.2f   (\$%.2f/MW-week)\n", l[1], l[2], b, cc, BRD[l].kappa)
    end
end
println("\nWrote results/curtailment_scenarios.csv, curtailment_hourly_baseline.csv,")
println("      curtailment_upgrades_baseline.csv, curtailment_kappa_up.csv, curtailment_inputs_hourly.csv")

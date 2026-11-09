# =============================================================================
#  dc_dispatch.jl
#  15.S09 project — Regime C: the ISO schedules a curtailable data center
#
#  Multi-period economic dispatch (DC-OPF, B-theta, LP) on a MATPOWER system.
#  A data center (DC) is injected at one bus as a load the ISO may curtail.
#  Scenario parameters (proposal notation):
#    sigma  maximum share of the hourly DC load the ISO may curtail onto BTM gas
#    beta   BTM gas capacity as a share of the peak DC load
#  Curtailed load is served by the DC's behind-the-meter (BTM) gas; load the
#  grid and BTM cannot serve is lost compute, priced at V. BTM solar offsets
#  DC load. BTM generation does not export to the grid, and there is no load
#  shedding of other customers (no VOLL): those cases are infeasible.
#
#  Usage: include("src/dc_dispatch.jl"); see setup.jl and the run_*.jl scripts.
# =============================================================================

using JuMP, HiGHS, CSV, DataFrames, Dates, Printf, XLSX

# -----------------------------------------------------------------------------
#  Data types
# -----------------------------------------------------------------------------

"""Transmission network in physical units (MW, \$/MWh). Indices are positions 1..N."""
struct Network
    baseMVA::Float64
    bus_ids::Vector{Int}        # original MATPOWER bus numbers
    ref::Int                    # index of the reference bus
    Pd::Vector{Float64}         # base active load at each bus (MW)
    gen_bus::Vector{Int}        # bus index of each grid generator
    pmax::Vector{Float64}       # (MW)
    cost::Vector{Float64}       # energy cost ($/MWh)
    ramp_up::Vector{Float64}    # MW/h (Inf = no limit)
    ramp_dn::Vector{Float64}    # MW/h (Inf = no limit)
    fr::Vector{Int}             # from-bus index of each line
    to::Vector{Int}             # to-bus index of each line
    b::Vector{Float64}          # series susceptance (p.u.)
    fmax::Vector{Float64}       # thermal limit (MW)
end

"""Data-center parameters (defaults = fully flexible case S12: sigma = 1, beta = 1)."""
Base.@kwdef struct DCParams
    dc_bus::Int          = 15        # MATPOWER bus number of the data center
    sigma::Float64       = 1.0       # σ: max share of hourly DC load the ISO may curtail onto BTM gas
    beta::Float64        = 1.0       # β: BTM gas capacity as a share of peak DC load
    dc_peak::Float64     = 100.0     # L̂: peak DC load (MW)
    btm_pv::Float64      = 50.0      # K^sol: BTM solar capacity (MW)
    imp_limit::Float64   = Inf       # K̄^imp: substation import limit (MW); Inf -> dc_peak
    heat_rate::Float64   = 9.717     # MMBtu/MWh   NREL ATB 2024, NG F-frame CT
    fuel_price::Float64  = 3.50      # $/MMBtu     assumption
    vom::Float64         = 6.94      # $/MWh       NREL ATB 2024
    co2_rate::Float64    = 0.05398   # tCO2/MMBtu  (119 lb/MMBtu)
    co2_price::Float64   = 0.0       # $/tCO2
    v_dc::Float64        = 3000.0    # V: $/MWh value of lost compute (GPU rental revenue per MWh)
    angle_max::Float64   = π / 6     # rad, bound on angle difference across a line
end

"""Marginal cost of BTM gas generation (\$/MWh)."""
gas_cost(par::DCParams) = par.heat_rate * (par.fuel_price + par.co2_price * par.co2_rate) + par.vom

# BTM capital + fixed O&M per MW per hour of the horizon: NREL ATB 2024, F-frame CT, Moderate, 2025
# (CAPEX $1,319.40/kW, FOM $25.70/kW-yr), annualized at the ATB real WACC 5.36% over 30 years.
const BTM_CRF = 0.0535848 / (1 - 1.0535848^-30)
const BTM_COST_PER_MW_HOUR = (BTM_CRF * 1319.40 + 25.70) * 1000 / 8760     # about $13.14/MW-h

"""BTM gas capacity (MW) = β L̂."""
btm_capacity(par::DCParams) = par.beta * par.dc_peak

"""Copy of `par` with some fields changed, e.g. update(par; sigma = 0.5, beta = 0.25)."""
function update(par::DCParams; kw...)
    vals = Dict{Symbol,Any}(f => getfield(par, f) for f in fieldnames(DCParams))
    for (key, v) in kw
        haskey(vals, key) || error("DCParams has no field $(key)")
        vals[key] = v
    end
    return DCParams(; vals...)
end

# -----------------------------------------------------------------------------
#  MATPOWER network
# -----------------------------------------------------------------------------

"""Read one numeric block `mpc.<name> = [ ... ];` of a MATPOWER file into a matrix."""
function _matpower_block(txt::AbstractString, name::AbstractString)
    mt = match(Regex("mpc\\.$(name)\\s*=\\s*\\[(.*?)\\]", "s"), txt)
    mt === nothing && error("block mpc.$(name) not found")
    rows = Vector{Vector{Float64}}()
    for line in split(mt.captures[1], '\n')
        line = strip(split(line, '%')[1])
        for row in split(line, ';')
            row = strip(row)
            isempty(row) && continue
            push!(rows, parse.(Float64, split(row)))
        end
    end
    return permutedims(reduce(hcat, rows))
end

"""Linear cost (\$/MWh) from a MATPOWER gencost row (quadratic -> average marginal cost)."""
function _linear_cost(row, pmin, pmax)
    model, n = Int(row[1]), Int(row[4])
    if model == 2
        n >= 2 || return 0.0
        c  = row[5:4+n]
        c2 = n >= 3 ? c[end-2] : 0.0
        return c[end-1] + c2 * (max(pmin, 0.0) + pmax)
    else
        p, fc = row[5:2:4+2n], row[6:2:4+2n]
        return (fc[end] - fc[1]) / max(p[end] - p[1], 1e-9)
    end
end

"""Mask of the buses in the largest connected island."""
function _largest_island(nbus::Int, fr::Vector{Int}, to::Vector{Int})
    adj = [Int[] for _ in 1:nbus]
    for (a, b) in zip(fr, to)
        push!(adj[a], b); push!(adj[b], a)
    end
    comp, ncomp = zeros(Int, nbus), 0
    for s in 1:nbus
        comp[s] == 0 || continue
        ncomp += 1; comp[s] = ncomp; stack = [s]
        while !isempty(stack)
            v = pop!(stack)
            for w in adj[v]
                if comp[w] == 0
                    comp[w] = ncomp; push!(stack, w)
                end
            end
        end
    end
    sizes = [count(==(c), comp) for c in 1:ncomp]
    return comp .== argmax(sizes)
end

"""
    load_matpower(path) -> Network

MATPOWER (.m) reader in MW and \$/MWh, mirroring the Homework data cleaning
(PowerModels make_basic_network): drops isolated / out-of-service elements and
small islands; fills missing thermal limits (rateA = 0) with the PowerModels rule;
b = x / (r^2 + x^2) (taps ignored, as in the homework). Pmin is not used: the
economic-dispatch LP has no commitment. Ramp limits default to Inf.
"""
function load_matpower(path::AbstractString)
    txt = read(path, String)
    baseMVA = parse(Float64, strip(match(r"mpc\.baseMVA\s*=\s*([^;]+);", txt).captures[1]))
    bus, gen, gencost, branch = (_matpower_block(txt, nm) for nm in ("bus", "gen", "gencost", "branch"))

    idx0  = Dict(Int(b) => i for (i, b) in enumerate(bus[:, 1]))
    busok = bus[:, 2] .!= 4
    brok  = [branch[l, 11] > 0 && busok[idx0[Int(branch[l, 1])]] && busok[idx0[Int(branch[l, 2])]]
             for l in 1:size(branch, 1)]
    br     = branch[brok, :]
    island = _largest_island(size(bus, 1), [idx0[Int(f)] for f in br[:, 1]], [idx0[Int(t)] for t in br[:, 2]])
    keepbus = busok .& island
    br  = br[[keepbus[idx0[Int(br[l, 1])]] && keepbus[idx0[Int(br[l, 2])]] for l in 1:size(br, 1)], :]
    bus = bus[keepbus, :]

    bus_ids = Int.(bus[:, 1])
    idx = Dict(b => i for (i, b) in enumerate(bus_ids))
    ref = findfirst(==(3), Int.(bus[:, 2]))
    ref === nothing && (ref = 1)

    keep = findall([gen[i, 8] > 0 && gen[i, 9] > 0 && haskey(idx, Int(gen[i, 1])) for i in 1:size(gen, 1)])
    pmax = gen[keep, 9]
    cost = [_linear_cost(gencost[i, :], gen[i, 10], gen[i, 9]) for i in keep]

    fr = [idx[Int(f)] for f in br[:, 1]]
    to = [idx[Int(t)] for t in br[:, 2]]
    r, x = br[:, 3], br[:, 4]
    z2 = r .^ 2 .+ x .^ 2
    b  = [z2[l] < 1e-12 ? 1e4 : x[l] / z2[l] for l in eachindex(x)]
    has_ang = size(br, 2) >= 13
    angmin = has_ang ? deg2rad.(br[:, 12]) : zeros(size(br, 1))
    angmax = has_ang ? deg2rad.(br[:, 13]) : zeros(size(br, 1))
    unset  = (angmin .== 0) .& (angmax .== 0)
    angmin = ifelse.((angmin .<= -pi / 2) .| unset, -1.0472, angmin)
    angmax = ifelse.((angmax .>= pi / 2) .| unset, 1.0472, angmax)
    thmax  = max.(abs.(angmin), abs.(angmax))
    vfr, vto = bus[fr, 12], bus[to, 12]
    cmax   = sqrt.(vfr .^ 2 .+ vto .^ 2 .- 2 .* vfr .* vto .* cos.(thmax))
    ymag   = [z2[l] < 1e-12 ? 1e4 : 1 / sqrt(z2[l]) for l in eachindex(z2)]
    fmax   = [br[l, 6] > 0 ? br[l, 6] : baseMVA * ymag[l] * max(vfr[l], vto[l]) * cmax[l] for l in eachindex(fr)]

    G = length(keep)
    return Network(baseMVA, bus_ids, ref, bus[:, 3], [idx[Int(g)] for g in gen[keep, 1]],
                   pmax, cost, fill(Inf, G), fill(Inf, G), fr, to, b, fmax)
end

"""Copy of `net` with new generator data (used by the HW1 / HW2 loaders)."""
with_gens(net::Network, gen_bus, pmax, cost, ru, rd) =
    Network(net.baseMVA, net.bus_ids, net.ref, net.Pd, gen_bus, pmax, cost, ru, rd,
            net.fr, net.to, net.b, net.fmax)

"""
    hw1_gen_costs(net, cost_path) -> Network

Homework 1 generator costs: picks a piecewise-linear offer curve for each generator
from gen_cost_curves.csv (same row-picking rule as HW1's load_opf_data) and uses its
average cost (HW1's get_gen_cost_vector).
"""
function hw1_gen_costs(net::Network, cost_path::AbstractString)
    cc = CSV.read(cost_path, DataFrame)
    n  = size(cc, 1)
    cost = zeros(length(net.pmax))
    j = 1
    for (i, pm) in enumerate(net.pmax)
        pg, v = Float64[], Float64[]
        while isempty(pg) || pg[end] < 1e-2
            j  = j % n
            pg = Float64.(collect(skipmissing(cc[j, 15:2:end-3])))
            v  = Float64.(collect(skipmissing(cc[j, 14:2:end-3])))
            j += 10
        end
        seg = pg .* (pm / sum(pg))
        cost[i] = sum(seg) > 0 ? sum(seg .* v) / sum(seg) : sum(v) / length(v)
    end
    return with_gens(net, net.gen_bus, net.pmax, cost, net.ramp_up, net.ramp_dn)
end

"""
    load_hw2_uc(net, uc_path) -> (net_uc, D24)

Homework 2 unit-commitment data (case_118_uc.xlsx): replaces the generators with the
54 UC units (bus, ramp up/down, variable cost c_pg, Pmax; Pmin and min up/down are not
used in the dispatch LP) and returns the 24-hour load at every bus, D24 (N x 24, MW).
"""
function load_hw2_uc(net::Network, uc_path::AbstractString)
    xf = XLSX.readxlsx(uc_path)
    g  = xf["Generators"]["B2:T55"][:, [1, 2, 3, 6, 7, 10, 14, 15, 19, 18]]
    # columns: bus, ru, rd, min_up, min_dn, c_su, c_on, c_pg, pmin, pmax
    bus_of = Dict(b => i for (i, b) in enumerate(net.bus_ids))
    gbus = [bus_of[Int(v)] for v in g[:, 1]]
    ru, rd = Float64.(g[:, 2]), Float64.(g[:, 3])
    cpg, pmax = Float64.(g[:, 8]), Float64.(g[:, 10])

    loadbus  = xf["Loads"]["A2:B92"]                  # load id ("Load1", ...), bus
    loadtime = xf["LoadTime"]["A2:D2185"][:, 1:3]     # load id, time, pd (dNominal, MW)
    bus_of_load = Dict(string(loadbus[r, 1]) => Int(loadbus[r, 2]) for r in 1:size(loadbus, 1))
    T24 = maximum(Int.(loadtime[:, 2]))
    D24 = zeros(length(net.bus_ids), T24)
    for r in 1:size(loadtime, 1)
        n = bus_of[bus_of_load[string(loadtime[r, 1])]]
        D24[n, Int(loadtime[r, 2])] += Float64(loadtime[r, 3])
    end
    # unit check: HW2 loads and capacities should be in MW
    if maximum(sum(D24, dims = 1)) < 0.1 * sum(net.Pd)
        @info "HW2 loads look like per-unit; multiplying by baseMVA"
        D24 .*= net.baseMVA
    end
    if sum(pmax) < 0.1 * sum(net.Pd)
        @info "HW2 capacities look like per-unit; multiplying by baseMVA"
        pmax .*= net.baseMVA; ru .*= net.baseMVA; rd .*= net.baseMVA
    end
    return with_gens(net, gbus, pmax, cpg, ru, rd), D24
end

# -----------------------------------------------------------------------------
#  Hourly profiles
# -----------------------------------------------------------------------------

const PJM_FMT = dateformat"m/d/yyyy I:M:S p"     # e.g. 7/14/2025 4:00:00 PM
const DC_FMT  = dateformat"yyyy-mm-ddTHH:MM"     # e.g. 2025-07-15T00:00

"""
    load_week(datadir) -> (time, L, CF, sigma)

168 hours, 15-21 July 2025 (EPT):
  L, CF  data-center load (MW) and solar capacity factor (single_dc_profile.csv)
  sigma  system load shape = PJM Dominion-zone load / its weekly maximum
"""
function load_week(datadir::AbstractString)
    dc  = CSV.read(joinpath(datadir, "single_dc_profile.csv"), DataFrame; types = Dict(:Timestamp => String))
    dom = CSV.read(joinpath(datadir, "DOM_hrl_load_metered_July14-22-2025.csv"), DataFrame;
                   types = Dict(:datetime_beginning_ept => String, :datetime_beginning_utc => String))
    time = DateTime.(dc.Timestamp, DC_FMT)
    ept  = DateTime.(dom.datetime_beginning_ept, PJM_FMT)
    row  = Dict(t => i for (i, t) in enumerate(ept))
    all(haskey(row, t) for t in time) || error("DOM load does not cover the DC profile hours")
    dom_mw = [Float64(dom.mw[row[t]]) for t in time]
    return time, Float64.(dc.DC_Load_MW), Float64.(dc.Solar_CF), dom_mw ./ maximum(dom_mw)
end

"""Daily scale factors (daily peak / largest daily peak) from an hourly shape."""
daily_factors(sigma::Vector{Float64}) =
    (pk = [maximum(sigma[(d - 1) * 24 + 1:d * 24]) for d in 1:length(sigma) ÷ 24]; pk ./ maximum(pk))

"""Repeat a 24-h load matrix over several days, scaled by the daily factors."""
multiday(D24::Matrix{Float64}, factors::Vector{Float64}) = reduce(hcat, [D24 .* f for f in factors])

# -----------------------------------------------------------------------------
#  Optimization model  (equation numbers follow Section 3 of the proposal)
# -----------------------------------------------------------------------------

"""
    solve_dispatch(net, D, L, CF, par; line_scale = 1.0, solver = HiGHS.Optimizer, silent = true)

Regime C economic dispatch over T = size(D, 2) hours.
  D      N x T non-DC load (MW)        L, CF  DC load (MW) and solar capacity factor (length T)
  par    DCParams (DC bus, sigma, beta, costs)
The ISO chooses utility generation, flows and the DC's grid withdrawal p^imp. Curtailed
withdrawal c = L - g^sol - p^imp is covered by BTM gas (at most sigma * L per hour) or lost
(u, at V). Returns a NamedTuple; errors if the LP is infeasible (no load shedding of other
customers is modeled).
"""
function solve_dispatch(net::Network, D::Matrix{Float64}, L::Vector{Float64}, CF::Vector{Float64},
                        par::DCParams; line_scale::Real = 1.0, solver = HiGHS.Optimizer, silent::Bool = true)
    tstart = time()
    N, G, NL, T = length(net.Pd), length(net.gen_bus), length(net.fr), size(D, 2)
    length(L) == T && length(CF) == T || error("L and CF must have length size(D, 2)")
    k = findfirst(==(par.dc_bus), net.bus_ids)
    k === nothing && error("DC bus $(par.dc_bus) not in the network")
    gens_at, lines_in, lines_out = [Int[] for _ in 1:N], [Int[] for _ in 1:N], [Int[] for _ in 1:N]
    for (i, n) in enumerate(net.gen_bus); push!(gens_at[n], i); end
    for l in 1:NL
        push!(lines_out[net.fr[l]], l); push!(lines_in[net.to[l]], l)
    end
    fmax  = line_scale .* net.fmax                       # F^max: 75% of rating in the base case
    gbar  = btm_capacity(par)                            # β L̂
    kimp  = isfinite(par.imp_limit) ? par.imp_limit : par.dc_peak
    cG    = gas_cost(par)

    m = Model(solver)
    silent && set_silent(m)
    N * T > 50_000 && set_string_names_on_creation(m, false)
    if solver == HiGHS.Optimizer
        set_optimizer_attribute(m, "dual_feasibility_tolerance", 1e-6)
        set_optimizer_attribute(m, "primal_feasibility_tolerance", 1e-6)
    end

    # ---- variables
    @variable(m, 0 <= pg[i = 1:G, t = 1:T] <= net.pmax[i])          # p_{g,t}          (12)
    @variable(m, theta[1:N, 1:T])                                   # θ_{i,t}
    @variable(m, -fmax[l] <= f[l = 1:NL, t = 1:T] <= fmax[l])       # f_{ij,t}         (10)
    @variable(m, 0 <= imp[t = 1:T] <= kimp)                         # p^imp_{d,t}      (8)
    @variable(m, 0 <= ggas[t = 1:T] <= gbar)                        # g^gas_{d,t}      (6)
    @variable(m, 0 <= gsol[t = 1:T] <= par.btm_pv * CF[t])          # g^sol_{d,t}      (7)
    @variable(m, u[t = 1:T] >= 0)                                   # u_{d,t}          (13)
    @expression(m, curt[t = 1:T], L[t] - gsol[t] - imp[t])          # c_{d,t}          (4)

    # ---- (2) nodal power balance, dual = LMP
    @expression(m, inj[n = 1:N, t = 1:T],
        sum(pg[i, t] for i in gens_at[n]) + sum(f[l, t] for l in lines_in[n])
        - sum(f[l, t] for l in lines_out[n]))
    for t in 1:T
        add_to_expression!(inj[k, t], -imp[t])
    end
    @constraint(m, balance[n = 1:N, t = 1:T], inj[n, t] == D[n, t])
    # ---- (3) data-center balance; (5) curtailment covered by BTM gas is at most σ L
    @constraint(m, dcbal[t = 1:T], ggas[t] + gsol[t] + imp[t] + u[t] == L[t])
    @constraint(m, curtcap[t = 1:T], curt[t] - u[t] <= par.sigma * L[t])
    # ---- (9) DC power flow; (11) angle bounds and reference bus
    @constraint(m, flowdef[l = 1:NL, t = 1:T],
        f[l, t] == net.baseMVA * net.b[l] * (theta[net.fr[l], t] - theta[net.to[l], t]))
    @constraint(m, [l = 1:NL, t = 1:T], -par.angle_max <= theta[net.fr[l], t] - theta[net.to[l], t] <= par.angle_max)
    @constraint(m, [t = 1:T], theta[net.ref, t] == 0)
    # ---- (14) generator ramp limits
    for i in 1:G, t in 2:T
        isfinite(net.ramp_up[i]) && @constraint(m, pg[i, t] - pg[i, t-1] <= net.ramp_up[i])
        isfinite(net.ramp_dn[i]) && @constraint(m, pg[i, t-1] - pg[i, t] <= net.ramp_dn[i])
    end

    # ---- (1) objective: total system operating cost ($)
    @objective(m, Min,
        sum(net.cost[i] * pg[i, t] for i in 1:G, t in 1:T)
        + sum(cG * ggas[t] + par.v_dc * u[t] for t in 1:T))

    tbuild = time()
    optimize!(m)
    tsolve = time()
    st = termination_status(m)
    st == MOI.OPTIMAL || error("solve failed: $(st) (sigma = $(par.sigma), beta = $(par.beta); " *
                              "an infeasible LP means the grid cannot serve the non-DC load)")

    pg_v, imp_v, g_v, y_v, u_v = value.(pg), value.(imp), value.(ggas), value.(gsol), value.(u)
    lmp  = dual.(balance)
    curt = L .- y_v .- imp_v                             # c_{d,t} = g^gas + u
    return (; objective = objective_value(m),
            grid_cost = sum(net.cost[i] * pg_v[i, t] for i in 1:G, t in 1:T),
            gas_cost = cG * sum(g_v), gas_mc = cG,
            pg = pg_v, flow = value.(f), lmp, lmp_dc = lmp[k, :],
            L, solar = y_v, imp = imp_v, curt, gas = g_v, lost = u_v,
            co2_t = par.heat_rate * par.co2_rate * sum(g_v),
            dc_bill = sum(lmp[k, :] .* imp_v) + cG * sum(g_v),
            par, line_scale, build_time = tbuild - tstart, solve_time = tsolve - tbuild, model = m)
end

# -----------------------------------------------------------------------------
#  Reporting
# -----------------------------------------------------------------------------

"""One-row summary of a solved case."""
function summarize(res)
    T = length(res.L)
    return (dc_bus = res.par.dc_bus, sigma = res.par.sigma, beta = res.par.beta,
            btm_gas_MW = btm_capacity(res.par), btm_pv_MW = res.par.btm_pv, gas_cost_MWh = res.gas_mc,
            total_cost_k = res.objective / 1e3, grid_gen_cost_k = res.grid_cost / 1e3,
            btm_fuel_cost_k = res.gas_cost / 1e3, lost_compute_cost_k = res.par.v_dc * sum(res.lost) / 1e3,
            curtail_hours = count(>(1e-6), res.curt), curtailed_MWh = sum(res.curt),
            curtail_pct = 100 * sum(res.curt) / sum(res.L),
            gas_MWh = sum(res.gas), dc_lost_MWh = sum(res.lost), co2_t = res.co2_t,
            lmp_dc_mean = sum(res.lmp_dc) / T, lmp_dc_max = maximum(res.lmp_dc),
            lmp_sys_mean = sum(res.lmp) / length(res.lmp),
            dc_bill_k = res.dc_bill / 1e3, solve_s = res.solve_time)
end

_clean(v) = round.(v; digits = 6) .+ 0.0        # drop solver noise and "-0.0"

"""Hourly data-center results as a DataFrame."""
hourly_table(res, time) =
    DataFrame(time = time, dc_load_MW = _clean(res.L), solar_MW = _clean(res.solar),
              grid_import_MW = _clean(res.imp), curtailed_MW = _clean(res.curt),
              btm_gas_MW = _clean(res.gas), dc_lost_MW = _clean(res.lost), lmp_dc = _clean(res.lmp_dc))

# -----------------------------------------------------------------------------
#  Scenario grid (sigma x beta) and the value of flexibility
# -----------------------------------------------------------------------------

"""
    scenario_grid(net, D, L, CF, base; sigmas, betas, kwargs...) -> DataFrame

Solve the baseline S0 (sigma = beta = 0) and every (sigma, beta) pair, numbered S1, S2, ...
with beta in the outer loop (S1-S4: beta = 25%, S5-S8: 50%, S9-S12: 100%). Adds the metrics of Section 3:
  VoF  = Z(S0) - Z(sigma, beta)                        value of flexibility (\$k)
  dVoF = VoF(sigma_k, beta) - VoF(sigma_{k-1}, beta)    marginal value of flexibility (\$k)
  v    = VoF / curtailed MWh                            value per curtailed MWh (\$/MWh)
  NSB  = VoF - C_btm * beta * peak                      net social benefit (\$k), with C_btm in \$/MW per hour
"""
function scenario_grid(net::Network, D::Matrix{Float64}, L::Vector{Float64}, CF::Vector{Float64},
                       base::DCParams; sigmas = [0.0, 0.25, 0.5, 1.0], betas = [0.25, 0.5, 1.0],
                       c_btm_hour::Real = BTM_COST_PER_MW_HOUR, kwargs...)
    cases = [(0.0, 0.0); [(s, b) for b in betas for s in sigmas]]
    df = DataFrame([summarize(solve_dispatch(net, D, L, CF, update(base; sigma = s, beta = b); kwargs...))
                    for (s, b) in cases])
    z0 = df.total_cost_k[1]
    df.VoF_k = z0 .- df.total_cost_k
    df.dVoF_k = [i == 1 ? 0.0 : (i == 2 || df.beta[i] != df.beta[i-1]) ? df.VoF_k[i] : df.VoF_k[i] - df.VoF_k[i-1]
                 for i in 1:nrow(df)]
    df.v_per_MWh = [c > 1e-6 ? 1e3 * x / c : 0.0 for (x, c) in zip(df.VoF_k, df.curtailed_MWh)]
    df.NSB_k = df.VoF_k .- c_btm_hour .* length(L) .* df.btm_gas_MW ./ 1e3
    insertcols!(df, 1, :scenario => ["S$(i)" for i in 0:nrow(df)-1])
    return df
end

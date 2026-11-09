# The Value of Data Center Flexibility

**15.S09 Optimization of Energy Systems, Fall 2026 — John Kwon and Diana Parnian**

How much does a curtailable data center with on-site (behind-the-meter, BTM) gas generation lower total system cost at a congested transmission bus, and is faster grid access worth the cost of flexibility for the data center?

The model is an hourly, multi-period DC optimal power flow (linear program) in Julia/JuMP, solved with HiGHS. The ISO dispatches utility generators, curtails up to a share σ of the data center's load, and the data center covers curtailed load with its own gas plant (size β × load). Unserved compute is priced at V.

## Repository contents

| Path | Contents |
|---|---|
| `toy_example.jl` | Proposal toy example: 3 buses × 2 hours, scenarios S0–S12, plus a V sweep |
| `setup.jl` | Settings and data loading for the full study (all assumptions at the top) |
| `src/dc_dispatch.jl` | The model (`solve_dispatch`), scenario grid (`scenario_grid`) and reporting |
| `run_sensitivity.jl` | Full study: scenario grid S0–S12, congestion sensitivity, V sensitivity |
| `run_day.jl`, `run_multiday.jl` | Hourly dispatch of S0 vs S12 for the peak day / the full week |
| `install_packages.jl` | Installs the required Julia packages |
| `data/` | Input data (sources below) |
| `results/` | Outputs used in the proposal (CSV and PNG) |
| `Data_Center_BTM_Project_Report v4.docx` | Project proposal |

## How to run

Requires Julia 1.9 or newer (https://julialang.org/downloads/).

```bash
git clone https://github.com/GitHubKwon/Data-Center-PJT.git
cd Data-Center-PJT
julia install_packages.jl        # once: JuMP, HiGHS, CSV, DataFrames, XLSX, Plots
julia toy_example.jl             # seconds  -> results/toy_example.csv, results/toy_V_sensitivity.csv
julia run_sensitivity.jl         # ~10 min  -> results/sensitivity_*.csv, results/sensitivity_sigma_beta.png
julia run_day.jl                 # optional -> results/day_*.csv, results/day_dispatch.png
julia run_multiday.jl            # optional -> results/multiday_*.csv, results/multiday_dispatch.png
```

In VS Code (with the Julia extension), open the repository folder, start a Julia REPL and run e.g. `include("toy_example.jl")`. Restart the REPL between scripts.

## Scenarios

S0 is the inflexible baseline (σ = β = 0). S1–S12 cross σ ∈ {0, 25, 50, 100%} with β ∈ {25, 50, 100%} of the 100 MW load:

| | σ = 0% | 25% | 50% | 100% |
|---|---|---|---|---|
| β = 25% (25 MW) | S1 | S2 | S3 | S4 |
| β = 50% (50 MW) | S5 | S6 | S7 | S8 |
| β = 100% (100 MW) | S9 | S10 | S11 | S12 |

## Base-case assumptions (`setup.jl`)

| Setting | Value | Meaning |
|---|---|---|
| `DC_BUS` | 15 | Data center location (most congested bus of the IEEE 118-bus HW2 system) |
| `DC_PEAK`, `FLAT_LOAD` | 100 MW, flat | Data center load |
| `V_LOST` | \$3,000/MWh | Value of lost compute (GPU rental revenue per MWh) |
| `LINE_SCALE` | 0.75 | Line limits = 75% of the network file's `rate_a` (HW2 template convention) |
| BTM gas cost | \$40.95/MWh | NREL ATB 2024 F-frame CT: 9.717 MMBtu/MWh × \$3.50/MMBtu (assumed) + \$6.94 VOM |
| BTM capital | ≈ \$115,000/MW-yr | NREL ATB 2024 CAPEX and fixed O&M, 5.36% real WACC, 30 years |
| BTM solar | 50 MW | Capacity factor from Open-Meteo, Ashburn, VA |

Sensitivities: V ∈ {\$600, \$1,000, \$3,000, \$6,600}/MWh (`V_SWEEP`) and line ratings at 75%, 80% and 90%.

## Data sources

| File | Source |
|---|---|
| `case_118_network.m`, `case_118_uc.xlsx` | IEEE 118-bus system, 15.S09 Homework 2 |
| `DOM_hrl_load_metered_July14-22-2025.csv` | PJM Data Miner 2, Dominion zone metered load: https://dataminer2.pjm.com/feed/hrl_load_metered |
| `single_dc_profile.csv` | Data center profile, July 15–21, 2025 (solar capacity factor from Open-Meteo, Ashburn, VA: https://open-meteo.com/en/docs/historical-weather-api) |
| `LMP_ashburn_July14-23-2025.csv` | PJM real-time hourly LMPs, Ashburn node (fallback analysis): https://dataminer2.pjm.com/feed/rt_hrl_lmps |

## Main results (7-day July week, bus 15, V = \$3,000/MWh)

| | S0 | S2 (25/25) | S7 (50/50) | S12 (100/100) |
|---|---|---|---|---|
| Weekly system cost | \$10.97M | \$9.94M | \$9.51M | \$9.28M |
| Value of flexibility (VoF) | – | \$1.03M | \$1.46M | \$1.69M |
| Net social benefit (after BTM capital) | – | \$0.98M | \$1.35M | \$1.47M |
| Lost compute | 517 MWh | 195 MWh | 71 MWh | 0 |

The first 25% of flexibility delivers about 60% of the value at every V tested. The benefit depends on congestion: with line ratings at 90%, VoF falls to \$0.11M and the net social benefit turns negative.

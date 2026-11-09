# The Value of Data Center Flexibility (15.S09 project)

John Kwon and Diana Parnian. An ISO dispatch model (DC-OPF linear program, Julia/JuMP/HiGHS) that measures how much a curtailable data center with on-site gas lowers system cost at a congested bus, and whether faster grid access makes that flexibility worthwhile for the data center.

## What is in this folder

| File | What it is |
|---|---|
| `Data_Center_BTM_Project_Report v4.docx` | The proposal (current version) |
| `toy_example.jl` | Toy example for the proposal: 3 buses x 2 hours, S0-S12, plus the V sweep |
| `setup.jl` | Shared settings and data loading for the full study (bus 15, flat 100 MW, V = $3,000/MWh) |
| `src/dc_dispatch.jl` | The model: `solve_dispatch`, `scenario_grid`, reporting |
| `run_sensitivity.jl` | Full study: scenario grid S0-S12 + the two sensitivities (congestion, V) |
| `run_day.jl` | One day, hourly: S0 vs S12 on the highest-load day |
| `run_multiday.jl` | All 7 days in one optimization: S0 vs S12 |
| `install_packages.jl` | Installs the Julia packages (run once) |
| `data/` | HW2 IEEE 118-bus files, data center profile, PJM load and LMP data |
| `results/` | CSV and PNG outputs (created by the scripts) |
| `_old/` | Earlier proposal draft, kept for reference |

## Running it in VS Code

1. **Install Julia** (1.9 or newer) from https://julialang.org/downloads/ and the **Julia extension** in VS Code (Extensions panel, search "Julia").
2. **Open this folder** in VS Code: File > Open Folder > `PJTv3`.
3. **Open a Julia terminal**: View > Command Palette > "Julia: Start REPL". The terminal starts in this folder.
4. **Install packages (once):**
   ```julia
   include("install_packages.jl")
   ```
5. **Run the toy example (seconds):**
   ```julia
   include("toy_example.jl")
   ```
   Writes `results/toy_example.csv` and `results/toy_V_sensitivity.csv`. These are the numbers in Section 5 of the proposal.
6. **Run the full study (about 8-10 minutes):**
   ```julia
   include("run_sensitivity.jl")
   ```
   Writes `results/sensitivity_sigma_beta.csv` (S0-S12), `results/sensitivity_congestion.csv`, `results/sensitivity_V.csv` and `results/sensitivity_sigma_beta.png`.
7. **Optional:** `include("run_day.jl")` and `include("run_multiday.jl")` for hourly dispatch plots of S0 vs S12.

You can also run any script from a normal terminal: `julia toy_example.jl`.

Tip: run each script in a fresh REPL (Command Palette > "Julia: Restart REPL") to avoid "redefinition of constant" warnings when you include several scripts in one session.

## Changing assumptions

All settings for the full study are at the top of `setup.jl`:

| Setting | Base value | Meaning |
|---|---|---|
| `DC_BUS` | 15 | Data center location |
| `DC_PEAK`, `FLAT_LOAD` | 100.0, true | Flat 100 MW load (set `FLAT_LOAD = false` for the hourly profile) |
| `V_LOST` | 3000.0 | Value of lost compute ($/MWh) |
| `V_SWEEP` | 600, 1000, 3000, 6600 | V sensitivity |
| `SIGMAS`, `BETAS` | 0/25/50/100%, 25/50/100% | Scenario grid |
| `LINE_SCALE` | 0.75 | Line limits as a share of the rating |

In `toy_example.jl` the inputs (network, loads, `V`, `V_SWEEP`) are at the top of the file.

## Scenarios

S0 is the inflexible baseline (sigma = beta = 0). S1-S12 cross sigma in {0, 25, 50, 100%} with beta in {25, 50, 100%} of the 100 MW load: S1-S4 beta = 25%, S5-S8 beta = 50%, S9-S12 beta = 100%. S12 is fully flexible.

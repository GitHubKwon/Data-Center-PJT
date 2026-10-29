using Downloads, JSON3, DataFrames, CSV, PowerModels

# ==============================================================================
# 1. Download IEEE 14-Bus Grid Topology (with offline fallback)
# ==============================================================================
grid_file = "pglib_opf_case14_ieee.m"
try
    Downloads.download(
        "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/master/pglib_opf_case14_ieee.m",
        grid_file
    )
    println("1. Downloaded \$grid_file from GitHub (", filesize(grid_file), " bytes)")
catch
    # Offline fallback: copy directly from your local PowerModels installation
    cp(joinpath(pkgdir(PowerModels), "test", "data", "matpower", "case14.m"), grid_file; force=true)
    println("1. Copied \$grid_file from local PowerModels package (", filesize(grid_file), " bytes)")
end

# ==============================================================================
# 2. Direct Download: Official NREL ATB Cost & Performance Workbook (.xlsx)
# ==============================================================================
atb_file = "NREL_ATB_Electricity_Workbook.xlsx"
if !isfile(atb_file)
    atb_url = "https://data.openei.org/files/6006/2024_v4_Workbook_7-22-2026.xlsx"
    Downloads.download(atb_url, atb_file)
end
println("2. Verified \$atb_file (", round(filesize(atb_file) / 1e6, digits=2), " MB)")

# ==============================================================================
# 3. Direct API Pull: Real Hourly Weather & Solar for 1 Data Center (Ashburn, VA)
# ==============================================================================
weather_file = "ashburn_weather_raw.json"
weather_url = "https://archive-api.open-meteo.com/v1/archive?latitude=39.0438&longitude=-77.4874&start_date=2025-07-15&end_date=2025-07-21&hourly=temperature_2m,shortwave_radiation&timezone=America%2FNew_York"
Downloads.download(weather_url, weather_file)

raw_json   = JSON3.read(read(weather_file, String))
timestamps = String.(raw_json.hourly.time)
temp_c     = Float64.(raw_json.hourly.temperature_2m)
solar_wm2  = Float64.(raw_json.hourly.shortwave_radiation)

# Construct 100 MW IT Single Data Center Load Profile + On-Site Solar CF
T_hours    = length(timestamps)
P_IT_rated = 100.0 # 100 MW single hyperscale campus IT capacity

# Server utilization (88% baseload + 6% daytime inference peak)
u_it = [0.88 + 0.06 * max(0.0, sin(pi * (mod(t-1, 24) - 6) / 14)) for t in 1:T_hours]

# Cooling PUE driven by real Ashburn, VA hourly temperature
pue = [1.10 + 0.008 * max(0.0, temp_c[t] - 15.0) for t in 1:T_hours]

# Total Data Center Power Requirement (MW) & Solar Capacity Factor [0, 1]
L_dc_mw  = [round(P_IT_rated * u_it[t] * pue[t], digits=2) for t in 1:T_hours]
cf_solar = [round(min(1.0, solar_wm2[t] / 1000.0), digits=3) for t in 1:T_hours]

df_dc = DataFrame(
    Hour = 1:T_hours,
    Timestamp = timestamps,
    Temp_C = temp_c,
    Solar_CF = cf_solar,
    IT_Utilization = round.(u_it, digits=3),
    Cooling_PUE = round.(pue, digits=3),
    DC_Load_MW = L_dc_mw
)

dc_csv_file = "single_dc_profile.csv"
CSV.write(dc_csv_file, df_dc)
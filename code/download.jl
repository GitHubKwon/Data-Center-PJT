using Downloads, JSON3, DataFrames, CSV, PowerModels, Dates

# Project root = parent of the code/ folder, so the script runs from any directory
const ROOT = normpath(joinpath(@__DIR__, ".."))

# ==============================================================================
# 1. Download IEEE 14-Bus Grid Topology (with offline fallback)
# ==============================================================================
grid_file = joinpath(ROOT, "data", "pglib_opf_case14_ieee.m")
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
atb_file = joinpath(ROOT, "data", "NREL_ATB_Electricity_Workbook.xlsx")
if !isfile(atb_file)
    atb_url = "https://data.openei.org/files/6006/2024_v4_Workbook_7-22-2026.xlsx"
    Downloads.download(atb_url, atb_file)
end
println("2. Verified \$atb_file (", round(filesize(atb_file) / 1e6, digits=2), " MB)")

# ==============================================================================
# 3. Direct API Pull: Real Hourly Weather & Solar for 1 Data Center (Ashburn, VA)
# ==============================================================================
weather_file = joinpath(ROOT, "data", "ashburn_weather_raw.json")
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

dc_csv_file = joinpath(ROOT, "data", "single_dc_profile.csv")
CSV.write(dc_csv_file, df_dc)

# ==============================================================================
# 4. Year-long weather (Nov 1, 2025 - Oct 31, 2026) to match the PJM LMP year file
#    Requested in UTC so it joins exactly with the LMP file's datetime_beginning_utc;
#    local (EPT) time is computed for the utilization profile.
# ==============================================================================
year_weather_file = joinpath(ROOT, "data", "ashburn_weather_year_raw.json")
year_url = "https://archive-api.open-meteo.com/v1/archive?latitude=39.0438&longitude=-77.4874&start_date=2025-11-01&end_date=2026-10-31&hourly=temperature_2m,shortwave_radiation&timezone=GMT"
Downloads.download(year_url, year_weather_file)

yj      = JSON3.read(read(year_weather_file, String))
tofloat(v) = v === nothing ? NaN : Float64(v)      # the archive returns null for hours not yet available
utc_y   = DateTime.(String.(yj.hourly.time), dateformat"yyyy-mm-ddTHH:MM")
temp_y  = tofloat.(yj.hourly.temperature_2m)
solar_y = tofloat.(yj.hourly.shortwave_radiation)

ok = .!(isnan.(temp_y) .| isnan.(solar_y))
n_missing = count(!, ok)
n_missing > 0 && println("   Warning: dropped $(n_missing) hours with missing weather data")
utc_y, temp_y, solar_y = utc_y[ok], temp_y[ok], solar_y[ok]

# US Eastern time: EDT (UTC-4) from 2:00 local on the 2nd Sunday of March
# to 2:00 local on the 1st Sunday of November; EST (UTC-5) otherwise.
nth_sunday(y, m, n) = (d = Date(y, m, 1); d + Day(mod(7 - dayofweek(d), 7)) + Week(n - 1))
function eastern_local(utc::DateTime)
    y = year(utc)
    dst_start = DateTime(nth_sunday(y, 3, 2)) + Hour(7)   # 2:00 EST = 07:00 UTC
    dst_end   = DateTime(nth_sunday(y, 11, 1)) + Hour(6)  # 2:00 EDT = 06:00 UTC
    return utc + Hour(dst_start <= utc < dst_end ? -4 : -5)
end
ept_y = eastern_local.(utc_y)

# Same load model as the July profile (Section 4.2)
u_y   = [0.88 + 0.06 * max(0.0, sin(pi * (hour(t) - 6) / 14)) for t in ept_y]
pue_y = [1.10 + 0.008 * max(0.0, T - 15.0) for T in temp_y]
L_y   = round.(P_IT_rated .* u_y .* pue_y, digits=2)
cf_y  = round.(min.(1.0, solar_y ./ 1000.0), digits=3)

df_year = DataFrame(
    Hour = 1:length(utc_y),
    Timestamp_UTC = Dates.format.(utc_y, "yyyy-mm-ddTHH:MM"),
    Timestamp_EPT = Dates.format.(ept_y, "yyyy-mm-ddTHH:MM"),
    Temp_C = temp_y,
    Solar_CF = cf_y,
    IT_Utilization = round.(u_y, digits=3),
    Cooling_PUE = round.(pue_y, digits=3),
    DC_Load_MW = L_y
)
year_csv_file = joinpath(ROOT, "data", "single_dc_profile_year.csv")
CSV.write(year_csv_file, df_year)
println("4. Saved year-long profile: ", nrow(df_year), " hours (",
        df_year.Timestamp_EPT[1], " to ", df_year.Timestamp_EPT[end], " EPT), peak load ",
        maximum(L_y), " MW")

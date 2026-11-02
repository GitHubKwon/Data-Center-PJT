using Downloads, JSON3, DataFrames, CSV, PowerModels, Dates

# Project root = parent of the code/ folder, so the script runs from any directory
const ROOT = normpath(joinpath(@__DIR__, ".."))

const P_IT_rated = 100.0   # 100 MW single hyperscale campus IT capacity

# US Eastern time: EDT (UTC-4) from 2:00 local on the 2nd Sunday of March
# to 2:00 local on the 1st Sunday of November; EST (UTC-5) otherwise.
# Weather is downloaded in UTC and converted here, because Open-Meteo applies the
# time zone's *current* UTC offset to a whole request (a local-time request made in
# winter would shift July data by one hour).
nth_sunday(y, m, n) = (d = Date(y, m, 1); d + Day(mod(7 - dayofweek(d), 7)) + Week(n - 1))
function eastern_local(utc::DateTime)
    y = year(utc)
    dst_start = DateTime(nth_sunday(y, 3, 2)) + Hour(7)   # 2:00 EST = 07:00 UTC
    dst_end   = DateTime(nth_sunday(y, 11, 1)) + Hour(6)  # 2:00 EDT = 06:00 UTC
    return utc + Hour(dst_start <= utc < dst_end ? -4 : -5)
end

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
#    July 15-21, 2025 in local time (EDT). Skipped if the profile already exists,
#    so re-running this script never changes the data used in the report.
# ==============================================================================
weather_file = joinpath(ROOT, "data", "ashburn_weather_raw.json")
dc_csv_file  = joinpath(ROOT, "data", "single_dc_profile.csv")
if isfile(dc_csv_file)
    println("3. Kept existing ", dc_csv_file, " (delete it to re-download and rebuild)")
else
    weather_url = "https://archive-api.open-meteo.com/v1/archive?latitude=39.0438&longitude=-77.4874&start_date=2025-07-14&end_date=2025-07-22&hourly=temperature_2m,shortwave_radiation&timezone=GMT"
    Downloads.download(weather_url, weather_file)

    raw_json  = JSON3.read(read(weather_file, String))
    utc_w     = DateTime.(String.(raw_json.hourly.time), dateformat"yyyy-mm-ddTHH:MM")
    ept_w     = eastern_local.(utc_w)
    keep      = DateTime(2025, 7, 15, 0) .<= ept_w .<= DateTime(2025, 7, 21, 23)
    ept_w     = ept_w[keep]
    temp_c    = Float64.(raw_json.hourly.temperature_2m)[keep]
    solar_wm2 = Float64.(raw_json.hourly.shortwave_radiation)[keep]
    T_hours   = length(ept_w)

    # Server utilization (88% baseload + 6% daytime inference peak), by local hour
    u_it = [0.88 + 0.06 * max(0.0, sin(pi * (hour(t) - 6) / 14)) for t in ept_w]
    # Cooling PUE driven by real Ashburn, VA hourly temperature
    pue  = [1.10 + 0.008 * max(0.0, T - 15.0) for T in temp_c]
    # Total Data Center Power Requirement (MW) & Solar Capacity Factor [0, 1]
    L_dc_mw  = round.(P_IT_rated .* u_it .* pue, digits=2)
    cf_solar = round.(min.(1.0, solar_wm2 ./ 1000.0), digits=3)

    df_dc = DataFrame(
        Hour = 1:T_hours,
        Timestamp = Dates.format.(ept_w, "yyyy-mm-ddTHH:MM"),
        Temp_C = temp_c,
        Solar_CF = cf_solar,
        IT_Utilization = round.(u_it, digits=3),
        Cooling_PUE = round.(pue, digits=3),
        DC_Load_MW = L_dc_mw
    )
    CSV.write(dc_csv_file, df_dc)
    println("3. Rebuilt ", dc_csv_file, ": ", T_hours, " hours")
end

# ==============================================================================
# 4. Year-long weather (Nov 1, 2025 - Nov 1, 2026) to match the PJM LMP year file
#    Requested in UTC so it joins exactly with the LMP file's datetime_beginning_utc;
#    local (EPT) time is computed for the utilization profile.
# ==============================================================================
year_weather_file = joinpath(ROOT, "data", "ashburn_weather_year_raw.json")
year_url = "https://archive-api.open-meteo.com/v1/archive?latitude=39.0438&longitude=-77.4874&start_date=2025-11-01&end_date=2026-11-01&hourly=temperature_2m,shortwave_radiation&timezone=GMT"
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

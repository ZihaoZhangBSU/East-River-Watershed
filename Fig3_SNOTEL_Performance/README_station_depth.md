# Supplemental station snow-depth performance

`Supp_Station_Depth_Performance.py` computes and plots the snow-depth counterpart
of the current station/year Figure 3. It has three station rows and three
metric columns: mean bias, RMSE, and depth peak-date error. It preserves the
current model colours, blue filled iSnobal circles and orange open Noah-MP
squares, paired horizontal connectors, year order, Upper Taylor asterisk and
WY2026 shading. It does not calculate or display a snow-disappearance score.

Python 3.10+, NumPy and Matplotlib are required; existing Fig. 3 dependencies
already cover these. Inputs, fonts and generated research assets are supplied
separately and are not included in the source-only repository.

```powershell
python Supp_Station_Depth_Performance.py `
  --snotel /data/snotel_daily.csv `
  --isnobal /data/isnobal_station_daily.csv `
  --noahmp /data/noahmp_station_daily.csv `
  --output-directory /path/to/figure-output `
  --font /fonts/NimbusSans-Regular.ttf
```

The optional font argument accepts absolute or relative paths. Without it the
script uses Matplotlib's bundled DejaVu Sans. The requested font must resolve
to the supplied file. No network access or model simulations are used.

Inputs require date (`YYYY-MM-DD`), water year, station ID/name, and depth in
millimetres: `snow_depth_mm` for SNOTEL/Noah-MP and `snow_depth_nearest_mm` for
iSnobal. Supported IDs/names are 380/Butte, 737/Schofield Pass and 1141/Upper
Taylor. Data must cover all daily dates for every station in WY2018–2025
(1 October–30 September) and WY2026 (1 October–30 June), with finite,
nonnegative depth. Duplicate/missing dates and nonfinite values are rejected;
the script never fills gaps, removes valid zeros, substitutes SWE, changes
dates, or pads July–September 2026. Other daily-table fields are ignored.

MB is mean(model minus SNOTEL) and RMSE is the square root of mean squared
differences on the identical complete dates shared by observations and both
models: 365/366 historical days and 273 WY2026 days. Each source's depth peak
is its own first chronological **exact maximum from 1 October through 30
June**. Later historical July–September values enter amount scores but cannot
replace that peak. Peak error is model depth-peak date minus SNOTEL depth-peak
date in calendar days. Positive errors indicate a later model peak. Repeated
peaks can correspond to different events; these dates do not measure a uniform
shift of the full season.

The workflow reads the existing native nearest-cell model series, preserving
their spatial sampling and archived dates. It verifies saved-data calendars,
values and scoring, and does not re-extract raw model fields, validate sensor
quality flags, reconcile coordinates or establish equal intraday measurement
clocks. Original model extraction converted depth from metres to millimetres;
the supplied daily CSVs are already in millimetres and are not converted again.

Outputs are `JHM_S_station_depth_performance.pdf`, a 600 dpi PNG, a native-size
96 dpi preview, a 54-row full-precision `station_depth_metrics.csv`, and an
execution JSON containing all source hashes, 81 peak/tie audits, 27 common-date
supports, 162 plotted marker values, definitions, font/runtime/output hashes.
Native size is **6.5 × 6.8 inches**, with the current Fig. 3 font hierarchy:
**8 pt labels, titles and panel letters; 7 pt ticks and legend**. Margins,
spacing, marker sizes/fills, line widths and year shading also match the
current Fig. 3. No caption is embedded; the manuscript supplement supplies
the caption and asterisk meaning.

The complete metric CSV can be replayed without daily inputs:

```powershell
python Supp_Station_Depth_Performance.py `
  --metrics /data/station_depth_metrics.csv `
  --output-directory /path/to/replay-output `
  --font /fonts/NimbusSans-Regular.ttf
```

Its schema requires all 54 station/year/model combinations, exact period/day
counts, observed and model depth peak dates/values/tie counts, and consistent
model-minus-observation errors. Replaying a table validates internal metadata;
it does not independently establish upstream data accuracy.

```powershell
python -m unittest discover -s . -p test_station_depth_performance.py -v
```

Synthetic checks cover depth rather than SWE selection, zero-day denominators,
exact ties and first dates, full-year versus partial/leap-year support, July
exclusion from peak timing, missing/nonfinite input rejection, score metadata,
and all nine panels' model symbols, marker counts, shading and print type.

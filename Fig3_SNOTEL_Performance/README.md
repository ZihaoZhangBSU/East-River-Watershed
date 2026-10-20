# Figure 3: SNOTEL performance

Generate the 12-panel Figure 3 comparing iSnobal and Noah-MP with SNOTEL at
Butte, Schofield Pass, and Upper Taylor for WY2018-WY2026. The existing station
labels and annotation convention are preserved.

This folder contains code, synthetic tests, requirements, and instructions.
Scientific inputs, figures, manuscripts, and fonts are supplied separately.
The scripts accept explicit paths and do not access the network.

## Install and generate the figure

Python 3.10 or newer is required. From this folder:

```shell
python -m pip install -r requirements.txt
python Fig3_SNOTEL_Performance.py --metrics /data/SNOTEL_figure_data.csv --output-dir figures
```

Replace the example input path with the saved metrics table. The command writes
`F3_snotel_performance.pdf`, a 600-dpi PNG, and `fig3_execution.json` containing
input/code/output hashes and runtime versions. The default font is Matplotlib's
bundled DejaVu Sans. To reproduce the publication typography, supply its original
Nimbus Sans font with `--font /fonts/NimbusSans-Regular.ttf`. A different font
changes typography, not plotted values.

## Recompute metrics from daily station inputs

The optional metric builder uses only Python's standard library:

```shell
python recompute_fig3_metrics.py --snotel /data/snotel_daily.csv --isnobal /data/isnobal_station_daily.csv --noahmp /data/noahmp_station_daily.csv --output-dir metrics
python Fig3_SNOTEL_Performance.py --metrics metrics/SNOTEL_figure_data.csv --output-dir figures
```

Add `--reference /data/SNOTEL_figure_data.csv` to the first command to compare all
216 plotted metric values and their event metadata with a saved reference.
The builder also writes `fig3_metrics_audit.json`, recording source hashes,
scoring support, observed zero counts, all 81 timing events, peak ties,
confirmation dates, and censor flags. Output files cannot replace inputs or
the supplied comparison reference.

Daily CSVs require `date` (YYYY-MM-DD), `station_id`, `station_name`, `water_year`,
and SWE in millimetres. SWE columns are `swe_mm` for SNOTEL and Noah-MP and
`swe_nearest_mm` for iSnobal. IDs/names are 380/Butte, 737/Schofield Pass, and
1141/Upper Taylor. Provide complete daily records for every station and year:
October 1-September 30 for WY2018-WY2025, October 1-June 30 for WY2026.
Additional columns are accepted but do not enter the scores.

## Metric definitions

- MB is the mean of model minus SNOTEL SWE; RMSE is the square root of the mean
  squared difference. Both use the same October-June dates in all three sources,
  including zero SWE and leap days. There is no interpolation or gap filling.
- Peak is each source's first chronological exact maximum in October-June.
  Exact ties use the first date; tie information is retained in the audit.
- Snow disappearance date (SDD) is the first day strictly after that source's
  own peak starting seven consecutive calendar days with SWE strictly below
  10 mm. It is the first day of the run, not its confirmation day.
- Historical SDD searches continue through September 30; WY2026 stops June 30.
  Missing/nonfinite days break the confirmation run. An incomplete endpoint run
  remains unresolved. A confirmed date is never replaced by a record endpoint.
- Timing errors are model date minus observed date in calendar days. SDD error
  is missing if either event is unconfirmed and is displayed as `C`, never zero.
  The original `sdd_censored` flag describes the model event; new tables also
  include `observed_sdd_censored`. When that optional field is absent, the
  plotter requires an uncensored observed event, as in the reviewed original.

The 54-row metrics table requires:

```text
station_id,station_name,wy,model,mb,rmse,peak_date,sdd_date,peak_error_days,sdd_error_days,sdd_censored,available_record_end
```

Rows must contain exactly three stations x nine years x two models (`iSnobal`
and `Noah-MP`). The plotter rejects duplicates, omissions, invalid/nonfinite
scores, inconsistent station metadata, impossible dates, invalid censor flags,
and disagreement between the two rows' implied SNOTEL timing dates. It preserves
the original colors, marker shapes, panel order, annotations, and WY2026 shading.

## Validation and scope

```shell
python -m unittest discover -s tests -v
```

Tests use synthetic data only. They cover daily calendars, zero inclusion,
leap days, tied peaks, post-own-peak SDD, the strict 10 mm threshold, seven-day
confirmation, missing days, censoring, invalid metric tables, and plot artists.

The publication check reproduced all 216 saved metric values (maximum amount
roundoff 2.84e-14 mm), with exact timing dates/errors, and checked the 216 plotted
positions. Tested with Python 3.12.14, NumPy 2.5.4, pandas 3.0.1, and Matplotlib
3.11.2. Lower dependency bounds describe supported APIs, not an exhaustive
version matrix. With the original publication font and the same runtime, the
600-dpi PNG is pixel-identical to replaying the archived Figure 3 plot function.

This validates the saved daily station inputs and Figure 3 calculations. It
does not re-extract raw model fields, check sensor quality flags, reconcile
station-coordinate provenance, or establish identical intraday timestamps.
Independent source checks remain necessary for those questions. Original
research inputs and manuscript files are unchanged by this code publication.

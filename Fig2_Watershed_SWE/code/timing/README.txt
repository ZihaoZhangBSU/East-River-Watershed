Dependent timing, F8(a), and S5 generators
========================================

Requires Python 3.10+, numpy, pandas, matplotlib, and an available font provided
separately. The local reviewed figures used Nimbus Sans; no font asset is
bundled with the public helper code. Inputs are supplied separately; research data and manuscripts are not
distributed with this code. These scripts only analyze existing model outputs.
They never run a snow model, alter a source cache, or remap one-dimensional means.

Run the edge-case tests from this directory:
  python timing_common_grid.py --test
  python dependent_F8_S5.py --test

1. F6, depletion-threshold figure, annual/near-maximum tables and old/new audit

  python timing_common_grid.py --test --daily INPUTS/Fig2_common_grid_daily.csv --old-daily INPUTS/figure04_native_watershed_daily.csv --old-thresholds INPUTS/timing_threshold_metrics.csv --upstream-manifest INPUTS/model_operator_provenance.json --output RESULTS/timing --figures RESULTS/figures --font INPUTS/NimbusSans-Regular.ttf

Required new daily columns: date, wy, swe_i_mm, swe_n_mm. Every actual daily
date from October 1 through June 30 is required for WY2018--2026, including leap
days. Explicit fixed-footprint/state-clock declarations, when present, are
validated. Peak is the first exact unsmoothed maximum; crossings are first
strictly later daily SWE <=40/50/60% of own model/year peak. No interpolation or
persistence is applied; an absent crossing remains missing. First/last
100/99/95% near-maximum anchors hold the original 50% crossing fixed. The old
native daily file and old threshold table verify that the established native
calculation is reproduced before comparing the changed spatial operator.

2. S5 selection and F8(a) curve/shading

  python dependent_F8_S5.py --test --isnobal-cache INPUTS/isnobal_wy2026_noahgrid_daily.npz --noah-cache INPUTS/noahmp_wy2026_noahgrid_daily.npz --grid INPUTS/grid_common.npz --daily INPUTS/Fig2_common_grid_daily.csv --timing RESULTS/timing/figure06_timing_metrics.csv --old-points INPUTS/wy2026_simultaneous_loss_points.csv --baseline-data INPUTS/F8_baseline --font INPUTS/NimbusSans-Regular.ttf --output RESULTS/timing --figures RESULTS/figures

Both NPZ caches require dates and swe_mm[day,row,col] on the common Noah-MP
grid. grid_common.npz requires the canonical fixed_domain_mask and fractional watershed
weights. The final WY2026 selection is February 21--May 10 inclusive, bounded
by the earlier own peak and later first half crossing of the new common-grid
means. Net loss is previous minus current D23 state, assigned to the ending
date. Both losses must be finite and strictly >1 mm/day; selected target
cell-days have equal weights. The script first reproduces all old February
21--May 9 point pairs, then records the added May 10 cell-days and new scores.

INPUTS/F8_baseline must contain the separately supplied archived files:
  F8_h1_recomputed.csv
  F8_wy2026_snow_state_footprints.csv
  F8_wy2026_precipitation_events.csv

F8(a) uses the final common-grid SWE curve and new date shading. F8(b) retains
the archived common-grid footprint/precipitation panel, independently checking
all footprint percentages against the reloaded spatial caches. F8(c) retains
the separately labeled native-grid solid-input/gain/retention calculation;
its gain is not inferred from the common-grid curve in panel (a).

3. Independent retained native-accounting audit

  python audit_F8_accounting.py --native-daily INPUTS/figure04_native_watershed_daily.csv --native-accounting INPUTS/F8_baseline/F8_h1_recomputed.csv --common-daily INPUTS/Fig2_common_grid_daily.csv --core-summary INPUTS/core_watershed_daily.csv --spatial-coverage INPUTS/spatial_cache_variable_coverage.csv --output RESULTS/timing

This verifies all nine original native SWE endpoint/gain/retention identities,
records the difference from new common-grid endpoints, and documents why the
separate native accounting terms remain unchanged. --spatial-coverage is
optional but records the independent fresh finite-coverage audit when supplied.

Each calculation writes full-precision local CSV/JSON records with input hashes.
The figures use exact 6.5-inch widths, the separately provided font, vector axes,
and 600-dpi PNG output (S5's dense hexbin layer is rasterized in its PDF).

4. Matching annual snowfall-input/retention supporting figure

  python S_Other_Year_Annual_Input_Retention.py --accounting INPUTS/F8_baseline/F8_h1_recomputed.csv --font INPUTS/NimbusSans-Regular.ttf --output RESULTS/figures

This standalone replay preserves the archived nine-year bars, net-gain line,
styles and native-grid accounting. Its legend and Fig. 8(c)'s first category now
both read "Snowfall input". This is a display-label change only: all supplied
values, ratios, initial-state correction, October 2-March 15 input window,
October 1-March 15 storage gain, and model-interface distinctions are unchanged.
The program writes JHM_S_other_year_annual_input_retention.pdf/.png and a local
execution record with code, input, font and output hashes. Data and fonts are
provided separately and are not included in this repository.

For the archived typography, the annual figure uses Matplotlib 3.10.9; the
refreshed Fig. 8 uses 3.11.2. Use separate environments with those versions
when exact rendering is required. The same accounting values are plotted in
either version, but a version change can shift text placement slightly.

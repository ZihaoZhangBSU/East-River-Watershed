# Fig2_Watershed_SWE

`Fig2_Watershed_SWE.py` generates Figure 2 from watershed SWE on the Noah-MP
analysis grid. This revision replaces the earlier native-grid display. It is
postprocessing of saved model states and original, unfiltered ASO products;
neither snow model is rerun. Research inputs and figures are distributed
separately from the public code.

## Prepare daily model means

Install NumPy and pandas, then run:

```sh
python code/models/extract_fixed_model_means.py --cache-dir /path/to/spatial_caches --grid /path/to/grid_common.npz --output-dir /path/to/prepared_data/models
```

The cache directory must contain the existing `isnobal_wyYYYY_noahgrid_daily.npz`
and `noahmp_wyYYYY_noahgrid_daily.npz` files for 2018 through 2026. These contain
`dates` and `swe_mm` arrays with shape (day, 66, 54). Historic years must contain
every October-September day; WY2026 ends June 30. The grid NPZ supplies
`weights`, `fixed_domain_mask`, `transform`, and projected metric `crs`.
The last-daily iSnobal states were area averaged from the native 50 m grid to
the 1 km Noah-MP grid in the verified production cache builder. This program
does not attempt to remap a one-dimensional basin mean.

All 842 positive watershed-area targets form the fixed model footprint. Both
models use this same footprint every day and fractional watershed-area
weights, accumulated in float64. Finite zeros remain. Missing or negative
retained model SWE aborts extraction, rather than changing daily masks. Input
hashes, independent scalar checks, and separate ancillary-variable coverage
are recorded. See `code/models/MODEL_OPERATOR_README.txt` for scope and limits.

## Prepare unfiltered ASO means

Install `requirements-common-grid.txt`, then run:

```sh
python extract_fig2_common_aso.py --aso-directory /path/to/raw_ASO --flights /path/to/flights.csv --grid-weights /path/to/grid_common.npz --output-directory /path/to/prepared_data/aso
```

The flight CSV contains ISO `date` and optional relative `source_file`. Source
filenames otherwise follow `ASO_50M_SWE_USCOGE_YYYYMMDD.tif`. Source masks,
declared NoData, and nonfinite values define validity; all other supplied
finite values, including zero, negative and extreme values, remain unfiltered.
The extractor checks source units before converting SWE to millimeters.

SWE is area averaged over valid contributing source area in each full target
cell. Separately, the source-validity raster is zero padded beyond its extent
before averaging, so the coverage denominator is the **whole target-cell
area**, including space outside the ASO raster. Retain targets with at least
50% valid source coverage, then compute `sum(a_j * SWE_j) / sum(a_j)`, where
`a_j` is target-cell area inside the watershed. Do not multiply these basin
weights by source coverage again. This two-stage operator averages full-cell
source fields before applying fractional watershed weights; it is not native
pixel clipping to the watershed polygon.

Reported ASO basin coverage is the watershed area represented by retained
targets divided by the fixed target watershed area. It is distinct from the
exact native observed footprint. Each flight can have a different retained
footprint. For direct numerical comparisons, sample both models on the same
flight-specific targets; the model curves themselves retain the fixed domain.

## Generate Figure 2

```sh
python -m pip install -r requirements.txt
python Fig2_Watershed_SWE.py --data-dir /path/to/prepared_data --output-dir /path/to/fig2_outputs
```

An optional `--font /path/to/font.ttf` selects an available font. Otherwise the
Matplotlib bundled sans-serif font is used. The prepared directory contains:

- `models/Fig2_common_grid_daily.csv`: `date`, `wy`, `swe_i_mm`, `swe_n_mm`,
  `fixed_model_area_fraction`, `daily_common_valid_area_fraction`,
  `fixed_model_cells`, `state_hour_utc`.
- `aso/aso_common_grid_flights.csv`: `date`, `wy`, `aso_watershed_mean_mm`,
  `aso_valid_area_fraction`, `retained_target_watershed_area_m2`,
  `fixed_model_watershed_area_m2`.

Model coverage declarations must equal 1, the target count must be 842, and
the state hour must be 23 UTC. ASO denominators must match the fixed analysis
area, 748679519.6193271 square meters. Dates, water years, complete calendars,
finite values, and all 14 documented flights are checked before plotting.
Outputs are a vector PDF, 600 dpi PNG, record-period statistics, a JSON file
defining the display statistics, and an execution record with hashes.
Each panel displays daily **intermodel agreement**, not ASO validation:

- `r`: Pearson correlation of the paired daily watershed means.
- `MB`: mean(Noah-MP minus iSnobal), in millimeters.
- `rMB`: 100 times MB divided by mean iSnobal SWE, in percent.

The iSnobal series is the reference for relative MB; the manuscript caption
identifies the sign and denominator. The original 14 ASO markers, daily curves,
axes and periods remain unchanged. No ASO value enters these statistics.
Undefined correlation (a constant series) or relative MB (zero iSnobal mean)
is displayed as `NA`, never as zero. Nonfinite pairs abort the calculation.

Primary amount summaries use full October-September
WY2018-WY2025 and the available October-June WY2026 record, retaining snow-free
days. `Fig2_model_statistics.csv` reports both model means, Noah-MP-minus-iSnobal
bias, Pearson correlation, and relative mean difference. The relative difference
is undefined when the iSnobal reference mean is zero. The extractor also supplies
separate common October-June statistics. These two period definitions must remain
distinct: differing record lengths do not provide equal-season absolute-mean
anomalies. Peak/depletion timing and forcing diagnostics keep their separately
declared October-June windows. The Figure 2 curves are unchanged.

## Validation and dependent figures

```sh
python -m pip install -r requirements-common-grid.txt -r requirements-native.txt
python -m unittest discover -s tests -v
python -m unittest discover -s code/models -p 'test*.py' -v
```

Run these commands from this figure folder. The native requirements are needed
only for the retained native audit utility and its tests; they are not required
to plot the prepared common-grid inputs.

Synthetic tests independently exercise partial coverage, source-edge
normalization, missing values, finite zeros, fractional weights, fixed masks,
and calendar/plot-input failures. Original sources are read only. Previously
published `extract_fig2_native_aso.py` remains an audit utility for the former
native operator and is not the current Figure 2 input pipeline.

Related generators retain figure-specific names:

- [Figure 5 ASO performance](code/f5/README.txt): corrected full-target ASO
  coverage, unfiltered SWE/depth scores, per-flight spatial correlation, and
  independent verification.
- [Supplemental SWE and depth elevation performance](code/f5/README_elevation.md):
  separate flight-by-elevation bias, RMSE and spatial correlation figures with
  variable-specific support, shared metric scales and explicit NA cells.
- [Figure 6 timing and dependent Figure 8/S5 updates](code/timing/README.txt):
  timing tables, revised conditional-loss dates, and separated native accounting.
- [Supplemental S2 three-source density](code/S2_Three_Source_Density.py):
  the same corrected ASO coverage on the established shared density support.

Spatial error scores use equal retained-cell weights;
Figure 2 basin means use fractional watershed-area weights. Storage/input
accounting remains on its original native support unless every term has been
verified on the same replacement support and clock. No new data or code
license is assigned by this update.

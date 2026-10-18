# Fig2_Watershed_SWE

`Fig2_Watershed_SWE.py` generates Figure 2, with native-grid model watershed
SWE curves and unfiltered native ASO flight-footprint means.

## Plot from prepared inputs

Install the plotting requirements, then supply the prepared data directory and
an available TrueType font (the local figure delivery includes Nimbus Sans):

```sh
python -m pip install -r requirements.txt
python Fig2_Watershed_SWE.py --data-dir /path/to/prepared_data --font /path/to/font.ttf --output-dir /path/to/fig2_outputs
```

The data directory must contain `watershed_swe_native_daily.csv` and
`aso_native_flight_footprint.csv`. Input data are distributed separately.
The first CSV needs date, wy, swe_i_mm, swe_n_mm and i_native_valid_area_fraction.
The second needs date, wy, aso_watershed_mean_mm, aso_valid_area_fraction,
valid_area_m2 and watershed_area_m2. Native area columns deliberately distinguish
the corrected input from the older remapped-cell observation table.

Outputs are a vector PDF, a 600-dpi PNG, model-comparison statistics and
an execution record with input/output hashes. Model summaries are preserved;
this is figure-generation code, not a model simulation.

## Recompute native ASO means

```sh
python -m pip install -r requirements-native.txt
python extract_fig2_native_aso.py --aso-directory /path/to/raw_ASO --watershed /path/to/WBDHU10.shp --watershed-field name --watershed-value "East River" --flights /path/to/flights.csv --output /path/to/prepared_data/aso_native_flight_footprint.csv --manifest /path/to/prepared_data/aso_native_provenance.json
python -m unittest discover -s tests -v
```

The flight CSV needs an ISO date column and may specify source_file filenames.
The default raw filename pattern is ASO_50M_SWE_USCOGE_YYYYMMDD.tif.

Native means are sum(SWE times exact valid source-pixel/watershed overlap area)
divided by total valid overlap area. Coverage divides that valid area by the
FULL projected watershed polygon area. Source-declared NoData/source masks
and nonfinite SWE are excluded; finite zero, negative and extreme values remain.
There is no Noah-MP remapping, magnitude filter or target-cell coverage threshold
for these markers. The code checks SWE units, uses projected metric areas,
rejects ambiguous units/geographic or rotated grids, and protects input paths
against accidental output overwriting. The supplied tests use synthetic data.

Model and ASO spatial supports differ; their overlay supplies seasonal context.
Other ASO comparisons can use separately defined spatial operators. Source
filenames, raw-data hashes and numerical validation records belong in the
separate scientific-data delivery. No new code/data license is assigned here.

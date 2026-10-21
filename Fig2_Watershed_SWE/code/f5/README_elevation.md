# Supplemental ASO SWE and depth elevation performance

`Supp_ASO_SWE_Elevation.py` plots an existing, validated 224-row score CSV.
It selects exactly 84 elevation-band rows for each requested variable and does
not reopen rasters, pool flights, recompute metrics, or run models. Keep it beside
`Fig5_ASO_Performance.py` in `code/f5/`; the helper uses that module's strict
score-schema validator and canonical fourteen flight dates. Dependencies are
the same as that figure generator (NumPy, Matplotlib, rasterio, affine and the
shared remapping module).

```powershell
python code/f5/Supp_ASO_SWE_Elevation.py `
  --scores-csv /path/to/F5_elevation_metrics_corrected.csv `
  --output-directory /path/to/elevation-output `
  --font /path/to/NimbusSans-Regular.ttf `
  --variable both
```

The font argument is optional; without it the helper uses Matplotlib's bundled
DejaVu Sans. Font and input files are not included in the code repository.
Absolute and relative font paths are accepted, including parent-directory
components. The resolved font registry entry must identify the supplied file.
`--variable SWE` is the compatibility default. `--variable Depth` generates only
the depth counterpart, while `--variable both` generates both. The helper name
and original SWE output basename are retained for existing workflows.

The PDF and 600 dpi PNG have native dimensions **6.5 × 7.2 inches**, with all
visible text **8.05 pt**. Placement at the manuscript template's 6.476-inch
width gives approximately 8.02 pt text and 7.17-inch height. Preserve a width
of at least 6.460 inches to retain text of at least 8 pt. An additional 96 dpi
PNG previews the native physical size.

Each figure's six panels are iSnobal (a–c) and Noah-MP (d–f), with MB in millimetres,
RMSE in millimetres, and spatial Pearson r across columns. Each panel has
fourteen chronological flight-date rows and Lower, Middle, Upper columns.
The unchanged terrain bands are Lower <2896 m, Middle 2896–<3353 m, and
Upper ≥3353 m on the canonical Noah grid. MB is model minus ASO. Equal-cell,
model-specific finite ASO–model pairs retain valid zeros and extremes; ASO
coverage and positive watershed eligibility are unchanged from Fig. 5. Depth
uses its own depth-valid pairs and counts; SWE support is never substituted.

Within each figure, both models and all dates/bands share fixed limits per
metric. SWE uses MB −300 to +300 mm and RMSE 0–400 mm. Depth uses MB −600 to
+600 mm and RMSE 0–800 mm. Both use the entire sequential `viridis` colour
map for r=0–1, since every defined correlation in these verified data is
positive. A future score outside these limits, including a negative r, causes
a validation error and requires an explicit scale revision. Values are never
silently clipped. MB retains a diverging colour map centred at zero; RMSE is
nonnegative. Undefined r is gray
hatched and marked NA; finite zero is retained as a numerical value and has
the ordinary zero colour. The two undefined r cells in the verified inputs
are Noah-MP Lower on 24 May 2018 and 10 June 2019 for each variable, where all
retained model values are zero. SWE has 187 pairs on each date; depth has
103 and 196 pairs, respectively. MB and RMSE remain defined there.

Outputs include `JHM_S_ASO_SWE_elevation_performance.pdf/png` and/or
`JHM_S_ASO_depth_elevation_performance.pdf/png`, variable-specific native-size
96 dpi previews, and the exact original scalar strings/support counts in a
separate 84-row CSV per variable. Each variable's JSON manifest contains full plotted arrays, NA masks,
shared scales, input/code/font hashes, software versions and output hashes.
The manifest records provenance; it does not independently validate the
upstream raster operator. A separate numerical audit against saved fields is
required before publication.

Run the presentation safety tests with:

```powershell
python -m unittest discover -s code/f5 -p test_supp_aso_swe_elevation.py
```

Research inputs and generated figure assets are deliberately excluded from
the source-only repository.

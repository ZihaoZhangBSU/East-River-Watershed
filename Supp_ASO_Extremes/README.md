# Supplement: ASO extreme locations

`Supp_ASO_Extremes.py` generates `JHM_S1_aso_extremes.pdf` and `.png` (currently
Figure S3 in the manuscript). It plots the previously verified original ASO
samples and overlays the actual East River HUC10 watershed, **1402000102**.
It does not rerun models, read raw model grids, filter or modify observations,
or derive a watershed from a raster edge.

Install `requirements.txt`, then provide separately supplied inputs:

```text
python Supp_ASO_Extremes.py --locations aso_extreme_locations_verified.csv --flights watershed_swe_aso_14_flights.csv --boundary WBDHU10.shp --output-dir figures
```

Keep the original `.shp`, `.dbf`, `.shx`, and `.prj` files together. An equivalent
vector file with a declared CRS and the same `huc10` attribute is accepted.
`--huc-field` overrides the attribute name. The source boundary used for the
manuscript is the original NAD83/EPSG:4269 WBDHU10 East River feature
(ObjectID 14, DBF load date 2012-10-24); the plotting program transforms it to
the observations' EPSG:32613 (WGS84 UTM zone 13N). Coordinates are displayed
in kilometers with equal map scaling. The entire vector boundary is shown.

The locations file must contain the existing 78 unique flight/row/column
records and columns `date,row,column,x_m,y_m,crs,longitude,latitude,swe_m,
snow_depth_m`. The flight file supplies the 14 distinct ISO dates in panel
order, including flights with no flagged samples. All 14 exact reference dates
must appear in chronological order. The program checks the
original strict criterion (SWE >5 m or depth >15 m), coordinate consistency,
and polygon placement. These checks never remove records.

The six flight colors, marker shapes, cell values and 14-flight counts are
retained. The black rings retain the original rule: another flight has a
flagged center within 125 m, including exactly 125 m. This is a spatial
proximity indicator, not independent evidence of the cause of extreme snow.
The black rings are fixed-size marker symbols, not drawn 125-m-radius buffers.

Use `--font path/to/font.ttf` to supply the same typeface as the manuscript.
The default uses Matplotlib's bundled DejaVu Sans. No font or research data
are distributed in this code-only folder. An execution JSON records input
hashes, boundary metadata, counts, output hashes, the generator's SHA256, and
Python and package versions.

```text
python -m unittest discover -s tests -v
```

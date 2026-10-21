# East River Watershed

Figure-generation code.

- [Fig2_Watershed_SWE](Fig2_Watershed_SWE/): Figure 2, fractional watershed-area SWE means on the Noah-MP analysis grid, with unfiltered ASO remapping and explicit flight coverage. Includes daily intermodel correlation and signed absolute/relative mean bias, validated extraction helpers, and synthetic tests.
- [Fig3_SNOTEL_Performance](Fig3_SNOTEL_Performance/): Figure 3, station SWE bias, RMSE, peak timing error, and snow-disappearance timing error. Includes a standalone plotter, daily metric calculator, input validation, and synthetic tests.
- [Fig4_SNOTEL_WY2026](Fig4_SNOTEL_WY2026/): Figure 4, WY2026 station SWE, model-cell ranges, and timing markers at a compact 5.85-inch print width with 8-point text. Replays validated cached tables.
- [Figure 5 ASO SWE and snow depth performance](Fig2_Watershed_SWE/code/f5/README.txt): SWE above snow depth, with per-flight bias, RMSE and spatial Pearson correlation. Validated score replay preserves the original unfiltered remapping and variable-specific scoring supports.
- [Supplemental ASO SWE and depth elevation performance](Fig2_Watershed_SWE/code/f5/README_elevation.md): Separate six-panel SWE and depth heatmaps by flight and elevation band, with shared metric scales, positive-range correlation colors and explicit undefined correlations.
- [Supp_ASO_Extremes](Supp_ASO_Extremes/): Supplemental ASO extreme-location map and flight counts, with the actual East River HUC10 watershed outline transformed to the sample coordinates.
- [Supp_Station_Bias_Storage_Density](Supp_Station_Bias_Storage_Density/): Supplemental station precipitation/SWE bias and peak-storage/density diagnostics, with explicit full historical and partial WY2026 periods, source-specific supports, and precipitation gap-fill records.

This repository upload contains source code, requirements, synthetic tests and usage instructions. Scientific input data and manuscript files are supplied separately by the author.

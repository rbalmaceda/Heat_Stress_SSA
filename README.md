
---
```rmarkdown
## Reproducible example for heat stress calculations using ERA5

This repository contains a reproducible R Markdown example illustrating part of the workflow
used in the manuscript:

Heat stress conditions in southern South America: characterization and long-term changes under
observational uncertainty

The notebook includes:
- collecting daily ERA5 data from the climate4R repository,
- estimation of specific humidity and relative humidity,
- computation of heat stress indices,
- percentile threshold calculation (P90/P99),
- a simple example of exceedance counts,
- linear trend estimation with Mann-Kendall significance testing.

## Files

- `reproduce_era5_example.Rmd`: main reproducible notebook
- `functions_humidity.R`: auxiliary humidity-related functions

## Notes

This is a lightweight reproducible example, not the full analysis pipeline used in the
manuscript.  
To reduce computational cost, the example is restricted to a small domain near Porto
Alegre and to DJF only.

## Requirements

Main R packages used:
- convertR
- HeatStress
- loadeR
- transformeR
- visualizeR
- Kendall
- trend
- lubridate

## Data access

ERA5 data are accessed through the climate4R repository:
https://data.meteo.unican.es/

## Citation

Balmaceda-Huarte R., Casanueva A.,Bettolli M.L., Heat stress conditions in southern South America: characterization and long-term changes under observational uncertainty, 2026.

```

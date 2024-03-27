# Trait-mapping-SBL

This repository is intented to support the manuscrit "Mapping canopy foliar functional traits in a mixed temperate forest using imaging spectroscopy".

The repository is maintained by Alice Gravel. Please feel free to contact me at alice.gravel@umontreal.ca for any questions.

## Scripts

The repository contains R scripts, numbered from 00 to 15. Scripts 00 are either functions or side analyses.

01. Preprocess and clean trait data
02. Extract spectra from tree crowns
03. Post-processing/cleaning of mean spectra per tree crown
04. Prepare data for PLSR modeling
05. Split data (calibration/validation)
06. Build models
07. Plot validations statistics of models
08. Plot models coefficients
09. Interpretation of models: Plot VIPs
10. Extract spectra from imaging spectroscopy of the entire study area
11. Post-processing/cleaning of all pixels
12. Inference : apply PLSR models to every pixel
13. Mapping traits: Creating rasters from predicted traits, uncertainties and relative uncertainties 
14. Confidence in predictions (pixel-wise)
15. Compare measured traits from field sampling with trait maps using delineated tree crowns

## References

The approach for model bulding (scripts 04-08) was adapted from the following article :
Kothari, S., Beauchamp‐Rioux, R., Blanchard, F., Crofts, A. L., Girard, A., Guilbeault‐Mayers, X., Hacker, P. W., Pardo, J., Schweiger, A. K., Demers‐Thibeault, S., Bruneau, A., Coops, N. C., Kalacska, M., Vellend, M., & Laliberté, E. (2023). Predicting leaf traits across functional groups using reflectance spectroscopy. _New Phytologist, 238_(2), 549–566. https://doi.org/10.1111/nph.18713
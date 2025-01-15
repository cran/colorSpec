# NEWS for **colorSpec** package


### Changes for version 1.6-0  [2025-01-15]
* added function `computeCRIdata()` and improved the capability of `computeCRI()`
* all logging done with package **logger**, which is imported
* improved the wording of some logged messages

### Changes for version 1.5-0  [2024-01-28]
* moved initialization of `colorSpec` options from `.onAttach()` to `.onLoad()`, so `colorSpec` can be used without attaching it; thanks to Pedro Aphalo
* fixed some internal warnings about `partial argument match`; thanks to Pedro Aphalo
* fixed `dimnames()` problem in `test-organization.R`
* removed `exportClasses` directive

### Changes for version 1.4-0  [2022-05-04]
* in function `calibrate()` added new option for `response` that is compliant with ASTM and CIE
* modified relevant vignettes to use the new calibration option
* in the man pages, changed mentions of vignettes to hyperlinks
* inactivated some (possibly) invalid URLs
* in User Guide vignette, fixed embedded table problem by replacing `cat()` with `knitr::raw_html()`
* increased width of User Guide vignette to match the width of embedded tables
* moved packages **MASS** and **spacesXYZ** from Imports to Suggests

### Changes for version 1.3-0  [2021-12-20]
* when reading .sp files, use value of `MEAS_TYPE` to assign the **colorSpec** `quantity`
* when reading .sp files, divide spectral values by `SPECTRAL_NORM`
* fix undefined color problem in some vignettes, caused by change in TeX package **xcolor**

### Changes for version 1.2-1  [2020-04-01]
* in function `invert()` added new `method='TLSS'`, and updated the vignette **Estimating a Spectrum from its Response - Inverse Colorimetry**
* in function `plot()` added new argument `type`, with custom option `type='step'`, and updated the vignette **Convexity and Transitions**
* added a README file

### Changes for version 1.1-1  [2019-12-07]
* added fix to `rotateOrganization()` in file `test-organization.R`; for upcoming change to `class(a matrix)` in R v 4.0, and when environment variable _R_CLASS_MATRIX_ARRAY_ is set to non-empty

### Changes for version 1.0-1  [2019-06-24]
* added new vignette **Convexity and Transitions - a strict examination of the CIE inverted-U**
* added new functions `bandMaterial()` and `bandRepresentation()`
* added new function `canonicalOptimalColors()`
* added new function `responsivityMetrics()`
* more efficient computation of zonohedra
* now Suggests package `quadprog`

### Changes for version 0.9-1  [2019-05-31]
* restored missing .R files in folder /inst/doc
* fixed NOTE:  "found 1 marked UTF-8 string"
* in `computeCCT()` suppress warning when spectrum is all 0s

### Changes for version 0.8-2  [2019-03-03]
* moved most CCT-related functions to package `spacesXYZ`, which is now imported
* in `probeOptimalColors()`, changed to zonohedral representation of the color solid.
* added function `sectionOptimalColors()`
* changed argument list for `plotOptimals3D()`
* added function `plotOptimals2D()`
* in `planckSpectra()` changed constant `c2` unit from nm*K to m*K to agree with the rest of the literature
* moved RGB-related functions to package `spacesRGB`, which is now Suggested
* package `minpack.lm` is no longer needed, or Imported
* added function `computeSSI()`, requested by Alex Forsythe
* add new theoretical camera `ACES.RGB` = ACES Reference Input Capture Device, from S-2008-001 Academy Color Encoding Specification.
* added 2 bonus spectra from EBU TECH 3355 - Method for the Assessment of the Colorimetric Properties of Luminaires
* add more keys to recognize CGATS files
* fixed documentation error regarding `readSpectraCGATS()`
* bug fix in `plot()`.  Spectra with NA values are now skipped.
* in all calls to `sprintf()`, changed `%d` to `%g`, unless obviously integral. Bug found by Dean Attali.

### Changes for version 0.7-5  [2018-11-19]
* add new function `invert()` plus new vignette **Estimating a Spectrum from its Response - Inverse Colorimetry**
* add new function `rectangularMaterial()`
* in `computeCCT()` etc., add new `method` `'mccamy'`
* now Suggests package `rootSolve`

### Changes for version 0.7-3  [2018-04-01]
* add new function `actinometric()`
* add new function `as.data.frame()`
* add new function `atmosTransmittance()`
* add new function `CCTfromuv()`
* add new functions `is.actinometric()` and `is.radiometric()`
* add new function `ptransform()`, and use it to compute `BT.709.RGB` and `Adobe.RGB`
* add new function `emulate()`
* add new function `as.colorSpec()`
* add new vignette **Emulation of one Camera by another Camera**
* add new vignette **Photon Counting**
* add built-in object `luminsivity.1nm`
* add new spectra files `moths.txt` and `sunglasses.txt`
* add new spectra files `Philips-HPS.txt`, `solar-exposure.txt`, `P4-phosphor-JEDEC.txt`, and `Cree-LED.txt`
* add "featured functions" to all vignettes
* spectral quantity `power` is deprecated, and replaced by `energy`.  `power` still works, but will eventually be removed.
* in `colorSpec()`, add argument `specnames`
* in `resample()`, add arguments `extrapolation` and `clamp`
* in `photometric()` add arguments `photopic`, `scotopic`, and `multiplier`
* in `cs.options()`, partial matching of the option name is enabled
* in `radiometric()`, add arguments `multiplier` and `warn`
* in `metadata()<-` add argument `add`
* in `extradata()<-` add argument `add`, and allow `value` to be `NULL`
* in `product()`, add argument `integration`, and added an ambiguity warning
* in `summary()`, the displayed Integral now works for irregular wavelengths
* in `summary()`, print attribute `ptransform` if present
* in `plotPatchesRGB()`, allow `background` to be linear RGB, fixed bug for `shape`
* in `planckSpectra()`, added new argument `c2`
* in `computeCCT()` and `CCTfromXYZ()` and `CCTfromuv()`, added new arguments `method`, `strict`, and `c2`
* in `computeCRI()`, CCT is now computed with `method='lm'`
* in `calibrate()`, fix special case when there is only 1 spectrum
* in `plot()`, fixed warning when `ylab` is an expression
* in `bind()`, fixed bug when binding the `extradata`
* **colorSpec** options are now stored in the global option list, and start with `'colorSpec.'`
* in vignette **Phenol Red - pH Indicator**, now display RGB matrix explicity, both before and after scaling
* additions and improvements to colorSpec User Guide
* now Imports package `minpack.lm`

### Changes for version 0.6-2  [2017-12-03]
* bug fix: `print.colorSpec()` and `summary.colorSpec()` were sending output to `stderr()`, instead of `stdout()`
* improved reading of CGATS files with standard whitespace convention
* in all sample CGATS files, changed spaces in field names to underscores

### Changes for version 0.6-1  [2017-11-16]
* renamed **proofs.Rmd** to **proofs.txt** - to avoid "Files named as vignettes but with no recognized vignette engine:"

### Changes for version 0.6  [2017-11-16]
* updated **colorSpec-guide.odt** to **colorSpec-guide.Rmd** (rmarkdown v 2)
* for physical models in **colorSpec-guide.Rmd**, replaced $L^2$ by $L^\infty$ and $L^1$, and added proofs
* added 4 new files with reflectance spectra
* added new function `interpolate()`
* added new vignette - an investigation of phenol red
* in `plot.colorSpec()`, compute automatic margin line for ylab, depending on the width of y-axis labels
* in `bind.colorSpec()`, add check for distinct `specnames()`
* bug fix to `"extradata<-.colorSpec"`; `specnames()` was not being preserved
* in logging mechanism, changed output from `stdout()` to `stderr()`, for better compatibility with RStudio
* export function `readCGATS()`, for easy access to non-spectral data in CGATS files

### Changes for version 0.5-3 [2016-05-16]
* Fixed filename case problem on non-Windows platforms and resubmitted

### Changes for version 0.5-2 [2016-05-15]
* Fixed 2 NOTEs and resubmitted

### Version: 0.5-1   [2016-05-14, the first submission]

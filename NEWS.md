# NEWS for **colorSpec** package

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

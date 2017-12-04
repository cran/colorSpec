# NEWS for colorSpec package

### Changes for version 0.6-2  [2017-12-xx]
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

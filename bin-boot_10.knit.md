---
title: "Composite curves and bootstrap C.I.'s via binning <br> (bin-boot.R)"
output:
  html_document:
    css: SI-md-08.css
    fig_caption: yes
    highlight: haddock
    number_sections: yes
    theme: united
    toc: yes
    toc_float: false
    collapsed: no
---



# Introduction #

This script (`binboot-curve.R`) fits step-like composite curves, with bootstrapping by site confidence intervals, via binning.  As is the case with other analysis scripts, the script has two parts, a set-up portion that changes from analysis to analysis, and a computation part that does not change.  This version reads "presampled" data.

# Set up #

The variables `queryname`,  `basename` and `binname` are used to compose the path and filenames for the presampled data.


```r
# age-bin averages of influx data
# bootstrap-by-site confidence intervals

# names
queryname <- "v3i" # 
basename <- "nt2kb" 

# paths for input and output .csv files -- modify as appropriate
datapath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/"
sitelistpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_sitelists/"
sitelist <- "v3i_nsa_globe"
outpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_curves/"

# prebinning bin width
pbw <- 10

pbw_char <- as.character(pbw)
if (pbw < 10) pbw_char <- paste("0", pbw_char, sep="")
```

Set the number of bootstrap samples or replications:


```r
# number of bootstrap samples/replications
nreps <- 200
```

Specifiy the age-bin centers, and spacing.  This version uses 20-year wide bins, with the first bin centered at -60 yr BP, or 2010 CE, and the last at 1950 yr BP, or 1 CE.


```r
# age bin centers
abw <- 20 
abinbeg <- -60 
abinend <- 1950
```

Set array sizes for saving bootstrap results.


```r
# array sizes
maxrecs <- 2000
maxreps <- 1000

# plotting set up
xmin <- 0; xmax <- 2020; ymin1 <- -1.0; ymax1 <- 1.0; ymin2 <- -0.5; ymax2 <- 1.0
xlim=c(xmin,xmax); ylim1=c(ymin1,ymax1); ylim2 <- c(ymin2,ymax2)
xlab="Year CE"; xminortick <- 50
ylab="Normalized Anomalies of Transformed Influx"

# plot output 
plotout <- "screen" # "pdf" # "png" # 
```

# Calculation preliminaries#

## Initial steps ##

Set various output paths and filenames.


```r
# no changes below here

# prebinning file name
binname <- paste("bw", pbw_char, sep="")

# presampled/binned files
csvpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_presamp_csv/"
csvname <- paste("_presamp_influx_",basename,"_",binname,".csv", sep="")

# site list file
sitelistfile <- paste(sitelistpath, sitelist, ".csv", sep="")
sitelistfile
```

```
## [1] "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_sitelists/v3i_nsa_globe.csv"
```

```r
# curve (output) path and file
curvecsvpath <- paste(datapath,queryname,"_curves/",sep="")

# if output folder does not exist, create it
dir.create(file.path(datapath, paste(queryname,"_curves/",sep="")), showWarnings=FALSE)
curvename <- paste(sitelist,"_binboot_",basename,"_",binname,"_abw",as.character(abw),"_",
  as.character(nreps), sep="")
curvefile <- paste(curvename, ".csv", sep="")
print(curvecsvpath)
```

```
## [1] "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_curves/"
```

```r
print(curvename)
```

```
## [1] "v3i_nsa_globe_binboot_nt2kb_bw10_abw20_200"
```

```r
print(curvefile)
```

```
## [1] "v3i_nsa_globe_binboot_nt2kb_bw10_abw20_200.csv"
```

Code block to implement writing .pdf to file:


```r
# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
pdffile <- paste(curvename, ".pdf", sep="")
print(pdffile)
}
# .png plot of bootstrap iterations
if (plotout == "png") {
  pngfile <- paste(curvename, ".png", sep="")
  print(pngfile)
}
```

Read the list of sites to be processed.  Note that this is the site list may contain sites that after transforming and presampling may have no useful data.  Those sites are ignored when the data are read in.






















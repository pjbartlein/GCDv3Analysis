# Introduction #

This is a set of web pages that describe the development of a composite curve of charcoal data drawn from the Global Charcoal Database version 3, (GCDv3) using a set of R scripts.  The intention here is to explicitly document the analysis steps that were developed originally as a set of Fortran programs, and which are now implemented in the R `paleofire` package.  The R scripts described here can also be used as a point of departure for the development of new analysis approaches.

The data source for these examples is a Microsoft Access (`.mdb`) database, downloaded from  [[http://www.gpwg.paleofire.org/]](http://www.gpwg.paleofire.org/), e.g., `GCDv03_Marlon_et_al_2015.mdb`.  In the example here, the scripts aim to reproduce the "Globe" curve in Fig. 6 of Marlon et al. (2016).

# Data #

The data used in this example are contained in two queries, saved as "database" views or internal tables in the Microsoft Access database named `GCDv03_Marlon_et_al_2015.mdb`, downloadable from the Global Charcoal Database [[http://paleofire.org/]](http://paleofire.org/) by exporting the full database.  The queries are named (for historical reasons) `ALL_BART_SITES` and `ALL_BART_DATA` and reside in the Access database as "view".  `ALL_BART_SITES` contains a list of sites, their names and locations, the depositional environment, and the units of measurement (i.e. influx, concentration, etc.) while `ALL_BART_DATA` contains the id, age, depth and quantity of charcoal in each sample.

For the examples here, the following folder structure for the data was used:

		/Projects/GPWG/GPWGv3/GCDv3Data/v3i/
			v3i_curves/
			v3i_debug/
			v3i_mdb/
			v3i_presamp_csv/
			v3i_query/
			v3i_sitelists/
			v3i_sites_csv/
			v3i_stats/
			v3i_trans_csv/

The the root folder and `v3i-mdb/` (into which the database should be copied) are created by the user, and the others are created during the analyses.

# Analysis steps #

There are four steps in the analysis:

1. reading the query results from the data base and making individual "site" .csv files
2. transforming and standardizing or normalizing (i.e. converting to anomalies) the individual records
3. implementing the presampling/prebinning step
4. composite-curve fitting using the R locfit() function, or via binning, and estimating uncertainties via bootstrapping

## 1 Query and site .csv files ##

(`mdb-to-csv.R`)

The first step in the analysis approach here involves examining the two query tables, checking for obvious issues in the data for individual sites, converting all charcoal data to both influx and concentration values, and finally extracting individual .csv files for each site.  All of the scripts here begin by setting appropriate path and folder names.  In the example script, the Access database files are read directly using the `RODBC` package.  (Note that on Windows, compiled versions of this package exist, and a the appropriate database drivers are built into the operating system.  On OS X, a third-party database driver must be used, and the `RODBC` package compiled from source.  There is a separate script, `mdb-to-csv_osx.R` that illustrates this.

The main part of the script loops over the individual sites that are specified in the site query, and does various checks and calculations, including 

- calculation of sedimentation rates and deposition times
- checking for age or depth reversals, or other data issues
- calculation of alternative quantities (e.g. influx, given concentrations)
- writing out a .csv file for each site

The calculation of sedimentation rates and deposition times is required for the conversion of concentration values to influx values and vice-versa.  Various checks for age reversals, zero sedimentation rates, missing data are done.  Typically, when a number of sites are added to the database, there will be issues, which are flagged by this step and resolved.  In the example, such is not the case, but the script illustrates those checks in any case.  The last part of this analysis step involves writing out one "site data" .csv file for each site (e.g. `0001_data.csv`), plus a single "sitelist" file (e.g. `v3i_all.csv`), which can be edited to control the particular selection of sites that are analyzed.

New data not included in the database can be added to the analysis by creating by hand a "site data" .csv file with the same format as those created by `mdb-to-csv.R` and adding a line to the sitelist file. 

## 2 Transformation and standardization/anomalization ##

(`trans-and-zscore.R` & `trans-and-norman.R`)

Charcoal data are reported in units that range over thirteen orders of magnitude, and charcoal records typically have "long-tailed" distributions (Power et al., 2010).  In order to compare or combine records, the data must therefore be transformed to approach normality (to reduce the impacts of non-constant variance) and rescaled to some common basis or range.  In one example here (`trans-and-zcore.R`), we apply the variance-stabilzing Box-Cox transformation, and rescale the data to "z-scores".  (For historical reasons, the transformed data are also rescaled using the "minimax" transformation so that all values lie between 0 and 1.  This reduces astonishment over transformed charcoal-influx or concentration values that may wind up being negative after transformation.)  This approach also requires the specification of a base period or time interval over which the transformation parameters are estimated and the mean and standard deviation used for calculating z-scores are calculated.  Further discussion of this approach can be found in Power et al. (2010) and Daniau et al. (2012).

A second example (`trans-and-norman.R`), illustrates the use of "normalized" anomalies, in which the deviations of the transformed charcoal influx values from a base period mean value are divided by that mean value, to produce a relative deviation, scaled by the overall level of the data.  This approach is useful for last-millennium type analyses, where records with few samples can produce standard deviations that are not robust, and hence z-scores that vary dramatically.

The specific tasks implemented by the script include, for each site: 

- censoring of samples with missing ages or ages after 2020 CE
- maximum likelihood estimation of of the Box-Cox transformation parameter `lambda`
- Box-Cox transformation of data
- minimax rescaling of the transformed data`tall`
- calculation z-scores `ztrans` or normalized anomalies `normans`
- writing out the transformed data for this site as a .csv file.

## 3 Presampling/prebinning ##

(`presample-bin.R`)

Charcoal data are available at all kinds of “native” resolutions, from samples that represent decades or centuries (or longer) to those that represent annual deposition. Further, some records have been interpolated to pseudo-annual time steps. In developing composite curves, those records with higher resolutions will contribute disproportionately to the curve. There are two general approaches for dealing with this: 1) weighting individual charcoal (influx or concentration) values according to their resolution, with lower-resolution records receiving higher weights, and vice-versa, or 2) reducing the sampling frequency of the records to some common interval (without interpolating or creating pseudo data).  We adopted the latter approach for its simplicity and transparency.

The binning is done by establishing a set of evenly spaced target points or bins, and then for each charcoal record, binning the individual observations.  If more that one observation falls in the same bin, the average (of he transformed and standardized data) is taken as the binned value.  No effort is made to interpolate between observations, to avoid pseudo-replication.  A .csv file is written out for each site.

## 4 Composite (i.e. smooth) curves and bootstrap C.I.'s ##

(`smooth-curve.R` & `bin-boot.R`)

The first script (`smooth-curve.R`) uses `locfit()` to get a (smooth) composite curve of the (presampled/binned) charcoal z-scores  or normans for a set of sites specified by an input "sitelist" (ultimately based on all or a subset of the sites listed in the `ALL_BART_SITES` query).  The smoothness of the curve is determined by the width of the smoothing window, customarily specified by the "half-width" (`hw`).  First, a "global" curve (in the sense of using all of the data from a particular list of sites is determined, and the number of sites with samples (`ndec_tot`, for historical reasons) and the number of samples that contribute to each fitted value (`ninwin_tot`) are also calculated.

Then, over `nreps` replications, the data are sampled by site (with replacement) to calculate bootstrap confidence intervals, and the upper and lower 95th-percentile confidence intervals are determined.  In the example here, the curves produced by individual bootstrap samples are plotted (in transparent gray), and the "global" curve is overplotted in red to give a visual indication of the uncertainty in the composite curve arising from the particular sample of sites.

The second script (`bin-boot.R`) creates a composite "curve" by directly binning each charcoal influx value in non-overlapping bins, and then calculating a simple average of the values in each bin.  This approach can be used over the past few millennia, where the median sample density is generally less than 20 years, with an appropriate bin width, also around 20 years.  This approach yields a composite curve that is more temporally variable than that provided by the local regression approach.

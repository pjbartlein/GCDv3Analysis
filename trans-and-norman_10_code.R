# trans-and-norman.R

# 1-parameter Box-Cox transformation of charcoal quanities for a single site
# (alpha (shift parameter) is specified, lambda (power transformation parameter) is estimated)

# This version calculates normalized anomalies ((tall-tmean)/tmean)

# input .csv files should contain at least the variables "est_age" and "quant",
# 	identified by labels in a header row  

# paths for input and output .csv files -- modify as appropriate
queryname <- "v3i"
datapath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/"
sitelistpath <- "/Projects/GPWG/GPWGv3/GCDv3Data/v3i/v3i_sitelists/"
sitelist <- "v3i_nsa_globe"

# set base period ages (or read a baseperiod info file)
basebeg <- 250 
baseend <- 1950 
basename <- "nt2kb"

# no changes below here

# various path and filenames
sitelistfile <- paste(sitelistpath, sitelist, ".csv", sep="")
sitelistfile

sitecsvpath <- paste(datapath,queryname,"_sites_csv/",sep="")
transcsvpath <- paste(datapath,queryname,"_trans_csv/",sep="")
# if output folder does not exist, create it
dir.create(file.path(datapath, paste(queryname,"_trans_csv/",sep="")), showWarnings=FALSE)
statscsvpath <- paste(datapath,queryname,"_stats/",sep="")
# if output folder does not exist, create it
dir.create(file.path(datapath, paste(queryname,"_stats/",sep="")), showWarnings=FALSE)
statsfile <- paste(statscsvpath,queryname,"_",basename,"_stats.csv", sep="")

# read list of sites
ptm <- proc.time()
sites <- read.csv(paste(sitelistfile), stringsAsFactors=FALSE)
nsites <- length(sites[,1])
print (nsites)

# storage for statistics
sn_save <- rep(0,nsites); lam_save <- rep(0,nsites); lik_save <- rep(0,nsites)
tmean_save <- rep(0,nsites); tstdev_save <- rep(0,nsites)

# main loop
for (j in seq(1,nsites)) {

  # 1. site .csv file name (input file)
  sitenum <- sites[j,1]
  sitename <- as.character(sites[j,5])
  siteidchar <- as.character(sitenum)
  if (sitenum >= 1) siteid <- paste("000", siteidchar, sep="")
  if (sitenum >= 10) siteid <- paste("00", siteidchar, sep="")
  if (sitenum >= 100) siteid <- paste("0", siteidchar, sep="")
  if (sitenum >= 1000) siteid <- paste(    siteidchar, sep="")
  inputfile <- paste(sitecsvpath, siteid, "_data.csv", sep="")
  print(j); print(sitenum); print(sitename); print(inputfile)

  # 2. read the input data
  sitedata <- read.csv(inputfile)
  
  # 3. discard samples with missing (-9999) ages
  sitedata <- sitedata[sitedata$est_age != -9999,]
  
  # 4. discard samples with ages > -60
  sitedata <- sitedata[sitedata$est_age > -70,]

  # 5. initial minimax rescaling of data
  minimax <- (sitedata$influx-min(sitedata$influx))/(max(sitedata$influx)-min(sitedata$influx))
  
  # 6. set `alpha` the Box-Cox transformation shift parameter
  alpha <- 0.01  # Box-Cox shift parameter
  # alternative alpha: 0.5 times the smallest nonzero value of influx
  # alpha <- 0.5*min(sitedata$influx[sitedata$influx != 0])  

  # 7. maximum likelihood estimation of lambda
  # derived from the boxcox.R function in the Venables and Ripley MASS library included in R 2.6.1

  npts <- 201 # number of estimates of lambda
  y <- minimax+alpha
  n <- length(y)
  logy <- log(y)
  ydot <- exp(mean(logy))
  lasave <- matrix(1:npts)
  liksave <- matrix(1:npts)
  for (i in 1:npts) {
    la <- -2.0+(i-1)*(4/(npts-1))
    if (la != 0.0) yt <- (y^la-1)/la else yt <- logy*(1+(la*logy)/2*(1+(la*logy)/3*(1+(la*logy)/4)))
    zt <- yt/ydot^(la-1)
    loglik <- -n/2*log(sum((zt - mean(zt))^2 ))
    lasave[i] <- la
    liksave[i] <- loglik
    }

  # save the maximum likelihood value and the associated lambda
  maxlh <- liksave[which.max(liksave)]
  lafit <- lasave[which.max(liksave)]
  print (c(sitenum, maxlh, lafit))

  # 8. Box-Cox transformation of data
  if (lafit == 0.0) tall <- log(y) else tall <- (y^lafit - 1)/lafit

  # 9. minimax rescaling
  tall <- (tall - min(tall))/(max(tall)-min(tall))
  
  # 10. calculate mean and standard deviation of data over base period
  tmean <- mean(tall[sitedata$est_age >= basebeg & sitedata$est_age <= baseend])
  tstdev <- sd(tall[sitedata$est_age >= basebeg & sitedata$est_age <= baseend])
  
  # 11. calculate "normans" normalized anomalies
  norman <- (tall-tmean)/tmean
  
  # 12. write out transformed data for this site
  siteout <- data.frame(cbind(sitedata[,1], sitedata$id_sample, sitedata$est_age,
    sitedata$depth, sitedata$influx, minimax, tall, norman))
  colnames(siteout) <- c("site_sample", "id_sample", "est_age", "depth", "influx", "influxmnx", "tall", "norman")

  outputfile <- paste(transcsvpath, siteid, "_trans_influx_", basename, ".csv", sep="")
  write.table(siteout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
  
  sn_save[j] <- sitenum
  lam_save[j] <- lafit
  lik_save[j] <- maxlh
  tmean_save[j] <- tmean
  tstdev_save[j] <- tstdev
  
}

# write out a file of statistics
stats <- data.frame(cbind(sn_save, lam_save, lik_save, tmean_save, tstdev_save))
colnames(stats) <- c("site", "lambda", "likelihood", "mean", "stdev")
write.table(stats, statsfile, col.names=TRUE, row.names=FALSE, sep=",")
proc.time() - ptm

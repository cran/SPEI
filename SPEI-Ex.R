pkgname <- "SPEI"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SPEI')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("PET")
### * PET

flush(stderr()); flush(stdout())

### Name: Potential evapotranspiration
### Title: Computation of potential evapotranspiration.
### Aliases: thornthwaite hargreaves penman

### ** Examples

# Load data for Tampa, lat=37.6475N, elevation=402.6 m. a.s.l.
# Data consists on monthly values since January 1980
data(wichita)
attach(wichita)
names(wichita)

# PET according to Thornthwaite
tho <- thornthwaite(TMED,37.6475)
# Hargreaves
har <- hargreaves(TMIN,TMAX,lat=37.6475)
# Penman, based on sun hours, ignore NAs
pen <- penman(TMIN,TMAX,AWND,tsun=TSUN,lat=37.6475,z=402.6,na.rm=TRUE)
# Penman, based on cloud cover
pen2 <- penman(TMIN,TMAX,AWND,CC=ACSH,lat=37.6475,z=402.6,na.rm=TRUE)
# Plot them together
plot(cbind(tho,har,pen,pen2))


# Now consider the data started in June 1900
thornthwaite(ts(TMED,start=c(1900,6),frequency=12),37.6475)



cleanEx()
nameEx("datasets")
### * datasets

flush(stderr()); flush(stdout())

### Name: Datasets
### Title: Datasets for illustrating the functions in the SPEI package.
### Aliases: wichita balance

### ** Examples

data(wichita)
names(wichita)
summary(wichita)

data(balance)
summary(balance)



cleanEx()
nameEx("kern")
### * kern

flush(stderr()); flush(stdout())

### Name: Kernel functions
### Title: Time kernel for computing the SPEI at different time scales.
### Aliases: kern kern.plot

### ** Examples

# A rectangular kernel with a time scale of 12 and no shift
kern(12)

# A gaussian kernel with a time scale of 12 and no shift
kern(12,'gaussian')

# Comparison of the four kernels, with and without shift
kern.plot(12)
kern.plot(12,2)



cleanEx()
nameEx("spei")
### * spei

flush(stderr()); flush(stdout())

### Name: Drought indices
### Title: Calculation of the Standardized Precipitation-Evapotranspiration
###   Index (SPEI) and the Standardized Precipitation Index (SPI).
### Aliases: spei spi

### ** Examples

# Load data
data(wichita)
attach(wichita)

# One and tvelwe-months SPI
spi_1 <- spi(PRCP,1)
spi_12 <- spi(PRCP,12)
plot(cbind(spi_1,spi_12))
# Notice that the first eleven values of spei_12 are NA

# One and tvelwe-months SPEI
wichita$PET <- thornthwaite(TMED,37.6475)
spei1 <- spei(PRCP-wichita$PET,1)
spei12 <- spei(PRCP-wichita$PET,12)
plot(cbind(spei1,spei12))

# Data series not starting in January: define the properties
# of the time series with ts()
plot(spei(ts(PRCP-wichita$PET,freq=12,start=c(1900,6)),12))

# Different kernels
spei24 <- spei(PRCP-wichita$PET,24)
spei24_gau <- spei(PRCP-wichita$PET,24,kernel=list(type='gaussian',shift=0))
plot(ts(cbind(spei24,spei24_gau),start=c(1900,1),freq=12))

# Several time series at a time
data(balance)
names(balance)
bal_spei12 <- spei(balance,12)
colnames(bal_spei12) <- names(balance)
plot(ts(bal_spei12[,1:6],start=c(1900,1),freq=12),main='12-month SPEI')
plot(ts(bal_spei12[,7:11],start=c(1900,1),freq=12),main='12-month SPEI')



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

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
###   Index
### Aliases: spei spi

### ** Examples

# Load data
data(wichita)

# One and tvelwe-months SPEI
wichita$PET <- thornthwaite(wichita$TMED,37.6475)
spei1 <- spei(wichita$PRCP-wichita$PET,1)
spei12 <- spei(wichita$PRCP-wichita$PET,12)
# not executed
#spei1
#summary(spei1)
#coefficients(spei1)
par(mfrow=c(2,1))
plot(spei1)
plot(spei12,'Wichita')

# One and tvelwe-months SPI
spi_1 <- spi(wichita$PRCP,1)
spi_12 <- spi(wichita$PRCP,12)
par(mfrow=c(2,1))
plot(spi_1)
plot(spi_12)

# Define the properties of the time series with ts()
plot(spei(ts(wichita$PRCP-wichita$PET,freq=12,start=c(1980,1)),12))

# Time series not starting in January
plot(spei(ts(wichita$PRCP-wichita$PET,freq=12,start=c(1980,6)),12))

# Using a particular reference period
plot(spei(ts(wichita$PRCP-wichita$PET,freq=12,start=c(1980,1)),12,
	ref.start=c(1980,1), ref.end=c(2000,1)))

# Different kernels
spei24 <- spei(wichita$PRCP-wichita$PET,24)
spei24_gau <- spei(wichita$PRCP-wichita$PET,24,kernel=list(type='gaussian',shift=0))
par(mfrow=c(2,1))
plot(spei24,'Rectangular kernel')
plot(spei24_gau,'Gaussian kernel')


# Several time series at a time
data(balance)
names(balance)
bal_spei12 <- spei(balance,12)
plot(bal_spei12)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
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

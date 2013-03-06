# Computation of the Standardized Precipitation Index (SPI).
#

# Generic function
spi <- function(x, y,...) UseMethod('spi')

# Fit SPI (previously spi() function). Default method.
spi <-
function(data, scale, kernel=list(type='rectangular',shift=0),
	distribution='Gamma', fit='ub-pwm', na.rm=FALSE,
	ref.start=NULL, ref.end=NULL, x=FALSE, params=NULL, ...) {
	return(spei(data, scale, kernel, distribution, fit, na.rm,
	ref.start, ref.end, x, params))
}

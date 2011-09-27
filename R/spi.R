# Computation of the Standardized Precipitation Index (SPI).
#
spi <-
function(data, scale, kernel=list(type='rectangular',shift=0),
	distribution='Gamma', fit='ub-pwm', na.rm=FALSE) {
	return(spei(data, scale, kernel, distribution, fit, na.rm))
}

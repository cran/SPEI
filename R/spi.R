# Computation of the Standardized Precipitation Index (SPI).
#

# Generic function
spi <- function(x, y,...) UseMethod('spi')

# Fit SPI (previously spi() function). Default method.
spi <-
function(data, scale, kernel=list(type='rectangular',shift=0),
	distribution='Gamma', fit='ub-pwm', na.rm=FALSE,
	ref.start=NULL, ref.end=NULL, x=FALSE, ...) {
	return(spei(data, scale, kernel, distribution, fit, na.rm,
	ref.start, ref.end, x))
}

# Print method
print.spi <- function (x, ...) {
	cat('Call:\n')
	print(x$call)
	cat('\nFitted:\n')
	print(x$fitted)
}

# Summary method
summary.spi <- function (object, ...) {
	x <- object
	cat('Call:\n')
	print(x$call)
	cat('\nCoefficients:\n')
	for (i in 1:dim(x$coeff)[2]) {
		cat('\t',dimnames(x$coeff)[[2]][i],':\n',sep='')
		tab <- cbind(t(x$coeff[,i,]))
		rownames(tab) <- 1:dim(x$coeff)[3]
		print(tab)
		cat('\n')
	}
}

# Plot method
plot.spi <- function (x, tit=NULL, ...) {
	ser <- ts(x$fitted[-c(1:x$scale),],end=end(x$fitted),fr=frequency(x$fitted))
	tit <- dimnames(x$coefficients)[2][[1]]
	#
	if (start(ser)[2]==1) {
		ns <- c(start(ser)[1]-1,12)
	} else {
		ns <- c(start(ser)[1],start(ser)[2]-1)	
	}
	if (end(ser)[2]==12) {
		ne <- c(end(ser)[1]+1,1)
	} else {
		ne <- c(end(ser)[1],end(ser)[2]+1)
	}
	#
	n <- ncol(ser)
	par(mar=c(4,4,2,1)+0.1)
	if (n>1 & n<5) par(mfrow=c(n,1))
	if (n>1 & n>=5) par(mfrow=c({n+1}%/%2,2))
	for (i in 1:n) {
		datt <- ts(c(0,ser[,i],0),freq=frequency(ser),start=ns,end=ne)
		datt.pos <- ifelse(datt>0,datt,0)
		datt.neg <- ifelse(datt<=0,datt,0)
		plot(datt,type='n',xlab='',ylab='SPI',main=tit[i])
		grid(col='black')
		polygon(datt.pos,col='blue',border=NA)
		polygon(datt.neg,col='red',border=NA)
		lines(datt,col='dark grey')
		abline(h=0)
		if (!is.null(x$ref.period)) {
			abline(v=x$ref.period[1,])
			abline(v=x$ref.period[2,])
		}
	}
}

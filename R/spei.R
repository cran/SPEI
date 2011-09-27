# Computation of the Standardized Precipitation-Evapotranspiration Index (SPEI).
#
spei <- function(data, scale, kernel=list(type='rectangular',shift=0),
	distribution='log-Logistic', fit='ub-pwm', na.rm=FALSE) {

	require(lmomco)

	#if (!exists("data",inherits=F) | !exists("scale",inherits=F)) {
	#	stop('Both data and scale must be provided')
	#}
	if (sum(is.na(data))>0 & na.rm==FALSE) {
		stop('Error: Data must not contain NAs')
	}
	if (distribution!='log-Logistic' & distribution!='Gamma' & distribution!='PearsonIII') {
		stop('Distrib must be one of "log-Logistic", "Gamma" or "PearsonIII"')
	}
	if (fit!='max-lik' & fit!='ub-pwm' & fit!='pp-pwm') {
		stop('Method must be one of "ub-pwm" (default), "pp-pwm" or "max-lik"')
	}

	if (!is.ts(data)) {
		data <- ts(as.matrix(data),freq=12)
	}
	m <- ncol(data)
	fr <- frequency(data)
	std <- data*NA
	
	# Loop through series (columns in data)
	for (s in 1:m) {
		# Cumulative series (acu)
		acu <- std[,s]
		if (scale>1) {
			wgt <- kern(scale,kernel$type,kernel$shift)
			for (t in scale:length(acu)) {
				acu[t] <- sum(data[t:{t-scale+1},s]*wgt)
			} # next t
		} else {
			acu <- data[,s]
		}

		# Loop through the seasons
		for (c in (1:fr)) {
			# Filter month m, excluding NAs
			f <- seq(c,length(acu),fr)
			f <- f[!is.na(acu[f])]
				
			# Monthly series, sorted
			month <- sort(acu[f])

			if (length(month)==0 | is.na(sd(month,na.rm=TRUE))) {
				std[f] <- NA
				next()
			}
		
			if (fit=='pp-pwm') {
				pwm <- pwm.pp(month,-0.35,0)
			} else {
				pwm <- pwm.ub(month)
			}
			lmom <- pwm2lmom(pwm)
			if (!are.lmom.valid(lmom) | is.nan(sum(lmom[[1]]))) {
				next()
			}
	
			if (distribution=='log-Logistic') {
				# Fit a generalized log-Logistic distribution
				llpar <- parglo(lmom)
				if (fit=='max-lik') {
					llpar <- parglo.maxlik(month,llpar$para)
				}
				# Compute standardized values
				std[f,s] <- qnorm(pglo(acu[f],llpar))
			} else {
				# Probability of monthly precipitation = 0 (pze)
				zeros <- sum(month==0)
				pze <- sum(month==0)/length(month)
# 				month <- sort(month)
				if (distribution =='Gamma') {
					# Fit a Gamma distribution
					gampar <- pargam(lmom.ub(month))
					# Compute standardized values
					std[f,s] <- qnorm(cdfgam(acu[f],gampar))
					std[f,s] <- qnorm(pze + (1-pze)*pnorm(std[f,s]))
				} else if (distribution =='PearsonIII') {
					# Fit a PearsonIII distribution
					p3par <- parpe3(lmom.ub(month))
					# Compute standardized values
					std[f,s] <- qnorm(cdfpe3(acu[f],p3par))
					std[f,s] <- qnorm(pze + (1-pze)*pnorm(std[f,s]))
				} # end if
			} # end if
		} # next c (month)

		#std[is.nan(std[,s]) | is.nan(std[,s]-std[,s]),s] <- NA
		#std[,s] <- std[,s]-mean(std[,s],na.rm=TRUE)
		#std[,s] <- std[,s]/sd(std[,s],na.rm=TRUE)
	} # next s (series)
	colnames(std) <- rep('SPEI',m)

	return(std)
}

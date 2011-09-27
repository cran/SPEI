\name{Drought indices}

\alias{spei}
\alias{spi}

\title{Calculation of the Standardized Precipitation-Evapotranspiration Index (SPEI) and the Standardized Precipitation Index (SPI).}

\description{Given a time series of the climatic water balance (precipitation minus potential evapotranspiration), gives a time series of the Standardized Precipitation-Evapotranspiration Index (SPEI).}

\usage{
spei(data, scale, kernel = list(type = 'rectangular', shift = 0),
	distribution = 'log-Logistic', fit = 'ub-pwm', na.rm = FALSE)

spi(data, scale, kernel = list(type = 'rectangular', shift = 0),
	distribution = 'Gamma', fit = 'ub-pwm', na.rm = FALSE)
}

\arguments{
  \item{data}{ a vector, matrix or data frame with time ordered values of precipitation (for the SPI) or of the climatic balance precipitation minus potential evapotranspiration (for the SPEI).}
  \item{scale}{ an integer, representing the time scale at which the SPEI / SPI will be computed.}
  \item{kernel}{ optional, a list defining the type of kernel used for computing the SPEI / SPI at scales higher than one.}
  \item{distribution}{ optional, name of the distribution function to be used for computing the SPEI / SPI (one of 'log-Logistic', 'Gamma' and 'PearsonIII').}
  \item{fit}{ optional, name of the method used for computing the distribution function parameters (one of 'ub-pwm', 'pp-pwm' and 'max-lik').}
  \item{na.rm}{ optional, a logical value indicating whether NA values should be stripped from the computations.}
}

\details{
The \code{spei} and \code{spi} functions allow computing the SPEI and the SPI indices. These are widely used for drought identification and monitoring. Both functions are identical (in fact, \code{spi} is just a wrapper for \code{spei}), but they are kept separated for clarity. Basically, the functions standardize a variable following a log-Logistic distribution function (i.e., they transform it to a standard Gaussian variate with zero mean and standard deviation of one). The input variable is a time ordered series of precipitation values for \code{spi}, and a series of the climate water balance (precipitation minus potential evapotranspiration) for \code{spei}. When used with the default options, they yield values of both indices as defined in the references given below.

Originally the SPEI and the SPI were defined for monthly data. Since the PDFs of the data are not homogenous from month to month, the data is split into twelve series (one for each month) and independent PDFs are fit to each series. If \code{data} is a vector or a matrix it will be treated as a sequence of monthly values starting in January. If it is a (univariate or multivariate) time series then the function \code{\link{cycle}} will be used to determine the position of each observation within the year (month), allowing the data to start in a month different than January.

An important advantage of the SPEI and the SPI is that they can be computed at different time scales. This means that it is possible to incorporate the influence of the past values of the variable in its computation, enabling the index to adapt to the memory of the system under study. The magnitude of this memory is controlled by parameter \code{scale}: for example, a value of six would imply that data from the current month and of the past five months will be used for computing the SPEI or SPI value for a given month. By default all past data will have the same weight in computing the index, as it was originally proposed by the authors. Other kernels, however, are available through parameter \code{kernel}. The parameter \code{kernel} is a list defining the shape of the kernel and a time shift. For more details refer to  the function \code{\link{kern}}.

Following the original definitions \code{spei} uses a log-Logistic distribution by default, and \code{spi} uses a Gamma distribution. This behaviour can be modified, however, through parameter \code{distribution}.

The default method for parameter fitting is based on unbiased Probability Weighted Moments ('ub-pwm'), but other methods can be used through parameter \code{fit}. Valid aleternatives are the plotting-position PWM ('pp-pwm') and the maximum likelihood method ('max-lik').
}

\value{
A time series with the values of the Standardized Precipitation-Evapotranspiration Index (SPEI) or the Standardized Precipitation Index (SPI). If data consists of several columns the function will treat each column as independent data, and the result will be a multivariate time series.

If a time scale greater than one has been used then the first elements will have NA value since the kernel can not be applied. The first element with valid data will be the one corresponding to the time scale chosen.
}

\references{
S.M. Vicente-Serrano, S. Beguería, J.I. López-Moreno. 2010. A Multi-scalar drought index sensitive to global warming: The Standardized Precipitation Evapotranspiration Index – SPEI. \emph{Journal of Climate} \bold{23}: 1696, DOI: 10.1175/2009JCLI2909.1.

\url{http://sac.csic.es/spei/}
}

\author{
Santiago Beguería and Sergio M. Vicente-Serrano. Maintainer: Santiago Beguería
}

\note{Dependencies: the spei function depends on the library \code{lmomco}.}

\seealso{\code{\link{kern}} for different kernel functions available. \code{\link{thornthwaite}}, \code{\link{hargreaves}} and \code{\link{penman}} for ways of calculating potential evapotranspiration.}

\examples{
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
}

%\keyword{}
%\keyword{}

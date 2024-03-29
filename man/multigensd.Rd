% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gensd.R
\name{multigensd}
\alias{multigensd}
\title{Multiple Attributes skewed and dependent}
\usage{
multigensd (n.attr, size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  method = "SAR", dst = "skewnormal")
}
\arguments{
\item{n.attr}{Number of attributes to be generated.}

\item{size}{a vector of two elements specifying respectively numbers of indivuals and time points of the attribute to be generated.}

\item{grid}{number of rows and columns for grids.}

\item{sk}{a scalar that represents the skewness coefficient to be considered.}

\item{dep}{a scalar that represents the dependence correlation between values within a given attribute.}

\item{mu}{the mean value}

\item{v}{the variance value}

\item{method}{a character string naming the the approach used to create spatial dependence between variables. The "SAR" method is the Kapoor, Kelejian, and Prucha (2007)-spatial autoregressive random effects (SAR-RE) model. The "SARMA" method represents the Lee and Yu (2012) spatial autoregressive moving average random effects (SARMA-RE) model.}

\item{dst}{a character string naming the distribution to be used for vectors generation. "skewnormal" the skew-normal distribution is limited to a range of [0,1[ for the coefficient of skewness. "lognormal" the lognormal distribution allows a large range of coefficients skewness values.}
}
\value{
Returns a list.
\item{mat.attributes}{Multiple spatial dependent and skewed attributes generated.}
\item{spatial.dependence}{Racall of the spatial dependence existing within values of a given attribute generated}
\item{skewness.coef}{Coefficient of skewness of spatial attributes generated.}
}
\description{
Generates multiple spatial dependent attributes with  fixed skewness.
}
\examples{

Consider a simulation study that requires 10 skewed (3.6) spatial attribbutes with a spatial autoregressive coefficient (spatial correlation) of 0.5.
Each attribute has 100 subjects and 7 observations (time points) per subject. We choose to design 10 x 10 square grids.

# Loading package
library(genspatattr)

# As the skewness required exceeds the interval [0,1], we use the lognormal distribution and SAR method to maintain null the spatial moving average coefficient.
k1 <- c(100, 7)
k2 <- c(10, 10)
dt <- multigensd(n.attr = 10, size = k1, grid = k2, sk = 3.6, dst = "lognormal")
> head(dt$mat.attributes)
       [,1]     [,2]      [,3]     [,4]      [,5]      [,6]      [,7]      [,8]      [,9]
[1,] 36.76639 29.65413 321.04957 54.14535  5.620346  38.27667 118.73450  44.93622  56.44374
[2,] 21.99757 27.46571  95.71541 37.93883  8.538977  82.19352  44.64520  39.25272  45.69996
[3,] 40.25177 32.81329  31.67804 34.86756 18.055213 384.21476  20.64373  42.84450  62.49243
[4,] 20.91775 50.49171  17.53125 35.79660 49.848041 171.73468  35.92800  44.00423 138.66885
[5,] 17.16706 37.16105  12.95207 50.99651 38.020209 159.57852  26.01867  75.07399  55.96092
[6,] 13.39583 35.46820  15.07685 36.63703 71.563533  43.86750  73.75986 173.56775  44.58139
      [,10]
[1,] 24.85834
[2,] 25.51857
[3,] 21.10237
[4,] 23.60568
[5,] 30.25069
[6,] 28.92886

> e1071::skewness(dt$mat.attributes[,1])
[1] 3.643574
> e1071::skewness(dt$mat.attributes[,10])
[1] 3.606651
> dt$spatial.dependence
[1] 0.5
> dt$skewness.coef
[1] 3.6  
}
\author{
Sewanou Honfo \url{<honfosewanou@gmail.com>}  and Brezesky Kotanmi \url{<kotangaebrezesky@gmail.com>} and Eunice Gnanvi \url{<eunysgnanvi@gmail.com>} and Chenangnon Tovissode \url{<chenangnon@gmail.com>} and Romain glèlè Kakaï \url{<glele.romain@gmail.com>}   / LABEF_08_2019
}

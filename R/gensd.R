#' @name gsd
#' @title Spatial dependent attribute with unknown skewness
#' @description Generates a spatial dependent attributes of which skewness is provided as output.
#'
#'
#' @usage gsd (size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  method = "SAR", dst = "skewnormal")
#'
#'
#' @param size a vector of two elements specifying respectively numbers of indivuals and time points of the attribute to be generated.
#' @param grid number of rows and columns for grids.
#' @param sk a scalar that represents the skewness coefficient to be considered.
#' @param dep a scalar that represents the dependence correlation between values within a given attribute.
#' @param mu the mean value
#' @param v the variance value
#' @param dst a character string naming the distribution to be used for vectors generation. "skewnormal" the skew-normal distribution is limited to a range of [0,1[ for the coefficient of skewness. "lognormal" the lognormal distribution allows a large range of coefficients skewness values.   
#' @param method a character string naming the the approach used to create spatial dependence between variables. The "SAR" method is the Kapoor, Kelejian, and Prucha (2007)-spatial autoregressive random effects (SAR-RE) model. The "SARMA" method represents the Lee and Yu (2012) spatial autoregressive moving average random effects (SARMA-RE) model.
#' @return Returns a list:
#' \item{spdata}{the generated attribute.}
#' \item{skew.coef}{effective skewness coefficient from generated data}
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>}  and Brezesky Kotanmi \url{<kotangaebrezesky@gmail.com>} and Eunice Gnanvi \url{<eunysgnanvi@gmail.com>} and Chenangnon Tovissode \url{<chenangnon@gmail.com>} and Romain glèlè Kakaï \url{<glele.romain@gmail.com>}   / LABEF_08_2019
#'
#' @references Mohadeseh Alsadat Farzammehr, Mohammad Reza Zadkarami & Geoffrey J. McLachlan (2019): Skew-normal generalized spatial panel data model, Communications in Statistics - Simulation and Computation, DOI: 10.1080/03610918.2019.1622718
#'
#'

gsd <- function (size, grid, sk = 0, dep = 0.5, mu = 0, v = 1, method = "SAR", 
                 dst= "skewnormal")
{
  
  mk <- match.call()
  miSiz <- missing(size)
  if(miSiz){
    stop("You must provide the size (number of subjects and number of observations per
         subject")
  }
  
  if(!miSiz && !is.vector(size)){
    stop("The size must be a vector")
  }
  
  if(!miSiz && is.vector(size) && length(size) != 2){
    stop("The size must be a vector with two elements: numbers of subjects
         and observations per subject, respectively.")
  }
  
  miGrid <- missing(grid)
  if(miGrid){
    stop("You must provide the grid (numbers of rows and columns")
  }
  
  if(!miGrid && !is.vector(grid)){
    stop("The grid must be a vector")
  }
  
  if(!miGrid && is.vector(grid) && length(grid) != 2){
    stop("The grid must be a vector with two elements: numbers of rows and columns,
          respectively.")
  }
  
  
  if(size[1] != grid[1] * grid[2]){
    stop("The number of observations must be equal to the product
         of numbers of rows and columns of the grid matrix")
  }
  
  if(sk < 0 || sk > 1 && dst == "skewnormal"){
    warning("The skewnormal distribution does not work on this skewness coefficient.
            You may choose the lognormal distribution instead.")
    return(gsd(size = size, grid = grid, sk = sk, dep = dep, mu = mu, v = v,  
                 method = method, dst = "lognormal"))
  }
  
  
  if(method == "SAR"){
    varrho <- numeric(2)
  }
  if(method == "SARMA"){
    varrho <- runif(n = 2, min = -1, max = 1)
  }
  
  k1 <- spdep::cell2nb(grid[1], grid[2])
  k2 <- spdep::nb2listw(k1, style="W")
  k3 <- spdep::invIrW(k2, rho = dep, method="solve", feasible=NULL)
  
  B1 <- B2 <- solve(k3)
  W <- (diag(dim(B1)[1]) - B1) / dep
  
  D1 <- diag(dim(B1)[1]) + varrho[1] * W
  D2 <- diag(dim(B2)[1]) + varrho[2] * W
  
  A1 <- solve(D1) %*% B1
  A2 <- solve(D2) %*% B2
  
  if(dst == "skewnormal"){
  
  ic <- ((2 * sk) / (4 - pi))^(1/3)
  theta_k <- ic / (sqrt(1 + ic^2))
  theta_t <- theta_k / sqrt(2 / pi)
  alph <- theta_t / sqrt(1 - theta_t^2)
  muz <- sqrt(2 / pi) * (alph / (sqrt(1 + alph^2)))
  
  sc <- sqrt(v / (1 - muz^2))
  loc <- mu - sc * muz
  mui <- VGAM::rskewnorm(n = size[1], location = loc, scale = sc, shape = alph)
  
  sc2 <- 4 * sqrt(v / (1 - muz^2))
  loc2 <- mu - sc2 * muz
  nuit <- VGAM::rskewnorm(n = size[1], location = loc2, scale = sc2, shape = alph)
  
  }
  
  if(dst == "lognormal"){
    
    sig <- function(gamma1){
      
      X <- polyroot(c(-gamma1^2 - 4, 0, 3, 1))
      X <- Re(X[Im(X)^2 <= min(Im(X)^2)])
      X <- log(X)
      X
    }
    
    sigma2 <- sig(sk)
    mui <- rlnorm(n = size[1], mean = mu, sd = sqrt(sigma2))
    nuit <- rlnorm(n = size[1], mean = mu, sd = 4 * sqrt(sigma2))
  }
  
  u1 <- solve(A1) %*% mui
  u2t <- solve(A2) %*% nuit
  ut <- u1 + u2t
  
  md <- function(N, t, u){
    
    x01 <-5 + 10 * runif(n = 1, min = -0.5, max = 0.5)
    x02 <-5 + 10 * runif(n = 1, min = -0.5, max = 0.5)
    x1 <- x2 <- yt <- matrix(numeric(N * t), nrow = N)
    
    
    x1[, 1] <- 0.1  + 0.5 * x01 + runif(n = 1, min = -0.5, max = 0.5)
    x2[, 1] <- 0.1  + 0.5 * x02 + runif(n = 1, min = -0.5, max = 0.5)
    yt[, 1] <- 0.5 * x1[, 1] + 7 * x2[, 1] + u
    
    for(i in 2 : t){
      
      x1[,i] <- 0.1 * i + 0.5 * x1[,i-1] + runif(n = 1, min = -0.5, max = 0.5)
      x2[,i] <- 0.1 * i + 0.5 * x2[,i-1] + runif(n = 1, min = -0.5, max = 0.5)
      
      yt[,i] <- 0.5 * x1[, 1] + 7 * x2[, 1] + u
      
    }
    
    return(c(yt))
  }
  
  y <- md(N = size[1], t = size[2], u = ut)
  
  
  return(list(spdata = y, skew.coef = round(e1071::skewness(y), 1)))

}




#' @name gensd
#' @title Spatial dependent attribute with fixed skewness
#' @description Generates a spatial dependent attribute with a fixed skewness.
#'
#'
#' @usage gensd (size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  method = "SAR", dst = "skewnormal")
#'
#'
#' @param size a vector of two elements specifying respectively numbers of indivuals and time points of the attribute to be generated.
#' @param grid number of rows and columns for grids.
#' @param sk a scalar that represents the skewness coefficient to be considered.
#' @param dep a scalar that represents the dependence correlation between values within a given attribute.
#' @param mu the mean value
#' @param v the variance value
#' @param dst a character string naming the distribution to be used for vectors generation. "skewnormal" the skew-normal distribution is limited to a range of [0,1[ for the coefficient of skewness. "lognormal" the lognormal distribution allows a large range of coefficients skewness values.   
#' @param method a character string naming the the approach used to create spatial dependence between variables. The "SAR" method is the Kapoor, Kelejian, and Prucha (2007)-spatial autoregressive random effects (SAR-RE) model. The "SARMA" method represents the Lee and Yu (2012) spatial autoregressive moving average random effects (SARMA-RE) model.
#' @return Returns a gsd object.
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>}  and Brezesky Kotanmi \url{<kotangaebrezesky@gmail.com>} and Eunice Gnanvi \url{<eunysgnanvi@gmail.com>} and Chenangnon Tovissode \url{<chenangnon@gmail.com>} and Romain glèlè Kakaï \url{<glele.romain@gmail.com>}   / LABEF_08_2019
#'
#'
#'


gensd <- function(size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  
                  method = "SAR", dst = "skewnormal")
{
  sk.calc <- sk + 100
  i <- 0
  n.iter <- 0
  
  while(sk.calc != sk){
    
    spd <- gsd(size = size, grid = grid, sk = sk, dep = dep, mu = mu, v = v, method = method, dst = dst)
    
    sk.calc <- spd$skew.coef
    
    if(sk.calc == sk) break
    if(n.iter == 10000) stop("Impossible to generate spatial attribute with the required skewness coefficient")
    
    i <- i + 1
    n.iter <- n.iter + 1
    
  }
  
  return(spd)
}



#' @name multigensd
#' @title Multiple Attributes skewed and dependent
#' @description Generates multiple spatial dependent attributes with  fixed skewness.
#'
#'
#' @usage multigensd (n.attr, size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  method = "SAR", dst = "skewnormal")
#'
#' @param n.attr Number of attributes to be generated.
#' @param size a vector of two elements specifying respectively numbers of indivuals and time points of the attribute to be generated.
#' @param grid number of rows and columns for grids.
#' @param sk a scalar that represents the skewness coefficient to be considered.
#' @param dep a scalar that represents the dependence correlation between values within a given attribute.
#' @param mu the mean value
#' @param v the variance value
#' @param dst a character string naming the distribution to be used for vectors generation. "skewnormal" the skew-normal distribution is limited to a range of [0,1[ for the coefficient of skewness. "lognormal" the lognormal distribution allows a large range of coefficients skewness values.   
#' @param method a character string naming the the approach used to create spatial dependence between variables. The "SAR" method is the Kapoor, Kelejian, and Prucha (2007)-spatial autoregressive random effects (SAR-RE) model. The "SARMA" method represents the Lee and Yu (2012) spatial autoregressive moving average random effects (SARMA-RE) model.
#' @return Returns a list.
#' \item{mat.attributes}{Multiple spatial dependent and skewed attributes generated.}
#' \item{spatial.dependence}{Racall of the spatial dependence existing within values of a given attribute generated}
#' \item{skewness.coef}{Coefficient of skewness of spatial attributes generated.}
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>}  and Brezesky Kotanmi \url{<kotangaebrezesky@gmail.com>} and Eunice Gnanvi \url{<eunysgnanvi@gmail.com>} and Chenangnon Tovissode \url{<chenangnon@gmail.com>} and Romain glèlè Kakaï \url{<glele.romain@gmail.com>}   / LABEF_08_2019
#' 
#'
#' @examples
#'
#' Consider a simulation study that requires 10 skewed (3.6) spatial attribbutes with a spatial autoregressive coefficient (spatial correlation) of 0.5.
#' Each attribute has 100 subjects and 7 observations (time points) per subject. We choose to design 10 x 10 square grids.
#'
#' # Loading package
#' library(genspatattr)
#' 
#' # As the skewness required exceeds the interval [0,1], we use the lognormal distribution and SAR method to maintain null the spatial moving average coefficient.
#' k1 <- c(100, 7)
#' k2 <- c(10, 10)
#' dt <- multigensd(n.attr = 10, size = k1, grid = k2, sk = 3.6, dst = "lognormal")
#' > head(dt$mat.attributes)
#'        [,1]     [,2]      [,3]     [,4]      [,5]      [,6]      [,7]      [,8]      [,9]
#' [1,] 36.76639 29.65413 321.04957 54.14535  5.620346  38.27667 118.73450  44.93622  56.44374
#' [2,] 21.99757 27.46571  95.71541 37.93883  8.538977  82.19352  44.64520  39.25272  45.69996
#' [3,] 40.25177 32.81329  31.67804 34.86756 18.055213 384.21476  20.64373  42.84450  62.49243
#' [4,] 20.91775 50.49171  17.53125 35.79660 49.848041 171.73468  35.92800  44.00423 138.66885
#' [5,] 17.16706 37.16105  12.95207 50.99651 38.020209 159.57852  26.01867  75.07399  55.96092
#' [6,] 13.39583 35.46820  15.07685 36.63703 71.563533  43.86750  73.75986 173.56775  44.58139
#'       [,10]
#' [1,] 24.85834
#' [2,] 25.51857
#' [3,] 21.10237
#' [4,] 23.60568
#' [5,] 30.25069
#' [6,] 28.92886
#' 
#' > e1071::skewness(dt$mat.attributes[,1])
#' [1] 3.643574
#' > e1071::skewness(dt$mat.attributes[,10])
#' [1] 3.606651
#' > dt$spatial.dependence
#' [1] 0.5
#' > dt$skewness.coef
#' [1] 3.6  


multigensd <- function(n.attr, size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  
                       method = "SAR", dst = "skewnormal")
{
  
  if(n.attr == 1) {
    warning("You may simply use the function gensd to generate one attribute.")
    return(gensd(size = size, grid = grid, sk = sk, dep = dep, mu = mu, v = v,  
                 method = method, dst = dst))
  }
  
  if(n.attr <= 0){
    stop("This function requires a non null positive integer to perform generation
         of attribute(s)")
  }
  
  at <- matrix(numeric(size[1] * size[2] * n.attr), ncol = n.attr)
  
  for(i in 1 : n.attr){
    
    k <- gensd(size = size, grid = grid, sk = sk, dep = dep, mu = mu, v = v,
                     method = method, dst = dst)
  at[, i] <- k$spdata
    
  }
  
  return(list(mat.attributes = at, spatial.dependence = dep, skewness.coef = sk))
}
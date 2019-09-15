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
#'
#' @examples
#'

gsd <- function (size, grid, sk = 0, dep = 0.5, mu = 0, v = 1, method = "SAR", dst= "skewnormal")
{
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
#' @examples
#'


gensd <- function(size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  method = "SAR", dst = "skewnormal")
{
  sk.calc <- sk + 100
  i <- 0
  
  
  while(sk.calc != sk){
    
    spd <- gsd(size = size, grid = grid, sk = sk, dep = dep, mu = mu, v = v, method = method, dst = dst)
    
    sk.calc <- spd$skew.coef
    
    if(sk.calc == sk) break
     
    i <- i + 1
    
    
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


multigensd <- function(n.attr, size, grid, sk = 1, dep = 0.5, mu = 0, v = 1,  method = "SAR", dst = "skewnormal")
{
  
  at <- matrix(numeric(size[1] * size[2] * n.attr), ncol = n.attr)
  
  for(i in 1 : n.attr){
    
    k <- gensd(size = size, grid = grid, sk = sk, dep = dep, mu = mu, v = v,
                     method = method, dst = dst)
  at[, i] <- k$spdata
    
  }
  
  return(list(mat.attributes = at, spatial.dependence = dep, skewness.coef = sk))
}
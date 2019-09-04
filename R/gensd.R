#' @name gensd
#' @title Attributes skewed and dependent
#' @description Generates spatial attributes that are skewed and dependent.
#'
#'
#' @usage     gensd (size, sk = 1, dep = 0.5, mu = 0, v = 1,  method = "standard", dst = "skewnormal", ...)
#'
#'
#' @param size a vector of two elements specifying respectively the number of rows and columns of the matrix of attributes to be generated.
#' @param sk a scalar that represents the skewness coefficient to be considered.
#' @param dep a scalar that represents the dependence correlation between values within a given attribute.
#' @param mu the mean value
#' @param v the variance value
#' @param method a character string naming the the approach used to create spatial dependence between variables. The "standard" method centers two vectors and residuals from a linear model that fits both vectors. The "orthogonal" method centers two vectors et computes their correlation based on their angle cosine.
#' @param dst a character string naming the distribution to be used for vectors generation. "skewnormal" the skew-normal distribution is limited to a range of [0,1[ for the coefficient of skewness. "lognormal" the lognormal distribution allows a large range of coefficients skewness values.   
#' 
#' @return Returns a list:
#' \item{spadata}{a matrix containing expected data.}
#' \item{init}{a recall of initial parameters}
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>}  and Brezesky Kotanmi \url{<kotangaebrezesky@gmail.com>} and Eunice Gnanvi \url{<eunysgnanvi@gmail.com>} and Chenangnon Tovissode \url{<chenangnon@gmail.com>} and Romain glèlè Kakaï \url{<glele.romain@gmail.com>}   / LABEF_08_2019
#'
#'
#' @examples
#'

gensd <- function (size, sk = 0, dep = 0.5, mu = 0, v = 1, method = "standard", 
                   dst = "skewnormal", ...)
{
  
  biv <- function(k, y1, y2, r, method = "standard"){
    
    if(method == "standard"){
      
      y <- scale(y2) * r  +  scale(residuals(lm(y1~y2))) * sqrt(1 - (r * r))
    }
    
    if(method == "orthogonal"){
      
      Y <- cbind(y2, y1)         
      Yctr <- scale(Y, center = TRUE, scale = FALSE)   
      Q <- tcrossprod(qr.Q(qr(Yctr[ , 1, drop=FALSE])))      
      Yo <- (diag(k) - Q) %*% Yctr[ , 2]                 
      Yc <- cbind(Yctr[ , 1], Yo)               
      Yf  <- Yc %*% diag(1 / sqrt(colSums(Yc^2)))  
      
      y <- Yf[ , 2] + (1 / tan(acos(r))) * Yf[ , 1]
    }
    
    return(y)
  }
  
    
  if(dst == "skewnormal"){
      
      ic <- ((2 * sk) / (4 - pi))^(1/3)
      theta_k <- ic / (sqrt(1 + ic^2))
      theta_t <- theta_k / sqrt(2 / pi)
      alph <- theta_t / sqrt(1 - theta_t^2)
      
      muz <- sqrt(2 / pi) * (alph / (sqrt(1 + alph^2)))
      sc <- sqrt(v / (1 - muz^2))
      loc <- mu - sc * muz
      
      x1 <- VGAM::rskewnorm(n = size[2], location = loc, scale = sc, shape = alph)
      x2 <- VGAM::rskewnorm(n = size[2], location = loc, scale = sc, shape = alph)
     
    }
    
    if(dst == "lognormal"){
      
      sig <- function(gamma1){
        
        X <- polyroot(c(-gamma1^2 - 4, 0, 3, 1))
        X <- Re(X[Im(X)^2 <= min(Im(X)^2)])
        X <- log(X)
        X
      }
      
      sigma2 <- sig(sk)
      x1 <- rlnorm(n = size[2], mean = mu, sd = sqrt(sigma2))
      x2 <- rlnorm(n = size[2], mean = mu, sd = sqrt(sigma2))
    }
  
  for(i in 1 : size[1]){
    
    res[,i] <- biv(k = size[2], y1 = x1, y2 = x2, r = dep, method = "standard")
    
  }
  
  gdt <- t(cbind(x2, res))
  depa <- list(skewness.coef = sk, moy = mu, variance = v)
  
  return(list(spadata = gdt, init = depa))
}
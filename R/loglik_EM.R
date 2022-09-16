#' Estimate parameters using the expectation-maximization algorithm
#'
#' \code{loglik_EM} Estimate the parameters using the expectation-maximization algorithm based on the complete log-likelihood of observed data in the presence of missing data.
#'
#' @param Dstar A numeric vector containing the observed outcome of interest, should be coded as 0/1.
#' @param X A numeric vector containing the predictor of interest.
#' @param V A numeric vector of covariate in sensitivity model.
#' @param L A numeric vector of covariate in specificity model.
#' @param start Numeric vectors of starting values for theta, beta and gamma. Theta is the parameter of
#'              binary regression model. Beta and gamma are regression coefficients in sensitivity and specificity model, respectively.
#' @param tol Numerical precision.
#' @param maxit Maximum number of interations to solve the estimating equations.
#'
#' @return A list
#' \describe{
#' \item{par}{the final estimates of the coeficients organized as (theta, beta and gamma).}
#' \item{param.seq}{the sequence of estimates at each itertaion.}
#' \item{loglik.seq}{the sequence of the log likelihood at each itertaion.}
#' \item{Dstar}{A numeric vector of the input observed outcome.}
#' \item{X}{A numeric vector of the input predictor variable.}
#' \item{V}{A numeric vector of the input covariate in sensitivity model.}
#' \item{L}{A numeric vector of the input covariate in specificity model.}
#' }
#' @importFrom stats predict
#' @export
#'
#' @examples
#' library(DMSS)
#' library(rje)
#'
#' # Generate the predictor variable and covariates
#' X = rnorm(10000)
#' L = rnorm(10000)
#' V = rnorm(10000)
#'
#' # Generate the covariate-related sensitivity and specificity
#' theta = c(-2,3)
#' p_Dtrue = expit(cbind(1,X)%*% theta)
#' beta = c(0.45,0.5)
#' sensitivity = expit(cbind(1,V)%*% beta)
#' gamma = c(4.5,1)
#' specificity = expit(cbind(1,L)%*% gamma)
#'
#' # Generate the observed outcome variable
#' p_Dstar = (1-specificity)+(specificity+sensitivity -1)*p_Dtrue
#' Dstar = rbinom(10000,1,p_Dstar)
#'
#' # initial value
#' start = c(-1,1,0,0,3,0)
#' result = loglik_EM (Dstar,X,V,L,start,tol = 1e-8, maxit = 1000)



loglik_EM <- function(Dstar, X, V, L,start, tol = 1e-8, maxit = 1000){

  n <- length(Dstar)
  #the dimension of X
  K = ncol(X)
  ####### initialise parameter for EM
  theta <- start[c(1:(K+1))]
  beta  <- start[c((K+2):(K+3))]
  gamma <- start[c((K+4):(K+5))]
  pred1 <- expit(cbind(1, X) %*% theta)
  pred2 <- expit(cbind(1, V) %*% beta)
  pred3 <- expit(cbind(1, L) %*% gamma)

  ####### E-step
  calculate.p <- function(pred1, pred2,pred3)
  {
    p1 = pred2*pred1/(pred2 * pred1 + (1-pred3)*(1-pred1))*Dstar
    p2 = (1-pred2)*pred1/((1-pred2)*pred1 + pred3 * (1-pred1))*(1-Dstar)
    p1 + p2
  }
  p <- calculate.p(pred1, pred2,pred3)

  it <- 1
  converged <- F

  param.seq  <- matrix(c(theta, beta, gamma), 1)
  loglik.seq <- -10 ^ 9

  ####### M step
  while (!converged && it < maxit) {
    suppressWarnings({
      fit.beta <- stats::glm(Dstar ~ 1 + V, weights = p,
                             family = stats::binomial())
    })
    suppressWarnings({
      fit.theta <- stats::glm(p ~ 1 + X, family = stats::binomial())
    })
    suppressWarnings({
      fit.gamma <- stats::glm(1-Dstar ~ 1 + L, family = stats::binomial(),
                              weights = 1-p)
    })
    pred1 <- predict(fit.theta, type = 'response')
    pred2 <- predict(fit.beta, type = 'response')
    pred3 <- predict(fit.gamma, type = 'response')
    p     <- calculate.p(pred1,pred2,pred3)

    loglik <- sum(p*log(pred1) + (1-p)*log(1-pred1) + Dstar * p*log(pred2) + (1-Dstar)*p*log(1-pred2) + (1-Dstar)*(1-p)*log(pred3) + Dstar * (1-p) *log(1-pred3))
    loglik.seq <- c(loglik.seq, loglik)

    it <- it + 1
    if (abs(loglik.seq[it] - loglik.seq[it - 1]) < tol)
      converged <- TRUE

    par <- c(stats::coef(fit.theta),stats::coef(fit.beta), stats::coef(fit.gamma))
    param.seq <- rbind(param.seq, par)
  }
  return(list(par = par, param.seq = param.seq,loglik.seq = loglik.seq, Dstar = Dstar, X = X, V = V, L=L))
}

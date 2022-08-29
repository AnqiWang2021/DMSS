#' Calculate the covariance matrix
#'
#' \code{loglik_EM_var} Consider the "Huber Sandwich Estimator" to calculate the covariance matrix using the estimates obtained via expectation-maximization algorithm.
#' @param Dstar A numeric vector containing the observed outcome of interest, should be coded as 0/1.
#' @param X A numeric vector containing the predictor of interest.
#' @param V A numeric vector of covariate in sensitivity model.
#' @param L A numeric vector of covariate in specificity model.
#' @param theta The estimates of parameters in general binary regression model.
#' @param beta The estimates of coefficients in sensitivity model.
#' @param gamma The estimates of coefficients in specificity model.
#' @param expected Whether or not to consider the expected information matrix.
#'
#' @return 'var', the covariance matrix of the final estimates.
#' @export
#'
#' @examples
loglik_EM_var <- function(Dstar,X,V,L, theta, beta, gamma,expected)
{
  if (!is.logical(expected) || length(expected) > 1)
    stop("'expected' must be a length one logical.")

  X1 <- cbind(1, X)
  V1 <- cbind(1, V)
  L1 <- cbind(1, L)

  XTheta <- X1 %*% theta
  VBeta <- V1 %*% beta
  LGamma <- L1 %*% gamma

  n <- length(Dstar)

  expit.xt <- expit(XTheta)
  expit.vb <- expit(VBeta)
  expit.lg <- expit(LGamma)

  exp.xt <- exp(XTheta)
  exp.vb <- exp(VBeta)
  exp.lg <- exp(LGamma)

  G1 <- as.vector(expit.vb * expit.xt + 1 - expit.xt - expit.lg + expit.lg * expit.xt)

  dG1.dT <- exp.xt / (1 + exp.xt) ^ 2 * (expit.vb + expit.lg -1)
  dG1.dG <- exp.lg / (1 + exp.lg) ^ 2 * (expit.xt -1)
  dG1.dB <- exp.vb / (1 + exp.vb) ^ 2 * (expit.xt)

  dG1.dTdT <- exp.xt*(1 - exp.xt) / (1 + exp.xt)^3 * (expit.vb + expit.lg-1)
  dG1.dTdB <- exp.xt / (1 + exp.xt) ^ 2 * exp.vb / (1 + exp.vb) ^ 2
  dG1.dTdG <-  exp.xt / (1 + exp.xt) ^ 2 * exp.lg / (1 + exp.lg) ^ 2
  dG1.dGdG <- exp.lg * (1 - exp.lg) / (1 + exp.lg)^3 * (expit.xt -1)
  dG1.dGdB <- 0
  dG1.dBdB <- exp.vb * (1 - exp.vb) / (1 + exp.vb)^3 * expit.xt

  # Calculate information matrix
  if (expected) {
    tmp <- 1 / (G1 * (1 - G1))
    bread.tt <- -dG1.dT * dG1.dT * tmp
    bread.tb <- -dG1.dT * dG1.dB * tmp
    bread.tg <- -dG1.dT * dG1.dG * tmp
    bread.bb <- -dG1.dB * dG1.dB * tmp
    bread.bg <- -dG1.dB * dG1.dG * tmp
    bread.gg <- -dG1.dG * dG1.dG * tmp
  } else {
    tmp <- 1 / G1 ^ 2
    bread.tt <- Ystar * (G1 * dG1.dTdT - dG1.dT * dG1.dT) * tmp
    bread.tb <- Ystar * (G1 * dG1.dTdB - dG1.dT * dG1.dB) * tmp
    bread.tg <- Ystar * (G1 * dG1.dTdG - dG1.dT * dG1.dG) * tmp
    bread.bb <- Ystar * (G1 * dG1.dBdB - dG1.dB * dG1.dB) * tmp
    bread.bg <- Ystar * (G1 * dG1.dGdB - dG1.dG * dG1.dB) * tmp
    bread.gg <- Ystar * (G1 * dG1.dGdG - dG1.dG * dG1.dG) * tmp


    tmp <- 1 / (1 - G1) ^ 2
    bread.tt <- bread.tt - (1 - Ystar) * ((1 - G1) * dG1.dTdT +
                                            dG1.dT * dG1.dT) * tmp
    bread.tb <- bread.tb - (1 - Ystar) * ((1 - G1) * dG1.dTdB +
                                            dG1.dT * dG1.dB) * tmp
    bread.tg <- bread.tg - (1 - Ystar) * ((1 - G1) * dG1.dTdG +
                                            dG1.dT * dG1.dG) * tmp
    bread.bb <- bread.bb - (1 - Ystar) * ((1 - G1) * dG1.dBdB +
                                            dG1.dB * dG1.dB) * tmp
    bread.bg <- bread.bg - (1 - Ystar) * ((1 - G1) * dG1.dGdB +
                                            dG1.dB * dG1.dG) * tmp
    bread.gg <- bread.gg - (1 - Ystar) * ((1 - G1) * dG1.dGdG +
                                            dG1.dG * dG1.dG) * tmp

  }

  # BREAD

  I.thetatheta <- t(apply(X1, 2, function(x) x* bread.tt)) %*% X1
  I.thetabeta <- t(apply(X1, 2, function(x) x* bread.tb)) %*% V1
  I.thetagamma <- t(apply(X1, 2, function(x) x * bread.tg)) %*% L1
  I.betabeta <- t(apply(V1, 2, function(x) x * bread.bb)) %*% V1
  I.betagamma <- t(apply(V1, 2, function(x) x * bread.bg)) %*% L1
  I.gammagamma <- t(apply(L1, 2, function(x) x * bread.gg)) %*% L1

  info <- rbind(cbind(I.thetatheta, I.thetabeta, I.thetagamma),
                cbind(t(I.thetabeta), I.betabeta, I.betagamma),
                cbind(t(I.thetagamma), t(I.betagamma), I.gammagamma))

  # MEAT
  tmp <- (Ystar - G1) / (G1 * (1 - G1))
  U.theta <- apply(X1, 2, function(x) tmp * dG1.dT * x)
  U.beta <- apply(V1, 2, function(v) tmp * dG1.dB * v)
  U.gamma <- apply(L1, 2, function(l) tmp * dG1.dG * l)
  U <- cbind(U.theta, U.beta, U.gamma)


  meat <- t(U) %*% U

  var <- solve(-info) %*% meat %*% solve(-info)
  return(var = var)
}

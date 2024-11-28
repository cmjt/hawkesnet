#' Calculate the negative log-likelihood of a Hawkes process
#' @param params A list of named parameters: \code{mu}, \code{abratio}, and \code{beta}
#' @param history 
hawkes_analytic <- function(params, times, marks = rep(1, length(times))) {
    mu <- params[["mu"]]
    K <- params[["abratio"]]
    beta <- params[["beta"]]
    eps <- 1e-7
    N <- length(times)
    if ((min(K, beta) < eps) | (K > (1 - eps))) {
        return(Inf)
    }
    sumlog <- log(mu)
    const <- K * beta
    for (j in 2:N) {
        sumterm <- 0
        for (i in 1:(j - 1)) {
            sumterm <- sumterm + exp(-beta * (times[j] - times[i])) * marks[i]
        }
        lamj <- mu + const * sumterm
        if (is.na(lamj) | (lamj < 0)) {
            return(Inf)
        }
        sumlog <- sumlog + log(lamj)
    }
    T_max <- max(times) + 1e-4 
    intlam <- mu * T_max + K * sum(marks) - K * sum(exp(-beta * (T_max - times)) * marks)
    loglik <- sumlog - intlam
    return(-1.0 * loglik)
}
#' Negative log-likelihood of a marked Hawkes process with 
#' both temporal and mark decay
#' @param params A list of named parameters: \code{mu}, \code{alpha},
#'  \code{beta01} (temporal decay), and \code{beta02} (marked decay)
#' @param data A named data frame with named variables
#' \code{times} and \code{marks}.
marked_hawkes_analytic <- function(params, data) {
  mu <- params[["mu"]]
  alpha <- params[["alpha"]]
  beta01 <- params[["beta01"]]
  beta02 <- params[["beta02"]]
  if (mu <= 0 || alpha <= 0 || beta01 <= 0 || beta02 <= 0) { 
    return(-Inf)  
  }
  N <- nrow(data)         
  T_max <- max(data$times)  
  logL <- 0
  for (i in 1:N) {
    if (i == 1) {
      lambda_i <- mu
    } else {
      time_diffs <- data$times[i] - data$times[1:(i-1)]
      excitations <- alpha * exp(-beta01 * time_diffs - beta02 * data$marks[1:(i-1)])
      lambda_i <- mu + sum(excitations)
    }
    logL <- logL + log(lambda_i)
  }
  integral_term <- mu * T_max + (alpha / beta01) * sum(exp(-beta02 * data$marks) * (1 - exp(-beta01 * (T_max - data$times))))
  logL <- logL - integral_term
  return(-logL)  
}


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

#' \code{RTMB} template for a univariate Hawkes model
#' see \url{https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html}
#' for more details.
#' @param params A list of named parameters: \code{log_mu}, \code{logit_abratio}, and \code{log_beta}
#' @details The function assumes (a \code{RTMB} thing) that there exists a data frame \code}
#' times with the variable \code{times} as the vector of ordered numeric times.
hawkes_rtmb <- function(params){
    getAll(times, params)
    mu <- exp(log_mu)
    beta <- exp(log_beta)
    alpha <- exp(logit_abratio) / (1 + exp(logit_abratio)) * beta
    n <- length(times)
    last <- times[n]
    nll <- 0
    A <- advector(numeric(n))
    for(i in 2:n){
        A[i] <- sum(exp(-beta * (times[i] - times[i - 1])) * (1 + A[i - 1]))
    }
    term_3vec <- log(mu + alpha * A)
    nll <- (mu * last) + ((alpha/beta) * (n - 1 - A[n])) - sum(term_3vec)
    ADREPORT(mu)
    ADREPORT(alpha)
    ADREPORT(beta)
    return(nll)
}
#' \code{RTMB} template for a marked Hawkes model
#' see \url{https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html}
#' for more details.
#' @param params A list of named parameters: \code{mu}, \code{alpha}, \code{beta01} and \code{beta02}
#' @details The function assumes (a \code{RTMB} thing) that there exists a data frame \code}
#' times with the variable \code{times} as the vector of ordered numeric times and the variable \code{marks}
#' as the vetor of marks.
mark_effect_hawkes <- function(params){
    getAll(times, params)
    mu <- exp(log_mu)
    beta01 <- exp(log_beta01)
    beta02 <- exp(log_beta02)
    alpha <- exp(logit_abratio) / (1 + exp(logit_abratio)) * beta01
    n <- length(times)
    last <- times[n]
    nll <- 0
    A <- advector(numeric(n)) 
    for(i in 2:n){
        A[i] <- sum(exp((-beta01 * (times[i] - times[1:i - 1])) - (beta02 * marks[1:i - 1])))
    }
    term_3vec <- log(mu + alpha * A) 
    nll <- (mu * last) + (alpha/beta01) * sum(exp(-beta02 * marks) * (1 - exp(-beta01 * (last - times)))) -
        sum(term_3vec)
    ADREPORT(mu)
    ADREPORT(alpha)
    ADREPORT(beta01)
    ADREPORT(beta02)
    return(nll)
}

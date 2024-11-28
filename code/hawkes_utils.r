mark_effect_hawkes_intensity <- function(times, mu, alpha, beta01, beta02, marks, f = 500) {
  p <- seq(min(times), max(times), length.out = f)
  n <- length(p)
  idx <- findInterval(p, times)
  ms <- marks[idx]
  lambda <- numeric(n)
  for (i in 1:n) {
    if (i == 1) {
      lambda[i] <- mu
    } else {
      time_diffs <- p[i] - times[times < p[i]]
      excitations <- alpha * exp(-beta01 * time_diffs - beta02 * ms[i])
      lambda[i] <- mu + sum(excitations)
    }
  }
  return(lambda)
}

#int <- mark_effect_hawkes_intensity(times = times$times, 
                                  #  mu = pars[1,1], alpha = pars[2,1],
                                   # beta01 = pars[3,1], beta02 = pars[4,1], marks = times$marks)

#int <- mark_effect_hawkes_intensity(times = data$times, mu = 1, alpha = 1, beta01 = 4, beta02 = 8, marks = data$marks)

#plot(seq(min(data$times), max(data$times), length.out = 500), int, type = "l")
#points(data$times, rep(1, nrow(data)), col = ifelse(data$marks == 0, "black", "red"), pch = 20) ## black causes bump

#         Estimate  Std. Error
#mu     2.141776e-01 0.003322477
#alpha  1.078124e-01 0.006133331
#beta01 5.067456e-01 0.033820412
#beta02 1.097612e-07 0.000059739

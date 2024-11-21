## Univariate Hawkes

    params <- list(mu = 0.5, alpha = 1, beta = 1.5)
    simulated_times <- stelfi::sim_hawkes(mu = params[["mu"]], alpha = params[["alpha"]], 
                                          beta = params[["beta"]], n = 40, seed = 654)

### Negative log-liklihood values

    ## benchmark
    require(hawkes)
    likelihoodHawkes(lambda0 = params[["mu"]], alpha = params[["alpha"]], 
                     beta = params[["beta"]], history = simulated_times)

    ## [1] 4.475382

    require(RTMB)
    ## check against repo RTMB template
    source("rtmb.r") 
    times <- list(times = simulated_times)
    abratio <- params[["alpha"]]/params[["beta"]]
    pars <- list(log_mu = log(params[["mu"]]), logit_abratio = log(abratio/(1 - abratio)), 
                 log_beta = log(params[["beta"]]))
    hawkes_rtmb(params = pars)

    ## class='advector'
    ## [1] 4.475382

    ## check against repo analytical fun
    source("analytical.r")
    pars <- list(mu = params[["mu"]], abratio = params[["alpha"]]/params[["beta"]], 
                 beta = params[["beta"]])
    hawkes_analytic(params = pars, times = simulated_times)

    ## [1] 4.476058

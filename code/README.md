
## Univariate Hawkes

``` r
params <- list(mu = 0.5, alpha = 1, beta = 1.5)
simulated_times <- stelfi::sim_hawkes(mu = params[["mu"]], alpha = params[["alpha"]], 
                                      beta = params[["beta"]], n = 40, seed = 654)
```

### Negative log-liklihood calculations

``` r
## benchmark
require(hawkes)
likelihoodHawkes(lambda0 = params[["mu"]], alpha = params[["alpha"]], 
                 beta = params[["beta"]], history = simulated_times)
```

    ## [1] 4.475382

``` r
require(RTMB)
## check against repo RTMB template
source("rtmb.r") 
times <- list(times = simulated_times)
abratio <- params[["alpha"]]/params[["beta"]]
pars <- list(log_mu = log(params[["mu"]]), logit_abratio = log(abratio/(1 - abratio)), 
             log_beta = log(params[["beta"]]))
hawkes_rtmb(params = pars)
```

    ## class='advector'
    ## [1] 4.475382

``` r
## check against repo analytical fun
source("analytical.r")
pars <- list(mu = params[["mu"]], abratio = params[["alpha"]]/params[["beta"]], 
             beta = params[["beta"]])
hawkes_analytic(params = pars, times = simulated_times)
```

    ## [1] 4.476058

### Parameter recovery

``` r
## benchmark
require(emhawkes)
h <- new("hspec", mu = params[["mu"]], alpha = params[["alpha"]], beta = params[["beta"]])
## emhawkes requires the inter arrival times to fit the model
inter <- diff(simulated_times)
fit_em <- hfit(object = h, inter_arrival = inter)
summary(fit_em)$estimate
```

    ##         Estimate Std. error  t value     Pr(> t)
    ## mu1    0.4115873  0.1909711 2.155234 0.031143510
    ## alpha1 1.5018005  0.5239613 2.866243 0.004153753
    ## beta1  1.8576517  0.6647470 2.794525 0.005197610

``` r
## check against repo RTMB template
## function from rtmb.r
abratio <- params[["alpha"]]/params[["beta"]]
pars <- list(log_mu = log(params[["mu"]]), logit_abratio = log(abratio/(1 - abratio)), 
             log_beta = log(params[["beta"]]))
obj <- MakeADFun(hawkes_rtmb, pars, silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj), "report")
```

    ##        Estimate Std. Error
    ## mu    0.4112733  0.1709051
    ## alpha 1.4453096  0.5029997
    ## beta  1.7225179  0.6324220

``` r
## check against repo analytical fun
## function from analytical.r
pars <- c(mu = params[["mu"]], abratio = params[["alpha"]]/params[["beta"]], 
             beta = params[["beta"]])
## using optim just for now
fit <- optim(par = pars, fn = hawkes_analytic, times = simulated_times, 
             method = "Nelder-Mead", hessian = TRUE)
## parameter ests
fit$par
```

    ##        mu   abratio      beta 
    ## 0.4113093 0.8390022 1.7227333

``` r
## alpha
as.numeric(fit$par[2]*fit$par[3])
```

    ## [1] 1.445377

``` r
## standard errors
fit$hessian |> solve() |> diag() |> sqrt()
```

    ##        mu   abratio      beta 
    ## 0.1709137 0.1443989 0.6325161

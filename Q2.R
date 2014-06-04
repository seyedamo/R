
rg <- function(alpha, n) {
    r_exp <- rexp(n, 0.5)
    r_beta <- rbeta(n, 5, 2)
    r_unif <- runif(n)
    ## here I assume the order in the result vector does not matter(iid)
    ## to simplify the method
    r_mix <- r_exp[r_unif < alpha]
    r_mix <- append(r_mix, r_beta[r_unif >= alpha])
    return(r_mix)
}

samp_300 <- rg(0.8, 300)

hrefnormal <- function(X){
    X <- sort(X)
    n <- length(X)
    InterquartileSigma <- (X[round(0.75 * n)] - X[round(0.25 * n)]) / 1.34
    RobustStd <- min(sd(X), InterquartileSigma)
    h <- 1.06 * RobustStd * n^(-1/5)
    return(h)
}

hrefnormal(samp_300)


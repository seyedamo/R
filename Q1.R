
## function rg genrate random values from g-alpha
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

set.seed(10)
samp_50<- rg(0.8, 50)

set.seed(10)
samp_150 <- rg(0.8, 150)

set.seed(10)
samp_300 <- rg(0.8, 300)



##hist(samp)

## Nine methods taught in lecture
BirgeRosenholz <- function(X) {
    n <- length(X)
    BR <- rep(0,n)
    p <- 2.5
    for(i in c(1:n)) {
        r <- range(X)
        breakSeq <- seq(r[1], r[2], by = (r[2] - r[1])/i)
        h <- hist(X, breaks = breakSeq, plot = FALSE)
        nn <- h$counts
        ## some bins contain zero value so log function will complain, here is to remove zero bins
        if ( sum(nn > 0) != length(nn)) {
            nn <- nn[nn > 0]
        }
        BR[i] <- t(nn) %*% log(nn) + n*log(i) - (i - 1 + (log(i)^p))
    }
    nbins <- which.max(BR)
    maxval4BR <- max(BR)
    return(nbins)
}

Rudemo <- function(X) {
    n <- length(X)
    BR <- rep(0,n)
    for(i in c(1:n)) {
        r <- range(X)
        breakSeq <- seq(r[1], r[2], by = (r[2] - r[1])/i)
        h <- hist(X, breaks = breakSeq, plot = FALSE)
        nn <- h$counts
        BR[i] <- (i * (n + 1) / n^2) * (t(nn) %*% nn) - 2 * i
    }
    nbins <- which.max(BR)
    maxval4BR <- max(BR)
    return(nbins)
}

KullbackCV <- function(X) {
    n <- length(X)
    BR <- rep(0,n)
    for(i in c(1:n)) {
        r <- range(X)
        breakSeq <- seq(r[1], r[2], by = (r[2] - r[1])/i)
        h <- hist(X, breaks = breakSeq, plot = FALSE)
        nn <- h$counts
        if ( sum(nn > 1) != length(nn)) {
            nn <- nn[nn > 0]
        }
        BR[i] <- t(nn) %*% log(nn-1) + n * log(i)
    }
    nbins <- which.max(BR)
    maxval4BR <- max(BR)
    return(nbins)
}

AkaikeCriterion <- function(X) {
    n <- length(X)
    BR <- rep(0,n)
    for(i in c(1:n)) {
        r <- range(X)
        breakSeq <- seq(r[1], r[2], by = (r[2] - r[1])/i)
        h <- hist(X, breaks = breakSeq, plot = FALSE)
        nn <- h$counts
        if ( sum(nn > 0) != length(nn)) {
            nn <- nn[nn > 0]
        }
        BR[i] <- t(nn) %*% log(nn) + n*log(i) - (i - 1)
    }
    nbins <- which.max(BR)
    maxval4BR <- max(BR)
    return(nbins)
}

## slide is different from code provided missing factorial
RisanneHallHannan <- function(X) {
    n <- length(X)
    BR <- rep(0,n)
    for(i in c(1:n)) {
        r <- range(X)
        breakSeq <- seq(r[1], r[2], by = (r[2] - r[1])/i)
        h <- hist(X, breaks = breakSeq, plot = FALSE)
        nn <- h$counts
        prod_list <- c(1)
        for(j in c(1:i)) {
            prod_list <- append(prod_list, lfactorial(nn[j]))
        }
		b = i;
		for(j in seq(from = i + 1, to = i + n - 1, by = 1)) {
			b = b * j
		}
		BR[i] <- (i^n * prod(prod_list))/b
    }
    nbins <- which.max(BR)
    maxval4BR <- max(BR)
    return(nbins)
}

## maximizing this formula
HallHannan <- function(X) {
    n <- length(X)
    BR <- rep(0,n)
    for(i in c(1:n)) {
        r <- range(X)
        breakSeq <- seq(r[1], r[2], by = (r[2] - r[1])/i)
        h <- hist(X, breaks = breakSeq, plot = FALSE)
        nn <- h$counts
        nn <- nn[nn > 0]
        BR[i] <-t(nn - 0.5) %*% log(nn - 0.5) - (n - i/2) * log(n - i/2) + n * log(i) - (i/2) * log(n)
    }
    nbins <- which.max(BR)
    maxval4BR <- max(BR)
    return(nbins)
}

## slide is different from code
TaylorKazanawa <- function(X) {
    s2 <- var(X)
    n <- length(X)
    r <- range(X)
    nbins <- floor((r[2] - r[1])*(n/s2)^(1/3)/2.29)
    return(nbins)
}

DevroyeGyorfi <- function(X) {
    s <- sd(X)
    r <- range(X)
    n <- length(X)
    nbins <- floor((r[2] - r[1])*(n^(1/3)/(s*2.72)))
    return(nbins)
}

Scott <- function(X) {
    s <- sd(X)
    r <- range(X)
    n <- length(X)
    nbins <- floor((r[2] - r[1])*(n^(1/3)/(s * 3.49)))
    return(nbins)
}

true_pdf <- function(samp) {
	x <- seq(from=min(samp), to=max(samp), length.out=300)
	alpha = 0.8
	d <- alpha * dexp(x, 0.5) + (1 - alpha) * dbeta(x, 5, 2)
	lines(x, d, col = 'red')
}

list_of_methods <- list(BirgeRosenholz, Rudemo, KullbackCV, AkaikeCriterion, RisanneHallHannan, HallHannan, TaylorKazanawa, DevroyeGyorfi, Scott)
list_of_methods_nm <- list("BirgeRosenholz","Rudemo","KullbackCV","AkaikeCriterion","RisanneHallHannan","HallHannan","TaylorKazanawa","DevroyeGyorfi","Scott")

## list is a list of functions, all have the same interface(takes X as argument)
## X is the random sample data
comp_methods <- function(l,ln, X) {
    n <- length(l)
    result_table <- matrix(rep(0, n*2), nrow = n)
    ## open 9 sub diagrams(n x 1)
    par(mfrow = c(3, 3))
    space <- range(X)
    X_min <- space[1]
    X_max <- space[2]
    ## for(method in l){
    ##     nbins <- method(X)
    ##     result_table[i,1] <- nbins
    ##     i = i + 1
    ##     hist(X, breaks = seq(X_min, X_max, by = (X_max - X_min)/(nbins)))
    ## }

    ## result_list <- lapply(l, FUN = function(method) {
    ##     nbins <- method(X);
    ##     hist(X, breaks = seq(X_min, X_max, by = (X_max - X_min)/nbins));
    ##     return(nbins)
    ## })

    result_list <- Map(function(method, title) {
        nbins <- method(X);
        hist(X, breaks = seq(X_min, X_max, by = (X_max - X_min)/nbins), 
			main = title, xlab = paste("number of bins",nbins), freq = FALSE, ylim =range(c(0,1)))
		true_pdf(X)
        return(nbins)
    }, l, ln)

    return(result_list)
}


result_list <- comp_methods(list_of_methods, list_of_methods_nm, samp_50)

dev.new()
result_list <- comp_methods(list_of_methods, list_of_methods_nm, samp_150)

dev.new()
result_list <- comp_methods(list_of_methods, list_of_methods_nm, samp_300)

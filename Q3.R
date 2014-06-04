library(ggplot2)

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

set.seed(100)
x <- rg(0.8, 300)
y <- 5 * x^(1/2) + rnorm(300, 0, 1)

## Q1

## function scatter.smooth
scatter.smooth(x,y)

## native loess in R
loess <- loess(y ~ x)
summary(loess)
plot(x, y)


## ggplot's loess smooth
p <- qplot(x = x, y = y)
p + geom_smooth(method = "loess", size = 1.5)

## put x,y in a data frame
samp <- data.frame(x, y)

## My own implementation of loess
##
my_loess <- function(degree, alpha, samp) {
    n = length(samp$x)
    k = floor(alpha * n)

    linear <- function(x0) {
        sort <- sort(abs(x0 - samp$x), index.return = TRUE)
        ind <- sort$ix[1:k]
        range <- range(sort$x[1:k])
        delx0 <- range[2] - range[1]
        window <- samp[ind,]
        u <- abs(x0 - window$x) / delx0
        w <- (1 - u^3)^3

        if(degree == 1) {
            lm.fit <- lm(y ~ x, data = window, weights = w)
            coef <- lm.fit$coefficients
            y0_hat = coef[1] + coef[2] * x0
        } else {
            lm.fit <- lm(y ~ x + I(x^2), data = window, weights = w)
            coef <- coef(lm.fit)
            y0_hat = coef[1] + coef[2] * x0 + coef[3] * x0^2
        }
        return(y0_hat)
    }

    y_hat <- sapply(samp$x, FUN = linear)

    x_sort <- sort(samp$x, index.return = TRUE)
    qplot(samp$x, samp$y) + geom_line(aes(x = x_sort$x, y = y_hat[x_sort$ix], colour = 'red', size = 1.5, alpha = 0.5))
}

my_loess(2,1/3,samp)

## Q3
plot(x,y)
lines(ksmooth(x, y, "normal", bandwidth = 2), col = 2)
lines(ksmooth(x, y, "normal", bandwidth = 3), col = 3)

?geom_smooth

my_NW <- function(h, samp) {
    x <- samp$x
    y <- samp$y
    y_hat <- sapply(x, FUN = function(x0) {
        w <- (exp(-0.5 * ((x - x0)/h)^2)/sqrt(2*pi)) / h
        y0 <- (t(w) %*% y) / sum(w)
        return(y0)
    })

    ## x_sort <- sort(samp$x, index.return = TRUE)
    ## plot(samp$x, samp$y)
    ## lines(x_sort$x, y_hat[x_sort$ix], col = 'red')

    df <- data.frame(samp$x, y_hat)
    names(df) <- c("x","y_hat")
    return(df)
}

df <- my_NW(1, samp)
df1 <- my_NW(0.5, samp)

qplot(samp$x, samp$y) + geom_line(data = df, aes(x = x, y = y_hat), color = "red", size = 1.5, alpha = 0.5) + geom_line(data = df1, aes(x = x, y = y_hat), color = "blue", size = 1.5, alpha = 0.5)

##K <- function(x0) {
##    return(exp(-0.5 * ((x - x0)/h)^2)/sqrt(2*pi))
##}

## Q4

my_local_linear <- function(h, samp) {
    x <- samp$x
    y <- samp$y
    n <- length(x)
    y_hat <- sapply(x, FUN = function(x0) {
        w <- (exp(-0.5 * ((x - x0)/h)^2)/sqrt(2*pi)) / h
        xc <- x - x0
        s0 <- sum(w) / n
        s1 <- sum(xc * w) / n
        s2 <- sum(xc^2 * w) / n
        y0 <- sum(((s2 - s1*xc)*w*y)/(s2 * s0 - s1^2)) / n
        return(y0)
    })

    ## x_sort <- sort(samp$x, index.return = TRUE)
    ## plot(samp$x, samp$y)
    ## lines(x_sort$x, y_hat[x_sort$ix], col = 'red')

    df <- data.frame(samp$x, y_hat)
    names(df) <- c("x","y_hat")
    return(df)
}

df <- my_local_linear(1, samp)
qplot(samp$x, samp$y) + geom_line(data = df, aes(x = x, y = y_hat), color = "red", size = 1.5, alpha = 0.5)

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
scatter.smooth(x,y)

loess <- loess(y ~ x)
summary(loess)
plot(x, y)

p <- qplot(x = x, y = y)
p + geom_smooth(method = "loess", size = 1.5)


samp <- data.frame(x, y)
x0 <- x[1]

sort <- sort(abs(x0 - samp$x), index.return = TRUE)
L <- sort$ix[1:5]
samp[L,]

## Q3
plot(x,y)
lines(ksmooth(x, y, "normal", bandwidth = 2), col = 2)
lines(ksmooth(x, y, "normal", bandwidth = 3), col = 3)

?geom_smooth


samp <- data.frame(x,y)

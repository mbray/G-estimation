
expit <- function(x) {exp(x)/(exp(x) + 1)}
logit <- function(x) {log(x/(1 - x))}

# For plotting residuals...
myplot <- function(x, ...) {
    require(lattice)

    mypanel <- function(x, y, ...) {
        if (is.factor(x)) {
            panel.bwplot(x, y, ..., horizontal=FALSE)
            xx <- 1:nlevels(x)
            yy <- tapply(y, x, median)
            panel.segments(
                xx[1:(nlevels(x)-1)],
                yy[1:(nlevels(x)-1)],
                xx[2:(nlevels(x))],
                yy[2:(nlevels(x))]
            )
        } else {
            panel.xyplot(x, y, col="gray60", ...)
            panel.loess(x, y, col="blue", ...)
        }
        panel.abline(h=0, col="red")
    }

    xyplot(x, panel=mypanel, ...)
}

### Simulation ###

simulate.dataset <- function(n=500) {
    l0 <- rnorm(n, 450, 150)
    l1 <- rnorm(n, 1.25*l0, 60)

    p1 <- expit(2 - 0.006*l0)
    p2 <- expit(0.8 - 0.004*l1)
    a1 <- rbinom(n, 1, p1)
    a2 <- pmax(a1, rbinom(n, 1, p2))

    psi10 <- 250
    psi11 <- -1
    psi20 <- 720
    psi21 <- -2
    gamma1 <- a1*(psi10 + psi11*l0)
    gamma2 <- a2*(psi20 + psi21*l1)
    mu1 <- pmax(psi10 + psi11*l0, 0) - gamma1
    mu2 <- pmax(psi20 + psi21*l1, 0) - gamma2

    y <- rnorm(n, 400 + 1.6*l0, 60) - mu1 - mu2

    df1 <- data.frame(l0, a1, l1, a2, y)

    return (df1)
}

### G-estimation ###

gest2 <- function(y.mod, data, blip.mod, treat.mod, ecfact.mod, keep=TRUE) {

    y <- model.response(model.frame(y.mod, data))
    n <- length(y)

    k <- length(blip.mod)
    if (length(treat.mod) != k || length(ecfact.mod) !=k) {
        stop("Lengths of blip.mod, treat.mod and ecfact.mod must all equal number of time intervals.")
    }

    kappa <- rep(0, n)

    psi.hat <- list()

    if (keep) {
        blip <- list()
        regret.sum <- list()
        residuals <- list()
        fitted <- list()
    }

    for (j in k:1) {
        z <- model.matrix(blip.mod[[j]], data)
        a <- model.response(model.frame(treat.mod[[j]], data))
        s <- a * z
        E.s <- fitted(glm(treat.mod[[j]], data=data, family=binomial)) * z
        w <- model.matrix(ecfact.mod[[j]], data)
        w.qr <- qr(w)
        temp1 <- t(s - E.s) %*% qr.resid(w.qr, y + kappa)
        temp2 <- t(s - E.s) %*% qr.resid(w.qr, s)
        psi.hat[[j]] <- drop(qr.solve(temp2, temp1))
        temp3 <- z %*% psi.hat[[j]]
        gamma.hat <- a * temp3

        if (keep) {
            blip[[j]] <- gamma.hat
            regret.sum[[j]] <- kappa
            residuals[[j]] <- qr.resid(w.qr, y - gamma.hat + kappa)
            fitted[[j]] <- qr.fitted(w.qr, y - gamma.hat + kappa) - kappa + gamma.hat
        }

        kappa <- kappa + pmax(temp3, 0) - gamma.hat
    }

    obj <- list()
    class(obj) <- "gest"
    obj$call <- match.call(expand.dots=TRUE)
    obj$samplesize <- n
    obj$intervals <- k
    obj$psi.hat <- psi.hat
    if (keep) {
        obj$blip <- blip
        obj$regret.sum <- regret.sum
        obj$residuals <- residuals
        obj$fitted <- fitted
    }
    return (obj)
}


# Together, gest6 and gest7 fit the true EC model.
# gest6 takes an initial guess at psi and produces a new guess.
# gest7 iterates this process until convergence.
gest6 <- function(y.mod, data, blip.mod, treat.mod, psi.initial=NULL, ...) {
    if(is.null(psi.initial)) {
        # Use a simple model for initial estimate
        ecfact.mod <- list(~ l0, ~ l0 + a1 + l1)
        gest2(y.mod, data, blip.mod, treat.mod, ecfact.mod, ...)
    } else {
        psi10 <- psi.initial[[1]][1]
        psi11 <- psi.initial[[1]][2]
        psi20 <- psi.initial[[2]][1]
        psi21 <- psi.initial[[2]][2]
        ecfact.mod <-
            #list(~ l0*I(psi10 + psi11*l0 > 0),
            #     ~ l0*I((a1 == 0 & psi10 + psi11*l0 > 0) - (a1 == 1 & psi10 + psi11*l0 < 0)) + l1*I(psi20 + psi21*l1 < 0))
            list(~ l0*I(psi10 + psi11*l0 > 0),
                 ~ l0*I(a1 == 0 & psi10 + psi11*l0 > 0) + l0*I(a1 == 1 & psi10 + psi11*l0 < 0) + l1*I(psi20 + psi21*l1 < 0))
        gest2(y.mod, data, blip.mod, treat.mod, ecfact.mod, ...)
    }
}

gest7 <- function(y.mod, data, blip.mod, treat.mod, psi.initial=NULL, verbose=FALSE, ...) {
    psi <- gest6(y.mod, data, blip.mod, treat.mod, psi.initial, ...)$psi.hat
    if (verbose) cat("Initial: ", unlist(psi), "\n")
    psi.new <- gest6(y.mod, data, blip.mod, treat.mod, psi, ...)$psi.hat
    iter <- 1
    while(any(abs((unlist(psi.new) - unlist(psi))/unlist(psi)) > 0.000001)) {
        if (verbose) cat("Iteration: ", unlist(psi.new), "\n")
        psi <- psi.new
        psi.new <- gest6(y.mod, data, blip.mod, treat.mod, psi, ...)$psi.hat

        # Note: this doesn't necessary converge.  It can oscillate between different values.  Usually, no more than 3 iterations are needed.
        if (iter >= 50) break
        iter <- iter + 1
    }
    gest6(y.mod, data, blip.mod, treat.mod, psi, ...)
}


# Works with gest2 and gest7
gest.boot <- function(data, gest.fn, B=200, verbose=FALSE, store=FALSE, ...) {
    psi.boot <- list()
    for (i in 1:B) {
        data.boot <- data[sample(1:nrow(data), replace=TRUE),]
        psi.boot[[i]] <- gest.fn(data=data.boot, ...)$psi.hat
        if (verbose) cat(i, "..")
    }
    if (verbose) cat("\n")

    psi.boot <- do.call(function(...) mapply(rbind, ..., SIMPLIFY=FALSE), psi.boot)
    psi.covmat <- lapply(psi.boot, var)

    obj <- gest.fn(data=data, ...)
    obj$psi.covmat <- psi.covmat
    if (store) {
        obj$psi.boot <- psi.boot
    }
    return (obj)
}


print.gest <- function(obj, ...) {
    cat("Linear g-estimation\n\n")
    cat("Sample size =", obj$samplesize, "\t\t Time intervals =", obj$intervals, "\n\n")
    cat("Estimates:\n")
    lapply(obj$psi.hat, print, ...)
    if (!is.null(obj$psi.covmat)) {
        cat("\nSE:\n")
        lapply(obj$psi.covmat, function(x) print(sqrt(diag(x)), ...), ...)
    }
    invisible(obj)
}

residuals.gest <- function(obj, j=NULL) {
    if (is.null(obj$residuals)) {
        warning("Residuals were not stored.  Use keep=TRUE to store residuals")
        return (NULL)
    }
    if (is.null(j)) {
        obj$residuals
    } else {
        obj$residuals[[j]]
    }
}

psimat.jk.gest <- function(obj, data) {
    cl <- obj$call
    cl[["data"]] <- expression(data[-i,])[[1]]
    psi <- list()
    for (i in 1:nrow(data)) {
        psi[[i]] <- eval(cl)$psi.hat
    }

    do.call(function(...) mapply(rbind, ..., SIMPLIFY=FALSE), psi)
}

qform <- function(vec, mat) {
    temp <- mat %*% vec
    if (is.matrix(vec)) {
        apply(vec * temp, 2, sum)
    } else {
        t(vec) %*% temp
    }
}

cook.d.gest <- function(psi.jackknife, psi.hat, psi.covmat) {
    cook <- list()
    for (j in 1:length(psi.hat)) {
        cook[[j]] <- qform(t(psi.jackknife[[j]]) - psi.hat[[j]], solve(psi.covmat[[j]]))
    }
    return (cook)
}

inflmeas.gest <- function(obj, data) {
    psi.jackknife <- psimat.jk.gest(obj, data)
    cook.d <- cook.d.gest(psi.jackknife, obj$psi.hat, obj$psi.covmat)
    list(psi.jackknife=psi.jackknife, cook.d=cook.d)
}

### Testing ###

df1 <- simulate.dataset(n=500)
y.mod <- y ~ 1
blip.mod <- list(~ l0, ~ l1)
treat.mod <- list(a1 ~ l0, a2 ~ l0 + a1 + l1)
ecfact.mod <- list(~ l0, ~ l0 + a1 + l1)
obj <- gest2(y.mod, df1, blip.mod, treat.mod, ecfact.mod)
print(obj)

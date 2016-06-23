est.reg <-
function (x, col_indep) 
{
    cols <- ncol(x)
    rows <- nrow(x)
    estimated_matrix <- matrix(ncol = cols, nrow = rows)
    if (col_indep == 1) {
        deps <- 2:cols
    }
    else {
        deps <- c(1:(col_indep - 1), (col_indep + 1):cols)
    }
    ndeps <- length(deps)
    indep <- x[, col_indep]
    rcoefs <- numeric()
    for (i in 1:ndeps) {
        vari <- deps[i]
        lm_fit <- lm(log10(x[, vari]) ~ log10(indep))
        lm_sum <- summary.lm(lm_fit)
        rcoefs[i] <- sqrt(lm_sum$r.squared)
    }
    ranks <- rank(rcoefs)
    newindep <- indep
    for (m in 1:length(deps)) {
        a <- m - 1
        whichvari <- ifelse(ranks == (length(rcoefs) - a), 1, 
            0)
        strongest <- sum(deps * whichvari)
        indeps_lm <- lm(log10(newindep) ~ log10(x[, strongest]))
        indeps_coef <- indeps_lm$coefficients
        logestimate_indep <- indeps_coef[1] + indeps_coef[2] * 
            log10(x[, strongest])
        estimate_indep <- 10^logestimate_indep
        missings <- ifelse(is.na(newindep), 1, 0)
        nonmissings <- ifelse(is.na(newindep), 0, 1)
        fillnonmissing <- ifelse(is.na(newindep), 0, newindep)
        fillmissing <- ifelse(is.na(newindep), estimate_indep, 
            0)
        newindep <- fillmissing + fillnonmissing
    }
    estimated_matrix[, col_indep] <- newindep
    for (k in 1:ndeps) {
        vari <- deps[k]
        variable <- x[, vari]
        lm_fit <- lm(log10(variable) ~ log10(newindep))
        lm_coef <- lm_fit$coefficients
        logest_vari <- lm_coef[1] + lm_coef[2] * log10(newindep)
        est_vari <- 10^logest_vari
        missings <- ifelse(is.na(variable), 1, 0)
        nonmissings <- ifelse(is.na(variable), 0, 1)
        fillnonmissing <- ifelse(is.na(variable), 0, variable)
        fillmissing <- ifelse(is.na(variable), est_vari, 0)
        estimated_matrix[, vari] <- fillmissing + fillnonmissing
    }
    return(estimated_matrix)
}

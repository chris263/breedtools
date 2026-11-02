# Reference adjusted means: average prediction of each genotype over all
# observed plots from a full lm() fit, with standard errors from vcov().
reference_adjusted_means <- function(fit, data, genotype = "genotype") {
  tt <- stats::delete.response(stats::terms(fit))
  ok <- !is.na(stats::coef(fit))
  b <- stats::coef(fit)[ok]
  V <- stats::vcov(fit, complete = FALSE)
  levels_g <- levels(data[[genotype]])
  L <- t(vapply(levels_g, function(g) {
    nd <- data
    nd[[genotype]] <- factor(g, levels = levels_g)
    mf <- stats::model.frame(tt, nd, xlev = fit$xlevels)
    colMeans(stats::model.matrix(tt, mf))[ok]
  }, numeric(sum(ok))))
  data.frame(
    genotype = levels_g,
    mean = drop(L %*% b),
    se = sqrt(rowSums((L %*% V) * L)),
    stringsAsFactors = FALSE
  )
}

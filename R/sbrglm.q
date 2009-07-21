glm <- function (formula, family = gaussian, data, weights, subset, 
  na.action, start = NULL, etastart, mustart, offset, control = glm.control(...), 
  model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...,
  separation = c("find", "test") ) 
{
  call <- match.call()
  separation <- match.arg(separation)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
    "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ", method))
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  if (!is.null(offset)) {
    if (length(offset) == 1) 
      offset <- rep(offset, NROW(Y))
    else if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
        length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")

  if(casefold(family$family) == "binomial" && length(unique(Y)) == 2) {
    if(separation == "test") {
      separation <- separator(X, Y, purpose = "test")$separation
      #separation <- separationTest(X, Y)
      if(separation)
        stop("Separation exists among the sample points.\n\tThis model cannot be fit by maximum likelihood.")
    }
    if(separation == "find") {
      separation <- separator(X, Y, purpose = "find")$beta
      #separation <- separationDirection(X, Y)
      separating.terms <- dimnames(X)[[2]][abs(separation) > 1e-09]
      if(length(separating.terms))
        stop(paste("The following terms are causing separation among the sample points:",
          paste(separating.terms, collapse = ", ")))
    }
  }

  fit <- glm.fit(x = X, y = Y, weights = weights, start = start, 
    etastart = etastart, mustart = mustart, offset = offset, 
    family = family, control = control, intercept = attr(mt, 
    "intercept") > 0)
  if (length(offset) && attr(mt, "intercept") > 0) {
    fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
      y = Y, weights = weights, offset = offset, family = family, 
      control = control, intercept = TRUE)$deviance
  }
  if (model) 
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x) 
    fit$x <- X
  if (!y) 
    fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula, terms = mt, 
    data = data, offset = offset, control = control, method = method, 
    contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))
  class(fit) <- c("glm", "lm")
  fit
}


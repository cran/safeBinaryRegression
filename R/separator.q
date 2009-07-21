## test for and/or find the direction of separation
## x a design matrix and y a 0-1 binary response vector

separator <- function(x, y, method = c("primal", "dual"), purpose = c("test", "find"),
                      tol = 1e-3)
{
  n <- dim(x)[1]
  p <- dim(x)[2]

  dimnames(x) <- NULL

  y.bar <- -sign(y - 0.5)
  x.bar <- y.bar * x

  ans <- list()

  if(method == "primal" && purpose == "test") {
    lp <- make.lp(n, p)
    for(j in 1:p)
      status <- set.column(lp, j, x.bar[, j])
    status <- set.rhs(lp, rep(0.0, n))
    status <- set.constr.type(lp, rep(1, n))
    status <- set.objfn(lp, -colSums(x.bar))
    status <- set.bounds(lp, lower = rep(-Inf, p), upper = rep(Inf, p))
    control <- lp.control(lp, pivoting = "firstindex", sense = "max",
                          simplextype = c("primal", "primal"))
    status <- solve(lp)

    if(status == 0)
      ans$separation <- FALSE
    else if(status == 3)
      ans$separation <- TRUE
    else
      stop("unexpected result from lpSolveAPI for primal test")
  }

  if(method == "primal" && purpose == "find") {
    lp <- make.lp(n, p)
    for(j in 1:p)
      status <- set.column(lp, j, x.bar[, j])
    status <- set.rhs(lp, rep(0.0, n))
    status <- set.constr.type(lp, rep(1, n))
    status <- set.objfn(lp, -colSums(x.bar))
    status <- set.bounds(lp, lower = rep(-1, p), upper = rep(1, p))
    control <- lp.control(lp, pivoting = "firstindex", sense = "max",
                          simplextype = c("primal", "primal"))
    status <- solve(lp)

    if(status != 0)
      stop("unexpected result from lpSolveAPI for primal test")

    beta <- get.variables(lp)

    if(sum(abs(beta)) > tol)
      ans$separation <- TRUE
    else
      ans$separation <- FALSE

    ans$beta <- beta
  }

  if(method == "dual" && purpose == "test") {
    lp <- make.lp(p, n)
    for(j in 1:n)
      status <- set.column(lp, j, x.bar[j, ])
    status <- set.rhs(lp, -colSums(x.bar))
    status <- set.constr.type(lp, rep(3, p))
    status <- set.objfn(lp, rep(0.0, n))
    status <- set.bounds(lp, lower = rep(0.0, n), upper = rep(Inf, n))
    control <- lp.control(lp, pivoting = "firstindex", sense = "min",
                          simplextype = c("primal", "primal"))
    status <- solve(lp)

    if(status == 0)
      ans$separation <- FALSE
    else if(status == 2)
      ans$separation <- TRUE
    else
      stop("unexpected result from lpSolveAPI for dual test")
  }

  if(method == "dual" && purpose == "find") {
    lp <- make.lp(p, n + 2*p)
    for(j in 1:n)
      status <- set.column(lp, j, x.bar[j, ])
    for(j in 1:p)
      status <- set.column(lp, n+j, -1.0, j)
    for(j in 1:n)
      status <- set.column(lp, n+p+j, 1.0, j)
    b <- -colSums(x.bar)
    status <- set.rhs(lp, b)
    status <- set.constr.type(lp, rep(3, p))
    status <- set.objfn(lp, rep(c(0.0, 1.0), c(n, 2*p)))
    status <- set.bounds(lp, lower = rep(0.0, n + 2*p), upper = rep(Inf, n + 2*p))
    control <- lp.control(lp, pivoting = "firstindex", sense = "min",
                          simplextype = c("primal", "primal"))
    basis <- 1:p
    basis[b >= 0.0] <- basis[b >= 0.0] + p
    status <- set.basis(lp, -(n + p + basis))
    status <- solve(lp)

    beta <- get.dual.solution(lp)[2:(p+1)]

    if(sum(abs(beta)) > tol)
      ans$separation <- TRUE
    else
      ans$separation <- FALSE

    ans$beta <- beta
  }

  ans
}



### Telescoping Sampling ####
# Arguments:
#       data = data in matrix format
#       K_init = initial value of K
#       alpha_init = initial value of alpha
#       s_alpha = std of the log normal alpha proposal
#       niter = number of posterior samples to draw
#       nburnin = number of samples to be used in burn-in before drawing as a posterior
#       pr_k = list of distribution name and corresponding parameters
#              so far accepted names are "pois" (poisson), "bnb" (bnb) and "geom" (geometric)
#              followed by appropriate number of parameters e.g. list("pois",1), list(bnb,3,4,3),...
#       pr_alpha = list of shape and rate parameters for gamma prior on alpha
#                  e.g. for vague gamma prior we have list(0.001,0.001)
#       K_max = max number that K can take
#       printprogress = TRUE/FALSE, if TRUE, estimated time until the finish will be printed (if run on console
#                       or rmd file, the printed value will stay in one line (won't create newline every time),
#                       but when knitting, it will print line by line on html output.
# Return values:
#       eta (dim 1*niter): log eta from Dirchlet posterior
#       Ks  (dim 1*niter): posterior K
#       Kps (dim 1*niter): posterior K_plus
#       alphas (dim 1*niter): posterior alpha
#       comp_means (dim K^(n)*r*niter): component means with varying 1st dimention depending on posterior K at iteration n
#                                       (so out of K^(n)*r*niter samples K_plus^(n)*r*niter from posterior and rest from prior)
#       comp_covs (dim K^(n)*r*r*niter): components covariances with varying 1st dim.. (same as above basically but cov this time)
#       bk (dim 1*niter): posterior bk
#       Bk (dim 1*niter): posterior Bk
#library(evd)
#library(Rcpp)
# modify rgig function in GIGrvg so that it's compatible with the new R version (.Call needs PACKAGE argument now)
# body(rgig)[[2]] = substitute(.Call("rgig", n, lambda, chi, psi,PACKAGE = "GIGrvg"))

# read cpp file used for multivariate normal pdf evaluation
#sourceCpp("hoge.cpp")

bmbclust <- function(data, K_init, alpha_init, s_alpha, niter, nburnin, pr_k, pr_alpha, K_max, printprogress = TRUE,v1,v2) {
  # set relevant parameters
  r <- ncol(data)
  K_0 <- K_init
  alpha <- alpha_init

  # initialize mean and cov
  meanmat <- array(runif(r * K_0), dim = c(K_0, r, 1))
  covmat <- array(runif(r * r * K_0), dim = c(K_0, r, r))
  covmat <- array(t(apply(covmat, 1, function(x) crossprod(x))), dim = c(K_0, r, r))

  lambdas <- rep(1, r)

  # prior for Wishart
  c0 <- 2.5 + (r - 1)/2
  g0 <- 0.5 + (r - 1)/2
  R <- apply(apply(data, 2, range), 2, function(x) x[2] - x[1])
  G0 <- 100 * g0/c0 * diag(1/R^2, nrow = r)
  tmp <- matrix(runif(r^2) * 2 - 1, ncol = r)
  C0_j <- crossprod(tmp)
  B0 <- diag(R^2, nrow = r)
  inv_B0 <- diag(1/(R^2), nrow = r)
  m0 <- matrix(apply(data, 2, median), ncol = 1)
  b0 <- m0

  etas <- lrdir(rep(1/K_0, K_0))

  K <- K_0


  res_lst <- list(eta = list(), Ks = numeric(niter), Kps = numeric(niter),
                  alphas = numeric(niter),comp_means = list(),
                  comp_covs = list(), bs = list(), Bs = list()
                  )
  # pb = txtProgressBar(min = 0, max = niter+nburnin, initial = 0,style = 3)
  start <- Sys.time()
  for (n in 1:(niter + nburnin)) {
    # setTxtProgressBar(pb, n)

    # for each comp_(means,covs), compute the log-likelihood of data + log eta, then add gumbel
    # distributed noise for each element and pick the max idx which equals the categorical sampling
    # from unnormalized log probabilities (see gumbel max trick for details)
    eta_prod_dens <- simplify2array(
      lapply(1:K,function(k)
        eval_pdf(data,list(dMvn,matrix(meanmat[k,,],ncol = r),matrix(covmat[k,,],ncol =r)),TRUE)+
          etas[k]
        )
      )
    gumbel_noise <- array(evd::rgumbel(K*dim(data)[1]),dim = c(dim(data)[1],1,K))
    S <- apply(eta_prod_dens+gumbel_noise,1,function(x) which.max(x))

    # number of data clusters
    K_plus <- length(unique(S))

    # reorder elements
    mapping <- 1:K_plus
    names(mapping) <- unique(S)

    not_filled <- seq(1, K)[!(1:K %in% unique(S))]
    if (length(not_filled) != 0) {
      mapping <- c(mapping, (K_plus + 1):K)
      names(mapping)[(K_plus + 1):K] <- not_filled
    }

    meanmat <- array(meanmat[as.numeric(names(mapping)), , ], dim = c(K, r, 1))
    covmat <- array(covmat[as.numeric(names(mapping)), , ], dim = c(K, r, r))
    S <- sapply(S, function(x) which(x == as.numeric(names(mapping))))
    N <- c(tabulate(S), rep(0, K - K_plus))

    # from here on it is specific to the Gaussian mixture
    ck <- c0 + N/2
    Ck <- lapply(1:K_plus, function(k) C0_j + 0.5 * crossprod(sweep(data[S == k, , drop = FALSE], 2, meanmat[k, , ], FUN = "-")))

    sigs <- lapply(1:K_plus, function(k) bayesm::rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[[k]]))))
    inv_sig <- lapply(1:K_plus, function(k) sigs[[k]]$W)
    sig <- lapply(1:K_plus, function(k) sigs[[k]]$IW)

    updated_vals <- lapply(1:K_plus, function(k) {
      Bk <- chol2inv(chol(inv_B0 + N[k] * inv_sig[[k]]))
      mean_yk <- apply(data, 2, function(col_i) mean(col_i[S == k]))
      bk <- Bk %*% ((inv_B0) %*% b0 + inv_sig[[k]] %*% mean_yk * N[k])
      mu_k <- t(chol(Bk)) %*% rnorm(r) + bk
      return(list(mu = mu_k, b = bk, B = Bk))
    })

    mu <- lapply(1:K_plus, function(k) updated_vals[[k]]$mu)
    bk <- lapply(1:K_plus, function(k) updated_vals[[k]]$b)
    Bk <- lapply(1:K_plus, function(k) updated_vals[[k]]$B)
    C0_j <- bayesm::rwishart(2 * (g0 + K_plus * c0), 0.5 * chol2inv(chol(G0 + Reduce("+", inv_sig))))$W

    pk <- v1 - K_plus/2
    aj <- 2 * v2
    bj <- matrix(apply(matrix(sapply(1:K_plus, function(k) (mu[[k]] - b0)^2), nrow = r), 1, sum)/(as.vector(R)^2), ncol = 1)
    #if(any(is.na(bj))){browser()}
    lambdas <- sapply(bj, function(x) rgig(1, pk, x, aj))
    B0 <- diag(lambdas * as.vector(R)^2/K_plus, nrow = r)
    inv_B0 <- diag(1/diag(B0), nrow = r)
    b0 <- crossprod(chol(B0),rnorm(r)) + Reduce("+", mu)/K_plus
    # when pk = 0 and one element of bj sufficiently small like 1e-15 ish rgig(1,pk,x,aj) will return Inf
    # causing b0 to be Inf as well
    if(any(is.infinite(b0))){browser()}
    if (pr_k[[1]] == "pois") {
      prior_k <- structure(list(dpois, pr_k[[2]]), class = "one_param")
    } else if (pr_k[[1]] == "geom") {
      prior_k <- structure(list(dgeom, pr_k[[2]]), class = "one_param")
    } else if (pr_k[[1]] == "bnb") {
      prior_k <- structure(list(dbnb, pr_k[[2]], pr_k[[3]], pr_k[[4]]), class = "three_param")
    }

    # compute unnormalized log prob to do gumbel max trick again to pick K
    unnorm_logprob <- sapply(K_plus:K_max, function(k) eval_pdf(k - 1, prior_k, TRUE) + lfactorial(k) - lfactorial(k - K_plus)
                             + sum(
                               log(k/(N * k + alpha)) + lgamma(1 + N + alpha/k) - log(k/alpha) - lgamma(1 + alpha/k)
                               )
                             )

    K <- which.max(evd::rgumbel(length(unnorm_logprob)) + unnorm_logprob) + K_plus - 1


    # sample alpha with metropolis random walk (on log space so with Jacobian corrections)
    alpha_prop <- exp(rnorm(1, log(alpha), s_alpha))
    tmp <- sapply(c(alpha_prop, alpha),
                  function(a) df(a, df1 = pr_alpha[[1]], df2 = pr_alpha[[2]], log = TRUE) + lgamma(a) - lgamma(sum(N) + a) +
                    sum(
                      log(K/(N * K + a)) + lgamma(1 + N + a/K) - log(K/a) - lgamma(1 + a/K)
                      )
                  )
    logr <- tmp[1] - tmp[2] + log(alpha_prop) - log(alpha)
    accept <- log(runif(1)) < logr
    if (accept) {
      alpha <- alpha_prop
    }

    # sample component means either from posterior (<=K_plus) or prior (>K_plus)
    meanmat <- simplify2array(lapply(1:K, function(k) if (k <= K_plus) {
      mu[[k]]
    } else {
      crossprod(chol(B0),rnorm(r)) + b0
    }))
    if (is.null(dim(meanmat))) {
      meanmat <- array(meanmat, dim = c(K, 1, 1))
    } else {
      meanmat <- aperm(meanmat, c(3, 1, 2))
    }

    # same thing for component covs
    covmat <- simplify2array(lapply(1:K, function(k) if (k <= K_plus) {
      sig[[k]]
    } else {
      bayesm::rwishart(2 * c0, 0.5 * chol2inv(chol(C0_j)))$IW
    }))
    if (is.null(dim(covmat))) {
      covmat <- array(covmat, dim = c(K, 1, 1))
    } else {
      covmat <- aperm(covmat, c(3, 1, 2))
    }

    N <- c(N[1:K_plus], rep(0, K - K_plus))
    etas <- lrdir(alpha/K + N)

    if (n > nburnin) {
      idx <- n - nburnin
      res_lst$eta[[idx]] <- etas
      res_lst$Ks[idx] <- K
      res_lst$Kps[idx] <- K_plus
      res_lst$alphas[idx] <- alpha
      res_lst$comp_means[[idx]] <- meanmat
      res_lst$comp_covs[[idx]] <- covmat
      res_lst$bs[[idx]] <- bk
      res_lst$Bs[[idx]] <- Bk
    }

    if (n%%100 == 0 && printprogress) {
      end <- Sys.time()
      cum_time <- as.numeric(difftime(end, start, units = "mins"))
      estimated_time <- cum_time/n * (niter + nburnin - n)
      cat("\r", "remaining time (in mins):\ ", as.character(format(round(estimated_time,digits = 3),nsmall = 3)))
    }

  }
  # close(pb)
  return(res_lst)
}
# Define an object for likelihood evaluation with any given distributions (so that we can expand for mixture of mixtures case etc.)
eval_pdf <- function(data, densfun, ...) {
  UseMethod("eval_pdf", densfun)
}

# Set default method as a 2 parameter family (like multivariate normal etc.)
eval_pdf.default <- function(data, densfun, ...) {
  params <- list(...)
  f <- densfun[[1]]
  mean <- densfun[[2]]
  cov <- densfun[[3]]
  Log <- params[[1]]
  f(data, mean, cov, Log)
}
eval_pdf.one_param <- function(data, densfun, ...) {
  params <- list(...)
  f <- densfun[[1]]
  mean <- densfun[[2]]
  Log <- params[[1]]
  f(data, mean, Log)
}
eval_pdf.three_param <- function(data, densfun, ...) {
  params <- list(...)
  f <- densfun[[1]]
  v1 <- densfun[[2]]
  v2 <- densfun[[3]]
  v3 <- densfun[[4]]
  Log <- params[[1]]
  f(data, v1, v2, v3, Log)
}

# BNB distri
dbnb <- function(K, alpha_l, a, b, Log) {
  val <- lgamma(alpha_l + K) + lbeta(alpha_l + a, K + b) - lgamma(alpha_l) - lgamma(K + 1) - lbeta(a, b)
  if (Log)
    return(val) else return(exp(val))
}
# logged Dir distri
lrdir <- function(alpha) {
  len <- length(alpha)
  y <- numeric(len)
  for (i in 1:len) {
    y[i] <- rgamma(1, alpha[i])
  }
  y <- log(y) - log(sum(y))
  return(y)
}


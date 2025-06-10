# All methods

setMethod("simulate", signature(bmbclust_obj = "bmbclust"), 
  function(bmbclust_obj, sampler = "Telescope", n.iter = 1e+05, 
    n.burnin = 1e+05, n.thin = 1, n.chains = 1, verbose = TRUE, 
    show.prior.Kplus = FALSE) {
    # This layer is useful once we implement another sampler
    # or optimaztion algorithm to do similar task
    if (sampler == "Telescope") {
      Telescope_sampler(bmbclust_obj, n.iter = n.iter, 
        n.burnin = n.burnin, n.thin = n.thin, n.chains = n.chains, 
        verbose = verbose, show.prior.Kplus = show.prior.Kplus)
    } else {
      stop(paste("no sampler named", sampler, "exists, please try from one of Telescope, INSERT SAMPLER NAME or INSERT SAMPLER NAME"))
    }
  })


setMethod("Telescope_sampler", signature(bmbclust_obj = "bmbclust"), 
  function(bmbclust_obj, n.iter, n.burnin, n.thin, n.chains, 
    verbose, show.prior.Kplus) {
    # prepare data and initial specifications for K and
    # alpha
    if (length(bmbclust_obj@posteriors) == 0) {
      idx <- 1
    } else {
      idx <- length(bmbclust_obj@posteriors$Ks) + 
        1
    }
    data <- data.matrix(bmbclust_obj@mixture.obj@data)
    bmbclust_obj@mixture.obj@data <- data
    r <- ncol(data)
    
    K_max <- bmbclust_obj@max.k
    # prepare other initial specifications not necessarily
    # same across various component densities
    bmbclust_obj <- Init_sampler(bmbclust_obj@mixture.obj, 
      bmbclust_obj)
    bmbclust_obj <- Init_params(bmbclust_obj@mixture.obj, 
      bmbclust_obj)
    
    if (length(bmbclust_obj@posteriors) == 0) {
      alpha <- bmbclust_obj@mixture.obj@init.alpha
      K <- bmbclust_obj@mixture.obj@init.k
      inits <- grep("init.*", names(attributes(bmbclust_obj@mixture.obj)), 
        value = TRUE)
      init_list <- lapply(setNames(inits, inits), 
        function(x) do.call("@", list(bmbclust_obj@mixture.obj, 
          x)))
      bmbclust_obj@inits <- init_list
    } else {
      last_idx <- length(bmbclust_obj@posteriors$Ks)
      alpha <- bmbclust_obj@posteriors$alphas[last_idx]
      K <- bmbclust_obj@posteriors$Ks[last_idx]
    }
    # sample weights
    etas <- lrdir(rep(1/K, K))
    
    if (verbose) 
      first_time <- TRUE
    
    bmbclust_obj <- Init_post_list(bmbclust_obj@mixture.obj, 
      bmbclust_obj, n.iter, n.thin)
    
    for (n in 1:(n.iter + n.burnin)) {
      tmp <- Index_update(bmbclust_obj@mixture.obj, 
        bmbclust_obj, etas, K)
      S <- tmp[[1]]
      alloc_probs <- tmp[[2]]
      # number of data clusters
      K_plus <- length(unique(S))
      
      # reorder elements
      mapping <- 1:K_plus
      # unique will return unique elements in the order of
      # appearance
      names(mapping) <- unique(S)
      
      not_filled <- seq(1, K)[!(1:K %in% unique(S))]
      if (length(not_filled) != 0) {
        mapping <- c(mapping, (K_plus + 1):K)
        names(mapping)[(K_plus + 1):K] <- not_filled
      }
      bmbclust_obj <- Reorder_components(bmbclust_obj@mixture.obj, 
        bmbclust_obj, mapping, K)
      S <- sapply(S, function(x) which(x == as.numeric(names(mapping))))
      N <- c(tabulate(S), rep(0, K - K_plus))
      
      if ((n > n.burnin) && ((n - n.burnin)%%n.thin == 
        0)) {
        loglik <- sum(alloc_probs)  #eval_partition_loglik(bmbclust_obj,K,S,etas)
        bmbclust_obj <- record_post_draws(bmbclust_obj@mixture.obj, 
          bmbclust_obj, idx, etas, K, K_plus, alpha, 
          N, loglik, S, alloc_probs)
      }
      
      bmbclust_obj <- Component_update(bmbclust_obj@mixture.obj, 
        bmbclust_obj, K_plus, N, S)
      if (!bmbclust_obj@fix.k) {
        unnorm_logprob <- sapply(K_plus:K_max, function(k) {
          do.call(bmbclust_obj@prior.k, append(list(x = k - 
          1), append(bmbclust_obj@prior.k.parameters, 
          list(log = TRUE)))) + lfactorial(k) - 
          lfactorial(k - K_plus) + sum(log(k/(N * 
          k + alpha)) + lgamma(1 + N + alpha/k) - 
          log(k/alpha) - lgamma(1 + alpha/k))
        })
        K <- which.max(evd::rgumbel(length(unnorm_logprob)) + 
          unnorm_logprob) + K_plus - 1
      }
      alpha_prop <- exp(rnorm(1, log(alpha), bmbclust_obj@mixture.obj@sig_alpha))
      tmp <- sapply(c(alpha_prop, alpha), function(a) {
        do.call(bmbclust_obj@prior.alpha, append(list(x = a), 
          append(bmbclust_obj@prior.alpha.parameters, 
          list(log = TRUE)))) + lgamma(a) - lgamma(sum(N) + 
          a) + sum(log(K/(N * K + a)) + lgamma(1 + 
          N + a/K) - log(K/a) - lgamma(1 + a/K))
      })
      logr <- tmp[1] - tmp[2] + log(alpha_prop) - 
        log(alpha)
      accept <- log(runif(1)) < logr
      if (accept) {
        alpha <- alpha_prop
      }
      bmbclust_obj <- Sample_component_params(bmbclust_obj@mixture.obj, 
        bmbclust_obj, K_plus, K)
      N <- c(N[1:K_plus], rep(0, K - K_plus))
      etas <- lrdir(alpha/K + N)
      
      if ((n > n.burnin) && ((n - n.burnin)%%n.thin == 
        0)) {
        # loglik <- eval_partition_loglik(bmbclust_obj,K,S,etas)
        # bmbclust_obj <-
        # record_post_draws(bmbclust_obj,idx,etas,K,K_plus,alpha,N,loglik,S,alloc_probs)
        if (verbose) {
          if (first_time) {
          x11(width = 7, height = 5)
          temporal_partitions <- matrix(prop.table(sort(N, 
            decreasing = TRUE)), nrow = 1)
          first_time <- FALSE
          par(mfrow = c(3, 3), mar = c(3, 3, 3, 
            0), mgp = c(1.75, 0.75, 0), oma = c(1, 
            0.1, 0.3, 0.3))
          layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 
            7, 7), 3, 3, byrow = TRUE), heights = c(3, 
            3, 3), widths = c(1.5, 0.9, 0.6))
          } else {
          if (length(N) > ncol(temporal_partitions)) {
            diff_ncols <- length(N) - ncol(temporal_partitions)
            add_col <- matrix(0, nrow = nrow(temporal_partitions), 
            ncol = diff_ncols)
            temporal_partitions <- cbind(temporal_partitions, 
            add_col)
          } else if (length(N) < ncol(temporal_partitions)) {
            diff_ncols <- ncol(temporal_partitions) - 
            length(N)
            N <- c(N, rep(0, diff_ncols))
          }
          temporal_partitions <- rbind(temporal_partitions, 
            prop.table(matrix(sort(N, decreasing = TRUE), 
            nrow = 1)))
          }
          dev.hold()
          
          # traceplot for alpha/K col = 100 is blue
          matplot(log(bmbclust_obj@posteriors$alphas[1:idx]/bmbclust_obj@posteriors$Ks[1:idx]), 
          type = "l", lty = 1, col = 100, ylab = "log(alpha/K)", 
          main = "Traceplot of log(alpha/K)")
          abline(h = 0, col = "black")
          # Traceplot and Histogram of Kplus/K (it may also be
          # useful to overlay the corresponding prior?)
          matplot(bmbclust_obj@posteriors$Kps[1:idx]/bmbclust_obj@posteriors$Ks[1:idx], 
          type = "l", lty = 1, col = 100, lwd = 0.3, 
          ylab = "Kplus/K", main = "Traceplot of Kplus/K")
          par(mar = c(3, 0, 3, 0))
          hist_Kp_K <- hist(bmbclust_obj@posteriors$Kps[1:idx]/bmbclust_obj@posteriors$Ks[1:idx], 
          plot = FALSE)
          barplot(hist_Kp_K$density, horiz = TRUE, 
          axes = FALSE, space = 0, main = "Histogram of Kplus/K")
          par(mar = c(3, 3, 3, 0))
          # barplot of partitions from niter-10 to niter (older
          # ones are discarded for storage purposes)
          rownames(temporal_partitions) <- if (nrow(temporal_partitions) == 
          1) {
          c("current")
          } else {
          c(paste("m-", (nrow(temporal_partitions) - 
            1):1), "current")
          }
          barplot(apply(t(temporal_partitions), 
          2, rev), main = "Current and past 10 recorded cluster allocations")
          # Traceplot and Histogram of Kplus (it may laso be
          # useful to overlay the corresponding prior?)
          matplot(bmbclust_obj@posteriors$Kps[1:idx], 
          type = "l", lty = 1, col = 100, ylab = "Kplus", 
          main = "Traceplot of Kplus")
          par(mar = c(3, 0, 3, 0))
          # hist_Kps <-
          # hist(bmbclust_obj@posteriors$Kps[1:idx],plot = FALSE)
          count_Kp <- table(bmbclust_obj@posteriors$Kps[1:idx])
          count_Kp <- count_Kp/sum(count_Kp)
          # barplot(hist_Kps$density,main= 'Histogram of
          # Kplus',horiz = TRUE,axes = FALSE,space = 0)
          barplot(count_Kp, main = "Histogram of Kplus", 
          horiz = TRUE, axes = FALSE, space = 0)
          if (show.prior.Kplus) {
          probs_Kplus <- compute_prior_Kplus(bmbclust_obj, 
            idx, sum(N))
          yaxis <- min(bmbclust_obj@posteriors$Kps[1:idx]):max(bmbclust_obj@posteriors$Kps[1:idx])
          lines(exp(probs_Kplus)[yaxis], yaxis - 
            min(bmbclust_obj@posteriors$Kps[1:idx]) + 
            0.5, col = "red", type = "o")
          }
          # curve(exp(probs_Kplus[x]),from =
          # max(bmbclust_obj@posteriors$Kps[1:idx]), to =
          # max(bmbclust_obj@posteriors$Kps[1:idx]),type = 'o',
          # col = 'red', n =
          # (max(bmbclust_obj@posteriors$Kps[1:idx])-min(bmbclust_obj@posteriors$Kps[1:idx])+1),
          # add = TRUE,lwd = 2)
          par(mar = c(3, 3, 3, 0))
          thresh <- quantile(bmbclust_obj@posteriors$loglik[1:idx], 
          probs = c(0.025))
          loglik_trunc <- bmbclust_obj@posteriors$loglik[1:idx]
          loglik_trunc <- sapply(loglik_trunc, function(x) ifelse(x >= 
          thresh, x, NA_real_))
          matplot(loglik_trunc, type = "l", lty = 1, 
          col = 100, main = "Mixture log-likelihood", 
          ylab = "log-likelihood")
          # histogram of K_plus
          # hist(bmbclust_obj@posteriors$Kps[1:idx],freq =
          # FALSE,main = 'Histogram of K_plus')
          if (((n - n.burnin)/n.thin) > 10) {
          temporal_partitions <- temporal_partitions[-1, 
            ]
          }
          # histogram of K_plus/K
          
          # if(!(n == (n.iter + n.burnin))){
          dev.flush()
          # }else{ dev.flush() dev.copy(device = x11) }
        }
        idx <- idx + 1
      }
      
    }
    # par(mfrow =c(1,1))
    bmbclust_obj@posteriors$columns <- colnames(data)
    return(bmbclust_obj)
  })

setMethod("compute_prior_Kplus", signature(bmbclust_obj = "bmbclust"), 
  function(bmbclust_obj, idx, N) {
    lprobs_Kplus <- sapply(1:max(bmbclust_obj@posteriors$Kps[1:idx]), 
      function(kp) {
        matrixStats::logSumExp(sapply(kp:bmbclust_obj@max.k, 
          function(KK) {
          gamma_k <- bmbclust_obj@posteriors$alphas[idx]/KK
          val <- lfactorial(KK) - lfactorial(KK - 
            kp) + lgamma(gamma_k * KK) - lgamma(N + 
            gamma_k * KK) - kp * lgamma(gamma_k) + 
            n_partitions_Nk(kp, N, gamma_k)[1] + 
            do.call(bmbclust_obj@prior.k, append(list(x = KK - 
            1), append(bmbclust_obj@prior.k.parameters, 
            list(log = TRUE))))
          val
          })) + lfactorial(N) - lfactorial(kp)
      })
    return(lprobs_Kplus)
  })


setMethod("Init_sampler", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj) {
    if (length(mixture_obj@component.priors_hyper) == 
      0) {
      data <- mixture_obj@data
      r <- ncol(data)
      
      c0 <- 2.5 + (r - 1)/2
      g0 <- 0.5 + (r - 1)/2
      R <- apply(apply(data, 2, range), 2, function(x) x[2] - 
        x[1])
      G0 <- 100 * g0/c0 * diag(1/R^2, nrow = r)
      tmp <- matrix(runif(r^2) * 2 - 1, ncol = r)
      C0_j <- crossprod(tmp)
      B0 <- diag(R^2, nrow = r)
      inv_B0 <- diag(1/(R^2), nrow = r)
      m0 <- matrix(apply(data, 2, median), ncol = 1)
      b0 <- m0
      mixture_obj@component.priors_hyper <- list(c0 = c0, 
        g0 = g0, G0 = G0, C0_j = C0_j, B0 = B0, 
        inv_B0 = inv_B0, b0 = b0, R = R)
    }
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Init_sampler", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj) {
    if (length(mixture_obj@component.priors_hyper) == 
      0) {
      data <- mixture_obj@data
      
      a0 <- 0.1
      g0 <- 0.5
      b0 <- a0/apply(data, 2, mean)
      G0 <- g0/b0
      mixture_obj@component.priors_hyper <- list(a0 = a0, 
        b0 = b0, g0 = g0, G0 = G0)
    }
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Init_sampler", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj) {
    if ((length(mixture_obj@component.priors_hyper) == 
      0) && (all(dim(mixture_obj@data_I) == 0))) {
      data <- mixture_obj@data
      cat <- apply(data, 2, function(x) max(as.numeric(x)) - 
        min(as.numeric(x)) + 1)
      a0 <- 1
      data_I <- Reduce(cbind, lapply(1:ncol(data), 
        function(col_i) {
          simplify2array(lapply(1:cat[col_i], function(cat_i) {
          as.numeric(data[, col_i]) == cat_i
          }))
        }))
      mixture_obj@component.priors_hyper <- list(cat = cat, 
        a0 = a0)
      mixture_obj@data_I <- data_I
    }
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Init_params", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj) {
    if ((length(mixture_obj@current.meanmat) == 0) && 
      (length(mixture_obj@current.covmat) == 0)) {
      data <- mixture_obj@data
      r <- ncol(data)
      K_0 <- mixture_obj@init.k
      if (K_0 > nrow(data)) {
        K_0 <- nrow(data)
        mixture_obj@init.k <- K_0
      }
      if (length(mixture_obj@init.meanmat) == 0) {
        clst <- kmeans(data, centers = K_0)
        init.meanmat <- array(clst$centers, dim = c(K_0, 
          r, 1))
        mixture_obj@init.meanmat <- init.meanmat
      }
      if (length(mixture_obj@init.covmat) == 0) {
        init.covmat <- array(runif(r * r * K_0), 
          dim = c(K_0, r, r))
        init.covmat <- array(t(apply(init.covmat, 
          1, function(x) crossprod(x))), dim = c(K_0, 
          r, r))
        mixture_obj@init.covmat <- init.covmat
      }
      mixture_obj@current.meanmat <- mixture_obj@init.meanmat
      mixture_obj@current.covmat <- mixture_obj@init.covmat
    } else {
      last_idx <- length(bmbclust_obj@posteriors$Ks)
      mixture_obj@current.meanmat <- bmbclust_obj@posteriors$comp_means[[last_idx]]
      mixture_obj@current.covmat <- bmbclust_obj@posteriors$comp_covs[[last_idx]]
    }
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Init_params", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj) {
    if ((length(mixture_obj@current.meanmat) == 0)) {
      if (length(mixture_obj@init.meanmat) == 0) {
        data <- mixture_obj@data
        K_0 <- mixture_obj@init.k
        if (K_0 > length(table(data))) {
          K_0 <- length(table(data))
          mixture_obj@init.k <- K_0
        }
        clst <- kmeans(data, centers = K_0)
        
        init.meanmat <- matrix(clst$centers, ncol = 1)
        mixture_obj@init.meanmat <- init.meanmat
      }
      mixture_obj@current.meanmat <- mixture_obj@init.meanmat
    } else {
      last_idx <- length(bmbclust_obj@posteriors$Ks)
      mixture_obj@current.meanmat <- bmbclust_obj@posteriors$comp_means[[last_idx]]
    }
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Init_params", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj) {
    if (length(mixture_obj@current.meanmat) == 0) {
      if (length(mixture_obj@init.meanmat) == 0) {
        K_0 <- mixture_obj@init.k
        data <- mixture_obj@data
        
        cat <- mixture_obj@component.priors_hyper$cat
        n_unique <- length(unique(apply(data, 1, 
          function(x) paste(x, collapse = ""))))
        if (K_0 > n_unique) {
          K_0 <- n_unique
          mixture_obj@init.k <- K_0
        }
        clst <- klaR::kmodes(data, modes = K_0, 
          iter.max = 20, weighted = FALSE)
        S_0 <- clst$cluster
        init.meanmat <- Reduce(rbind, lapply(1:length(cat), 
          function(col_i) {
          # matrix of number of categories in col_i times number
          # of clusters filled with empirical probabilities
          simplify2array(lapply(1:K_0, function(k) {
            tabulate(data[S_0 == k, col_i], cat[col_i])/sum(S_0 == 
            k)
          }))
          }))
        mixture_obj@init.meanmat <- init.meanmat
      }
      mixture_obj@current.meanmat <- mixture_obj@init.meanmat
    } else {
      last_idx <- length(bmbclust_obj@posteriors$Ks)
      mixture_obj@current.meanmat <- bmbclust_obj@posteriors$comp_means[[last_idx]]
    }
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Init_post_list", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, n.iter, n.thin) {
    if (length(bmbclust_obj@posteriors) == 0) {
      bmbclust_obj@posteriors <- list(eta = list(), 
        Ks = numeric(n.iter/n.thin), Kps = numeric(n.iter/n.thin), 
        alphas = numeric(n.iter/n.thin), comp_means = list(), 
        alloc_vectors = list(), comp_covs = list(), 
        bs = list(), Bs = list(), max_kplus = 0, 
        partition = list(), MAP = list(), columns = NA_character_, 
        loglik = numeric(n.iter/n.thin), alloc_probs = list())
    } else {
      bmbclust_obj@posteriors <- list(eta = bmbclust_obj@posteriors$eta, 
        alphas = c(bmbclust_obj@posteriors$alphas, 
          numeric(n.iter/n.thin)), Ks = c(bmbclust_obj@posteriors$Ks, 
          numeric(n.iter/n.thin)), Kps = c(bmbclust_obj@posteriors$Kps, 
          numeric(n.iter/n.thin)), comp_means = bmbclust_obj@posteriors$comp_means, 
        alloc_vectors = bmbclust_obj@posteriors$alloc_vectors, 
        comp_covs = bmbclust_obj@posteriors$comp_covs, 
        bs = bmbclust_obj@posteriors$bs, Bs = bmbclust_obj@posteriors$Bs, 
        max_kplus = bmbclust_obj@posteriors$max_kplus, 
        partition = bmbclust_obj@posteriors$partition, 
        MAP = bmbclust_obj@posteriors$MAP, columns = bmbclust_obj@posteriors$columns, 
        loglik = c(bmbclust_obj@posteriors$loglik, 
          numeric(n.iter/n.thin)), alloc_probs = bmbclust_obj@posteriors$alloc_probs)
    }
    return(bmbclust_obj)
  })

setMethod("Init_post_list", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, n.iter, n.thin) {
    if (length(bmbclust_obj@posteriors) == 0) {
      bmbclust_obj@posteriors <- list(eta = list(), 
        Ks = numeric(n.iter/n.thin), Kps = numeric(n.iter/n.thin), 
        alphas = numeric(n.iter/n.thin), comp_means = list(), 
        alloc_vectors = list(), b0 = list(), max_kplus = 0, 
        partition = list(), MAP = list(), columns = NA_character_, 
        loglik = numeric(n.iter/n.thin), alloc_probs = list())
    } else {
      bmbclust_obj@posteriors <- list(eta = bmbclust_obj@posteriors$eta, 
        alphas = c(bmbclust_obj@posteriors$alphas, 
          numeric(n.iter/n.thin)), Ks = c(bmbclust_obj@posteriors$Ks, 
          numeric(n.iter/n.thin)), Kps = c(bmbclust_obj@posteriors$Kps, 
          numeric(n.iter/n.thin)), comp_means = bmbclust_obj@posteriors$comp_means, 
        alloc_vectors = bmbclust_obj@posteriors$alloc_vectors, 
        bs = bmbclust_obj@posteriors$bs, max_kplus = bmbclust_obj@posteriors$max_kplus, 
        partition = bmbclust_obj@posteriors$partition, 
        MAP = bmbclust_obj@posteriors$MAP, columns = bmbclust_obj@posteriors$columns, 
        loglik = c(bmbclust_obj@posteriors$loglik, 
          numeric(n.iter/n.thin)), alloc_probs = bmbclust_obj@posteriors$alloc_probs)
      
    }
    return(bmbclust_obj)
  })

setMethod("Init_post_list", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, n.iter, n.thin) {
    if (length(bmbclust_obj@posteriors) == 0) {
      bmbclust_obj@posteriors <- list(eta = list(), 
        Ks = numeric(n.iter/n.thin), Kps = numeric(n.iter/n.thin), 
        alphas = numeric(n.iter/n.thin), comp_means = list(), 
        alloc_vectors = list(), max_kplus = 0, partition = list(), 
        MAP = list(), columns = NA_character_, loglik = numeric(n.iter/n.thin), 
        alloc_probs = list())
    } else {
      bmbclust_obj@posteriors <- list(eta = bmbclust_obj@posteriors$eta, 
        alphas = c(bmbclust_obj@posteriors$alphas, 
          numeric(n.iter/n.thin)), Ks = c(bmbclust_obj@posteriors$Ks, 
          numeric(n.iter/n.thin)), Kps = c(bmbclust_obj@posteriors$Kps, 
          numeric(n.iter/n.thin)), comp_means = bmbclust_obj@posteriors$comp_means, 
        alloc_vectors = bmbclust_obj@posteriors$alloc_vectors, 
        max_kplus = bmbclust_obj@posteriors$max_kplus, 
        partition = bmbclust_obj@posteriors$partition, 
        MAP = bmbclust_obj@posteriors$MAP, columns = bmbclust_obj@posteriors$columns, 
        loglik = c(bmbclust_obj@posteriors$loglik, 
          numeric(n.iter/n.thin)), alloc_probs = bmbclust_obj@posteriors$alloc_probs)
      
    }
    return(bmbclust_obj)
  })

setMethod("Index_update", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, etas, K) {
    data <- mixture_obj@data
    r <- ncol(data)
    
    meanmat <- mixture_obj@current.meanmat
    covmat <- mixture_obj@current.covmat
    # long_meanmat <- Reduce(function(l,r)
    # cbind(l,r),lapply(1:K,function(x)
    # matrix(meanmat[x,,],ncol = r))) block_covmat <-
    # array2blkdiag(covmat) lik <-
    # dMvn_multi(data,long_meanmat,block_covmat,K,log=TRUE)
    # pos_prob <- sweep(lik,2,etas,'+')
    
    eta_prod_dens <- simplify2array(lapply(1:K, function(k) {
      dMvn(data, matrix(meanmat[k, , ], ncol = r), 
        matrix(covmat[k, , ], ncol = r), log = TRUE) + 
        etas[k]
    }))
    gumbel_noise <- array(evd::rgumbel(K * dim(data)[1]), 
      dim = c(dim(data)[1], 1, K))
    S <- apply(eta_prod_dens + gumbel_noise, 1, function(x) which.max(x))
    # gumbel_noise <-
    # matrix(evd::rgumbel(K*dim(data)[1]),ncol = K) S <-
    # apply(pos_prob+gumbel_noise,1,function(x)
    # which.max(x))
    return(list(S, eta_prod_dens))
  })

setMethod("Index_update", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, etas, K) {
    data <- mixture_obj@data
    
    meanmat <- mixture_obj@current.meanmat
    eta_prod_dens <- simplify2array(lapply(1:K, function(k) {
      dpois(data, meanmat[k, ], log = TRUE) + etas[k]
    }))
    
    gumbel_noise <- array(evd::rgumbel(K * dim(data)[1]), 
      dim = c(dim(data)[1], 1, K))
    S <- apply(eta_prod_dens + gumbel_noise, 1, function(x) which.max(x))
    return(list(S, eta_prod_dens))
  })

setMethod("Index_update", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, etas, K) {
    data_I <- mixture_obj@data_I
    meanmat <- mixture_obj@current.meanmat
    
    eta_prod_dens <- simplify2array(lapply(1:K, function(k) {
      apply(data_I, 1, function(x) sum(log(meanmat[x, 
        k]))) + etas[k]
    }))
    gumbel_noise <- matrix(evd::rgumbel(K * nrow(data_I)), 
      nrow = nrow(data_I))
    S <- apply(eta_prod_dens + gumbel_noise, 1, function(x) which.max(x))
    return(list(S, eta_prod_dens))
  })

setMethod("Reorder_components", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, mapping, K) {
    r <- ncol(mixture_obj@data)
    
    mixture_obj@current.meanmat <- array(mixture_obj@current.meanmat[as.numeric(names(mapping)), 
      , ], dim = c(K, r, 1))
    mixture_obj@current.covmat <- array(mixture_obj@current.covmat[as.numeric(names(mapping)), 
      , ], dim = c(K, r, r))
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Reorder_components", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, mapping, K) {
    r <- ncol(mixture_obj@data)
    
    mixture_obj@current.meanmat <- mixture_obj@current.meanmat[as.numeric(names(mapping)), 
      , drop = FALSE]
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Reorder_components", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, mapping, K) {
    mixture_obj@current.meanmat <- mixture_obj@current.meanmat[, 
      as.numeric(names(mapping)), drop = FALSE]
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Component_update", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, K_plus, N, S) {
    data <- mixture_obj@data
    r <- ncol(data)
    
    c0 <- mixture_obj@component.priors_hyper$c0
    C0_j <- mixture_obj@component.priors_hyper$C0_j
    g0 <- mixture_obj@component.priors_hyper$g0
    G0 <- mixture_obj@component.priors_hyper$G0
    B0 <- mixture_obj@component.priors_hyper$B0
    R <- mixture_obj@component.priors_hyper$R
    inv_B0 <- mixture_obj@component.priors_hyper$inv_B0
    b0 <- mixture_obj@component.priors_hyper$b0
    v1 <- mixture_obj@hyperprior.params$nu1
    v2 <- mixture_obj@hyperprior.params$nu2
    
    meanmat <- mixture_obj@current.meanmat
    covmat <- mixture_obj@current.covmat
    
    ck <- c0 + N/2
    Ck <- lapply(1:K_plus, function(k) C0_j + 0.5 * 
      crossprod(sweep(data[S == k, , drop = FALSE], 
        2, meanmat[k, , ], FUN = "-")))
    
    sigs <- tryCatch(expr = {
      lapply(1:K_plus, function(k) bayesm::rwishart(2 * 
        ck[k], 0.5 * chol2inv(chol(Ck[[k]]))))
    }, error = function(e) {
      lapply(1:K_plus, function(k) bayesm::rwishart(2 * 
        ck[k], 0.5 * chol2inv(chol(Ck[[k]] + diag(runif(dim(Ck[[k]])[1]) * 
        0.01)))))
    })
    # sigs <- lapply(1:K_plus, function(k)
    # bayesm::rwishart(2 * ck[k], 0.5 *
    # chol2inv(chol(Ck[[k]]))))
    inv_sig <- lapply(1:K_plus, function(k) sigs[[k]]$W)
    sig <- lapply(1:K_plus, function(k) sigs[[k]]$IW)
    
    updated_vals <- lapply(1:K_plus, function(k) {
      Bk <- tryCatch(expr = {
        chol2inv(chol(inv_B0 + N[k] * inv_sig[[k]]))
      }, error = function(e) {
        chol2inv(chol(inv_B0 + N[k] * inv_sig[[k]] + 
          diag(runif(dim(inv_sig[[k]])[1]) * 0.01)))
      })
      # Bk <- chol2inv(chol(inv_B0 + N[k] * inv_sig[[k]]))
      mean_yk <- apply(data, 2, function(col_i) mean(col_i[S == 
        k]))
      bk <- Bk %*% ((inv_B0) %*% b0 + inv_sig[[k]] %*% 
        mean_yk * N[k])
      mu_k <- t(chol(Bk)) %*% rnorm(r) + bk
      return(list(mu = mu_k, b = bk, B = Bk))
    })
    
    mu <- lapply(1:K_plus, function(k) updated_vals[[k]]$mu)
    bk <- lapply(1:K_plus, function(k) updated_vals[[k]]$b)
    Bk <- lapply(1:K_plus, function(k) updated_vals[[k]]$B)
    C0_j <- bayesm::rwishart(2 * (g0 + K_plus * c0), 
      0.5 * chol2inv(chol(G0 + Reduce("+", inv_sig))))$W
    
    pk <- v1 - K_plus/2
    aj <- 2 * v2
    bj <- matrix(apply(matrix(sapply(1:K_plus, function(k) (mu[[k]] - 
      b0)^2), nrow = r), 1, sum)/(as.vector(R)^2), 
      ncol = 1)
    lambdas <- sapply(bj, function(x) GIGrvg::rgig(1, 
      pk, x, aj))
    # browser() if(!is.numeric(R)) browser()
    B0 <- diag(lambdas * as.vector(R)^2/K_plus, nrow = r)
    inv_B0 <- diag(1/diag(B0), nrow = r)
    b0 <- crossprod(chol(B0), rnorm(r)) + Reduce("+", 
      mu)/K_plus
    lambdas <- sapply(bj, function(x) GIGrvg::rgig(1, 
      pk, x, aj))
    # browser() if(!is.numeric(R)) browser()
    B0 <- diag(lambdas * as.vector(R)^2/K_plus, nrow = r)
    inv_B0 <- diag(1/diag(B0), nrow = r)
    b0 <- crossprod(chol(B0), rnorm(r)) + Reduce("+", 
      mu)/K_plus
    mixture_obj@component.priors_hyper <- list(c0 = c0, 
      g0 = g0, G0 = G0, C0_j = C0_j, B0 = B0, inv_B0 = inv_B0, 
      b0 = b0, R = R)
    mixture_obj@component.params <- list(mu = mu, sig = sig, 
      bk = bk, Bk = Bk)
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Component_update", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, K_plus, N, S) {
    a0 <- mixture_obj@component.priors_hyper$a0
    b0 <- mixture_obj@component.priors_hyper$b0
    g0 <- mixture_obj@component.priors_hyper$g0
    G0 <- mixture_obj@component.priors_hyper$G0
    data <- mixture_obj@data
    
    mus <- simplify2array(lapply(1:K_plus, function(k) {
      rgamma(1, shape = a0 + N[k] * mean(data[S == 
        k, ]), rate = b0 + N[k])
    }))
    mean_yk <- simplify2array(lapply(1:K_plus, function(k) {
      colMeans(data[S == k, , drop = FALSE])
    }))
    
    ak <- a0 + mean_yk * N[1:K_plus]
    bk <- b0 + N
    
    b0 <- rgamma(1, shape = g0 + K_plus * a0, rate = G0 + 
      sum(mus))
    
    mixture_obj@component.priors_hyper$b0 <- b0
    mixture_obj@component.params <- list(mu = mus, ak = ak, 
      bk = bk)
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Component_update", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, K_plus, N, S) {
    cat <- mixture_obj@component.priors_hyper$cat
    a0 <- mixture_obj@component.priors_hyper$a0
    data <- mixture_obj@data
    mu <- Reduce(cbind, lapply(1:K_plus, function(k) {
      matrix(Reduce(function(l, r) c(l, r), lapply(1:length(cat), 
        function(x) {
          Nk_jd <- tabulate(data[S == k, x], cat[x])
          a_kj <- a0 + Nk_jd
          logprobs <- lrdir(a_kj)
          probs <- exp(logprobs)/sum(exp(logprobs))
        })), ncol = 1)
    }))
    mixture_obj@component.params <- list(mu = mu)
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Sample_component_params", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, K_plus, K) {
    r <- ncol(mixture_obj@data)
    
    b0 <- mixture_obj@component.priors_hyper$b0
    B0 <- mixture_obj@component.priors_hyper$B0
    c0 <- mixture_obj@component.priors_hyper$c0
    C0_j <- mixture_obj@component.priors_hyper$C0_j
    
    mu <- mixture_obj@component.params$mu
    sig <- mixture_obj@component.params$sig
    
    meanmat <- mixture_obj@current.meanmat
    covmat <- mixture_obj@current.covmat
    
    # sample component means either from posterior
    # (<=K_plus) or prior (>K_plus)
    meanmat <- simplify2array(lapply(1:K, function(k) if (k <= 
      K_plus) {
      mu[[k]]
    } else {
      crossprod(chol(B0), rnorm(r)) + b0
    }))
    if (is.null(dim(meanmat))) {
      meanmat <- array(meanmat, dim = c(K, 1, 1))
    } else {
      meanmat <- aperm(meanmat, c(3, 1, 2))
    }
    
    # same thing for component covs
    covmat <- simplify2array(lapply(1:K, function(k) if (k <= 
      K_plus) {
      sig[[k]]
    } else {
      bayesm::rwishart(2 * c0, 0.5 * chol2inv(chol(C0_j)))$IW
    }))
    if (is.null(dim(covmat))) {
      covmat <- array(covmat, dim = c(K, 1, 1))
    } else {
      covmat <- aperm(covmat, c(3, 1, 2))
    }
    
    mixture_obj@current.meanmat <- meanmat
    mixture_obj@current.covmat <- covmat
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Sample_component_params", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, K_plus, K) {
    mu <- mixture_obj@component.params$mu
    a0 <- mixture_obj@component.priors_hyper$a0
    b0 <- mixture_obj@component.priors_hyper$b0
    
    meanmat <- simplify2array(lapply(1:K, function(k) {
      if (k <= K_plus) {
        mu[[k]]
      } else {
        rgamma(1, shape = a0, rate = b0)
      }
    }))
    mixture_obj@current.meanmat <- matrix(meanmat)
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("Sample_component_params", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, K_plus, K) {
    mu <- mixture_obj@component.params$mu
    a0 <- mixture_obj@component.priors_hyper$a0
    cat <- mixture_obj@component.priors_hyper$cat
    
    meanmat <- Reduce(cbind, lapply(1:K, function(k) {
      if (k < K_plus) {
        mu[, k]
      } else {
        matrix(Reduce(function(l, r) c(l, r), lapply(1:length(cat), 
          function(x) {
          exp(lrdir(rep(a0, cat[x])))
          })), ncol = 1)
      }
    }), init = NULL)
    mixture_obj@current.meanmat <- meanmat
    bmbclust_obj@mixture.obj <- mixture_obj
    return(bmbclust_obj)
  })

setMethod("eval_partition_loglik", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, K, S, etas) {
    data <- mixture_obj@data
    r <- ncol(data)
    
    meanmat <- mixture_obj@current.meanmat
    covmat <- mixture_obj@current.covmat
    
    loglik <- simplify2array(lapply(1:nrow(data), function(idx) {
      dMvn(matrix(data[idx, ], ncol = r), matrix(meanmat[S[idx], 
        , ], ncol = r), matrix(covmat[S[idx], , 
        ], ncol = r), log = TRUE) + etas[S[idx]]
    }))
    return(sum(loglik))
  })

setMethod("eval_partition_loglik", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, K, S, etas) {
    data <- mixture_obj@data
    r <- ncol(data)
    
    meanmat <- mixture_obj@current.meanmat
    
    loglik <- simplify2array(lapply(1:nrow(data), function(idx) {
      dpois(matrix(data[idx, ], ncol = r), matrix(meanmat[S[idx], 
        ], ncol = r), log = TRUE) + etas[S[idx]]
    }))
    return(sum(loglik))
  })

setMethod("eval_partition_loglik", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, K, S, etas) {
    data_I <- mixture_obj@data_I
    meanmat <- mixture_obj@current.meanmat
    eta_prod_dens <- simplify2array(lapply(1:K, function(k) {
      apply(data_I, 1, function(x) sum(log(meanmat[x, 
        k]))) + etas[k]
    }))
    # loglik <- simplify2array(lapply(1:nrow(data_I),
    # function(idx){
    # sum(log(meanmat[data_I[idx,],S[idx]]))+etas[S[idx]]
    # }))
    return(sum(eta_prod_dens))
  })

setMethod("record_post_draws", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, idx, etas, K, K_plus, 
    alpha, N, loglik, S, alloc_probs) {
    bmbclust_obj@posteriors$eta[[idx]] <- etas
    bmbclust_obj@posteriors$Ks[idx] <- K
    bmbclust_obj@posteriors$Kps[idx] <- K_plus
    bmbclust_obj@posteriors$alphas[idx] <- alpha
    bmbclust_obj@posteriors$comp_means[[idx]] <- mixture_obj@current.meanmat
    bmbclust_obj@posteriors$comp_covs[[idx]] <- mixture_obj@current.covmat
    bmbclust_obj@posteriors$bs[[idx]] <- mixture_obj@component.params$bk
    bmbclust_obj@posteriors$Bs[[idx]] <- mixture_obj@component.params$Bk
    bmbclust_obj@posteriors$partition[[idx]] <- N
    bmbclust_obj@posteriors$loglik[idx] <- loglik
    bmbclust_obj@posteriors$alloc_vectors[[idx]] <- S
    bmbclust_obj@posteriors$alloc_probs[[idx]] <- alloc_probs
    return(bmbclust_obj)
  })

setMethod("record_post_draws", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, idx, etas, K, K_plus, 
    alpha, N, loglik, S, alloc_probs) {
    bmbclust_obj@posteriors$eta[[idx]] <- etas
    bmbclust_obj@posteriors$Ks[idx] <- K
    bmbclust_obj@posteriors$Kps[idx] <- K_plus
    bmbclust_obj@posteriors$alphas[idx] <- alpha
    bmbclust_obj@posteriors$comp_means[[idx]] <- mixture_obj@current.meanmat
    bmbclust_obj@posteriors$b0[[idx]] <- mixture_obj@component.params$b0
    bmbclust_obj@posteriors$partition[[idx]] <- N
    bmbclust_obj@posteriors$loglik[idx] <- loglik
    bmbclust_obj@posteriors$alloc_vectors[[idx]] <- S
    bmbclust_obj@posteriors$alloc_probs[[idx]] <- alloc_probs
    return(bmbclust_obj)
  })

setMethod("record_post_draws", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, idx, etas, K, K_plus, 
    alpha, N, loglik, S, alloc_probs) {
    bmbclust_obj@posteriors$eta[[idx]] <- etas
    bmbclust_obj@posteriors$Ks[idx] <- K
    bmbclust_obj@posteriors$Kps[idx] <- K_plus
    bmbclust_obj@posteriors$alphas[idx] <- alpha
    cat <- cumsum(mixture_obj@component.priors_hyper$cat)
    cumcat_and_dif <- rbind(cat, c(cat[1], diff(cat)))
    bmbclust_obj@posteriors$comp_means[[idx]] <- lapply(1:ncol(cumcat_and_dif), 
      function(idx) {
        mixture_obj@current.meanmat[(cumcat_and_dif[1, 
          idx] - cumcat_and_dif[2, idx] + 1):cumcat_and_dif[1, 
          idx], ]
      })
    bmbclust_obj@posteriors$partition[[idx]] <- N
    bmbclust_obj@posteriors$loglik[idx] <- loglik
    bmbclust_obj@posteriors$alloc_vectors[[idx]] <- S
    bmbclust_obj@posteriors$alloc_probs[[idx]] <- alloc_probs
    return(bmbclust_obj)
  })

setMethod("tolabel.switching", signature(bmbclust_obj = "bmbclust"), 
  function(bmbclust_obj, K_MAP = NULL, methods = NULL) {
    if (!require(label.switching)) {
      stop("To run this method, you need to have label.switching package installed")
    }
    if (is.null(methods)) {
      methods <- c("ECR", "ECR-ITERATIVE-1", "PRA", 
        "ECR-ITERATIVE-2", "STEPHENS", "DATA-BASED")
    }
    if (any(c("SJW", "AIC") %in% methods)) {
      stop("Methods 'SJW' and 'AIC' cannot be run, please exclude those and try again")
    }
    if (is.null(K_MAP)) {
      K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors$Kps))))
    }
    MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
    return_list <- return_label_switching_args(bmbclust_obj@mixture.obj, 
      bmbclust_obj, MAP_idx, K_MAP)
    
    post_lblsw <- do.call("label.switching", args = append(list(method = methods), 
      return_list))
    # turn perm matrix to numeric indicators and compare
    # across methods to obtain a vector of permutations that
    # all methods agreed on
    perm_indicator_mat <- simplify2array(lapply(post_lblsw$permutations, 
      function(method_i) {
        as.numeric(as.factor(apply(method_i, 1, 
          function(x) paste(x, collapse = ""))))
      }))
    # logical vector whether all methods agreed on
    # permutaitons or not
    nonperm <- apply(perm_indicator_mat, 1, function(x) length(unique(x)) != 
      1)
    perms <- post_lblsw$permutations[[1]]
    similarity <- post_lblsw$similarity
    return(list(perms = perms, nonperm = nonperm, similarity = similarity))
  })

setMethod("return_label_switching_args", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, MAP_idx, K_MAP) {
    mcmc.mean <- simplify2array(lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(x) {
        if (is.null(dim(x[, , 1]))) {
          matrix(x[, , 1][1:K_MAP])
        } else {
          x[, , 1][1:K_MAP, ]
        }
      }))
    mcmc.mean <- aperm(mcmc.mean, c(3, 1, 2))
    mcmc.cov <- simplify2array(lapply(bmbclust_obj@posteriors$comp_covs[MAP_idx], 
      function(x) t(simplify2array(lapply(1:K_MAP, 
        function(k) as.vector(x[k, , ]))))))
    if (ncol(bmbclust_obj@mixture.obj@data) == 1) {
      mcmc.cov <- aperm(mcmc.cov, c(3, 2, 1))
    } else {
      mcmc.cov <- aperm(mcmc.cov, c(3, 1, 2))
    }
    mcmc.weights <- t(simplify2array(lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) exp(x)[1:K_MAP])))
    mcmc.params <- abind::abind(mcmc.mean, mcmc.cov, 
      mcmc.weights, along = 3)
    p <- simplify2array(lapply(bmbclust_obj@posteriors$alloc_probs[MAP_idx], 
      function(x) apply(x[, , 1:K_MAP], 1, function(row) exp(row)/sum(exp(row)))))
    p <- aperm(p, c(3, 2, 1))
    z <- t(simplify2array(bmbclust_obj@posteriors$alloc_vectors[MAP_idx]))
    loglik_max_idx <- which.max(bmbclust_obj@posteriors$loglik[MAP_idx])
    zpivot <- bmbclust_obj@posteriors$alloc_vectors[MAP_idx][[loglik_max_idx]]
    prapivot <- mcmc.params[loglik_max_idx, , ]
    return(list(zpivot = zpivot, z = z, K = K_MAP, prapivot = prapivot, 
      p = p, mcmc = mcmc.params, data = bmbclust_obj@mixture.obj@data))
  })

setMethod("return_label_switching_args", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, MAP_idx, K_MAP) {
    mcmc.mean <- simplify2array(lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(x) x[1:K_MAP, , drop = FALSE]))
    mcmc.mean <- aperm(mcmc.mean, c(3, 1, 2))
    mcmc.weights <- t(simplify2array(lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) exp(x)[1:K_MAP])))
    mcmc.params <- abind::abind(mcmc.mean, mcmc.weights, 
      along = 3)
    p <- simplify2array(lapply(bmbclust_obj@posteriors$alloc_probs[MAP_idx], 
      function(x) apply(x[, , 1:K_MAP], 1, function(row) exp(row)/sum(exp(row)))))
    p <- aperm(p, c(3, 2, 1))
    z <- t(simplify2array(bmbclust_obj@posteriors$alloc_vectors[MAP_idx]))
    loglik_max_idx <- which.max(bmbclust_obj@posteriors$loglik[MAP_idx])
    zpivot <- bmbclust_obj@posteriors$alloc_vectors[MAP_idx][[loglik_max_idx]]
    prapivot <- mcmc.params[loglik_max_idx, , ]
    return(list(zpivot = zpivot, z = z, K = K_MAP, prapivot = prapivot, 
      p = p, mcmc = mcmc.params, data = bmbclust_obj@mixture.obj@data))
  })

setMethod("return_label_switching_args", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, MAP_idx, K_MAP) {
    if (is.null(K_MAP)) {
      K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors$Kps))))
    }
    MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
    mcmc.mean <- simplify2array(lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(x) Reduce(rbind, x)[, 1:K_MAP]))
    mcmc.mean <- aperm(mcmc.mean, c(3, 2, 1))
    mcmc.weights <- t(simplify2array(lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) exp(x)[1:K_MAP])))
    mcmc.params <- abind::abind(mcmc.mean, mcmc.weights, 
      along = 3)
    p <- simplify2array(lapply(bmbclust_obj@posteriors$alloc_probs[MAP_idx], 
      function(x) apply(x[, 1:K_MAP], 1, function(row) exp(row)/sum(exp(row)))))
    p <- aperm(p, c(3, 2, 1))
    z <- t(simplify2array(bmbclust_obj@posteriors$alloc_vectors[MAP_idx]))
    loglik_max_idx <- which.max(bmbclust_obj@posteriors$loglik[MAP_idx])
    zpivot <- bmbclust_obj@posteriors$alloc_vectors[MAP_idx][[loglik_max_idx]]
    prapivot <- mcmc.params[loglik_max_idx, , ]
    return(list(zpivot = zpivot, z = z, K = K_MAP, prapivot = prapivot, 
      p = p, mcmc = mcmc.params, data = bmbclust_obj@mixture.obj@data_I))
  })

setMethod("identify", signature(bmbclust_obj = "bmbclust"), 
  function(bmbclust_obj, K_MAP = NULL, perms = NULL, nonperm = NULL) {
    if (is.null(perms)) {
      if (is.null(K_MAP)) {
        K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors$Kps))))
      }
      tmp <- identify_w_kmeans(bmbclust_obj@mixture.obj, 
        bmbclust_obj, K_MAP)
      nonperm_rate <- tmp$nonperm_rate
      perm_idx <- tmp$perm_idx
      clst_perm_idx <- tmp$clst_perm_idx
    } else {
      if (!is.logical(nonperm)) {
        stop("'nonperm' should be a logical vector")
      }
      if (length(nonperm) != nrow(perms)) {
        stop("the number of rows in the 'perms' should equal the length of the 'nonperm'")
      }
      if (is.null(K_MAP)) {
        K_MAP <- ncol(perms)
      }
      MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
      perm_idx <- which(MAP_idx)[!nonperm]
      perm_idx <- cbind(which(!nonperm), perm_idx)
      clst_perm_idx <- perms[!nonperm, ]
      nonperm_rate <- mean(nonperm)
    }
    bmbclust_obj <- reorder_and_record(bmbclust_obj@mixture.obj, 
      bmbclust_obj, K_MAP, perm_idx, clst_perm_idx)
    bmbclust_obj@post.identification$K_MAP <- K_MAP
    bmbclust_obj@post.identification$nonperm_rate <- nonperm_rate
    return(bmbclust_obj)
  })

setMethod("reorder_and_record", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, K_MAP, perm_idx, 
    clst_perm_idx) {
    res_ordered_mean <- lapply(1:nrow(perm_idx), function(x) {
      if (ncol(bmbclust_obj@mixture.obj@data) == 1) {
        bmbclust_obj@posteriors$comp_means[[perm_idx[x, 
          2]]][, , 1][clst_perm_idx[x, ]]
      } else {
        bmbclust_obj@posteriors$comp_means[[perm_idx[x, 
          2]]][, , 1][clst_perm_idx[x, ], ]
      }
    })
    
    res_ordered_cov <- lapply(1:nrow(perm_idx), function(x) {
      if (ncol(bmbclust_obj@mixture.obj@data) == 1) {
        bmbclust_obj@posteriors$comp_covs[[perm_idx[x, 
          2]]][clst_perm_idx[x, ], , , drop = FALSE]
      } else {
        bmbclust_obj@posteriors$comp_covs[[perm_idx[x, 
          2]]][clst_perm_idx[x, ], , ]
      }
    })
    
    res_ordered_etas <- lapply(1:nrow(perm_idx), function(x) {
      bmbclust_obj@posteriors$etas[[perm_idx[x, 2]]][clst_perm_idx[x, 
        ]]
    })
    res_ordered_mean <- lapply(res_ordered_mean, function(x) {
      if (is.null(dim(x))) {
        tmp <- matrix(x)
        colnames(tmp) <- bmbclust_obj@posteriors$columns
        tmp
      } else {
        colnames(x) <- bmbclust_obj@posteriors$columns
        x
      }
    })
    
    res_ordered_cov <- lapply(res_ordered_cov, function(x) {
      lapply(1:K_MAP, function(k) {
        if (ncol(bmbclust_obj@mixture.obj@data) == 
          1) {
          tmp <- x[k, , , drop = FALSE]
        } else {
          tmp <- x[k, , ]
        }
        colnames(tmp) <- bmbclust_obj@posteriors$columns
        rownames(tmp) <- bmbclust_obj@posteriors$columns
        tmp
      })
    })
    bmbclust_obj@post.identification$mean <- res_ordered_mean
    bmbclust_obj@post.identification$cov <- res_ordered_cov
    bmbclust_obj@post.identification$etas <- res_ordered_etas
    return(bmbclust_obj)
  })

setMethod("reorder_and_record", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, K_MAP, perm_idx, 
    clst_perm_idx) {
    res_ordered_mean <- lapply(1:nrow(perm_idx), function(x) {
      if (ncol(bmbclust_obj@mixture.obj@data) == 1) {
        bmbclust_obj@posteriors$comp_means[[perm_idx[x, 
          2]]][, , 1][clst_perm_idx[x, ]]
      } else {
        bmbclust_obj@posteriors$comp_means[[perm_idx[x, 
          2]]][, , 1][clst_perm_idx[x, ], ]
      }
    })
    
    res_ordered_etas <- lapply(1:nrow(perm_idx), function(x) {
      bmbclust_obj@posteriors$etas[[perm_idx[x, 2]]][clst_perm_idx[x, 
        ]]
    })
    res_ordered_mean <- lapply(res_ordered_mean, function(x) {
      if (is.null(dim(x))) {
        tmp <- matrix(x)
        colnames(tmp) <- bmbclust_obj@posteriors$columns
        tmp
      } else {
        colnames(x) <- bmbclust_obj@posteriors$columns
        x
      }
    })
    
    bmbclust_obj@post.identification$mean <- res_ordered_mean
    bmbclust_obj@post.identification$etas <- res_ordered_etas
    return(bmbclust_obj)
  })

setMethod("reorder_and_record", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, K_MAP, perm_idx, 
    clst_perm_idx) {
    res_ordered_mean <- lapply(1:nrow(perm_idx), function(x) {
      lapply(bmbclust_obj@posteriors$comp_means[[perm_idx[x, 
        2]]], function(M) {
        M[, clst_perm_idx[x, ]]
      })
    })
    cat <- cumsum(mixture_obj@component.priors_hyper$cat)
    res_ordered_mean <- lapply(res_ordered_mean, function(x) {
      names(x) <- bmbclust_obj@posteriors$columns
      x
    })
    res_ordered_etas <- lapply(1:nrow(perm_idx), function(x) {
      bmbclust_obj@posteriors$etas[[perm_idx[x, 2]]][clst_perm_idx[x, 
        ]]
    })
    bmbclust_obj@post.identification$mean <- res_ordered_mean
    bmbclust_obj@post.identification$etas <- res_ordered_etas
    return(bmbclust_obj)
  })

setMethod("identify_w_kmeans", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, K_MAP) {
    MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
    res_df_mean <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(x) {
        if (is.null(dim(x[, , 1]))) {
          matrix(x[, , 1][1:K_MAP])
        } else {
          x[, , 1][1:K_MAP, ]
        }
      }))
    res_df_covdet <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_covs[MAP_idx], 
      function(x) {
        if (is.null(dim(x[, , 1]))) {
          matrix(log(x[, , 1][1:K_MAP]))
        } else {
          matrix(apply(x[1:K_MAP, , ], 1, function(M) log(det(M))))
        }
      }))
    res_df_eta <- Reduce(rbind, lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) {
        matrix(x[1:K_MAP])
      }))
    res_df <- cbind(res_df_mean, res_df_mean, res_df_eta)
    res_df <- scale(res_df)
    idx_center <- K_MAP * (which.max(bmbclust_obj@posteriors$loglik[MAP_idx]) - 
      1)
    res_center <- res_df[(idx_center + 1):(idx_center + 
      K_MAP), ]
    clst_idx <- kmeans(res_df, res_center)$cluster
    permornot <- sapply(1:sum(MAP_idx), function(x) {
      length(unique(clst_idx[((x - 1) * K_MAP + 1):(x * 
        K_MAP)])) == K_MAP
    })
    nonperm_rate <- 1 - mean(permornot)
    perm_idx <- which(MAP_idx)[permornot]
    perm_idx <- cbind(which(permornot), perm_idx)
    clst_perm_idx <- matrix(clst_idx, ncol = K_MAP, 
      byrow = TRUE)
    clst_perm_idx <- clst_perm_idx[permornot, ]
    return(list(clst_perm_idx = clst_perm_idx, perm_idx = perm_idx, 
      nonperm_rate = nonperm_rate))
  })

setMethod("identify_w_kmeans", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, K_MAP) {
    MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
    res_df_mean <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(x) {
        x[1:K_MAP, ]
      }))
    res_df_eta <- Reduce(rbind, lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) {
        matrix(x[1:K_MAP])
      }))
    res_df <- cbind(res_df_mean, res_df_eta)
    res_df <- scale(res_df)
    idx_center <- K_MAP * (which.max(bmbclust_obj@posteriors$loglik[MAP_idx]) - 
      1)
    res_center <- res_df[(idx_center + 1):(idx_center + 
      K_MAP), ]
    clst_idx <- kmeans(res_df, res_center)$cluster
    permornot <- sapply(1:sum(MAP_idx), function(x) {
      length(unique(clst_idx[((x - 1) * K_MAP + 1):(x * 
        K_MAP)])) == K_MAP
    })
    nonperm_rate <- 1 - mean(permornot)
    perm_idx <- which(MAP_idx)[permornot]
    perm_idx <- cbind(which(permornot), perm_idx)
    clst_perm_idx <- matrix(clst_idx, ncol = K_MAP, 
      byrow = TRUE)
    clst_perm_idx <- clst_perm_idx[permornot, ]
    return(list(clst_perm_idx = clst_perm_idx, perm_idx = perm_idx, 
      nonperm_rate = nonperm_rate))
  })

setMethod("identify_w_kmeans", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, K_MAP) {
    MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
    # cat <- cumsum(mixture_obj@component.priors_hyper$cat)
    res_df_mean <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(M) {
        t(Reduce(rbind, lapply(M, function(x) x[-nrow(x), 
          1:K_MAP]), init = NULL))
      }))
    res_df_eta <- Reduce(rbind, lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) {
        matrix(x[1:K_MAP])
      }))
    res_df <- cbind(res_df_mean, res_df_eta)
    res_df <- scale(res_df)
    idx_center <- K_MAP * (which.max(bmbclust_obj@posteriors$loglik[MAP_idx]) - 
      1)
    res_center <- res_df[(idx_center + 1):(idx_center + 
      K_MAP), ]
    clst_idx <- kmeans(res_df, res_center)$cluster
    permornot <- sapply(1:sum(MAP_idx), function(x) {
      length(unique(clst_idx[((x - 1) * K_MAP + 1):(x * 
        K_MAP)])) == K_MAP
    })
    nonperm_rate <- 1 - mean(permornot)
    perm_idx <- which(MAP_idx)[permornot]
    perm_idx <- cbind(which(permornot), perm_idx)
    clst_perm_idx <- matrix(clst_idx, ncol = K_MAP, 
      byrow = TRUE)
    clst_perm_idx <- clst_perm_idx[permornot, ]
    return(list(clst_perm_idx = clst_perm_idx, perm_idx = perm_idx, 
      nonperm_rate = nonperm_rate))
  })

setMethod("compute_CI", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, meanmat, K_MAP, quantiles) {
    CIs <- lapply(1:K_MAP, function(k) {
      if (is.null(dim(simplify2array(meanmat)[k, , 
        ]))) {
        tmp <- Reduce(rbind, simplify2array(meanmat)[k, 
          , , drop = FALSE], init = NULL)
        qt <- quantile(tmp, probs = quantiles)
        qt <- t(t(qt))
        colnames(qt) <- colnames(meanmat[[1]])
        qt
      } else {
        apply(t(simplify2array(meanmat)[k, , ]), 
          2, function(x) {
          
          quantile(x, probs = quantiles)
          })
      }
    })
    return(CIs)
  })

setMethod("compute_CI", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, meanmat, K_MAP, quantiles) {
    CIs <- lapply(1:K_MAP, function(k) {
      if (is.null(dim(simplify2array(meanmat)[k, , 
        ]))) {
        tmp <- Reduce(rbind, simplify2array(meanmat)[k, 
          , , drop = FALSE], init = NULL)
        qt <- quantile(tmp, probs = quantiles)
        qt <- t(t(qt))
        colnames(qt) <- colnames(meanmat[[1]])
        qt
      } else {
        apply(t(simplify2array(meanmat)[k, , ]), 
          2, function(x) {
          
          quantile(x, probs = quantiles)
          })
      }
    })
    return(CIs)
  })


setMethod("compute_CI", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, meanmat, K_MAP, quantiles) {
    cat <- cumsum(mixture_obj@component.priors_hyper$cat)
    CIs <- lapply(1:K_MAP, function(k) {
      lapply(setNames(1:length(mixture_obj@col_levels), 
        names(mixture_obj@col_levels)), function(idx) {
        x <- mixture_obj@col_levels[[idx]]
        lapply(setNames(1:length(x), x), function(y) {
          quantile(simplify2array(lapply(meanmat, 
          function(M) do.call("$", args = list(M, 
            names(cat[idx])))[y, k])), probs = quantiles)
        })
      })
    })
    return(CIs)
  })

setMethod("summary", signature(object = "bmbclust"), function(object, 
  quantiles = NULL) {
  if (length(object@post.identification) == 0) {
    stop("parameters are unindentified, please run identify() function prior to calling summary()")
  }
  if (is.null(quantiles)) {
    quantiles <- c(0.025, 0.5, 0.975)
  }
  K_MAP <- object@post.identification$K_MAP
  meanmat <- object@post.identification$mean
  CIs <- compute_CI(object@mixture.obj, meanmat, K_MAP, 
    quantiles)
  invisible(sapply(1:length(CIs), function(cluster_idx) {
    cat("Cluster ", cluster_idx, "\n")
    print(CIs[[cluster_idx]])
    cat("\n")
  }))
})


setMethod("plot", signature(x = "bmbclust"), function(x, 
  K_MAP = NULL) {
  bmbclust_obj <- x
  if (is.null(K_MAP)) {
    K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors$Kps))))
  }
  op <- par(ask = TRUE)
  if (length(bmbclust_obj@post.identification) == 0) {
    MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
    df_pointprocess <- prepare_df(bmbclust_obj@mixture.obj, 
      bmbclust_obj, K_MAP, MAP_idx)
    if (ncol(df_pointprocess) > 1) {
      combs <- combn(ncol(df_pointprocess), 2)
      for (i in 1:ncol(combs)) {
        x_name <- colnames(df_pointprocess)[combs[, 
          i][1]]
        y_name <- colnames(df_pointprocess)[combs[, 
          i][2]]
        plot_df <- df_pointprocess[, combs[, i]]
        colnames(plot_df) <- c(x_name, y_name)
        plot(plot_df, cex = 0.1, col = 100, main = "Point process representation of parameters")
      }
    } else {
      colname <- colnames(df_pointprocess_i)[1]
      stripchart(do.call("$", args = list(df_pointprocess, 
        colname)), main = paste("Strip chart of ", 
        colname), xlab = colname, col = 100, pch = "+")
    }
  } else {
    ## still in development
    means_by_clust <- lapply(1:K_MAP, function(k) {
      t(simplify2array(lapply(bmbclust_obj@post.identification$mean, 
        function(meanmat_i) {
          meanmat_i[k, ]
        })))
    })
    covs_by_clust <- lapply(1:K_MAP, function(k) {
      lapply(bmbclust_obj@post.identification$cov, 
        function(covmat_i) {
          covmat_i[[k]]
        })
    })
    min_sample <- 1000
    n_post_sample <- nrow(means_by_clust[[1]])
    sample_per_draw <- min_sample%/%n_post_sample + 
      1
    post_sample_per_clust <- lapply(1:K_MAP, function(k) {
      Reduce(function(l, r) rbind(l, r), lapply(1:n_post_sample, 
        function(i) {
          mvtnorm::rmvnorm(sample_per_draw, mean = means_by_clust[[k]][i, 
          ], sigma = covs_by_clust[[k]][[i]])
        }))
    })
    combs <- combn(ncol(post_sample_per_clust[[1]]), 
      2)
    for (i in 1:ncol(combs)) {
      plot_df <- bmbclust_obj@mixture.obj@data[, combs[, 
        i]]
      plot(plot_df, cex = 0.1, col = 100, main = "Data clusters with probabilistic cluster membership labels")
      for (k in 1:K_MAP) {
        post_sample_per_clust[[k]][, combs[, i]]
      }
    }
  }
  par(op)
})

setMethod("prepare_df", signature(mixture_obj = "BMBCMx_gauss"), 
  function(mixture_obj, bmbclust_obj, K_MAP, MAP_idx) {
    res_df_mean <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(x) {
        if (is.null(dim(x[, , 1]))) {
          matrix(x[, , 1][1:K_MAP])
        } else {
          x[, , 1][1:K_MAP, ]
        }
      }))
    res_df_covdet <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_covs[MAP_idx], 
      function(x) {
        if (is.null(dim(x[, , 1]))) {
          matrix(log(x[, , 1][1:K_MAP]))
        } else {
          matrix(apply(x[1:K_MAP, , ], 1, function(M) log(det(M))))
        }
      }))
    res_df_eta <- Reduce(rbind, lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) {
        matrix(exp(x[1:K_MAP]))
      }))
    if (ncol(res_df_mean) > 4) {
      df_pca <- prcomp(res_df_mean)
      return_df_mean <- as.data.frame(df_pca$x[, 1:4])
    } else {
      if (is.null(dim(res_df_mean))) {
        df_col_named <- matrix(res_df_mean)
        colnames(df_col_named) <- bmbclust_obj@posteriors$columns
        tmp
      } else {
        df_col_named <- res_df_mean
        colnames(df_col_named) <- bmbclust_obj@posteriors$columns
      }
      return_df_mean <- as.data.frame(df_col_named)
    }
    colnames(res_df_covdet) <- "Log determinant of the Covariance"
    return_df_covdet <- as.data.frame(res_df_covdet)
    colnames(res_df_eta) <- "Weights"
    return_df_etas <- as.data.frame(res_df_eta)
    return_dfs <- cbind(return_df_mean, return_df_covdet, 
      return_df_etas)
    return(return_dfs)
  })

setMethod("prepare_df", signature(mixture_obj = "BMBCMx_pois"), 
  function(mixture_obj, bmbclust_obj, K_MAP, MAP_idx) {
    res_df_mean <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(x) {
        if (is.null(dim(x[, , 1]))) {
          matrix(x[, , 1][1:K_MAP])
        } else {
          x[, , 1][1:K_MAP, ]
        }
      }))
    res_df_eta <- Reduce(rbind, lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) {
        matrix(exp(x[1:K_MAP]))
      }))
    if (ncol(res_df_mean) > 4) {
      df_pca <- prcomp(res_df_mean)
      return_df_mean <- as.data.frame(df_pca$x[, 1:4])
    } else {
      if (is.null(dim(res_df_mean))) {
        df_col_named <- matrix(res_df_mean)
        colnames(df_col_named) <- bmbclust_obj@posteriors$columns
        tmp
      } else {
        df_col_named <- res_df_mean
        colnames(df_col_named) <- bmbclust_obj@posteriors$columns
      }
      return_df_mean <- as.data.frame(df_col_named)
    }
    colnames(res_df_eta) <- "Weights"
    return_df_etas <- as.data.frame(res_df_eta)
    return_dfs <- cbind(return_df_mean, return_df_etas)
    return(return_dfs)
  })

setMethod("prepare_df", signature(mixture_obj = "BMBCMx_categ"), 
  function(mixture_obj, bmbclust_obj, K_MAP, MAP_idx) {
    res_df_mean <- Reduce(rbind, lapply(bmbclust_obj@posteriors$comp_means[MAP_idx], 
      function(M) {
        t(Reduce(rbind, lapply(M, function(x) x[-nrow(x), 
          1:K_MAP]), init = NULL))
      }))
    data_col_levels <- unlist(lapply(1:length(mixture_obj@col_levels), 
      function(x) {
        paste(names(mixture_obj@col_levels)[x], 
          mixture_obj@col_levels[[x]], sep = "_")
      }), use.names = FALSE)
    cat <- cumsum(mixture_obj@component.priors_hyper$cat)
    res_df_eta <- Reduce(rbind, lapply(bmbclust_obj@posteriors$eta[MAP_idx], 
      function(x) {
        matrix(exp(x[1:K_MAP]))
      }))
    if (ncol(res_df_mean) > 4) {
      df_pca <- prcomp(res_df_mean)
      return_df_mean <- as.data.frame(df_pca$x[, 1:4])
    } else {
      if (is.null(dim(res_df_mean))) {
        df_col_named <- matrix(res_df_mean)
        colnames(df_col_named) <- data_col_levels[-cat]
        tmp
      } else {
        df_col_named <- res_df_mean
        colnames(df_col_named) <- data_col_levels[-cat]
      }
      return_df_mean <- as.data.frame(df_col_named)
    }
    colnames(res_df_eta) <- "Weights"
    return_df_etas <- as.data.frame(res_df_eta)
    return_dfs <- cbind(return_df_mean, return_df_etas)
    return(return_dfs)
  })


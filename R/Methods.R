# All methods

setMethod("compute_prior_Kplus",
          signature(bmbclust_obj = "bmbclust"),
          function(bmbclust_obj, idx, N){
              lprobs_Kplus <- sapply(1:max(bmbclust_obj@posteriors$Kps[1:idx]), 
                                     function(kp){ 
                                       matrixStats::logSumExp(sapply(kp:bmbclust_obj@max.k, 
                                                                     function(KK){
                                                                       gamma_k <- bmbclust_obj@posteriors$alphas[idx]/KK
                                                                       val <- lfactorial(KK) - lfactorial(KK-kp) + lgamma(gamma_k*KK) - 
                                                                         lgamma(N+gamma_k*KK) - kp*lgamma(gamma_k) + n_partitions_Nk(kp,N,gamma_k)[1] +
                                                                         do.call(bmbclust_obj@prior.k,
                                                                                 append(list(x=KK-1),
                                                                                        append(bmbclust_obj@prior.k.parameters,
                                                                                               list(log = TRUE)
                                                                                        )
                                                                                 )
                                                                         )
                                                                       return(val)  
                                                                     }
                                       )
                                       ) + lfactorial(N) - lfactorial(kp)
                                     }
              )
            return(lprobs_Kplus)
          }
)

setMethod("preprocess_data",
          signature(component_obj = "BMBCMx"),
          function(component_obj,data){
            columns <- as.list(colnames(data))
            if(is.data.frame(data)){
              data <- data.matrix(data)
            }
            numericornot <- is.numeric(data)
            if(!numericornot){
              stop("data is not numeric, numeric data is required for gaussian mixtures")
            }
            if(component_obj@component_density == "poisson"){
              # is.true instead of all.equal
              integerornot <- isTRUE(all.equal(floor(data),data))
              if(!integerornot){
                stop("not a count data, all values should be integers")
              }
            }
            return(list(data,columns))
          })


setMethod("preprocess_data",
          signature(component_obj = "BMBCMx_categ"),
          function(component_obj,data){
            #as.data.frame(x, ...,
            #stringsAsFactors = default.stringsAsFactors())
            # what purpose did as.data.frame serve?
            factorornot <- all(simplify2array(lapply(as.data.frame(data), is.factor)))
            if(!factorornot){
              stop("not a categorical data, columns should be factors")
            }
            columns <- lapply(data,levels)
            data <- Reduce(cbind,lapply(1:length(columns),function(col_i){
              simplify2array(lapply(columns[[col_i]],function(level_i){
                data[,col_i] == level_i
              }))
            }))
            return(list(data,columns))
          })

setMethod("simulate", signature(bmbclust_obj = "bmbclust"), 
          function(bmbclust_obj, sampler = "Telescope", n_iter = 1e5, n_burnin = 1e5, n_thin = 1,
                   n_chains = 1, verbose = TRUE){
            # This layer is useful once we implement another sampler or optimaztion algorithm to do similar task
            # sampler could also be a function that user specified
            if(sampler == "Telescope"){
              Telescope_sampler(bmbclust_obj, n_iter = n_iter, 
                                n_burnin = n_burnin, n_thin = n_thin, n_chains = n_chains, 
                                verbose = verbose)
            }else{
              stop(paste("no sampler named",sampler,"exists, please try from one of Telescope, INSERT SAMPLER NAME or INSERT SAMPLER NAME"))
            }
          })


setMethod("Telescope_sampler",
          signature(bmbclust_obj = "bmbclust"),
          function(bmbclust_obj, n_iter, n_burnin, n_thin, n_chains, verbose){
            # prepare data and initial specifications for K and alpha
            if(length(bmbclust_obj@posteriors_unlabeled) == 0){
              idx <- 1
            }else{
              # check this part
              idx <- length(bmbclust_obj@posteriors_unlabeled$K) + 1
            }
            
            # prepare other initial specifications not necessarily same across various component densities
            bmbclust_obj@component_obj <- Init_component(bmbclust_obj@component_obj,bmbclust_obj)

            bmbclust_obj <- Init_latent_params(bmbclust_obj@component_obj, bmbclust_obj@weight_obj,bmbclust_obj)
            
            K <- bmbclust_obj@weight_obj@K

            # sample weights
            bmbclust_obj@weight_obj@weights <- lrdir(rep(1/K, K))
            if(verbose) first_time <- TRUE
            
            bmbclust_obj <- Init_post_list(bmbclust_obj, n_iter, n_thin)
            
            for(n in 1:(n_iter + n_burnin)){
              if(!verbose){
                if(n %% 1000 == 0) print(n)
              }
              bmbclust_obj <- Index_update(bmbclust_obj)
              K_plus <- bmbclust_obj@weight_obj@K_plus
              S <- bmbclust_obj@weight_obj@cluster_allocations

              # reorder elements
              mapping <- 1:K_plus
              # unique will return unique elements in the order of appearance
              names(mapping) <- unique(S)
              # combine this with Reorder_components
              # generally make it more modular so that profiling is easy
              not_filled <- seq(1, K)[!(1:K %in% unique(S))]
              if(length(not_filled) != 0){
                mapping <- c(mapping, (K_plus+1):K)
                names(mapping)[(K_plus+1):K] <- not_filled
              }
              bmbclust_obj <- Reorder_components(bmbclust_obj@component_obj, bmbclust_obj, mapping)
              S <- sapply(S, function(x) which(x == as.numeric(names(mapping))))
              N <- c(tabulate(S), rep(0, K-K_plus))
              bmbclust_obj@weight_obj@cluster_allocations <- S
              bmbclust_obj@weight_obj@cluster_sizes <- N
              
              if((n > n_burnin) && ((n - n_burnin)%% n_thin == 0)){
                loglik <- sum(bmbclust_obj@weight_obj@allocation_probs)
                bmbclust_obj <- record_post_draws(bmbclust_obj@component_obj,bmbclust_obj@weight_obj, bmbclust_obj, idx, loglik,accept)
              }
              
              bmbclust_obj@component_obj <- Component_update(bmbclust_obj)
              
              bmbclust_obj@weight_obj <- update_K(bmbclust_obj@weight_obj,N)
              K <- bmbclust_obj@weight_obj@K
              bmbclust_obj@weight_obj <- update_gamma(bmbclust_obj@weight_obj,N)
              
              bmbclust_obj@component_obj <- Sample_component_params(bmbclust_obj@component_obj, bmbclust_obj)
              N <- c(N[1:K_plus], rep(0, K-K_plus))
              etas <- lrdir(bmbclust_obj@weight_obj@gamma+N)
              bmbclust_obj@weight_obj@cluster_sizes <- N
              bmbclust_obj@weight_obj@weights <- etas
              
              if((n > n_burnin) && ((n - n_burnin) %% n_thin == 0)){
                #loglik <- eval_partition_loglik(bmbclust_obj,K,S,etas)
                #bmbclust_obj <- record_post_draws(bmbclust_obj,idx,etas,K,K_plus,alpha,N,loglik,S,alloc_probs)
                if(verbose){
                  #browser()
                  if(first_time){
                    x11(width = 7, height = 5)
                    temporal_partitions <- matrix(prop.table(sort(N, decreasing = TRUE)), nrow = 1)
                    first_time <- FALSE
                    par(mfrow=c(3,3), mar=c(3,3,3,0), mgp=c(1.75,0.75,0), oma = c(1.0,0.1,0.3,0.3))
                    layout(matrix(c(1,2,3,4,5,6,7,7,7), 3, 3, byrow = TRUE), heights = c(3.0,3.0,3.0), widths = c(1.5,0.9,0.6))
                  }else{
                    if(length(N) > ncol(temporal_partitions)){
                      diff_ncols <- length(N) - ncol(temporal_partitions)
                      add_col <- matrix(0.0, nrow = nrow(temporal_partitions), ncol = diff_ncols)
                      temporal_partitions<- cbind(temporal_partitions,add_col)
                    }else if(length(N) < ncol(temporal_partitions)){
                      diff_ncols <- ncol(temporal_partitions) - length(N)
                      N <- c(N, rep(0, diff_ncols))
                    }
                    temporal_partitions <- rbind(temporal_partitions, prop.table(matrix(sort(N, decreasing = TRUE), nrow = 1)))
                  }
                  dev.hold()
                  
                  # traceplot for alpha/K
                  # col = 100 is blue
                  matplot(log(bmbclust_obj@posteriors_unlabeled$gamma[1:idx]),#/bmbclust_obj@posteriors_unlabeled$K[1:idx]), 
                              type = "l", lty = 1, col = 100, ylab = "log(gamma)", main = "Traceplot of log(gamma)")
                  abline(h = 0.0, col = "black")
                  
                  # Traceplot and Histogram of Kplus/K (it may also be useful to overlay the corresponding prior?)
                  matplot(bmbclust_obj@posteriors_unlabeled$K_plus[1:idx]/bmbclust_obj@posteriors_unlabeled$K[1:idx], type = "l",
                          lty = 1, col = 100, lwd = 0.3, ylab = "Kplus/K", main = "Traceplot of Kplus/K")
                  par(mar = c(3,0,3,0))
                  hist_Kp_K <- hist(bmbclust_obj@posteriors_unlabeled$K_plus[1:idx]/bmbclust_obj@posteriors_unlabeled$K[1:idx], plot = FALSE)
                  barplot(hist_Kp_K$density, horiz = TRUE, axes = FALSE, space = 0, main = "Histogram of Kplus/K")

                  # barplot of partitions from niter-10 to niter (older ones are discarded for storage purposes)
                  par(mar = c(3,3,3,0))
                  rownames(temporal_partitions) <- if(nrow(temporal_partitions) == 1){
                      c("current")
                  }else{
                      c(paste("m-",(nrow(temporal_partitions)-1):1), "current")
                    }
                  barplot(apply(t(temporal_partitions), 2, rev), main = "Current and past 10 recorded cluster allocations")
                  
                  # Traceplot and Histogram of Kplus (it may laso be useful to overlay the corresponding prior?)
                  matplot(bmbclust_obj@posteriors_unlabeled$K_plus[1:idx], type = "l", lty = 1, col = 100, ylab = "Kplus", main = "Traceplot of Kplus")
                  par(mar = c(3,0,3,0))
                  #hist_Kps <- hist(bmbclust_obj@posteriors$Kps[1:idx],plot = FALSE)
                  count_Kp <- table(bmbclust_obj@posteriors_unlabeled$K_plus[1:idx])
                  count_Kp <- count_Kp/sum(count_Kp)
                  #barplot(hist_Kps$density,main= "Histogram of Kplus",horiz = TRUE,axes = FALSE,space = 0)
                  barplot(count_Kp, main= "Histogram of Kplus", horiz = TRUE, axes = FALSE, space = 0)
                  #if(show.prior.Kplus){
                  #  probs_Kplus <- compute_prior_Kplus(bmbclust_obj, idx, sum(N))
                  #  yaxis <- min(bmbclust_obj@posteriors_unlabeled$K_plus[1:idx]):max(bmbclust_obj@posteriors_unlabeled$K_plus[1:idx])
                  #  lines(exp(probs_Kplus)[yaxis], yaxis - min(bmbclust_obj@posteriors_unlabeled$K_plus[1:idx]) + 0.5, col="red", type = "o")
                  #}
                  
                  # Traceplot of mixture log-likelihood
                  par(mar = c(3,3,3,0))
                  thresh <- quantile(bmbclust_obj@posteriors_unlabeled$loglikelihood[1:idx], probs = c(0.025))
                  loglik_trunc <- bmbclust_obj@posteriors_unlabeled$loglikelihood[1:idx]
                  loglik_trunc <- sapply(loglik_trunc, function(x) ifelse(x >= thresh, x, NA_real_))
                  matplot(loglik_trunc, type = "l", lty = 1, col = 100, main = "Mixture log-likelihood", ylab = "log-likelihood")

                  
                  if(((n-n_burnin)/n_thin) > 10){
                    temporal_partitions <- temporal_partitions[-1,]
                  }
                  
                  #if(!(n == (n_iter + n_burnin))){
                  dev.flush()
                  #  }else{
                  #    dev.flush()
                  #    dev.copy(device = x11)
                  #  }
                }
                idx <- idx + 1
              }
              
            }
            par(mfrow =c(1,1))
            #bmbclust_obj@posteriors$columns <- colnames(data)
            return(bmbclust_obj)
          }
)


setMethod("Init_component",
          signature(component_obj = "BMBCMx_gauss",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj){
            # no argument for update_rule function -> assume that user did not provide
            # update rule -> use our default specification
            empty_update_func <- is.null(names(as.list(args(component_obj@update_rule))))
            empty_compdens_func <- is.null(names(as.list(args(component_obj@component_density_func))))
            if(empty_compdens_func){
              component_obj@component_density_func <- dMvn
            }
            if(empty_update_func){
              # define our own default update based on conjugate update with GIG shrinkage
              # make it possible to move between hierarchical/non-hierarchical setting
              component_obj@update_rule <- function(data,N,K_plus,S,# data and cluster sizes
                                                    mus,sigs,b0,B0, # parameters that are updated
                                                    C0_j, #parameters that carry over from last iteration
                                                    c0,g0,G0,R,v1,v2, #hyperparameters
                                                    shrinkage # whether we add hierarchical shrinkage priors or not
                                                    ){
                inv_B0 <- chol2inv(chol(B0))
                
                # Sample posterior component precision matrix
                cks <- c0 + N/2
                #browser()
                sigs_tmp <- lapply(1:K_plus,function(k){
                  Ck <- C0_j + 0.5 * crossprod(sweep(data[S == k, , drop = FALSE], 2, mus[k, , ], FUN = "-"))
                  sig_k <- bayesm::rwishart(2 * cks[k], 0.5*chol2inv(chol(Ck)))
                  return(sig_k)
                })
                # Sample posterior component precision
                #sigs <- #tryCatch(
                #  expr = {
                #    lapply(1:K_plus, function(k) bayesm::rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[[k]]))))
                #  },
                #  error = function(e){
                #    lapply(1:K_plus, function(k) bayesm::rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[[k]]+diag(runif(dim(Ck[[k]])[1])*0.01)))))
                #  }
                #)
                # for computational conveninence
                inv_sigs <- lapply(1:K_plus, function(k) sigs_tmp[[k]]$W)
                sigs <- lapply(1:K_plus, function(k) sigs_tmp[[k]]$IW)
                
                # Sample posterior component mean
                mus <- lapply(1:K_plus, function(k) {
                  Bk <- chol2inv(chol(inv_B0 + N[k] * inv_sigs[[k]]))
                  mean_yk <- apply(data, 2, function(col_i) mean(col_i[S == k]))
                  bk <- Bk %*% ((inv_B0) %*% b0 + inv_sigs[[k]] %*% mean_yk * N[k])
                  mu_k <- crossprod(chol(Bk),rnorm(r)) + bk
                  return(mu_k)
                })
               
                # Sample hyperparameters
                C0_j <- bayesm::rwishart(2 * (g0 + K_plus * c0), 0.5 * chol2inv(chol(G0 + Reduce("+", inv_sigs))))$W
                
                # Extra steps for normal gamma shrinkage to posterior mean
                if(shrinkage){
                  pk <- v1 - K_plus/2 
                  aj <- 2 * v2
                  bj <- matrix(apply(matrix(sapply(1:K_plus, function(k) (mus[[k]] - b0)^2), nrow = r), 1, sum)/(as.vector(R)^2), ncol = 1)
                  lambdas <- sapply(bj, function(x) GIGrvg::rgig(1, pk, x, aj))
                  B0 <- diag(lambdas * as.vector(R)^2/K_plus, nrow = r)
                  b0 <- crossprod(chol(B0),rnorm(r)) + Reduce("+", mus)/K_plus
                }
                
                
                # only return parameters that will be used & updated in next iteration
                return(list(mus = mus,sigs = sigs, b0 = b0, B0 = B0, C0_j = C0_j))
              }
              
              # it first checks whether some of these parameters are specified by user
              # or not as user may not have their own update rule, but may have different 
              # set of values for hyperparameters for our default update rule
              user_hyperparams <- component_obj@update_rule_hyperparams
              r <- ncol(bmbclust_obj@data)
              c0 <- ifelse(is.null(user_hyperparams$c0), 2.5 + (r - 1)/2, user_hyperparams$c0)
              g0 <- ifelse(is.null(user_hyperparams$g0), 0.5 + (r - 1)/2, user_hyperparams$g0)
              v1 <- ifelse(is.null(user_hyperparams$v1), 0.5, user_hyperparams$v1)
              v2 <- ifelse(is.null(user_hyperparams$v2), 0.5, user_hyperparams$v2)
              shrinkage <- ifelse(is.null(user_hyperparams$shrinkage),TRUE,user_hyperparams$shrinkage)
              R <- if(is.null(user_hyperparams$R)){
                apply(apply(bmbclust_obj@data, 2, range), 2, function(x) x[2] - x[1])
              }else{
                user_hyperparams$R
              } 
              G0 <- if(is.null(user_hyperparams$G0)){
                100 * g0/c0 * diag(1/R^2, nrow = r)
              }else{
                user_hyperparams$G0
              }
              C0_j <- if(is.null(user_hyperparams$C0_j)){
                tmp <- matrix(runif(r^2) * 2 - 1, ncol = r)
                C0_j <- crossprod(tmp)
                C0_j
              }else{
                user_hyperparams$C0_j
              }
              B0 <- if(is.null(user_hyperparams$B0)){
                diag(R^2, nrow = r)
              }else{
                user_hyperparams$B0
              }
              inv_B0 <- if(is.null(user_hyperparams$B0)){
                diag(1/(R^2), nrow = r)
              }else{
                chol2inv(chol(B0))
              }
              m0 <- if(is.null(user_hyperparams$b0)){
                m0 <- matrix(apply(bmbclust_obj@data, 2, median), ncol = 1)
                b0 <- m0
                b0
              }else{
                user_hyperparams$b0
              }
              component_obj@update_rule_hyperparams <- list(c0=c0,g0=g0,G0=G0,R=R,v1 = v1, v2 = v2,shrinkage = shrinkage)
              component_obj@update_rule_params <- list(b0 = b0, B0 = B0, C0_j = C0_j)
            }
            
            return(component_obj)
          }
)

setMethod("Init_component",
          signature(component_obj = "BMBCMx_pois",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj){
            empty_update_func <- is.null(names(as.list(args(component_obj@update_rule))))
            empty_compdens_func <- is.null(names(as.list(args(component_obj@component_density_func))))
            if(empty_compdens_func){
              component_obj@component_density_func <- stats::dpois
            }
            if(empty_update_func){
              component_obj@update_rule <- function(data,N,K_plus,S,# data and cluster sizes
                                                    mus, b0, # parameters to be updated
                                                    a0, g0, G0 # fixed hyperparameters
                                                    ){
                mus <- simplify2array(lapply(1:K_plus,function(k){
                  rgamma(1, shape = a0 + N[k] * mean(data[S == k,]), rate = b0 + N[k])
                }))
                mean_yk <- simplify2array(lapply(1:K_plus,function(k){
                  colMeans(data[S == k,,drop = FALSE])
                }))
                
                # maybe useful as initial values for kmeans
                #ak <- a0 + mean_yk*N[1:K_plus]
                #bk <- b0 + N
                
                #update b0
                b0 <- rgamma(1,shape = g0+K_plus*a0,rate = G0 + sum(mus))
                return(list(mus = mus, b0 = b0))
              }
              data <- bmbclust_obj@data
              user_hyperparams <- component_obj@update_rule_hyperparams
              a0 <- ifelse(is.null(user_hyperparams$a0), 0.1, user_hyperparams$a0)
              g0 <- ifelse(is.null(user_hyperparams$g0), 0.5, user_hyperparams$g0)
              b0 <- if(is.null(user_hyperparams$b0)){
                a0/apply(bmbclust_obj@data,2,mean)
                }else{
                  user_hyperparams$b0
                } 
              G0 <- if(is.null(user_hyperparams$G0)){
                g0/b0
              }else{
                  user_hyperparams$G0
                }
              
              component_obj@update_rule_hyperparams <- list(a0 = a0, g0 = g0, G0 = G0)
              component_obj@update_rule_params <- list(b0 = b0)
            }
            return(component_obj)
          }
)

setMethod("Init_component",
          signature(component_obj = "BMBCMx_categ",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj){
            empty_update_func <- is.null(names(as.list(args(component_obj@update_rule))))
            empty_compdens_func <- is.null(names(as.list(args(component_obj@component_density_func))))
            if(empty_compdens_func){
              component_obj@component_density_func <- function(x,mu){
                apply(x,1,function(indicator) sum(log(mu[indicator])))
              }
            }
            if(empty_update_func){
              component_obj@update_rule <- function(data,N,K_plus,S,# data and cluster sizes
                                                    mus, # parameters to be updated
                                                    a0, columns # fixed hyperparameters
                                                    ){
                col_index_last <- as.vector(cumsum(sapply(columns,length)))
                col_index_first <- as.vector(c(0,col_index_last[-length(columns)])+1)
                
                mus <- Reduce(cbind,lapply(1:K_plus, function(k){
                  matrix(Reduce(function(l,r) c(l,r), lapply(1:length(columns),function(x){
                    Nk_col_i <- colSums(data[S == k,col_index_first[x]:col_index_last[x],drop = FALSE])
                    ak_col_i <- a0 + Nk_col_i
                    logprobs <- lrdir(ak_col_i)
                    probs <- exp(logprobs)/sum(exp(logprobs))
                  })),ncol = 1)
                }))
                return(list(mus = mus))
              }
              a0 <- ifelse(is.null(component_obj@update_rule_hyperparams$a0), 1, component_obj@update_rule_hyperparams$a0)
              component_obj@update_rule_hyperparams <- list(a0 = a0,columns = bmbclust_obj@columns)
            }
            return(component_obj)
          })

setMethod("Init_latent_params",
          signature(component_obj = "BMBCMx_gauss",bmbclust_obj = "bmbclust"),
          function(component_obj, weight_obj,bmbclust_obj){
            data <- bmbclust_obj@data
            K_0 <- ifelse(length(weight_obj@K)==1,weight_obj@K,nrow(data)-1)
            if(K_0>=nrow(data)){
              K_0 <-nrow(data)-1
            }
            weight_obj@K <- K_0
            if(length(weight_obj@gamma)==0){
              weight_obj@gamma <- 1.0
            }
            if((is.null(component_obj@update_rule_params$mus))||
               (is.null(component_obj@update_rule_params$sigs))){
              r <- ncol(data)
              # initialize mus
              if(is.null(component_obj@update_rule_params$mus)){
                clst <- kmeans(data, centers = K_0)
                mus <- array(clst$centers, dim = c(K_0,r,1))
                component_obj@update_rule_params$mus <- mus  
              }
              
              # initialize sigs
              if(is.null(component_obj@update_rule_params$sigs)){
                # sigs should also be defined in data driven way
                # make it scale invariant (read arxiv paper)
                #sigs <- array(runif(r*r*K_0),dim = c(K_0,r,r))
                #sigs <- array(t(apply(sigs,1,function(x) crossprod(x))), dim = c(K_0,r,r))
                if(ncol(data) == 1){
                  tmp_rng <- range(data)
                  sig <- matrix((tmp_rng[2] - tmp_rng[1])^2)
                  sigs <- array(replicate(K_0,sig),dim = c(K_0,1,1))
                }else{
                  sig <- diag(apply(apply(data, 2, range), 2, function(x) (x[2] - x[1])^2))
                  sigs <- aperm(replicate(K_0,sig),c(3,2,1))
                }
                component_obj@update_rule_params$sigs <- sigs
              }
              bmbclust_obj@component_obj <- component_obj
              bmbclust_obj@weight_obj <- weight_obj
            }
            return(bmbclust_obj)
          }
)

setMethod("Init_latent_params",
          signature(component_obj = "BMBCMx_pois",bmbclust_obj = "bmbclust"),
          function(component_obj,weight_obj,bmbclust_obj){
            data <- bmbclust_obj@data
            K_0 <- weight_obj@max_K
            if(K_0>=length(unique(data))){
              K_0 <- length(unique(data))-1
            }
            weight_obj@K <- K_0
            if(is.null(component_obj@update_rule_params$mus)){
              clst <- kmeans(data,centers = K_0)
              mus <- matrix(clst$centers,ncol =1)
              component_obj@update_rule_params$mus <- mus
            }
            bmbclust_obj@component_obj <- component_obj
            bmbclust_obj@weight_obj <- weight_obj
            return(bmbclust_obj)
          })

setMethod("Init_latent_params",
          signature(component_obj = "BMBCMx_categ",bmbclust_obj = "bmbclust"),
          function(component_obj,weight_obj,bmbclust_obj){
            data <- bmbclust_obj@data
            # unique is applicable also matrix unique(nrow(matrix))
            n_unique <- length(unique(apply(data,1,function(x) paste(x,collapse = ""))))
            K_0 <- floor(n_unique/2)
            #if(K_0> n_unique){
            #  K_0 <- n_unique
            #}
            weight_obj@K <- K_0
            if(is.null(component_obj@update_rule_params$mus)){
              columns <- bmbclust_obj@columns
              clst <- klaR::kmodes(data,modes = K_0,iter.max = 20, weighted = FALSE)
              S_0 <- clst$cluster
              
              col_index_last <- as.vector(cumsum(sapply(columns,length)))
              col_index_first <- as.vector(c(0,col_index_last[-length(columns)])+1)
          
              mus <- Reduce(cbind,lapply(1:K_0, function(k){
                matrix(Reduce(function(l,r) c(l,r), lapply(1:length(columns),function(x){
                  probs_col_i <- colSums(data[S_0 == k,col_index_first[x]:col_index_last[x],drop = FALSE])/sum(S_0 ==k)
                })),ncol = 1)
              }))
              
              component_obj@update_rule_params$mus <- mus
            }
            bmbclust_obj@component_obj <- component_obj
            bmbclust_obj@weight_obj <- weight_obj
            return(bmbclust_obj)
          })


setMethod("Init_post_list",
          signature(bmbclust_obj = "bmbclust"),
          function(bmbclust_obj,n_iter,n_thin){
            if(length(bmbclust_obj@posteriors_unlabeled) == 0){
              component_class <-bmbclust_obj@component_obj@component_density
              bmbclust_obj@posteriors_unlabeled <- 
                switch(component_class,
                       gaussian = list(
                         # across clusters parameters
                         weights = list(),
                         K = numeric(n_iter/n_thin),
                         K_plus = numeric(n_iter/n_thin),
                         gamma = numeric(n_iter/n_thin),
                         gamma_accepted = rep(NA,n_iter/n_thin),
                         loglikelihood = numeric(n_iter/n_thin),
                         cluster_sizes = list(), # N
                         cluster_allocations = list(), # S
                         allocation_probs = list(),
                         # within cluster parameters
                         component_means = list(),
                         component_covs = list()
                         ),
                       poisson = list(
                         # across clusters parameters
                         weights = list(),
                         K = numeric(n_iter/n_thin),
                         K_plus = numeric(n_iter/n_thin),
                         gamma = numeric(n_iter/n_thin),
                         gamma_accepted = rep(NA,n_iter/n_thin),
                         loglikelihood = numeric(n_iter/n_thin),
                         cluster_sizes = list(), # N
                         cluster_allocations = list(), # S
                         allocation_probs = list(),
                         # within cluster parameters
                         component_means = list()
                       ),
                       categorical = list(
                         # across clusters parameters
                         weights = list(),
                         K = numeric(n_iter/n_thin),
                         K_plus = numeric(n_iter/n_thin),
                         gamma = numeric(n_iter/n_thin),
                         gamma_accepted = rep(NA,n_iter/n_thin),
                         loglikelihood = numeric(n_iter/n_thin),
                         cluster_sizes = list(), # N
                         cluster_allocations = list(), # S
                         allocation_probs = list(),
                         # within cluster parameters
                         component_means = list()
                       ),stop(paste("No class named",component_class," found"))
                )
            }else{
              bmbclust_obj@posteriors_unlabeled$K <- c(bmbclust_obj@posteriors_unlabeled$K,numeric(n_iter/n_thin))
              bmbclust_obj@posteriors_unlabeled$K_plus <- c(bmbclust_obj@posteriors_unlabeled$K_plus,numeric(n_iter/n_thin))
              bmbclust_obj@posteriors_unlabeled$gamma <- c(bmbclust_obj@posteriors_unlabeled$gamma,numeric(n_iter/n_thin))
            }
            return(bmbclust_obj)
          })


setMethod("Index_update",
          signature(bmbclust_obj = "bmbclust"),
          function(bmbclust_obj){
            data <- bmbclust_obj@data
            etas <- bmbclust_obj@weight_obj@weights
            K <- bmbclust_obj@weight_obj@K
            r <- ncol(data)
            
            mus <- bmbclust_obj@component_obj@update_rule_params$mus
            sigs <- bmbclust_obj@component_obj@update_rule_params$sigs
            
            component_density <- bmbclust_obj@component_obj@component_density
            component_density_func <- bmbclust_obj@component_obj@component_density_func

            eta_prod_dens <- simplify2array(
              lapply(1:K,function(k){
                switch(component_density,
                       gaussian = {
                         do.call(component_density_func,args = list(x = data,
                                                                    mean = matrix(mus[k,,],ncol = r),
                                                                    sigma = matrix(sigs[k,,],ncol =r),
                                                                    log = TRUE))
                       },
                       poisson = {
                         do.call(component_density_func,args = list(x = data,
                                                                    lambda =mus[k,],
                                                                    log = TRUE))
                       },
                       categorical = {
                         matrix(do.call(component_density_func,args = list(x = data,
                                                                    mu =mus[,k]
                                                                    )),ncol = 1)
                       },
                       stop(paste("no density named",component_density)))+ etas[k]
                }
              )
            )
            
            gumbel_noise <- array(evd::rgumbel(K*dim(data)[1]),dim = c(dim(data)[1],1,K))
            S <- apply(eta_prod_dens+gumbel_noise,1,function(x) which.max(x))
            K_plus <- length(unique(S))
            bmbclust_obj@weight_obj@allocation_probs <- eta_prod_dens
            bmbclust_obj@weight_obj@cluster_allocations <- S
            bmbclust_obj@weight_obj@K_plus <- K_plus
            return(bmbclust_obj)
          }
)


setMethod("Reorder_components",
          signature(component_obj = "BMBCMx_gauss",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj,mapping){
            r <- ncol(bmbclust_obj@data)
            K <- bmbclust_obj@weight_obj@K
            component_obj@update_rule_params$mus <- array(component_obj@update_rule_params$mus[as.numeric(names(mapping)), , ], dim = c(K, r, 1))
            component_obj@update_rule_params$sigs<- array(component_obj@update_rule_params$sigs[as.numeric(names(mapping)), , ], dim = c(K, r, r))
            bmbclust_obj@component_obj <- component_obj
            return(bmbclust_obj)
          }
)



setMethod("Reorder_components",
          signature(component_obj = "BMBCMx_pois",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj,mapping){
            r <- ncol(bmbclust_obj@data)
            component_obj@update_rule_params$mus <- component_obj@update_rule_params$mus[as.numeric(names(mapping)), ,drop=FALSE]
            bmbclust_obj@component_obj <- component_obj
            return(bmbclust_obj)
          }
)

setMethod("Reorder_components",
          signature(component_obj = "BMBCMx_categ",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj,mapping){
            component_obj@update_rule_params$mus <- component_obj@update_rule_params$mus[,as.numeric(names(mapping)) ,drop=FALSE]
            bmbclust_obj@component_obj <- component_obj
            return(bmbclust_obj)
          }
)

setMethod("Component_update",
          signature(bmbclust_obj = "bmbclust"),
          function(bmbclust_obj){
            weight_obj <- bmbclust_obj@weight_obj
            component_obj <- bmbclust_obj@component_obj
            component_obj@update_rule_params <- do.call(component_obj@update_rule, 
                                                        # c() 
                                                        #vector("list", 4) 
                                                        args = append(list(data = bmbclust_obj@data, N = weight_obj@cluster_sizes,
                                                                           K_plus = weight_obj@K_plus, S = weight_obj@cluster_allocations),append(
              component_obj@update_rule_params, component_obj@update_rule_hyperparams
            )))
            return(component_obj)
          }
)


setMethod("Sample_component_params",
          signature(component_obj = "BMBCMx_gauss",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj){
            r <- length(component_obj@update_rule_params$b0)
            K_plus <- bmbclust_obj@weight_obj@K_plus
            K <- bmbclust_obj@weight_obj@K
            
            b0 <- component_obj@update_rule_params$b0
            B0 <- component_obj@update_rule_params$B0
            c0 <- component_obj@update_rule_hyperparams$c0
            C0_j <- component_obj@update_rule_params$C0_j
            
            mus <-   component_obj@update_rule_params$mus
            sigs <-  component_obj@update_rule_params$sigs
            
            
            # sample component means either from posterior (<=K_plus) or prior (>K_plus)
            mus <- simplify2array(lapply(1:K, function(k) if (k <= K_plus) {
              mus[[k]]
            } else {
              crossprod(chol(B0),rnorm(r)) + b0
            }))
            if (is.null(dim(mus))) {
              mus <- array(mus, dim = c(K, 1, 1))
            } else {
              mus <- aperm(mus, c(3, 1, 2))
            }
            # same thing for component covs
            sigs <- simplify2array(lapply(1:K, function(k) if (k <= K_plus) {
              sigs[[k]]
            } else {
              bayesm::rwishart(2 * c0, 0.5 * chol2inv(chol(C0_j)))$IW
            }))
            if (is.null(dim(sigs))) {
              sigs <- array(sigs, dim = c(K, 1, 1))
            } else {
              sigs <- aperm(sigs, c(3, 1, 2))
            }
            
            component_obj@update_rule_params$mus  <- mus
            component_obj@update_rule_params$sigs <- sigs
            return(component_obj)
          })

setMethod("Sample_component_params",
          signature(component_obj = "BMBCMx_pois",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj){
            r <- length(component_obj@update_rule_params$b0)
            K_plus <- bmbclust_obj@weight_obj@K_plus
            K <- bmbclust_obj@weight_obj@K
            
            b0 <- component_obj@update_rule_params$b0
            a0 <- component_obj@update_rule_hyperparams$a0
            
            mus <-   component_obj@update_rule_params$mus
            
            mus <- simplify2array(lapply(1:K, function(k){ if (k <= K_plus) {
              mus[[k]]
            } else {
              rgamma(1,shape = a0, rate = b0)
            }}))
            
            component_obj@update_rule_params$mus  <- matrix(mus)
            return(component_obj)
          })

setMethod("Sample_component_params",
          signature(component_obj = "BMBCMx_categ",bmbclust_obj = "bmbclust"),
          function(component_obj,bmbclust_obj){
            mus <- component_obj@update_rule_params$mus
            a0 <- component_obj@update_rule_hyperparams$a0
            columns <- bmbclust_obj@columns
            K_plus <- bmbclust_obj@weight_obj@K_plus
            K <- bmbclust_obj@weight_obj@K
            
            cat <- sapply(1:length(columns),function(x) length(columns[[x]]))

            mus <- Reduce(cbind,lapply(1:K,function(k){
              if(k<K_plus){mus[,k]}
              else{
                matrix(Reduce(function(l,r) c(l,r),lapply(1:length(cat),function(x){
                  exp(lrdir(rep(a0,cat[x])))
                  })),ncol =1)
              }
            }),init = NULL)
            component_obj@update_rule_params$mus  <- mus
            return(component_obj)
          })

setMethod("update_K",
          signature(weight_obj = "Weights"),
          function(weight_obj,N){
            return(weight_obj)
          })

setMethod("update_K",
          signature(weight_obj = "Weights_dynamic"),
          function(weight_obj,N){
            K_max <- weight_obj@max_K
            K_plus <- weight_obj@K_plus
            K <- weight_obj@K
            alpha <- weight_obj@gamma*K
            unnorm_logprob <- sapply(K_plus:K_max, function(k){
              do.call(weight_obj@prior_K,
                      append(list(x=k-1), append(weight_obj@prior_K_params, list(log = TRUE))))+
                lfactorial(k) - lfactorial(k-K_plus)+
                sum(log(k/(N*k+alpha)) + lgamma(1+N+alpha/k) - log(k/alpha) - lgamma(1+alpha/k))
            })
            K <- which.max(evd::rgumbel(length(unnorm_logprob)) + unnorm_logprob) + K_plus - 1
            weight_obj@K <- K
            weight_obj@alpha <- alpha
            return(weight_obj)
          })

setMethod("update_gamma",
          signature(weight_obj = "Weights"),
          function(weight_obj,N){
            if(!weight_obj@fix_gamma){
              K <- weight_obj@K
              gamma <- weight_obj@gamma
              sig_gamma <- weight_obj@sig_gamma
              gamma_prop <- exp(rnorm(1,log(gamma), sig_gamma))
              tmp <- sapply(c(gamma_prop,gamma),function(g){
                do.call(weight_obj@prior_gamma,c(list(x=g),weight_obj@prior_gamma_params,list(log=TRUE)))+
                  lgamma(K*g)-K*lgamma(g)+(gamma-1)*sum(weight_obj@weights)
              })
            logr <- tmp[1] - tmp[2] + log(gamma_prop) - log(gamma)
            accept <- log(runif(1)) < logr
            if(accept){
              gamma <- gamma_prop
            }
            weight_obj@gamma <- gamma
            weight_obj@accepted <- accept
            }
            return(weight_obj)
          })

setMethod("update_gamma",
          signature(weight_obj = "Weights_dynamic"),
          function(weight_obj,N){
            K <- weight_obj@K
            alpha <- weight_obj@gamma*K
            sig_alpha <- weight_obj@sig_alpha
            alpha_prop <- exp(rnorm(1, log(alpha), sig_alpha))
            tmp <- sapply(c(alpha_prop, alpha),function(a){
              do.call(weight_obj@prior_alpha,
                      append(list(x=a), append(weight_obj@prior_alpha_params,list(log=TRUE))))+
                lgamma(a) - lgamma(sum(N)+a)+
                sum(log(K/(N*K+a)) + lgamma(1+N+a/K) - log(K/a) - lgamma(1+a/K))
            })
            logr <- tmp[1] - tmp[2] + log(alpha_prop) - log(alpha)
            accept <- log(runif(1)) < logr
            if(accept){
              alpha <- alpha_prop
            }
            weight_obj@accepted <- accept
            weight_obj@alpha <- alpha
            weight_obj@gamma <- alpha/K
            return(weight_obj)
          })

setMethod("record_post_draws",
          signature(component_obj = "BMBCMx_gauss",bmbclust_obj = "bmbclust"),
          function(component_obj,weight_obj,bmbclust_obj,idx, loglik,accept){
            # Across components
            bmbclust_obj@posteriors_unlabeled$weights[[idx]] <- weight_obj@weights
            bmbclust_obj@posteriors_unlabeled$K[idx] <- weight_obj@K
            bmbclust_obj@posteriors_unlabeled$K_plus[idx] <- weight_obj@K_plus
            bmbclust_obj@posteriors_unlabeled$gamma[idx] <- weight_obj@gamma
            bmbclust_obj@posteriors_unlabeled$cluster_sizes[[idx]] <- weight_obj@cluster_sizes
            bmbclust_obj@posteriors_unlabeled$cluster_allocations[[idx]] <- weight_obj@cluster_allocations
            bmbclust_obj@posteriors_unlabeled$loglikelihood[idx] <- loglik
            bmbclust_obj@posteriors_unlabeled$allocation_probs[[idx]] <- weight_obj@allocation_probs
            bmbclust_obj@posteriors_unlabeled$accepted[idx] <- weight_obj@accepted
            # Component Specific
            bmbclust_obj@posteriors_unlabeled$component_means[[idx]] <- component_obj@update_rule_params$mus
            bmbclust_obj@posteriors_unlabeled$component_covs[[idx]] <- component_obj@update_rule_params$sigs
            return(bmbclust_obj)
          })

setMethod("record_post_draws",
          signature(component_obj = "BMBCMx_pois",bmbclust_obj = "bmbclust"),
          function(component_obj,weight_obj,bmbclust_obj,idx, loglik,accept){
            # Across components
            bmbclust_obj@posteriors_unlabeled$weights[[idx]] <- weight_obj@weights
            bmbclust_obj@posteriors_unlabeled$K[idx] <- weight_obj@K
            bmbclust_obj@posteriors_unlabeled$K_plus[idx] <- weight_obj@K_plus
            bmbclust_obj@posteriors_unlabeled$gamma[idx] <- weight_obj@gamma
            bmbclust_obj@posteriors_unlabeled$cluster_sizes[[idx]] <- weight_obj@cluster_sizes
            bmbclust_obj@posteriors_unlabeled$cluster_allocations[[idx]] <- weight_obj@cluster_allocations
            bmbclust_obj@posteriors_unlabeled$loglikelihood[idx] <- loglik
            bmbclust_obj@posteriors_unlabeled$allocation_probs[[idx]] <- weight_obj@allocation_probs
            bmbclust_obj@posteriors_unlabeled$accepted[idx] <- weight_obj@accept
            #Component_specific
            bmbclust_obj@posteriors_unlabeled$component_means[[idx]] <- component_obj@update_rule_params$mus
            return(bmbclust_obj)
          })

setMethod("record_post_draws",
          signature(component_obj = "BMBCMx_categ",bmbclust_obj = "bmbclust"),
          function(component_obj,weight_obj,bmbclust_obj,idx, loglik,accept){
            # Across components
            bmbclust_obj@posteriors_unlabeled$weights[[idx]] <- weight_obj@weights
            bmbclust_obj@posteriors_unlabeled$K[idx] <- weight_obj@K
            bmbclust_obj@posteriors_unlabeled$K_plus[idx] <- weight_obj@K_plus
            bmbclust_obj@posteriors_unlabeled$gamma[idx] <- weight_obj@gamma
            bmbclust_obj@posteriors_unlabeled$cluster_sizes[[idx]] <- weight_obj@cluster_sizes
            bmbclust_obj@posteriors_unlabeled$cluster_allocations[[idx]] <- weight_obj@cluster_allocations
            bmbclust_obj@posteriors_unlabeled$loglikelihood[idx] <- loglik
            bmbclust_obj@posteriors_unlabeled$allocation_probs[[idx]] <- weight_obj@allocation_probs
            bmbclust_obj@posteriors_unlabeled$accepted[idx] <- weight_obj@accepted
            #Component_specific
            bmbclust_obj@posteriors_unlabeled$component_means[[idx]] <- component_obj@update_rule_params$mus
            return(bmbclust_obj)
          })

setMethod("tolabel.switching",
          signature(bmbclust_obj = "bmbclust"),
          function(bmbclust_obj,K_MAP = NULL,methods = NULL){
            if(!require(label.switching)){
              stop("To run this method, you need to have label.switching package installed")
            }
            if(is.null(methods)){
              methods = c("ECR","ECR-ITERATIVE-1","PRA","ECR-ITERATIVE-2","STEPHENS","DATA-BASED")
            }
            if(any(c("SJW","AIC") %in% methods)){
              stop("Methods 'SJW' or 'AIC' cannot be run, please exclude those and try again")
            }
            if(is.null(K_MAP)){
              K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors_unlabeled$K_plus))))
            }
            MAP_idx <- bmbclust_obj@posteriors_unlabeled$K_plus == K_MAP
            return_list <- return_label_switching_args(bmbclust_obj@component_obj,bmbclust_obj,MAP_idx,K_MAP)
            
            post_lblsw <- do.call("label.switching",args = append(list(method = methods),return_list))
            # turn perm matrix to numeric indicators and compare across methods to obtain a vector of 
            # permutations that all methods agreed on
            perm_indicator_mat <- simplify2array(
              lapply(post_lblsw$permutations,function(method_i){
                as.numeric(as.factor(
                  apply(method_i,1,function(x) paste(x,collapse = ""))
                ))
              })
            )
            # logical vector whether all methods agreed on permutaitons or not
            nonperm <- apply(perm_indicator_mat,1,function(x) length(unique(x))!=1)
            perms <- post_lblsw$permutations[[1]]
            similarity <- post_lblsw$similarity
            return(list(perms = perms,nonperm = nonperm,similarity = similarity))
          })

setMethod("return_label_switching_args",
          signature(component_obj = "BMBCMx_gauss"),
          function(component_obj,bmbclust_obj,MAP_idx,K_MAP){
            mcmc.mean <- simplify2array(lapply(bmbclust_obj@posteriors_unlabeled$component_means[MAP_idx],function(x){
              if(is.null(dim(x[,,1]))){matrix(x[,,1][1:K_MAP])}else{x[,,1][1:K_MAP,]}
            }))
            mcmc.mean <- aperm(mcmc.mean,c(3,1,2))
            mcmc.cov <- simplify2array(lapply(bmbclust_obj@posteriors_unlabeled$component_covs[MAP_idx],function(x){
              t(simplify2array(lapply(1:K_MAP,function(k)as.vector(x[k,,]))))
            }))
            if(ncol(bmbclust_obj@data)==1){
              mcmc.cov <- aperm(mcmc.cov,c(3,2,1))
            }else{
              mcmc.cov <- aperm(mcmc.cov,c(3,1,2))
            }
            mcmc.weights <- t(simplify2array(lapply(bmbclust_obj@posteriors_unlabeled$weights[MAP_idx],function(x) exp(x)[1:K_MAP])))
            mcmc.params <- abind::abind(mcmc.mean,mcmc.cov,mcmc.weights,along = 3)
            p <- simplify2array(lapply(bmbclust_obj@posteriors_unlabeled$allocation_probs[MAP_idx],function(x){
              apply(x[,,1:K_MAP],1,function(row) exp(row)/sum(exp(row)))
              }))
            p <- aperm(p,c(3,2,1))
            z <- t(simplify2array(bmbclust_obj@posteriors_unlabeled$cluster_allocations[MAP_idx]))
            loglik_max_idx <- which.max(bmbclust_obj@posteriors_unlabeled$loglikelihood[MAP_idx])
            zpivot <- bmbclust_obj@posteriors_unlabeled$cluster_allocations[MAP_idx][[loglik_max_idx]]
            prapivot <- mcmc.params[loglik_max_idx,,]
            return(list(zpivot = zpivot,z = z, K = K_MAP,prapivot = prapivot,
                        p = p, mcmc = mcmc.params,data = bmbclust_obj@data))
          })

#etMethod("return_label_switching_args",
#         signature(mixture_obj = "BMBCMx_pois"),
#         function(mixture_obj,bmbclust_obj,MAP_idx,K_MAP){
#           mcmc.mean <- simplify2array(lapply(bmbclust_obj@posteriors$comp_means[MAP_idx],function(x)x[1:K_MAP,,drop=FALSE]))
#           mcmc.mean <- aperm(mcmc.mean,c(3,1,2))
#           mcmc.weights <- t(simplify2array(lapply(bmbclust_obj@posteriors$eta[MAP_idx],function(x) exp(x)[1:K_MAP])))
#           mcmc.params <- abind::abind(mcmc.mean,mcmc.weights,along =3)
#           p <- simplify2array(lapply(bmbclust_obj@posteriors$alloc_probs[MAP_idx],function(x) apply(x[,,1:K_MAP],1,function(row) exp(row)/sum(exp(row)))))
#           p <- aperm(p,c(3,2,1))
#           z <- t(simplify2array(bmbclust_obj@posteriors$alloc_vectors[MAP_idx]))
#           loglik_max_idx <- which.max(bmbclust_obj@posteriors$loglik[MAP_idx])
#           zpivot <- bmbclust_obj@posteriors$alloc_vectors[MAP_idx][[loglik_max_idx]]
#           prapivot <- mcmc.params[loglik_max_idx,,]
#           return(list(zpivot = zpivot,z = z, K = K_MAP,prapivot = prapivot,
#                       p = p, mcmc = mcmc.params,data = bmbclust_obj@mixture.obj@data))
#         })

#etMethod("return_label_switching_args",
#         signature(mixture_obj = "BMBCMx_categ"),
#         function(mixture_obj,bmbclust_obj,MAP_idx,K_MAP){
#           if(is.null(K_MAP)){
#             K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors$Kps))))
#           }
#           MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
#           mcmc.mean <- simplify2array(lapply(bmbclust_obj@posteriors$comp_means[MAP_idx],function(x)Reduce(rbind,x)[,1:K_MAP]))
#           mcmc.mean <- aperm(mcmc.mean,c(3,2,1))
#           mcmc.weights <- t(simplify2array(lapply(bmbclust_obj@posteriors$eta[MAP_idx],function(x) exp(x)[1:K_MAP])))
#           mcmc.params <- abind::abind(mcmc.mean,mcmc.weights,along =3)
#           p <- simplify2array(lapply(bmbclust_obj@posteriors$alloc_probs[MAP_idx],function(x) apply(x[,1:K_MAP],1,function(row) exp(row)/sum(exp(row)))))
#           p <- aperm(p,c(3,2,1))
#           z <- t(simplify2array(bmbclust_obj@posteriors$alloc_vectors[MAP_idx]))
#           loglik_max_idx <- which.max(bmbclust_obj@posteriors$loglik[MAP_idx])
#           zpivot <- bmbclust_obj@posteriors$alloc_vectors[MAP_idx][[loglik_max_idx]]
#           prapivot <- mcmc.params[loglik_max_idx,,]
#           return(list(zpivot = zpivot,z = z, K = K_MAP,prapivot = prapivot,
#                       p = p, mcmc = mcmc.params,data = bmbclust_obj@mixture.obj@data_I))
#         })

setMethod("identify",
          signature(bmbclust_obj="bmbclust"),
          function(bmbclust_obj,K_MAP=NULL,perms=NULL,nonperm=NULL){
            if(is.null(perms)){
              if(is.null(K_MAP)){
                K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors_unlabeled$K_plus)))) 
              }
              tmp <- identify_w_kmeans(bmbclust_obj@component_obj,bmbclust_obj,K_MAP)
              nonperm_rate <- tmp$nonperm_rate
              perm_idx <- tmp$perm_idx
              clst_perm_idx <- tmp$clst_perm_idx
            }else{
              if(!is.logical(nonperm)){
                stop("'nonperm' should be a logical vector")
              }
              if(length(nonperm)!= nrow(perms)){
                stop("the number of rows in the 'perms' should equal the length of the 'nonperm'")
              }
              if(is.null(K_MAP)){
                K_MAP <- ncol(perms)
              }
              MAP_idx <- bmbclust_obj@posteriors_unlabeled$K_plus == K_MAP
              perm_idx <- which(MAP_idx)[!nonperm]
              perm_idx <- cbind(which(!nonperm),perm_idx)
              clst_perm_idx <- perms[!nonperm,]
              nonperm_rate <- mean(nonperm)
            }
            bmbclust_obj <- reorder_and_record(bmbclust_obj@component_obj,bmbclust_obj,K_MAP,perm_idx,clst_perm_idx)
            bmbclust_obj@posteriors_unlabeled$K_MAP <- K_MAP
            bmbclust_obj@posteriors_unlabeled$nonperm_rate <- nonperm_rate
            return(bmbclust_obj)
          })

setMethod("reorder_and_record",
          signature(component_obj="BMBCMx_gauss"),
          function(component_obj,bmbclust_obj,K_MAP,perm_idx,clst_perm_idx){
            res_ordered_mean <- lapply(1:nrow(perm_idx),function(x){
              if(ncol(bmbclust_obj@data)==1){
                bmbclust_obj@posteriors_unlabeled$component_means[[perm_idx[x,2]]][,,1][clst_perm_idx[x,]] 
              }else{
                bmbclust_obj@posteriors_unlabeled$component_means[[perm_idx[x,2]]][,,1][clst_perm_idx[x,],] 
              }
            })
            
            res_ordered_cov <- lapply(1:nrow(perm_idx),function(x){
              if(ncol(bmbclust_obj@data)==1){
                bmbclust_obj@posteriors_unlabeled$component_covs[[perm_idx[x,2]]][clst_perm_idx[x,],,,drop=FALSE]
              }else{
                bmbclust_obj@posteriors_unlabeled$component_covs[[perm_idx[x,2]]][clst_perm_idx[x,],,]
              }})
            
            res_ordered_etas <- lapply(1:nrow(perm_idx),function(x){
              bmbclust_obj@posteriors_unlabeled$weights[[perm_idx[x,2]]][clst_perm_idx[x,]]
            })
            res_ordered_mean <- lapply(res_ordered_mean,function(x){
              if(is.null(dim(x))){
                tmp <- matrix(x)
                colnames(tmp) <- bmbclust_obj@columns
                tmp
              }else{
                colnames(x) <- bmbclust_obj@columns
                x
              }
            })
            
            res_ordered_cov<- lapply(res_ordered_cov,function(x){
              lapply(1:K_MAP,function(k){
                if(ncol(bmbclust_obj@data)==1){
                  tmp <- x[k,,,drop = FALSE]
                }else{
                  tmp <- x[k,,]
                }
                colnames(tmp) <- bmbclust_obj@columns
                rownames(tmp) <- bmbclust_obj@columns
                tmp
              })
            })
            bmbclust_obj@posteriors_labeled$component_means <- res_ordered_mean
            bmbclust_obj@posteriors_labeled$component_covs <- res_ordered_cov
            bmbclust_obj@posteriors_labeled$weights <- res_ordered_etas
            return(bmbclust_obj)
          })

setMethod("reorder_and_record",
          signature(component_obj="BMBCMx_pois"),
          function(component_obj,bmbclust_obj,K_MAP,perm_idx,clst_perm_idx){
            res_ordered_mean <- lapply(1:nrow(perm_idx),function(x){
              if(ncol(bmclust_obj@data)==1){
                bmbclust_obj@posteriors_unlabeled$component_means[[perm_idx[x,2]]][,,1][clst_perm_idx[x,]] 
              }else{
                bmbclust_obj@posteriors_unlabeled$component_means[[perm_idx[x,2]]][,,1][clst_perm_idx[x,],] 
              }
            })
            
            res_ordered_etas <- lapply(1:nrow(perm_idx),function(x){
              bmbclust_obj@posteriors_unlabeled$weights[[perm_idx[x,2]]][clst_perm_idx[x,]]
            })
            res_ordered_mean <- lapply(res_ordered_mean,function(x){
              if(is.null(dim(x))){
                tmp <- matrix(x)
                colnames(tmp) <- bmbclust_obj@columns
                tmp
              }else{
                colnames(x) <- bmbclust_obj@columns
                x
              }
            })
            
            bmbclust_obj@posterior_labeled$component_means <- res_ordered_mean
            bmbclust_obj@posterior_labeled$weights <- res_ordered_etas
            return(bmbclust_obj)
          })

#setMethod("reorder_and_record",
#          signature(mixture_obj="BMBCMx_categ"),
#          function(mixture_obj,bmbclust_obj,K_MAP,perm_idx,clst_perm_idx){
#            res_ordered_mean <- lapply(1:nrow(perm_idx),function(x){
#              lapply(bmbclust_obj@posteriors$comp_means[[perm_idx[x,2]]],function(M){ 
#                M[,clst_perm_idx[x,]]
#              })
#            })
#            cat <- cumsum(mixture_obj@component.priors_hyper$cat)
#            res_ordered_mean <- lapply(res_ordered_mean,function(x){
#              names(x) <- bmbclust_obj@posteriors$columns
#              x
#            })
#            res_ordered_etas <- lapply(1:nrow(perm_idx),function(x){
#              bmbclust_obj@posteriors$etas[[perm_idx[x,2]]][clst_perm_idx[x,]]
#            })
#            bmbclust_obj@post.identification$mean <- res_ordered_mean
#            bmbclust_obj@post.identification$etas <- res_ordered_etas
#            return(bmbclust_obj)
#          })
#
setMethod("identify_w_kmeans",
          signature(component_obj="BMBCMx_gauss"),
          function(component_obj,bmbclust_obj,K_MAP){
            MAP_idx <- bmbclust_obj@posteriors_unlabeled$K_plus == K_MAP
            res_df_mean <- Reduce(rbind,
                                  lapply(bmbclust_obj@posteriors_unlabeled$component_means[MAP_idx],function(x){
                                    if(is.null(dim(x[,,1]))){
                                      matrix(x[,,1][1:K_MAP])
                                    }else{
                                      x[,,1][1:K_MAP,] 
                                    }
                                  })
            )
            res_df_covdet <- Reduce(rbind,
                                    lapply(bmbclust_obj@posteriors_unlabeled$component_covs[MAP_idx],function(x){
                                      if(is.null(dim(x[,,1]))){
                                        matrix(log(x[,,1][1:K_MAP]))
                                      }else{
                                        matrix(apply(x[1:K_MAP,,],1,function(M) log(det(M))))
                                      }
                                    })
            )
            res_df_eta <- Reduce(rbind,lapply(bmbclust_obj@posteriors_unlabeled$weights[MAP_idx],function(x){
              matrix(x[1:K_MAP])
            }))
            res_df <- cbind(res_df_mean,res_df_mean,res_df_eta)
            res_df <- scale(res_df)
            idx_center <- K_MAP * (which.max(bmbclust_obj@posteriors_unlabeled$loglikelihood[MAP_idx])-1)
            res_center <- res_df[(idx_center+1):(idx_center+K_MAP),]
            clst_idx <- kmeans(res_df,res_center)$cluster
            permornot <- sapply(1:sum(MAP_idx),function(x){
              length(unique(clst_idx[((x-1)*K_MAP+1):(x*K_MAP)])) == K_MAP
            })
            nonperm_rate <- 1-mean(permornot)
            perm_idx <- which(MAP_idx)[permornot]
            perm_idx <- cbind(which(permornot),perm_idx)
            clst_perm_idx <- matrix(clst_idx,ncol = K_MAP,byrow = TRUE)
            clst_perm_idx <- clst_perm_idx[permornot,]
            return(list(clst_perm_idx=clst_perm_idx,perm_idx = perm_idx,nonperm_rate = nonperm_rate))
          })

#setMethod("identify_w_kmeans",
#          signature(mixture_obj="BMBCMx_pois"),
#          function(mixture_obj,bmbclust_obj,K_MAP){
#            MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
#            res_df_mean <- Reduce(rbind,
#                                  lapply(bmbclust_obj@posteriors$comp_means[MAP_idx],function(x){
#                                    x[1:K_MAP,]
#                                  })
#            )
#            res_df_eta <- Reduce(rbind,lapply(bmbclust_obj@posteriors$eta[MAP_idx],function(x){
#              matrix(x[1:K_MAP])
#            }))
#            res_df <- cbind(res_df_mean,res_df_eta)
#            res_df <- scale(res_df)
#            idx_center <- K_MAP * (which.max(bmbclust_obj@posteriors$loglik[MAP_idx])-1)
#            res_center <- res_df[(idx_center+1):(idx_center+K_MAP),]
#            clst_idx <- kmeans(res_df,res_center)$cluster
#            permornot <- sapply(1:sum(MAP_idx),function(x){
#              length(unique(clst_idx[((x-1)*K_MAP+1):(x*K_MAP)])) == K_MAP
#            })
#            nonperm_rate <- 1-mean(permornot)
#            perm_idx <- which(MAP_idx)[permornot]
#            perm_idx <- cbind(which(permornot),perm_idx)
#            clst_perm_idx <- matrix(clst_idx,ncol = K_MAP,byrow = TRUE)
#            clst_perm_idx <- clst_perm_idx[permornot,]
#            return(list(clst_perm_idx=clst_perm_idx,perm_idx = perm_idx,nonperm_rate = nonperm_rate))
#          })
#
#setMethod("identify_w_kmeans",
#          signature(mixture_obj="BMBCMx_categ"),
#          function(mixture_obj,bmbclust_obj,K_MAP){
#            MAP_idx <- bmbclust_obj@posteriors$Kps == K_MAP
#            #cat <- cumsum(mixture_obj@component.priors_hyper$cat)
#            res_df_mean <- Reduce(rbind,lapply(bmbclust_obj@posteriors$comp_means[MAP_idx],function(M){ 
#              t(Reduce(rbind,lapply(M,function(x)x[-nrow(x),1:K_MAP]),init=NULL))
#            })
#            )
#            res_df_eta <- Reduce(rbind,lapply(bmbclust_obj@posteriors$eta[MAP_idx],function(x){
#              matrix(x[1:K_MAP])
#            }))
#            res_df <- cbind(res_df_mean,res_df_eta)
#            res_df <- scale(res_df)
#            idx_center <- K_MAP * (which.max(bmbclust_obj@posteriors$loglik[MAP_idx])-1)
#            res_center <- res_df[(idx_center+1):(idx_center+K_MAP),]
#            clst_idx <- kmeans(res_df,res_center)$cluster
#            permornot <- sapply(1:sum(MAP_idx),function(x){
#              length(unique(clst_idx[((x-1)*K_MAP+1):(x*K_MAP)])) == K_MAP
#            })
#            nonperm_rate <- 1-mean(permornot)
#            perm_idx <- which(MAP_idx)[permornot]
#            perm_idx <- cbind(which(permornot),perm_idx)
#            clst_perm_idx <- matrix(clst_idx,ncol = K_MAP,byrow = TRUE)
#            clst_perm_idx <- clst_perm_idx[permornot,]
#            return(list(clst_perm_idx=clst_perm_idx,perm_idx = perm_idx,nonperm_rate = nonperm_rate))
#          })

setMethod("compute_CI",
          signature(component_obj ="BMBCMx_gauss"),
          function(component_obj,meanmat,K_MAP,quantiles){
            CIs <- lapply(1:K_MAP,function(k){
              if(is.null(dim(simplify2array(meanmat)[k,,]))){
                tmp <- Reduce(rbind,simplify2array(meanmat)[k,,,drop = FALSE],init = NULL)
                qt <-quantile(tmp,probs = quantiles)
                qt <- t(t(qt))
                colnames(qt) <- colnames(meanmat[[1]])
                qt
              }else{
                apply(t(simplify2array(meanmat)[k,,]),2,function(x){
                  
                  quantile(x,probs = quantiles)
                })
              }})
            return(CIs)
          })

#setMethod("compute_CI",
#          signature(mixture_obj ="BMBCMx_pois"),
#          function(mixture_obj,meanmat,K_MAP,quantiles){
#            CIs <- lapply(1:K_MAP,function(k){
#              if(is.null(dim(simplify2array(meanmat)[k,,]))){
#                tmp <- Reduce(rbind,simplify2array(meanmat)[k,,,drop = FALSE],init = NULL)
#                qt <-quantile(tmp,probs = quantiles)
#                qt <- t(t(qt))
#                colnames(qt) <- colnames(meanmat[[1]])
#                qt
#              }else{
#                apply(t(simplify2array(meanmat)[k,,]),2,function(x){
#                  
#                  quantile(x,probs = quantiles)
#                })
#              }})
#            return(CIs)
#          })
#
#
#setMethod("compute_CI",
#          signature(mixture_obj ="BMBCMx_categ"),
#          function(mixture_obj,meanmat,K_MAP,quantiles){
#            cat <- cumsum(mixture_obj@component.priors_hyper$cat)
#            CIs <- lapply(1:K_MAP,function(k) {
#              lapply(setNames(1:length(mixture_obj@col_levels),names(mixture_obj@col_levels)),function(idx){ 
#                x<-mixture_obj@col_levels[[idx]] 
#                lapply(setNames(1:length(x),x),function(y){ 
#                  quantile(simplify2array(lapply(meanmat,function(M) do.call("$",args = list(M,names(cat[idx])))[y,k])),probs = quantiles)
#                })
#              })
#            })
#            return(CIs)
#          })

setMethod("summary",
          signature(object="bmbclust"),
          function(object,quantiles=NULL){
            if(length(object@posteriors_labeled) ==0){
              stop("parameters are unindentified, please run identify() function prior to calling summary()")
            }
            if(is.null(quantiles)){quantiles = c(.025,.5,.975)}
            K_MAP <- object@posteriors_unlabeled$K_MAP
            meanmat <- object@posteriors_labeled$component_means
            CIs <- compute_CI(object@component_obj,meanmat,K_MAP,quantiles)
            invisible(
              sapply(1:length(CIs),function(cluster_idx){
                cat("Cluster ",cluster_idx,"\n")
                print(CIs[[cluster_idx]])
                cat("\n")
              })
            )
          })


setMethod("plot",
          signature(x = "bmbclust"),
          function(x,K_MAP = NULL){
            bmbclust_obj <- x
            if(is.null(K_MAP)){
              K_MAP <- as.numeric(names(which.max(table(bmbclust_obj@posteriors_unlabeled$K_plus)))) 
            }
            op <- par(ask = TRUE)
            if(length(bmbclust_obj@posteriors_labeled) == 0){
              MAP_idx <- bmbclust_obj@posteriors_unlabeled$K_plus == K_MAP
              df_pointprocess <- prepare_df(bmbclust_obj@component_obj,bmbclust_obj,K_MAP,MAP_idx)
              if(ncol(df_pointprocess)>1){
                combs <- combn(ncol(df_pointprocess),2)
                for(i in 1:ncol(combs)){
                  x_name <- colnames(df_pointprocess)[combs[,i][1]]
                  y_name <- colnames(df_pointprocess)[combs[,i][2]]
                  plot_df <- df_pointprocess[,combs[,i]]
                  colnames(plot_df) <- c(x_name,y_name)
                  plot(plot_df,cex = 0.1,col = 100,main = "Point process representation of parameters")
                }}else{
                  colname <- colnames(df_pointprocess_i)[1]
                  stripchart(do.call("$",args = list(df_pointprocess,colname)),main = paste("Strip chart of ",colname),
                             xlab = colname, col = 100,pch = "+")
                }
            }else{
              ## still in development
              means_by_clust<- lapply(1:K_MAP,function(k){
                t(
                  simplify2array(
                    lapply(bmbclust_obj@post.identification$mean,function(meanmat_i){
                      meanmat_i[k,]
                    })
                  )
                )
              })
              covs_by_clust<- lapply(1:K_MAP,function(k){
                lapply(bmbclust_obj@post.identification$cov,function(covmat_i){
                  covmat_i[[k]]
                })
              })
              min_sample <- 1000
              n_post_sample <- nrow(means_by_clust[[1]])
              sample_per_draw <- min_sample %/% n_post_sample + 1
              post_sample_per_clust <-
                lapply(1:K_MAP,function(k){
                  Reduce(function(l,r) rbind(l,r),lapply(1:n_post_sample,function(i){
                    mvtnorm::rmvnorm(sample_per_draw,mean = means_by_clust[[k]][i,],sigma = covs_by_clust[[k]][[i]])
                  }))
                })
              combs <- combn(ncol(post_sample_per_clust[[1]]),2)
              for(i in 1:ncol(combs)){
                plot_df <- bmbclust_obj@mixture.obj@data[,combs[,i]]
                plot(plot_df,cex = 0.1,col = 100, main ="Data clusters with probabilistic cluster membership labels")
                for(k in 1:K_MAP){
                  post_sample_per_clust[[k]][,combs[,i]]
                }
              }
            }
            par(op)
          })

setMethod("prepare_df",
          signature(component_obj = "BMBCMx_gauss"),
          function(component_obj,bmbclust_obj,K_MAP,MAP_idx){
            res_df_mean <- Reduce(rbind,
                                  lapply(bmbclust_obj@posteriors_unlabeled$component_means[MAP_idx],function(x){
                                    if(is.null(dim(x[,,1]))){
                                      matrix(x[,,1][1:K_MAP])
                                    }else{
                                      x[,,1][1:K_MAP,] 
                                    }
                                  })
            )
            res_df_covdet <- Reduce(rbind,
                                    lapply(bmbclust_obj@posteriors_unlabeled$component_covs[MAP_idx],function(x){
                                      if(is.null(dim(x[,,1]))){
                                        matrix(log(x[,,1][1:K_MAP]))
                                      }else{
                                        matrix(apply(x[1:K_MAP,,],1,function(M) log(det(M))))
                                      }
                                    })
            )
            res_df_eta <- Reduce(rbind,lapply(bmbclust_obj@posteriors_unlabeled$weights[MAP_idx],function(x){
              matrix(exp(x[1:K_MAP]))
            }))
            if(ncol(res_df_mean)>4){
              df_pca <- prcomp(res_df_mean)
              return_df_mean <- as.data.frame(df_pca$x[,1:4])
            }else{
              if(is.null(dim(res_df_mean))){
                df_col_named <- matrix(res_df_mean)
                colnames(df_col_named) <- bmbclust_obj@columns
                tmp
              }else{
                df_col_named <- res_df_mean
                colnames(df_col_named) <- bmbclust_obj@columns
              }
              return_df_mean <- as.data.frame(df_col_named)
            }
            colnames(res_df_covdet) <- "Log determinant of the Covariance"
            return_df_covdet <- as.data.frame(res_df_covdet)
            colnames(res_df_eta) <- "Weights"
            return_df_etas <- as.data.frame(res_df_eta)
            return_dfs <- cbind(return_df_mean,return_df_covdet,return_df_etas)
            return(return_dfs)
          })

#setMethod("prepare_df",
#          signature(mixture_obj = "BMBCMx_pois"),
#          function(mixture_obj,bmbclust_obj,K_MAP,MAP_idx){
#            res_df_mean <- Reduce(rbind,
#                                  lapply(bmbclust_obj@posteriors$comp_means[MAP_idx],function(x){
#                                    if(is.null(dim(x[,,1]))){
#                                      matrix(x[,,1][1:K_MAP])
#                                    }else{
#                                      x[,,1][1:K_MAP,] 
#                                    }
#                                  })
#            )
#            res_df_eta <- Reduce(rbind,lapply(bmbclust_obj@posteriors$eta[MAP_idx],function(x){
#              matrix(exp(x[1:K_MAP]))
#            }))
#            if(ncol(res_df_mean)>4){
#              df_pca <- prcomp(res_df_mean)
#              return_df_mean <- as.data.frame(df_pca$x[,1:4])
#            }else{
#              if(is.null(dim(res_df_mean))){
#                df_col_named <- matrix(res_df_mean)
#                colnames(df_col_named) <- bmbclust_obj@posteriors$columns
#                tmp
#              }else{
#                df_col_named <- res_df_mean
#                colnames(df_col_named) <- bmbclust_obj@posteriors$columns
#              }
#              return_df_mean <- as.data.frame(df_col_named)
#            }
#            colnames(res_df_eta) <- "Weights"
#            return_df_etas <- as.data.frame(res_df_eta)
#            return_dfs <- cbind(return_df_mean,return_df_etas)
#            return(return_dfs)
#          })
#
#setMethod("prepare_df",
#          signature(mixture_obj = "BMBCMx_categ"),
#          function(mixture_obj,bmbclust_obj,K_MAP,MAP_idx){
#            res_df_mean <- Reduce(rbind,lapply(bmbclust_obj@posteriors$comp_means[MAP_idx],function(M){ 
#              t(Reduce(rbind,lapply(M,function(x)x[-nrow(x),1:K_MAP]),init=NULL))
#            })
#            )
#            data_col_levels <- unlist(lapply(1:length(mixture_obj@col_levels),function(x){
#              paste(names(mixture_obj@col_levels)[x],mixture_obj@col_levels[[x]],sep = "_")
#            }),use.names = FALSE)
#            cat <- cumsum(mixture_obj@component.priors_hyper$cat)
#            res_df_eta <- Reduce(rbind,lapply(bmbclust_obj@posteriors$eta[MAP_idx],function(x){
#              matrix(exp(x[1:K_MAP]))
#            }))
#            if(ncol(res_df_mean)>4){
#              df_pca <- prcomp(res_df_mean)
#              return_df_mean <- as.data.frame(df_pca$x[,1:4])
#            }else{
#              if(is.null(dim(res_df_mean))){
#                df_col_named <- matrix(res_df_mean)
#                colnames(df_col_named) <- data_col_levels[-cat]
#                tmp
#              }else{
#                df_col_named <- res_df_mean
#                colnames(df_col_named) <- data_col_levels[-cat]
#              }
#              return_df_mean <- as.data.frame(df_col_named)
#            }
#            colnames(res_df_eta) <- "Weights"
#            return_df_etas <- as.data.frame(res_df_eta)
#            return_dfs <- cbind(return_df_mean,return_df_etas)
#            return(return_dfs)
#          })




#setMethod("Component_update",
#          signature(mixture_obj = "BMBCMx_categ"),
#          function(mixture_obj,bmbclust_obj,K_plus,N,S){
#            cat <-  mixture_obj@component.priors_hyper$cat
#            a0 <-   mixture_obj@component.priors_hyper$a0
#            data <- mixture_obj@data
#            mu <- Reduce(cbind,lapply(1:K_plus,function(k){
#              matrix(Reduce(function(l,r) c(l,r),lapply(1:length(cat),function(x){
#                Nk_jd <- tabulate(data[S == k,x],cat[x])
#                a_kj = a0+ Nk_jd
#                logprobs <- lrdir(a_kj)
#                probs <- exp(logprobs)/sum(exp(logprobs))
#              })),ncol = 1)
#            }))
#            mixture_obj@component.params <- list(mu = mu)
#            bmbclust_obj@mixture.obj <- mixture_obj
#            return(bmbclust_obj)
#          })
#
#setMethod("Sample_component_params",
#          signature(mixture_obj = "BMBCMx_pois"),
#          function(mixture_obj,bmbclust_obj,K_plus,K){
#            mu <- mixture_obj@component.params$mu
#            a0 <- mixture_obj@component.priors_hyper$a0
#            b0 <- mixture_obj@component.priors_hyper$b0
#            
#            meanmat <- simplify2array(lapply(1:K, function(k){ if (k <= K_plus) {
#              mu[[k]]
#            } else {
#              rgamma(1,shape = a0, rate = b0)
#            }}))
#            mixture_obj@current.meanmat <- matrix(meanmat)
#            bmbclust_obj@mixture.obj <- mixture_obj
#            return(bmbclust_obj)
#          })
#
#setMethod("Sample_component_params",
#          signature(mixture_obj = "BMBCMx_categ"),
#          function(mixture_obj,bmbclust_obj,K_plus,K){
#            mu <- mixture_obj@component.params$mu
#            a0 <- mixture_obj@component.priors_hyper$a0
#            cat <- mixture_obj@component.priors_hyper$cat
#            
#            meanmat <- Reduce(cbind,lapply(1:K,function(k){
#              if(k<K_plus){mu[,k]}
#              else{
#                matrix(Reduce(function(l,r) c(l,r),lapply(1:length(cat),function(x){exp(lrdir(rep(a0,cat[x])))})),ncol =1)
#              }
#            }),init = NULL)
#            mixture_obj@current.meanmat <- meanmat
#            bmbclust_obj@mixture.obj <- mixture_obj
#            return(bmbclust_obj)
#          })
#
#
#setMethod("record_post_draws",
#          signature(mixture_obj = "BMBCMx_categ"),
#          function(mixture_obj,bmbclust_obj,idx,etas,K,K_plus,alpha,N,loglik,S,alloc_probs){
#            bmbclust_obj@posteriors_unlabeled$eta[[idx]] <- etas
#            bmbclust_obj@posteriors_unlabeled$Ks[idx] <- K
#            bmbclust_obj@posteriors_unlabeled$Kps[idx] <- K_plus
#            bmbclust_obj@posteriors_unlabeled$alphas[idx] <- alpha
#            cat <- cumsum(mixture_obj@component.priors_hyper$cat)
#            cumcat_and_dif <- rbind(cat,c(cat[1],diff(cat)))
#            bmbclust_obj@posteriors_unlabeled$comp_means[[idx]] <-lapply(1:ncol(cumcat_and_dif),function(idx){
#              mixture_obj@current.meanmat[(cumcat_and_dif[1,idx]-cumcat_and_dif[2,idx]+1):cumcat_and_dif[1,idx],]
#            })
#            bmbclust_obj@posteriors$partition[[idx]] <- N
#            bmbclust_obj@posteriors$loglik[idx] <- loglik
#            bmbclust_obj@posteriors$alloc_vectors[[idx]] <- S
#            bmbclust_obj@posteriors$alloc_probs[[idx]] <- alloc_probs
#            return(bmbclust_obj)
#          })
#
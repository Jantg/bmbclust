# useful functions
# BNB distri
dbnb <- function(x, alpha_l, a, b, log=FALSE) {
  val <- lgamma(alpha_l + x) + lbeta(alpha_l + a, x + b) - lgamma(alpha_l) - lgamma(x + 1) - lbeta(a, b)
  if (log) return(val)
  else return(exp(val))
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

array2blkdiag<- function(array){
  if(length(dim(array))!=3){stop("wrong dimension")}
  k<- dim(array)[1]
  l<- dim(array)[2]
  mat <- matrix(0,nrow = k*l,ncol = k*l)
  for(i in 1:k){
    mat[((i-1)*l+1):(i*l),((i-1)*l+1):(i*l)] = array[i,,]
  }
  return(mat)
}

distriname2func <- function(name,quantity = "d"){
  if(quantity == "d"){
    if(name == "gaussian") return(expression(stats::dnorm))
    else if(name == "Gamma") return(expression(stats::dgamma))
    else if(name == "poisson") return(expression(stats::dpois))
    else if(name == "F") return(expression(stats::df))
    else if(name == "bnb") return(expression(dbnb))
    else if(name == "geom") return(expression(stats::dgeom))
    else{
      stop("no such distribution is currently available")
    }
    
  }else if(quantity == "r"){
   if(name == "gaussian") return(expression(stats::rnorm))
   else if(name == "Gamma") return(expression(stats::rgamma))
   else if(name == "poisson") return(expression(stats::rpois))
   else if(name == "F") return(expression(stats::rf))
   else if(name == "geom") return(expression(stats::rgeom))
    else{
      stop("no such distribution is currently available")
    }
  }else{
    stop("quatity must be either 'd' or 'r'")
  }
}
#All Generics

setGeneric("bmbclust",
           function(data,model,control=NULL,init=NULL,parallel=FALSE)
             standardGeneric("bmbclust"))
setGeneric("Init_Sampler",
           function(model,data) standardGeneric("Init_Sampler"))
setGeneric("Component_update",
           function(data,model,...) standardGeneric("Component_update"))
setGeneric("Index_update",
           function(data,model,...) standardGeneric("Index_update"))
setGeneric("Tel_Sampler",
           function(model,data,control,init) standardGeneric("Tel_Sampler"))
setGeneric("eval_partition_loglik",
           function(data,model,...) standardGeneric("eval_partition_loglik"))
setGeneric("eval_pdf",
           function(x,d,p1,p2,p3,...) standardGeneric("eval_pdf"))
setGeneric("sensitivity_analysis",
           function(result,...) standardGeneric("sensitivity_analysis"))
setGeneric("identification",
           function(result,...) standarGeneric("identification"))

# Virtual Class
setClass("BMBC",
         slots = c(
           n.iter = "numeric",
           n.burnin = "numeric",
           n.thin = "numeric",
           progress.bar = "character",
           verbose = "logical"
         ),
         prototype = list(
           n.iter = 100000,
           n.burnin = 10000,
           n.thin = 1,
           progress.bar = "both",
           verbose = TRUE
         ),
         validity = function(object){
           flag = all(c(object@n.iter,object@n.burnin,object@n.thin)>0)
           if(!flag) return ("Positive numeric values have to be specified for n.* arguments")
           else return(flag)
         },
         contains = "VIRTUAL"
)
# Mixture class (Mixtures of Mixtures will be a different one)
setClass("BMBCMx",
         slots = c(component.density = "character",
                   prior.component = "list",
                   hyperprior = "list",
                   prior.alpha = "list",
                   prior.k = "list",
                   init.k = "numeric",
                   init.alpha = "numeric",
                   init.compmean = "array",
                   init.compcov = "array",
                   max.k = "numeric",
                   s_alpha = "numeric"
         ),
         prototype = list(
           hyperprior = list("GIG",nu1=0.4,nu2=0.4),
           prior.alpha = list("F",6,5),
           prior.k = list("bnb",6,4,3),
           init.k = 10,
           init.alpha = 1.0,
           max.k = 150,
           s_alpha = 2.0
         ),
         validity = function(object)object@component.density %in% c("gaussian","poisson")
         ,
         contains = "BMBC"
)

# Mixture of gaussian class
setClass("BMBCMx_gauss",
         slots = c(priors="list"),
         contains = "BMBCMx")

setClass("post_identification",
         slots = c(parameters = "list",num_samples = "numeric",K_mode = "numeric"))

setClass("bmbclust",
         slots = c(posteriors = "list",
                   priors = "list",
                   control.params = "list"),
)
# The main function bmbclust
setMethod("bmbclust",
          signature(model = "list"),
          function(data = list(),model = list(),
                   control = NULL,init = NULL,parallel = FALSE){
            across_model_args = append(control,init)
            if(!is.null(across_model_args)){
              for(i in names(across_model_args)){
                assign(i,do.call("$",list(across_model_args,i)))
              }
              model <- lapply(model,
                              function(m){
                                for(i in setdiff(slotNames(m),names(across_model_args))){
                                  assign(i,do.call("@",list(m,i)))
                                }
                                new("BMBCMx_gauss",component.density = component.density,
                                    prior.component = prior.component,
                                    hyperprior = hyperprior,prior.alpha = prior.alpha,
                                    prior.k = prior.k,init.k = init.k,init.alpha = init.alpha,
                                    init.compmean = init.compmean,init.compcov = init.compcov,
                                    n.iter = n.iter, n.burnin =n.burnin, n.thin = n.thin,
                                    max.k = max.k, s_alpha = s_alpha,
                                    progress.bar = progress.bar,verbose = verbose)
                              })
            }else{
              model <- lapply(model,
                              function(m){
                                new("BMBCMx_gauss",component.density = m@component.density,
                                    prior.component =  m@prior.component,
                                    hyperprior = m@hyperprior,prior.alpha = m@prior.alpha,
                                    prior.k = m@prior.k,init.k = m@init.k,init.alpha = m@init.alpha,
                                    init.compmean = m@init.compmean,init.compcov = m@init.compcov,
                                    n.iter = m@n.iter, n.burnin =m@n.burnin, n.thin = m@n.thin,
                                    max.k = m@max.k, s_alpha = m@s_alpha,
                                    progress.bar = m@progress.bar,verbose = m@verbose)
                              }
              )
            }
            params_list <- lapply(model,function(m) sapply(setdiff(slotNames(m),c("prior.k","prior.alpha")),function(x) do.call("@",list(m,x))))
            prior_list <-lapply(model,function(m) sapply(c("prior.k","prior.alpha"),function(x) do.call("@",list(m,x))))
            if(parallel){res_list <- parallel::mclapply(model,function(m) Tel_Sampler(m,data),mc.cores = min(c(parallel::detectCores(),length(model))))}
            else{res_list <- lapply(model,function(m) Tel_Sampler(m,data))}
            return_obj <- new("bmbclust",posteriors = res_list,priors = prior_list,control.params = params_list)
            return(return_obj)
          }
)

setMethod("Init_Sampler",
          signature(model = "BMBCMx_gauss"),
          function(model,data){
            r <- ncol(data)
            K_0 <- model@init.k
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
            priors <- list(c0=c0,g0=g0,G0=G0,C0_j=C0_j,B0=B0,inv_B0=inv_B0,b0=b0,R=R)
            return(priors)
          })

setMethod("Component_update",
          signature(model = "BMBCMx_gauss"),
          function(data,model,r,K_plus,N,S,meanmat,covmat){
            c0 <- model@priors$c0
            C0_j <- model@priors$C0_j
            g0 <- model@priors$g0
            G0 <- model@priors$G0
            B0 <- model@priors$B0
            R <- model@priors$R
            inv_B0 <- model@priors$inv_B0
            b0 <- model@priors$b0
            v1 <- model@hyperprior$nu1
            v2 <- model@hyperprior$nu2
            ck <- c0 + N/2
            Ck <- lapply(1:K_plus, function(k)
              C0_j + 0.5 * crossprod(sweep(data[S == k, , drop = FALSE], 2, meanmat[k, , ], FUN = "-")))
            
            sigs <- tryCatch(
              expr = {
                lapply(1:K_plus, function(k) bayesm::rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[[k]]))))
              },
              error = function(e){
                lapply(1:K_plus, function(k) bayesm::rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[[k]]+diag(runif(dim(Ck[[k]])[1])*0.01)))))
              }
            )
            #sigs <- lapply(1:K_plus, function(k) bayesm::rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[[k]]))))
            inv_sig <- lapply(1:K_plus, function(k) sigs[[k]]$W)
            sig <- lapply(1:K_plus, function(k) sigs[[k]]$IW)

            updated_vals <- lapply(1:K_plus, function(k) {
              Bk <- tryCatch(
                expr = {
                  chol2inv(chol(inv_B0 + N[k] * inv_sig[[k]]))
                },
                error = function(e){
                  chol2inv(chol(inv_B0 + N[k] * inv_sig[[k]]+diag(runif(dim(inv_sig[[k]])[1])*0.01)))
                }
              )    
              #Bk <- chol2inv(chol(inv_B0 + N[k] * inv_sig[[k]]))
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
            lambdas <- sapply(bj, function(x) GIGrvg::rgig(1, pk, x, aj))
            #browser()
            #if(!is.numeric(R)) browser()
            B0 <- diag(lambdas * as.vector(R)^2/K_plus, nrow = r)
            inv_B0 <- diag(1/diag(B0), nrow = r)
            b0 <- crossprod(chol(B0),rnorm(r)) + Reduce("+", mu)/K_plus

            priors <- list(c0=c0,g0=g0,G0=G0,C0_j=C0_j,B0=B0,inv_B0=inv_B0,b0=b0,lambdas = lambdas,R=R)
            model@priors <- priors
            return_params <- list(mu = mu,sig = sig,bk = bk,Bk=Bk,model = model)
            return(return_params)
          })

setMethod("Index_update",
          signature(model="BMBCMx_gauss"),
          function(data,model,r,K,meanmat,covmat,etas){
            #model is for mixture of mixture case (not implemented yet)
            #long_meanmat <- Reduce(function(l,r) cbind(l,r),lapply(1:K,function(x) matrix(meanmat[x,,],ncol = r)))
            #block_covmat <- array2blkdiag(covmat)
            #lik <- dMvn_multi(data,long_meanmat,block_covmat,K,log=TRUE)
            #pos_prob <- sweep(lik,2,etas,"+")
        
           eta_prod_dens <- simplify2array(
             lapply(1:K,function(k){
               dMvn(data,matrix(meanmat[k,,],ncol = r),matrix(covmat[k,,],ncol =r),log = TRUE)+
                 etas[k]}#,mc.cores = 8
             )
           )
           gumbel_noise <- array(evd::rgumbel(K*dim(data)[1]),dim = c(dim(data)[1],1,K))
           S <- apply(eta_prod_dens+gumbel_noise,1,function(x) which.max(x))
           #gumbel_noise <- matrix(evd::rgumbel(K*dim(data)[1]),ncol = K)
           #S <- apply(pos_prob+gumbel_noise,1,function(x) which.max(x))
            return(S)
          })

setMethod("eval_partition_loglik",
          signature(model = "BMBCMx_gauss"),
          function(data,model,r,K,meanmat,covmat,S){
            loglik <- simplify2array(lapply(1:nrow(data),function(idx){
              dMvn(matrix(data[idx,],ncol = r),matrix(meanmat[S[idx],,],ncol = r),matrix(covmat[S[idx],,],ncol = r),log = TRUE)
            }))
            return(sum(loglik))
          })

setMethod("Tel_Sampler",signature(model = "BMBCMx"),
          function(model,data,control,init){
            data <- data.matrix(data)
            niter <- model@n.iter
            nburnin <- model@n.burnin
            nthin <- model@n.thin
            r <- ncol(data)
            K_0 <- model@init.k
            alpha <- model@init.alpha
            meanmat <- model@init.compmean
            covmat <- model@init.compcov
            priors <- Init_Sampler(model,data)
            model@priors <- priors
            etas <- lrdir(rep(1/K_0, K_0))
            K <- K_0
            K_max <- model@max.k
            s_alpha <- model@s_alpha
            verbose <- model@verbose
            if(verbose) first_time <- TRUE
            max_kplus <- 0
            v1 = 0.4;v2 = 0.4
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

            res_lst <- list(eta = list(), Ks = numeric(niter), Kps = numeric(niter),
                            alphas = numeric(niter),comp_means = list(),
                            comp_covs = list(), bs = list(), Bs = list(),max_kplus = 0,
                            partition = list(), MAP = list(),columns = "character",loglik = numeric(niter)
            )
            for (n in 1:(niter + nburnin)) {
              #if(n>nburnin+11){browser()}
              if(n %% 2 == 0){print(c(n,K_plus))}
              S <- Index_update(data,model,r,K,meanmat,covmat,etas)

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

              ret_vals <- Component_update(data,model,r,K_plus,N,S,meanmat,covmat)
              mu <- ret_vals$mu
              sig <- ret_vals$sig
              bk <- ret_vals$bk
              Bk <- ret_vals$Bk
              model <- ret_vals$model
              C0_j <- model@priors$C0_j
              B0 <- model@priors$B0
              b0 <- model@priors$b0
              lambdas <- model@priors$lambdas
              c0 <- model@priors$c0
              # compute unnormalized log prob to do gumbel max trick again to pick K
              unnorm_logprob <- sapply(K_plus:K_max, function(k){
                do.call("eval_pdf",
                        append(append(append(list(x=k-1),list(distriname2func(model@prior.k[[1]]))),
                                      model@prior.k[2:length(model@prior.k)]),list(log=TRUE)))+
                  lfactorial(k) - lfactorial(k - K_plus)+sum(log(k/(N * k + alpha)) + lgamma(1 + N + alpha/k) - log(k/alpha) - lgamma(1 + alpha/k))}
              )

              K <- which.max(evd::rgumbel(length(unnorm_logprob)) + unnorm_logprob) + K_plus - 1

              # sample alpha with metropolis random walk (on log space so with Jacobian corrections)
              alpha_prop <- exp(rnorm(1, log(alpha), s_alpha))
              tmp <- sapply(c(alpha_prop, alpha),function(a){
                do.call("eval_pdf",
                        append(append(append(list(x=a),list(distriname2func(model@prior.alpha[[1]]))),
                                      model@prior.alpha[2:length(model@prior.alpha)]),list(log=TRUE)))+
                  lgamma(a) - lgamma(sum(N) + a)+
                  sum(log(K/(N * K + a)) + lgamma(1 + N + a/K) - log(K/a) - lgamma(1 + a/K))}
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
              loglik <- eval_partition_loglik(data,model,r,K,meanmat,covmat,S)
              if (n > nburnin) {
                if(K_plus>max_kplus) max_kplus <- K_plus
                idx <- n - nburnin
                res_lst$eta[[idx]] <- etas
                res_lst$Ks[idx] <- K
                res_lst$Kps[idx] <- K_plus
                res_lst$alphas[idx] <- alpha
                res_lst$comp_means[[idx]] <- meanmat
                res_lst$comp_covs[[idx]] <- covmat
                res_lst$bs[[idx]] <- bk
                res_lst$Bs[[idx]] <- Bk
                res_lst$partition[[idx]] <- N
                res_lst$loglik[idx] <- loglik
                if(verbose){
                  if(first_time){
                    x11()
                    temporal_partitions <- matrix(prop.table(sort(N,decreasing = TRUE)),nrow = 1)
                    first_time <- FALSE
                    par(mfrow=c(3,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
                    layout(matrix(c(1,5,3,4,2,2),3,2,byrow = TRUE),heights = c(2.5,2.5,3.0))
                  }else{
                    #par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
                    #par(mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
                    #layout(matrix(c(2,1,3,4,5,5),3,2,byrow = TRUE),heights = c(2.5,2.5,3.0))
                    if(length(N) > ncol(temporal_partitions)){
                      diff_ncols <- length(N) - ncol(temporal_partitions)
                      add_col <- matrix(0.0,nrow=nrow(temporal_partitions),ncol = diff_ncols)
                      temporal_partitions<- cbind(temporal_partitions,add_col)
                    }else if(length(N) < ncol(temporal_partitions)){
                      diff_ncols <- ncol(temporal_partitions) - length(N)
                      N <- c(N,rep(0,diff_ncols))
                    }
                    temporal_partitions <- rbind(temporal_partitions,prop.table(matrix(sort(N,decreasing = TRUE),nrow = 1)))
                  }
                  dev.hold()
                  # traceplot for alpha/K
                  # col = 100 is blue
                  matplot(log(res_lst$alphas[1:idx]/res_lst$Ks[1:idx]), type = "l", lty = 1,col = 100, ylab = "log(alpha/K)")
                  abline(h = 0.0,col = "black")
                  matplot(res_lst$loglik[1:idx],type = "l",lty = 1, col = 100,ylab = "log-likelihood")
                  abline(h = -50.0,col = "black")
                  # histogram of K_plus
                  hist(res_lst$Kps[1:idx],freq = FALSE,main = "Histogram of K_plus")
                  # barplot of partitions from niter-10 to niter (older ones are discarded for storage purposes)
                  barplot(apply(t(temporal_partitions),2,rev))
                  if(n-nburnin>10){
                    temporal_partitions <- temporal_partitions[-1,]
                  }
                  # histogram of K_plus/K
                  hist(res_lst$Kps[1:idx]/res_lst$Ks[1:idx],main = "Histogram of K_plus/K",freq = FALSE)
                  dev.flush()
                }
              }

            }
            par(mfrow =c(1,1))
            res_lst$max_kplus <- temporal_partitions
            res_lst$columns <- colnames(data)
      
            #K_mode <- as.numeric(names(which.max(table(res_lst$Kps))))
            #res_df <- Reduce(function(l,r)
            #  rbind(l,r) ,lapply(which(res_lst$Kps == K_mode),function(x)
            #    res_lst$comp_means[[x]][,,1][1:K_mode,]))
            #res_df <- scale(res_df) 
            #res_center <-simplify2array(lapply(which(res_lst$Kps == K_mode),function(x)
            #  Reduce(function(l,r) rbind(l,r),lapply(res_lst$bs[[x]][1:K_mode],t))))
            #clst_idx <- kmeans(res_df,res_center[,,dim(res_center)[3]])$cluster
            #perm_or_not <- sapply(1:sum(res_lst$Kps==K_mode),function(x)
            #  all(sort(clst_idx[((x-1)*K_mode+1):(x*K_mode)])== 1:K_mode))
            #res_ordered = lapply(idx,function(i)
            #  res_df[((i-1)*K_mode+1):(K_mode*i),][clst_idx[((i-1)*K_mode+1):(K_mode*i)],])
            #res_lst$MAP = res_ordered
            return(res_lst)
          })




BMBCMx <- function(data,component.density = c("gaussian","poisson"),
                   hyperprior = NULL,prior.component=NULL,prior.alpha=NULL,prior.k=NULL,
                   init.k= NULL, init.alpha = NULL,n.iter = NULL
                   ,n.burnin = NULL,n.thin = NULL,init.compmean=NULL,
                   init.compcov=NULL,max.k=NULL,s_alpha=NULL,progress.bar = NULL){
  r <- ncol(data)
  component.density <- match.arg(component.density)
  mx_obj <- new("BMBCMx",component.density = component.density)
  if (is.null(prior.component)){
    if(component.density == "gaussian"){
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
      if(is.null(hyperprior)){
        lambdas <- rep(1,r)
        prior.component<- list(c0=c0,g0=g0,G0=G0,C0_j=C0_j,B0=B0,
                               b0=b0,lambdas =lambdas)

      }else{
        prior.component<- list(c0=c0,g0=g0,G0=G0,C0_j=C0_j,B0=B0,
                               b0=b0)
      }
    }
    else if(component.density == "poisson") prior.component <- list("poisson",1)
  }
  for(i in names(which(sapply(formals(BMBCMx),is.null)))){
    if(do.call("is.null",list(as.name(i)))){
      do.call("<-",list(i,do.call("@",list(mx_obj,i))))
    }
  }
  if(length(init.compmean)==0) init.compmean <- array(runif(r* init.k), dim = c(init.k, r, 1))
  if(length(init.compcov)==0){
    init.compcov <- array(runif(r * r * init.k), dim = c(init.k, r, r))
    init.compcov <- array(t(apply(init.compcov, 1, function(x) crossprod(x))), dim = c(init.k, r, r))
  }
  obj_return <- new("BMBCMx",component.density = component.density,
                    prior.component = prior.component, prior.alpha = prior.alpha,
                    prior.k = prior.k,n.iter = n.iter, n.burnin = n.burnin,
                    n.thin = n.thin,init.compmean = init.compmean,
                    init.compcov = init.compcov,max.k=max.k,
                    s_alpha=s_alpha,progress.bar = progress.bar)
  return(obj_return)
}


setMethod("eval_pdf",
          signature(d="expression",p1 = "numeric",p2 = "missing",p3 = "missing"),
          function(x,d,p1,log=FALSE){
            print(class(d))
            do.call(eval(d),list(x=x,p1,log=log))
          })
setMethod("eval_pdf",
          signature(d="expression",p1 = "numeric",p2 = "numeric",p3 = "missing"),
          function(x,d,p1,p2,log=FALSE){
            do.call(eval(d),list(x=x,p1,p2,log=log))
          })
setMethod("eval_pdf",
          signature(d="expression",p1="numeric",p2="numeric",p3="numeric"),
          function(x,d,p1,p2,p3,log=FALSE){
            do.call(eval(d),list(x=x,p1,p2,p3,log=log))
          })
setMethod("eval_pdf",
          signature(d="expression",p1="matrix",p2="matrix"),
          function(x,d,p1,p2,log=FALSE){
            do.call(eval(d),list(x=x,p1,p2,log=log))
          })

setMethod("sensitivity_analysis",
          signature(result="bmbclust"),
          function(result,model.idx=NULL){
            if(is.null(model.idx)){
              model.idx = 1:length(result@posteriors)
            }
            max_kplus = max(sapply(model.idx,function(m)result@posteriors[[m]]$max_kplus))
            n_models <- length(model.idx)
            bins <- seq(0.5, max_kplus+0.5, by=1)
            if(n_models>8){stop("Only 8>= models can be compared")}
            col_idx <- sapply(0:7,function(x){ as.integer(intToBits(x))})[1:3,]
            #col_idx <- col_idx[,c(7,4,1,6,3,0,5,2)]
            col_idx <- col_idx[,7:1]
            cols = vector(mode="character", length=n_models)
            for(i in 1:n_models){
              if(i == 1){
                h<- hist(result@posteriors[[model.idx[i]]]$Kps,plot = FALSE,breaks = bins)
                h$counts <- h$counts/sum(h$counts)
                col <- rgb(col_idx[1,i],col_idx[2,i],col_idx[3,i],1/4,max = 1)
                cols[i] <- col
                #browser()
                plot(h, col=col,xlim = c(1,max_kplus),
                     ylim=c(0,1),ylab = "prob",main ="Histogram of posterior K_plus across models",xlab = "K_plus")
              }else{
                h<- hist(result@posteriors[[model.idx[i]]]$Kps,plot = FALSE,breaks = bins)
                h$counts <- h$counts/sum(h$counts)
                col <- rgb(col_idx[1,i],col_idx[2,i],col_idx[3,i],1/4,max = 1)
                cols[[i]] <- col
                plot(h,col=col,add= TRUE)
              }
            }
            legend('topright',sapply(model.idx,function(m) paste(gsub("list","",paste("K~",result@priors[[m]])[1]),
                                                                 gsub("list","",paste("Alpha~",result@priors[[m]])[2]))),
                   fill = cols, bty = 'n',
                   border = NA)
          })

setMethod("identification",
          signature(result="bmbclust"),
          function(result){
            ret_vals <- lapply(result@posteriors,function(res_lst){
              #MAP <- lapply(result,function(res_lst){
              K_mode <- as.numeric(names(which.max(table(res_lst$Kps))))
              res_df <- Reduce(function(l,r)
                rbind(l,r) ,lapply(which(res_lst$Kps == K_mode),function(x){
                  matrix(res_lst$comp_means[[x]][,,1],nrow = res_lst$Ks[[x]])[1:K_mode,]}))
              res_df <- scale(res_df)
              res_center <-simplify2array(lapply(which(res_lst$Kps == K_mode),function(x)
                Reduce(function(l,r) rbind(l,r),lapply(res_lst$bs[[x]][1:K_mode],t))))
              res_center_last <- sweep(sweep(matrix(res_center[,,dim(res_center)[3]],nrow = K_mode),2,attributes(res_df)$'scaled:center',"-"),
                                       2,attributes(res_df)$'scaled:scale',"/")
              clst_idx <- kmeans(res_df,res_center_last)$cluster
              perm_or_not <- sapply(1:sum(res_lst$Kps==K_mode),function(x)
                all(sort(clst_idx[((x-1)*K_mode+1):(x*K_mode)])== 1:K_mode))
              idx = which(perm_or_not)
              res_ordered = lapply(idx,function(i){
                sweep(sweep(matrix(res_df[((i-1)*K_mode+1):(K_mode*i),,drop=FALSE],nrow=K_mode)[clst_idx[((i-1)*K_mode+1):(K_mode*i)],,drop=FALSE]
                ,2,attributes(res_df)$'scaled:scale',"*"),2,attributes(res_df)$'scaled:center',"+")})
              #print(res_ordered)
              MAP = lapply(res_ordered,function(x){
                colnames(x) <- res_lst$columns
                x
              })
              #browser()
              return(MAP)
              
            })
            ret_obj <- new("post_identification",parameters = ret_vals,num_samples = sapply(ret_vals,length),
                           K_mode = sapply(ret_vals,function(x) nrow(x[[1]])))
            return(ret_obj)
          })
setMethod("summary",
          signature(object="post_identification"),
          function(object){
            ret_output <-lapply(object@parameters,function(model){
              K_mode <- nrow(model[[1]])
              ret <- lapply(1:K_mode,function(i){
                pos_mat <- simplify2array(model)[i,,,drop=FALSE]
                pos_mat <- t(as.matrix(pos_mat[1,,]))
                apply(pos_mat,2,function(x){
                  quantile(x,probs = c(.025,.5,.975))
                })
              })
            })
            invisible(sapply(1:length(ret_output),function(model_idx){
              cat("Model ",model_idx,"\n")
              cat("\n")
              sapply(1:length(ret_output[[model_idx]]),function(cluster_idx){
                cat("Cluster ",cluster_idx,"\n")
                print(ret_output[[model_idx]][[cluster_idx]])
                cat("\n")
              })
            }))
          })



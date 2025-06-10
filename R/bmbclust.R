# Useful functions and the main bmbclust function


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


array2blkdiag <- function(array) {
  if (length(dim(array)) != 3) {
    stop("wrong dimension")
  }
  k <- dim(array)[1]
  l <- dim(array)[2]
  mat <- matrix(0, nrow = k*l, ncol = k*l)
  for(i in 1:k){
    mat[((i-1)*l+1):(i*l),((i-1)*l+1):(i*l)] = array[i,,]
  }
  return(mat)
}

logTopdotexp <- function(ln_x, ln_y) {
  ln_x_max <- max(ln_x)
  ln_y_max <- max(ln_y)
  Top <- toeplitz(exp(ln_x - ln_x_max))
  Top[lower.tri(Top)] <- 0
  Top <- cbind(matrix(0, nrow = nrow(Top)), Top)
  logTop <- log(Top %*% exp(ln_y - ln_y_max)) + ln_y_max + ln_x_max
  return(logTop)
}

n_partitions_Nk <- function(k, N, gamma_k) {
  init_vec <- sapply(N:1, function(n) lgamma(n + gamma_k) - lgamma(n + 1))
  res <- list(matrix(init_vec, ncol = 1))
  c_K_k <- function(k, N, gamma_k){
    if(k == 1){
      return(matrix(init_vec))
    }
    if(length(res) < k){ 
      res <<- "length<-" (res,k)
    }
    if(!is.null(res[[k]])){ 
      return(res[[k]])
    }
    ws <- sapply(1:(N-k+1),function(n) lgamma(n+gamma_k) - lgamma(n+1))
    res[[k]]<<- logTopdotexp(ws, c_K_k(k-1, N, gamma_k))
    res[[k]]
  }
  c_K_k(k, N, gamma_k)
}


# component_density -> either string or user provided function (fitdistr from MASS)
defineComponentObj <- function(type, component_density = NULL,component_density_func = NULL,update_rule = NULL,
                               update_rule_params = NULL,update_rule_hyperparams = NULL){
  if(is.null(component_density)){
    component_density <- switch(type,
                                real_valued = "gaussian",
                                count = "poisson",
                                categorical = "categorical",
                                stop(
                                  paste("data type", type,"is not supported")
                                  )
                                )
  }
  # need to check the validity of component density if user provided it
  # perhaps inside bmbclust() function
  call_list <- as.list(match.call())
  names_matched <- !(names(call_list) %in% c("","type"))
  objargs_list <- call_list[names_matched]
  component_obj <- switch(component_density,
                          gaussian = {
                            tmp_args <- append(append(list("BMBCMx_gauss"),objargs_list),
                                               list(component_density = component_density))
                            do.call("new",tmp_args)
                          },
                          poisson = {
                            tmp_args <- append(append(list("BMBCMx_pois"),objargs_list),
                                               list(component_density = component_density))
                            do.call("new",tmp_args)
                          },
                          categorical = {
                              tmp_args <- append(append(list("BMBCMx_categ"),objargs_list),
                                                 list(component_density = component_density))
                              do.call("new",tmp_args)
                          },
                          stop(paste("No component density named", component_density," found"))
                          )
  return(component_obj)
}

## 
defineWeightObj <- function(type,prior_K = NULL, prior_K_params = NULL, prior_alpha = NULL,
                            prior_alpha_params = NULL, sig_alpha = NULL ,prior_gamma = NULL,
                            prior_gamma_params = NULL,init_K = NULL,init_alpha = NULL
                            ,init_gamma = NULL, max_K = NULL,fix_gamma = NULL){
  # three speicifications for type: unknown_K, known_K,fixed_K
  # should have some sanity check function that gives warning message when 
  # combination of user-specified priors and that of type does not match
  # e.g. fixed_K but alpha not small enough to be sparse mixtures
  call_list <- as.list(match.call()) 
  names_matched <- !(names(call_list) %in% c("","type","init_K","init_alpha","init_gamma"))
  objargs_list <- call_list[names_matched]
  if(!is.null(init_K)){
    objargs_list$K <- init_K
  }
  if(!is.null(init_alpha)){
    objargs_list$alpha <- init_alpha
  }
  if(!is.null(init_gamma)){
    objargs_list$gamma <- init_gamma
  }
  weight_obj <- switch(type,
                       unknown_K = {
                         if(!is.null(init_alpha)){
                           objargs_list$alpha <- init_alpha
                         }
                         tmp_args <- append(list("Weights_dynamic"),objargs_list)
                         do.call("new",tmp_args)
                       },
                       fixed_K ={
                         tmp_args <- append(list("Weights_SparseFM"),objargs_list)
                         do.call("new",tmp_args)
                       },
                       known_K = {
                         tmp_args <- append(list("Mixture_Weights"),objargs_list)
                         do.call("new",tmp_args)
                       })
  return(weight_obj)
}


bmbclust <- function(data,weight_obj,component_obj) {
  if(!is.matrix(data)&&!is.data.frame(data)) stop("data should be of type matrix or data.frame")
  # rather condition on the type of the data (matrix or)
  if(is.data.frame(data)){
    data_types <- unique(sapply(data, class))
    is_mixed <- length(data_types)!=1
    if(is_mixed){
      stop(paste0("mixed type data (",
           paste(data_types,collapse = ", "),") are not supported"))
    }
  }
  tmp <- preprocess_data(component_obj,data)
  data_mat <- tmp[[1]]
  columns <- tmp[[2]]
  return_obj <- new("bmbclust",weight_obj = weight_obj, component_obj = component_obj,
                    data = data_mat,columns = columns)

  return(return_obj)
}


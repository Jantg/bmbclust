# All classes

# BNB distri
dbnb <- function(x, mean, shape1, shape2, log = FALSE) {
  val <- lgamma(mean + x) + lbeta(mean + shape1, x + shape2) - lgamma(mean) - lgamma(x + 1) - lbeta(shape1, shape2)
  if (log){
    return(val)
  }else{
      return(exp(val))
  }
}

# Virtual class for the mixture object within the bmbclust class
setClass("BMBC",
         slots = c(
           init.k = "numeric",
           init.alpha = "numeric",
           sig_alpha = "numeric"
         ),
         prototype = list(
           init.k = 30,
           init.alpha = 1.0,
           sig_alpha = 2.0
         ),
         validity = function(object){
         },
         contains = "VIRTUAL"
)

setClass("Weights",
         slots = c(
           K = "numeric",
           gamma = "numeric",
           K_plus = "numeric",
           cluster_sizes = "vector",# N
           cluster_allocations = "vector", # S
           allocation_probs = "array",
           weights = "vector",
           accepted = "logical"
         ),prototype = list(
           accepted = TRUE
         ),
         contains = "VIRTUAL")

setClass("Weights_dynamic",
         slots = c(
           alpha = "numeric",
           prior_alpha = "function",
           prior_alpha_params = "list",
           sig_alpha = "numeric",
           prior_K = "function",
           prior_K_params = "list",
           max_K = "numeric"
         ),prototype = list(
           prior_alpha = stats::df,
           prior_alpha_params = list(df1 =6.0, df2 = 3.0),
           sig_alpha = 2.0,
           prior_K = dbnb,
           prior_K_params = list(mean =4.0, shape1 = 4.0, shape2 = 3.0),
           max_K = 150
         ),contains = "Weights")

setClass("Weights_SparseFM",
         slots = c(
           prior_gamma = "function",
           prior_gamma_params = "list",
           fix_gamma = "logical",
           sig_gamma = "numeric"
         ),prototype = list(
           sig_gamma = 2.0,
           prior_gamma = dgamma,
           prior_gamma_params = list(shape = 1.0, rate= 2.0),
           fix_gamma = FALSE
         ),contains = "Weights")

# The main class 
setClass("bmbclust",
         slots = c(
           data = "matrix",
           columns = "list",
           weight_obj = "Weights",
           component_obj = "BMBC",# later rename this as well
           posteriors_unlabeled = "list",
           posteriors_labeled = "list"
           #data.type = "character", # keep it with mixture.obj
           #component.density = "character",
           #mixture.obj = "BMBC", # component object
           #posteriors = "list",
           #post.identification = "list",
           #inits = "list",
           # combine those to a new class?
           # for the weight model
           #  constructor function()?
           #max.k = "numeric",
           #fix.k = "logical",
           #prior.k = "function",
           #prior.alpha = "function",
           #prior.k.parameters = "list",
           #prior.alpha.parameters = "list"
         ),
         prototype = list(
           #data.type = "real-valued",
           #component.density = "gaussian",
           #prior.k = dbnb,
           #prior.alpha = stats::df,
           #prior.k.parameters = list(mean =4.0, shape1 = 4.0, shape2 = 3.0),
           #prior.alpha.parameters = list(df1 =6.0, df2 = 3.0),
           #max.k = 150,
           #fix.k = FALSE
         ))

# Subclass of BMBC for mixtures (another subclass for mixtures of mixtures will be implemented soonish)
setClass("BMBCMx",
         slots = c(
           hyperprior = "logical"
         ),
         prototype = list(
           hyperprior = TRUE
         ),
         #validity = function(object){
         #},
         contains = "BMBC"
)

setClassUnion("matrixOrdata.frame", c("matrix", "data.frame"))


# Subclass of BMBCMx with Gaussian component distribution for real-valued data
setClass("BMBCMx_gauss",
         slots = c(
           # separate parameters that are fixed and those that change
           # also separate filled and unfilled parameters
           component_density = "character",
           component_density_func = "function",
           update_rule = "function",
           update_rule_params = "list",
           update_rule_hyperparams = "list"
         ),
         prototype = list(
         ),
         validity = function(object){
           #flag = all(simplify2array(lapply(as.data.frame(object@data),is.numeric)))
           #if(!flag) return("columns of data must be all numeric type")
           #else return(flag)
         },
         contains = "BMBCMx"
)

# Subclass of BMBCMx with Poisson component distribution for count data
setClass("BMBCMx_pois",
         slots = c(
           # separate parameters that are fixed and those that change
           # also separate filled and unfilled parameters
           component_density = "character",
           component_density_func = "function",
           update_rule = "function",
           update_rule_params = "list",
           update_rule_hyperparams = "list"
         ),
         prototype = list(),
         validity = function(object){
           #flag = all(simplify2array(lapply(as.data.frame(object@data),is.numeric)))
           #if(!flag) return("columns of data must be all numeric type")
           #else{
           #  flag2 = all(simplify2array(lapply(as.data.frame(object@data),function(x) all(x == floor(x)))))
           #  if(!flag2) return("not a count data, all values should be integers")
           #  else return(flag2)
           #}
         },
         contains = "BMBCMx"
)

# Subclass of BMBCMx with categorical component distribution for categorical data
setClass("BMBCMx_categ",
         slots = c(
           # mat, perhaps name it as data_matrix, w/ some function we can 
           # return the original data
           # separate parameters that are fixed and those that change
           # also separate filled and unfilled parameters
           component_density = "character",
           component_density_func = "function",
           update_rule = "function",
           update_rule_params = "list",
           update_rule_hyperparams = "list"
         ),
         prototype = list(),
         validity = function(object){
           # length of the collevels to be equal to ncol of y
           #flag = all(simplify2array(lapply(as.data.frame(object@data), is.factor)))
           #if(!flag) return("not a categorical data, columns must be factors")
           #else return(flag)
         },
         contains = "BMBCMx"
)


# All generics

setGeneric("preprocess_data",
           function(component_obj,data) standardGeneric("preprocess_data"))
setGeneric("simulate",
           function(bmbclust_obj,...) standardGeneric("simulate"))

setGeneric("Telescope_sampler",
           function(bmbclust_obj,...) standardGeneric("Telescope_sampler"))

setGeneric("compute_prior_Kplus",
           function(bmbclust_obj,...) standardGeneric("compute_prior_Kplus"))

setGeneric("Init_component",
           function(component_obj,bmbclust_obj,...) standardGeneric("Init_component"), signature = c("component_obj","bmbclust_obj"))

setGeneric("Init_latent_params",
           function(component_obj,weight_obj,bmbclust_obj,...) standardGeneric("Init_latent_params"),signature = c("component_obj","bmbclust_obj"))

setGeneric("Init_post_list",
           function(bmbclust_obj,...) standardGeneric("Init_post_list"))

setGeneric("Index_update",
           function(bmbclust_obj,...) standardGeneric("Index_update"))

setGeneric("Reorder_components",
           function(component_obj,bmbclust_obj,mapping,...) standardGeneric("Reorder_components"),signature = c("component_obj","bmbclust_obj"))

setGeneric("Component_update",
           function(bmbclust_obj,...) standardGeneric("Component_update"))

setGeneric("Sample_component_params",
           function(component_obj,bmbclust_obj,...) standardGeneric("Sample_component_params"),signature = c("component_obj","bmbclust_obj"))

setGeneric("update_K",
           function(weight_obj,...) standardGeneric("update_K"))

setGeneric("update_gamma",
           function(weight_obj,...) standardGeneric("update_gamma"))

setGeneric("record_post_draws",
           function(component_obj,weight_obj,bmbclust_obj,...) standardGeneric("record_post_draws"),signature = c("component_obj","bmbclust_obj"))

setGeneric("tolabel.switching",
           function(bmbclust_obj,...) standardGeneric("tolabel.switching"))

setGeneric("return_label_switching_args",
           function(component_obj,bmbclust_obj,...) standardGeneric("return_label_switching_args"))

setGeneric("identify",
           function(bmbclust_obj,...) standardGeneric("identify"))

setGeneric("reorder_and_record",
           function(component_obj,bmbclust_obj,...) standardGeneric("reorder_and_record"))

setGeneric("identify_w_kmeans",
           function(component_obj,bmbclust_obj,...) standardGeneric("identify_w_kmeans"))

setGeneric("compute_CI",
           function(component_obj,...) standardGeneric("compute_CI"))

setGeneric("prepare_df",
           function(component_obj,...) standardGeneric("prepare_df"))

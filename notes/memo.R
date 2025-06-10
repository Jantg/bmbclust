Reduce(rbind,lapply(res_lca@posteriors$comp_means[MAP_idx],function(M) Reduce(cbind,lapply(M,function(x) apply(x[,1:K_MAP],2,sort,decreasing = FALSE)[-1,]),init=NULL)))


lapply(names(object@col_levels),function(col_i) lapply(meanmat,function(M) lapply(do.call("$",args= list(object@col_levels,col_i,),function(x) do.call("$",args = list(M,col_i))))[[1]]
                                                       
lapply(object@col_levels,function(levels_col) lapply(setNames(1:length(levels_col),levels_col),function(x) x))

lapply(object@col_levels,function(x) lapply(setNames(1:length(x),x),function(y) ))

lapply(setNames(1:length(object@col_levels),names(object@col_levels)),function(idx){ x<-object@col_levels[[idx]]; lapply(setNames(1:length(x),x),function(y){ quantile(simplify2array(lapply(meanmat_new,function(M) do.call("$",args = list(M,names(cat[idx])))[y,1])))})})

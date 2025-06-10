Kp_star <- 3
MAP_idx <- res@posteriors$Kps == Kp_star
mcmc.mean <- simplify2array(lapply(res@posteriors$comp_means[MAP_idx],function(x)x[,,1][1:Kp_star,]))
mcmc.mean <- aperm(mcmc.mean,c(3,1,2))
mcmc.cov <- simplify2array(lapply(res@posteriors$comp_covs[MAP_idx],function(x) t(simplify2array(lapply(1:3,function(k)as.vector(x[k,,]))))))
mcmc.cov <- aperm(mcmc.cov,c(3,1,2))
library(abind)
mcmc.params <- abind::abind(mcmc.mean,mcmc.cov,along = 3)
# is stephens method even appropriate?
p <- simplify2array(lapply(res@posteriors$alloc_probs[MAP_idx],function(x) apply(x[,1,1:Kp_star],1,function(row) exp(row)/sum(exp(row)))))
p <- aperm(p,c(3,2,1))
z <- t(simplify2array(res@posteriors$alloc_vectors[MAP_idx]))

ret_ste <- label.switching::label.switching("STEPHENS",z=z,p = p)
ret_ecriter1 <- label.switching::label.switching("ECR-ITERATIVE-1",z = z,K = Kp_star) 
ret_ecriter2 <- label.switching::label.switching("ECR-ITERATIVE-2",z = z,p=p)
ret_pra <- label.switching::label.switching("PRA",z=z,mcmc = mcmc.params,prapivot = mcmc.params[nrow(mcmc.params),,])
perms <-  apply(ret_ste$permutations$STEPHENS,1,function(x) paste(x,collapse = ""))
perms_ecriter2 <-  apply(ret_ecriter2$permutations$`ECR-ITERATIVE-2`,1,function(x) paste(x,collapse = ""))
perms_pra <-  apply(ret_pra$permutations$PRA,1,function(x) paste(x,collapse = ""))
perms_ecriter1 <-  apply(ret_ecriter1$permutations$`ECR-ITERATIVE-1`,1,function(x) paste(x,collapse = ""))
  
table(perms)

mcmc.mean <- simplify2array(lapply(res@posteriors$comp_means[MAP_idx],function(x)x[1:Kp_star,]))
# for categ
mcmc.mean <- simplify2array(lapply(res@posteriors$comp_means[MAP_idx],function(x)Reduce(rbind,x)[,1:Kp_star]))
mcmc.params <- aperm(mcmc.mean,c(3,2,1))
p <- simplify2array(lapply(res@posteriors$alloc_probs[MAP_idx],function(x) apply(x[,1:Kp_star],1,function(row) exp(row)/sum(exp(row)))))
p <- aperm(p,c(3,2,1))
z <- t(simplify2array(res@posteriors$alloc_vectors[MAP_idx]))

# for count
mcmc.mean <- simplify2array(lapply(res@posteriors$comp_means[MAP_idx],function(x)x[1:Kp_star,]))
mcmc.mean <- t(mcmc.mean)


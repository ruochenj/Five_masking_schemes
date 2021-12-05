# five ways of introducing zeros to a matrix:
set.seed(1234)
# can change this value to set zero_prop of the no zero values to zero.
zero_prop = 0.5
# a toy example matrix
complete_mat <- matrix(rnorm(10000,2,2), 100, 100)
complete_mat[complete_mat < 0.1] = 0

# random mask (all genes)
sce_ct_nzi <- complete_mat
zi_idx <- matrix(rbinom(dim(sce_ct_nzi)[1] * dim(sce_ct_nzi)[2], size = 1, prob = (1-zero_prop)), nrow =  dim(sce_ct_nzi)[1], ncol = dim(sce_ct_nzi)[2])
sce_ct_zi <- sce_ct_nzi * zi_idx
sce_ct_zi1 <- sce_ct_zi
sum(sce_ct_zi1 == 0) / (dim(sce_ct_zi)[1] * dim(sce_ct_zi)[2])

# quantile mask (all genes)
# introduce zero by truncation
sce_ct_nzi2 <- complete_mat
sce_ct_zi2 <- sce_ct_nzi2
idx_nz1 <- which(sce_ct_nzi2 > 0)
cutoff <- quantile(sce_ct_nzi2[idx_nz1], zero_prop)
idx_zero2 <- which(sce_ct_nzi2[idx_nz1] <= cutoff)
idx_zero2 <- sample(idx_nz1[idx_zero2], floor(zero_prop * length(idx_nz1)))
sce_ct_zi2[idx_zero2] <- 0
sum(sce_ct_zi2 == 0) / (dim(sce_ct_zi2)[1] * dim(sce_ct_zi2)[2])
sum(complete_mat == 0) / (dim(sce_ct_nzi)[1] * dim(sce_ct_nzi)[2])
zp_rec <- sum(sce_ct_zi2 == 0) / (dim(sce_ct_nzi)[1] * dim(sce_ct_zi)[2])

# random mask (gene specific)
sce_ct_nzi3 <- complete_mat
sce_ct_zi3 <- sce_ct_nzi3
ZC_avg <- apply(sce_ct_nzi3, 1, FUN = function(x){
  if(sum(x == 0) == length(x)){
    return(c(0,0))
  }else{
    return(c(sum(x > 0), mean(x[x>0])))
  }
})
f <- function(lambda){
  sce_ct_zi3 <- apply(sce_ct_nzi3, 1, FUN = function(x){
    if(sum(x == 0) == length(x)){
      return(x)
    }else{
      nz_idx <- which(x > 0)
      x[nz_idx] <- x[nz_idx] * (1-rbinom(length(nz_idx), size = 1, prob = exp(-lambda*log(mean(x[nz_idx]) + 1.01)^2)))
      return(x)
    }
  })
  return(sum(sce_ct_zi3 == 0) / (dim(sce_ct_nzi3)[1] * dim(sce_ct_nzi3)[2]) - zp_rec)
}
lambda = uniroot(f, c(0,20))$root
sce_ct_zi3 <- apply(sce_ct_nzi3, 1, FUN = function(x){
  if(sum(x == 0) == length(x)){
    return(x)
  }else{
    nz_idx <- which(x > 0)
    x[nz_idx] <- x[nz_idx] * (1-rbinom(length(nz_idx), size = 1, prob = exp(-lambda*log(mean(x[nz_idx]) + 1.01)^2)))
    return(x)
  }
})
sum(sce_ct_zi3 == 0) / (dim(sce_ct_nzi3)[1] * dim(sce_ct_nzi3)[2])


# quantile mask (same percentage)
sce_ct_nzi4 <- complete_mat
sum(sce_ct_nzi4 == 0) / (dim(sce_ct_nzi4)[1] * dim(sce_ct_nzi4)[2])
sce_ct_zi4 <- apply(sce_ct_nzi4, 1, FUN = function(x){
  idx_nz1 <- which(x > 0)
  cutoff <- quantile(x[idx_nz1], zero_prop)
  idx_zero2 <- which(x[idx_nz1] <= cutoff)
  idx_zero2 <- sample(idx_nz1[idx_zero2], floor(zero_prop * length(idx_nz1)))
  x[idx_zero2] <- 0
  return(x)
})
sum(sce_ct_zi4 == 0) / (dim(sce_ct_zi4)[1] * dim(sce_ct_zi4)[2])


# quantile mask (gene specific)
sce_ct_nzi5 <- complete_mat
g <- function(lambda){
  sce_ct_zi5 <- apply(sce_ct_nzi5, 1, FUN = function(x){
    if(sum(x == 0) == length(x)){
      return(x)
    }else{
      x[rank(x) <= (length(x) * exp(-lambda*log(mean(x[x > 0]))^2))] = 0
      return(x)
    }
  })
  return(sum(sce_ct_zi5 == 0) / (dim(sce_ct_nzi5)[1] * dim(sce_ct_nzi5)[2]) - zp_rec)
}
lambda = uniroot(g, c(0,200))$root
sce_ct_zi5 <- apply(sce_ct_nzi5, 1, FUN = function(x){
  if(sum(x == 0) == length(x)){
    return(x)
  }else{
    x[rank(x) <= (length(x) * exp(-lambda*log(mean(x[x > 0]))^2))] = 0
    return(x)
  }
})
sum(sce_ct_zi5 == 0) / (dim(sce_ct_nzi5)[1] * dim(sce_ct_nzi5)[2])

sum(sce_ct_zi1 == 0) / (dim(sce_ct_nzi)[1] * dim(sce_ct_nzi)[2])
sum(sce_ct_zi2 == 0) / (dim(sce_ct_nzi)[1] * dim(sce_ct_nzi)[2])
sum(sce_ct_zi3 == 0) / (dim(sce_ct_nzi)[1] * dim(sce_ct_nzi)[2])
sum(sce_ct_zi4 == 0) / (dim(sce_ct_nzi)[1] * dim(sce_ct_nzi)[2])
sum(sce_ct_zi5 == 0) / (dim(sce_ct_nzi)[1] * dim(sce_ct_nzi)[2])

# Next steps
# Run imputation methods on the scRNA-seq data with increased zero values (sce_ct_zi1, sce_ct_zi2, sce_ct_zi3, sce_ct_zi4, sce_ct_zi5)







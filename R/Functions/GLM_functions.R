# Function from publically available code from  Zuur et al. (2009) - 
## A protocol for data exploration to avoid common statistical problems 
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

## Assessing multicollinearity between predictors when running the dredge function (MuMIn package)
### From user r.jaffe (https://stats.stackexchange.com/questions/82984/how-to-test-and-avoid-multicollinearity-in-mixed-linear-model)
max.r <- function(x){
  corm <- cov2cor(vcov(x))
  corm <- as.matrix(corm)
  if (length(corm)==1){
    corm <- 0
    max(abs(corm))
  } else if (length(corm)==4){
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    cormf <- 0
    max(abs(cormf))
  } else {
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    diag(cormf) <- 0
    max(abs(cormf))
  }
}

## Modified max.r function 
mod.max.r <- function(x) {
  corm <- cov2cor(vcov(x))
  if (ncol(corm) <= 1) {
    return(0) # No predictors or only intercept present
  }
  
  cormf <- corm[-1, -1] # Remove the intercept row and column
  diag(cormf) <- 0 # Set the diagonal to 0 (to ignore self-correlation)
  
  return(max(abs(cormf)))
}

#################################################
# custom functions used in main analysis
#################################################

# function to calculate heterogeneity - squared version, including CVH2 and M2
h.calc2 <- function(mod){
  # I2
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  # s^2_t = total variance
  I2_total <- 100 * (sum(mod$sigma2) / (sum(mod$sigma2) + sigma2_v))
  I2_each <- 100 * (mod$sigma2 / (sum(mod$sigma2) + sigma2_v))
  #names(I2_each) <- paste0("I2_", model$s.names)
  #names(I2_total) <- "I2_Total"
  I2s_Shinichi <- c(I2_total, I2_each)
  
# matrix version  
  W <- solve(mod$V)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total2 <- 100* (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  I2_each2 <- 100* (mod$sigma2 / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  #names(I2_each2) <- paste0("I2_", model$s.names)
  #names(I2_total2) <- "I2_Total2"
  I2s_Wolfgang <- c(I2_total2, I2_each2)
  
  
  # CVH2
  CV_total <- (sum(mod$sigma2) / (mod$beta[1])^2)
  CV_each <- (mod$sigma2 / (mod$beta[1])^2)

  #names(CVB_each) <- paste0("CVB_", mod$s.names)
  #names(CVB_total) <- "CVB_total"
  CVHs <- c(CV_total, CV_each)
  
  # M2
  M_total <- (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$beta[1])^2))
  M_each <- (mod$sigma2) / (sum(mod$sigma2) + (mod$beta[1])^2)
  Ms <- c(M_total, M_each)

  hs <- data.frame(I2s_Shinichi,CVHs,Ms)
  rownames(hs) <- c("Total", mod$s.names)
  return(hs)

}


# function to calculate heterogeneity - original version, including CVH1 and M1
h.calc1 <- function(mod){
  # I2
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  # s^2_t = total variance
  I2_total <- 100 * (sum(mod$sigma2) / (sum(mod$sigma2) + sigma2_v))
  I2_each <- 100 * (mod$sigma2 / (sum(mod$sigma2) + sigma2_v))
  I2s_Shinichi <- c(I2_total, I2_each)
  
  # CVH1
  CV_total <- ( sqrt(sum(mod$sigma2)) / abs(mod$beta[1]) )
  CV_each <- ( sqrt(mod$sigma2) / abs(mod$beta[1]) )
  CVHs <- c(CV_total, CV_each)
  
  # M1
  M_total <- ( sqrt(sum(mod$sigma2)) / (abs(mod$beta[1]) + sqrt(sum(mod$sigma2))) )
  M_each <- ( sqrt(mod$sigma2) / (sqrt(sum(mod$sigma2)) + abs(mod$beta[1])) )
  Ms <- c(M_total, M_each)
 
  
  hs <- data.frame(I2s_Shinichi, CVHs, Ms)
  rownames(hs) <- c("Total", mod$s.names)
  return(hs)
  
}

# function to estimate typical sampling error variance
sigma2_v <- function(mod){
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  return(sigma2_v)
}

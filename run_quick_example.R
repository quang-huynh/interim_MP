
# Code to reduce the number of simulations for testing demonstration

myOM <- readRDS("OM/POP_OM.rds")
source('fn/iMP.R')

nsim <- 3
myOM@cpars <- lapply(myOM@cpars, function(x) {
  if(is.matrix(x)) {
    x[1:nsim, ]
  } else if(is.array(x)) {
    x[1:nsim, , ]
  } else {
    x[1:nsim]
  }
})
myOM@nsim <- nsim

myOM@interval <- 1

MSE_batch_1 <- runMSE(myOM, MPs = c("iMP_avg_5", "iMP_avg_10"), parallel = FALSE, PPD = TRUE, ntrials = 200)
diagnostic_AM(MSE_batch_1)

MSE_batch_3 <- runMSE(myOM, MPs = c("iMP_buffer_5", "iMP_buffer_10"), parallel = FALSE, PPD = TRUE, ntrials = 200)
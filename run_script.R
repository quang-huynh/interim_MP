
######## R version 3.5.3

######## Install packages
# Note that this version of DLMtool is different than the version on CRAN.
# It has been modified so that the index follows vulnerable biomass instead of total biomass
install.packages("DLMtool_5.3.1.tar.gz", repos = NULL) 
install.packages("MSEtool_1.4.1.tar.gz", repos = NULL)

library(MSEtool)

# Setup loops
scenario <- c("base", "hs", "dep", "hd", "lf", "epiM")
OM_name <- c("capelin", "POP", "vs")
seeds <- c(54, 86, 24)

DLMtool::setup(12) # Run in parallel over 12 cores


for(k in seq_along(OM_name)) { ######### Loop over operating model
  
  myOM <- readRDS(paste0("OM/", OM_name[k], "_OM.rds"))
  
  source('fn/iMP.R')
  if(k == 3) source("fn/iMP_vs.R") ######## For Vermilion Snapper, selectivity is dome-shaped in SCA

  
  for(i in seq_along(scenario)) { ######### Loop over scenario
    
    if(i == 2) myOM@beta <- c(1/3, 2/3)
    if(i == 3) myOM@cpars$D <- 0.3 * myOM@cpars$D
    if(i == 4) myOM@beta <- c(1.5, 3)
    if(i == 5) myOM@cpars$D <- 2 * myOM@cpars$D
    if(i == 6) {
      set.seed(seed[k])
      M_mult <- rbinom(myOM@proyears * myOM@nsim, 1, 0.1) * pmin(exp(rnorm(myOM@proyears * myOM@nsim, 0, 2)), 4)
      M_y <- myOM@M[1] * (1 + M_mult)
      M_array_hist <- array(myOM@M[1], dim = c(myOM@nsim, myOM@maxage, myOM@nyears))
      M_array_future <- aperm(array(M_y, dim = c(myOM@nsim, myOM@proyears, myOM@maxage)), perm = c(1, 3, 2))
      myOM@cpars$M_ageArray <- abind::abind(M_array_hist, M_array_future, along = 3)
    }
    
    ######## Fixed TAC MPs and Averaged Index MPs
    myOM@interval <- c(5, 10, 1, 1)
    message(paste("Fixed TAC and Averaged MPs for scenario:", scenario[i]))
    
    MSE_batch_1 <- runMSE(myOM, MPs = c("SCA_5", "SCA_10", "iMP_avg_5", "iMP_avg_10"), 
                          parallel = TRUE, PPD = TRUE, ntrials = 200)
    
    ######## Annual Assessment MP
    myOM@interval <- 1
    MSE_batch_2 <- runMSE(myOM, MPs = "SCA_1", parallel = TRUE, PPD = TRUE, ntrials = 200)
    
    ######## Buffered interim MPs
    myOM@interval <- 1
    MSE_batch_3 <- runMSE(myOM, MPs = c("buffer_iMP_5", "buffer_iMP_10"), parallel = TRUE, PPD = TRUE, ntrials = 200)
    
    ######## Projection MPs
    myOM@interval <- 1
    MSE_batch_4 <- runMSE(myOM, MPs = c("pMP_5", "pMP_10"), parallel = TRUE, PPD = TRUE, ntrials = 200)
    
    ######## Merge MSE and save output
    res <- merge_MSE(MSE_batch_1, MSE_batch_2, MSE_batch_3, MSE_batch_4)
    saveRDS(res, file = paste0("MSE_obj/", OM_name[k], "_", scenario[i], ".rds"))
    
  }
}


sfStop()

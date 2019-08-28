

source('iMP.R')

scenario <- c("base", "hs", "dep")
OM_name <- "POP"
OM_base <- readRDS(paste0("OM/", OM_name, "_OM.rds"))
OM_base <- Replace(OM_base, Precise_Unbiased)

max_ind <- OM_base@maxage - 1 + OM_base@proyears + OM_base@nyears
start_proj <- OM_base@maxage + OM_base@nyears 
#plot(apply(OM_base@cpars$Perr_y[, start_proj:max_ind], 2, mean), typ = 'o')

#corrections <- matrix(apply(OM_base@cpars$Perr_y[, start_proj:max_ind], 2, mean), 250, 50, byrow = TRUE)
corrections <- mean(OM_base@cpars$Perr_y[, start_proj:max_ind])

OM_base@cpars$Perr_y[, start_proj:max_ind] <- 
  OM_base@cpars$Perr_y[, start_proj:max_ind] / corrections

for(i in 3:3) {
  myOM <- OM_base
  if(i == 2) myOM@beta <- c(0.3, 0.6)
  if(i == 3) myOM@cpars$D <- 0.3 * myOM@cpars$D
  
  ######## Reference MPs
  myOM@interval <- 1
  
  setup_parallel(12)
  sfExport(list = "SCA_1")
  
  tim <- proc.time()
  message(paste("RefMPs for scenario:", scenario[i]))
  
  res <- runMSE(myOM, MPs = c("NFref", "FMSYref75", "FMSYref", "curE", "SQ", "AvC"), parallel = TRUE, PPD = TRUE)
  #res <- runMSE(myOM, MPs = c("NFref", "FMSYref75", "FMSYref"), parallel = TRUE, PPD = TRUE)
  
  message("Elapsed time: ", round((proc.time()-tim)[3], 2), " seconds.")
  
  saveRDS(res, file = paste0(scenario[i], "/", OM_name, "_refMP.rds"))
  
  sfStop()
  rm(res)
  
  ######## SCA interval
  myOM@interval <- c(5, 10, 1, 1)
  
  setup_parallel(12)
  tim <- proc.time()
  message(paste("Interval SCAs for scenario:", scenario[i]))
  
  res <- runMSE(myOM, MPs = c("SCA_5", "SCA_10", "iMP_5", "iMP_10"), parallel = TRUE, PPD = TRUE)
  
  message("Elapsed time: ", round((proc.time()-tim)[3], 2), " seconds.")
  
  saveRDS(res, file = paste0(scenario[i], "/", OM_name, "_intervalSCA.rds"))
  
  sfStop()
  rm(res)
  
  ######## SCA annual
  myOM@interval <- 1
  
  setup_parallel(12)
  tim <- proc.time()
  message(paste("Annual SCAs for scenario:", scenario[i]))
  
  res <- runMSE(myOM, MPs = "SCA_1", parallel = TRUE, PPD = TRUE)
  
  message("Elapsed time: ", round((proc.time()-tim)[3], 2), " seconds.")
  
  saveRDS(res, file = paste0(scenario[i], "/", OM_name, "_annualSCA.rds"))
  
  sfStop()
  rm(res)
  
  ######## Buffered interim MPs
  myOM@interval <- 1
  
  setup_parallel(12)
  tim <- proc.time()
  message(paste("Buffered interim MPs for scenario:", scenario[i]))
  
  res <- runMSE(myOM, MPs = c("buffer_iMP_5", "buffer_iMP_10"), parallel = TRUE, PPD = TRUE)
  
  message("Elapsed time: ", round((proc.time()-tim)[3], 2), " seconds.")
  
  saveRDS(res, file = paste0(scenario[i], "/", OM_name, "_buffer_iMP.rds"))
  
  sfStop()
  rm(res)

}




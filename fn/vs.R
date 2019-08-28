

source('iMP_vs.R')

scenario <- c("base", "hs", "dep")
OM_name <- "vs"
OM_base <- readRDS(paste0("OM/", OM_name, "_OM.rds"))
OM_base <- Replace(OM_base, Precise_Unbiased)

for(i in 1:length(scenario)) {
  myOM <- OM_base
  if(i == 2) myOM@beta <- c(1/3, 2/3)
  if(i == 3) myOM@cpars$D <- 0.3 * myOM@cpars$D
  
  ######## Reference MPs
  myOM@interval <- 1
  
  setup_parallel(12)
  sfExport(list = "SCA_1")
  
  tim <- proc.time()
  message(paste("RefMPs for scenario:", scenario[i]))
  
  res <- runMSE(myOM, MPs = c("NFref", "FMSYref75", "FMSYref", "curE", "SQ", "AvC"), parallel = TRUE, PPD = TRUE)
  
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




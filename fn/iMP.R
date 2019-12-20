

interim_MP <- eval(bquote(function(x, Data, reps = 1, assessment_interval, HCR = HCR_MSY, SCA_arg = list(I_type = "VB", CAA_multiplier = 20), 
                                   I_smooth = c("none", "loess", "mean", "buffer"), smooth_par = NULL) {
  dependencies <- .(MSEtool:::get_dependencies("SCA"))
  I_smooth <- match.arg(I_smooth)
  
  current_yr <- Data@Year[length(Data@Year)]
  run_SCA <- current_yr == Data@LHYear
  if (current_yr > Data@LHYear) {
    run_SCA <- current_yr == Data@Misc[[x]]$next_assess_yr
  }
  if (!run_SCA) run_interim_MP <- TRUE else run_interim_MP <- FALSE

  if (run_SCA) {
    # Return TAC = UMSY * VB_current when run_SCA = TRUE
    SCA_formals <- list(x = x, Data = Data)
    do_Assessment <- do.call(SCA_Pope, c(SCA_formals, SCA_arg))
    Assess_output <- Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
    
	if (do_Assessment@conv) {
      Rec <- do.call(match.fun(HCR), list(Assessment = do_Assessment, reps = reps))
      # Set-up references for interim MP
      q <- as.numeric(do_Assessment@SD$value[names(do_Assessment@SD$value) == "q"])
      sigma_buffer <- sd(log(do_Assessment@Obs_Index/do_Assessment@Index), na.rm = TRUE)
      buffer_quantities <- structure(c(Rec@TAC, do_Assessment@Index[length(do_Assessment@Index)], sigma_buffer),
                                     names = c("C_ref", "I_ref", "sigma_buffer"))
      Rec@Misc <- c(list(I = Data@Ind[x, ], q = q, buffer_quantities = buffer_quantities, Last_Assess = current_yr,
                       UMSY = do_Assessment@UMSY, MSY = do_Assessment@MSY, SSB_dep = do_Assessment@SSB_SSB0[length(do_Assessment@SSB_SSB0)],
                       next_assess_yr = current_yr + assessment_interval), Assess_output)
      
      if (cap_TAC) {
        Rec@TAC <- max(Rec@TAC, Rec@Misc$MSY)
      }
      run_interim_MP <- FALSE
      
    } else {
      
      if (current_yr == Data@LHYear || length(Data@Misc[[x]]) == 2) {
        Rec <- new("Rec")
        Rec@TAC <- TACfilter(rep(NA, reps))
        Rec@Misc <- c(list(next_assess_yr = current_yr + 1), Assess_output)
        
        run_interim_MP <- FALSE
      } else {
        run_interim_MP <- TRUE
        next_assess_yr <- current_yr + 1
      }
    }
  }

  if (run_interim_MP) {
    # Estimate new_VB as new_I/q, then new TAC = UMSY * new_VB - also equivalent: TAC = MSY * I_y / I_MSY
    q_ratio <- Data@Ind[x, 1]/Data@Misc[[x]]$I[1]
    q_update <- Data@Misc[[x]]$q * q_ratio
    
    if (I_smooth == "none") new_Index <- Data@Ind[x, length(Data@Ind[x, ])]
    if (I_smooth == "loess") {
      I_df <- data.frame(Year = Data@Year, Ind = Data@Ind[x, ])
      fit <- loess(Ind ~ Year, I_df)
      new_Index <- fit$fitted[length(fit$fitted)]
    }
    
    if (I_smooth == "mean") {
      nyr <- smooth_par[1]
      new_Index <- mean(Data@Ind[x, (length(Data@Ind[x, ])-nyr+1):length(Data@Ind[x, ])], na.rm = TRUE)
    }
    
    if (I_smooth == "buffer") {
      new_Index <- Data@Ind[x, length(Data@Ind[x, ])]
      new_Iref <- q_ratio * Data@Misc[[x]]$buffer_quantities["I_ref"]
      sd_buffer <- Data@Misc[[x]]$buffer_quantities["sigma_buffer"]
      new_TAC <- Data@Misc[[x]]$buffer_quantities["C_ref"] * (new_Index + smooth_par * sd_buffer)/
        (new_Iref + smooth_par * sd_buffer)
    } else {
      VB_curr <- as.numeric(new_Index/q_update)
      new_TAC <- Data@Misc[[x]]$UMSY * VB_curr
    }
    
    Rec <- new("Rec")
    if(is.infinite(TAC_used) || is.na(TAC_used)) stop("Error in TAC during interim")
    Rec@TAC <- TACfilter(TAC_used)
    if (run_SCA) {
      Rec@Misc <- c(list(I = Data@Ind[x, ], q = q_update, buffer_quantities = Data@Misc[[x]]$buffer_quantities, 
                       Last_Assess = Data@Misc[[x]]$Last_Assess,
                       UMSY = Data@Misc[[x]]$UMSY, MSY = Data@Misc[[x]]$MSY,
                       SSB_dep = Data@Misc[[x]]$SSB_dep,
                       next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr)),
                       Assess_output)
      
    } else {
      
      Rec@Misc <- list(I = Data@Ind[x, ], q = q_update, buffer_quantities = Data@Misc[[x]]$buffer_quantities, 
                       Last_Assess = Data@Misc[[x]]$Last_Assess,
                       UMSY = Data@Misc[[x]]$UMSY, MSY = Data@Misc[[x]]$MSY,
                       SSB_dep = Data@Misc[[x]]$SSB_dep,
                       next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr),
                       diagnostic = Data@Misc[[x]]$diagnostic)
    }
  }
  return(Rec)
}))
class(interim_MP) <- "MP"
environment(interim_MP) <- asNamespace("MSEtool")


# Function to make MP with variable options
make_iMP <- function(...) {
  fn <- interim_MP
  dots <- list(...)
  arg_ind <- pmatch(names(dots), names(formals(fn)))
  formals(fn)[arg_ind] <- dots
  class(fn) <- "MP"
  return(fn)
}

# Averaged Index and Buffered Index MPs with (numbers indicate assessment interval)
iMP_avg_5 <- make_iMP(assessment_interval = 5, I_smooth = "mean", smooth_par = 3)
iMP_avg_10 <- make_iMP(assessment_interval = 10, I_smooth = "mean", smooth_par = 3)

iMP_buffer_5 <- make_iMP(assessment_interval = 5, I_smooth = "buffer", smooth_par = 1)
iMP_buffer_10 <- make_iMP(assessment_interval = 10, I_smooth = "buffer", smooth_par = 1)

# Fixed TAC MPs (numbers indicate assessment interval) and annual assessment MPs
SCA_5 <- SCA_10 <- SCA_1 <- make_MP(SCA_Pope, HCR_MSY, diagnostic = "min", I_type = "VB", CAA_multiplier = 20)


# Merge batches of MSE output
merge_MSE <- function(...) {
  dots <- list(...)
  
  slots_identical <- function(slotname, x = dots, is_logical = FALSE) {
    res <- lapply(x, getElement, slotname)
    is_identical <- all(vapply(res[-1], identical, logical(1), res[[1]]))
    if(is_logical) {
      return(is_identical)
    } else return(unique(do.call(c, res)))
  }
  
  slots_identical("Name")
  slots_identical("nyears")
  slots_identical("proyears")
  slots_identical("nsim")
  
  stopifnot(slots_identical("OM", is_logical = TRUE))
  stopifnot(slots_identical("Obs", is_logical = TRUE))
  stopifnot(slots_identical("SSB_hist", is_logical = TRUE))
  stopifnot(slots_identical("CB_hist", is_logical = TRUE))
  stopifnot(slots_identical("FM_hist", is_logical = TRUE)) 
  
  nMPs <- vapply(dots, getElement, numeric(1), "nMPs")
  
  slotvec <- c("B_BMSY", "F_FMSY", "B", "SSB", "VB", "FM", "C", "TAC", "Effort")
  res <- list()
  for(i in 1:length(slotvec)) {
    new_mat <- array(NA, dim = c(slots_identical("nsim"), sum(nMPs), slots_identical("proyears")))
    for(j in 1:length(dots)) {
      if(j == 1) new_mat[, 1:nMPs[1], ] <- getElement(dots[[j]], slotvec[i])
      if(j > 1) new_mat[, (sum(nMPs[1:(j-1)]) + 1):(sum(nMPs[1:(j-1)]) + nMPs[j]), ] <- getElement(dots[[j]], slotvec[i])
    }
    res[[i]] <- new_mat
  }
  
  slotvec2 <- c("PAA", "CAA", "CAL")
  res2 <- list()
  for(i in 1:length(slotvec2)) {
    maxage <- dim(dots[[1]]@PAA)[3]
    if(i < 3) new_mat <- array(NA, dim = c(slots_identical("nsim"), sum(nMPs), maxage))
    if(i == 3) new_mat <- array(NA, dim = c(slots_identical("nsim"), sum(nMPs), length(slots_identical("CALbins"))))
    for(j in 1:length(dots)) {
      if(j == 1) new_mat[, 1:nMPs[1], ] <- getElement(dots[[j]], slotvec2[i])
      if(j > 1) new_mat[, (sum(nMPs[1:(j-1)]) + 1):(sum(nMPs[1:(j-1)]) + nMPs[j]), ] <- getElement(dots[[j]], slotvec2[i])
    }
    res2[[i]] <- new_mat
  }
  
  ## Create MSE Object ####
  MSEout <- new("MSE", Name = slots_identical("Name"), nyears = slots_identical("nyears"), 
                proyears = slots_identical("proyears"), nMPs = length(slots_identical("MPs")), 
                MPs = slots_identical("MPs"), nsim = slots_identical("nsim"), 
                OM = dots[[1]]@OM, Obs = dots[[1]]@Obs, B_BMSY = res[[1]], F_FMSY = res[[2]], B = res[[3]], SSB = res[[4]], 
                VB = res[[5]], FM = res[[6]], res[[7]], TAC = res[[8]], SSB_hist = dots[[1]]@SSB_hist, CB_hist = dots[[1]]@CB_hist, 
                FM_hist = dots[[1]]@FM_hist, Effort = res[[9]], PAA = res2[[1]], CAA = res2[[2]], CAL = res2[[2]], 
                CALbins = slots_identical("CALbins"), Misc = list(Data = do.call(c, lapply(dots, function(x) x@Misc$Data))))
  
  # Store MSE info
  attr(MSEout, "version") <- packageVersion("DLMtool")
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version	
  
  MSEout
}


projection_MP <- eval(bquote(function(x, Data, reps = 1, assessment_interval, HCR = HCR_MSY, SCA_arg = list(I_type = "VB", CAA_multiplier = 20)) {
  dependencies <- .(MSEtool:::get_dependencies("SCA_Pope"))
  
  current_yr <- Data@Year[length(Data@Year)]
  run_SCA <- current_yr == Data@LHYear
  if (current_yr > Data@LHYear) {
    run_SCA <- current_yr == Data@Misc[[x]]$next_assess_yr
  }
  if (!run_SCA) get_projected_TAC <- TRUE else get_projected_TAC <- FALSE
  
  if (run_SCA) {
    # Return TAC = UMSY * VB_current when run_SCA = TRUE
    SCA_formals <- list(x = x, Data = Data)
    do_Assessment <- do.call(SCA_Pope, c(SCA_formals, SCA_arg))
    Assess_output <- Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
    
    if (do_Assessment@conv) {
      Rec <- do.call(match.fun(HCR), list(Assessment = do_Assessment, reps = reps))
      pro <- projection(do_Assessment, FMort = do_Assessment@UMSY, p_years = 52, p_sim = 1, obs_error = c(0, 0), process_error = 0)
      Rec@Misc <- c(list(project_catch = pro@Catch[1, ], last_assess_yr = current_yr, next_assess_yr = current_yr + assessment_interval), Assess_output)
      
      get_projected_TAC <- FALSE
      
    } else {
      
      if (current_yr == Data@LHYear || length(Data@Misc[[x]]) == 2) {
        Rec <- new("Rec")
        Rec@TAC <- TACfilter(rep(NA, reps))
        Rec@Misc <- c(list(next_assess_yr = current_yr + 1), Assess_output)
        
        get_projected_TAC <- FALSE
      } else {
        get_projected_TAC <- TRUE
        next_assess_yr <- current_yr + 1
      }
    }
  }
  
  if (get_projected_TAC) {
    Rec <- new("Rec")
    TAC_used <- Data@Misc[[x]]$project_catch[current_yr - Data@Misc[[x]]$last_assess_yr + 1]
    if(is.infinite(TAC_used) || is.na(TAC_used)) stop("Error in TAC during interim")
    Rec@TAC <- TACfilter(TAC_used)
    if (run_SCA) {
      Rec@Misc <- c(list(project_catch = Data@Misc[[x]]$project_catch, last_assess_yr = Data@Misc[[x]]$last_assess_yr,
                         next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr)),
                    Assess_output)
      
    } else {
      Rec@Misc <- list(project_catch = Data@Misc[[x]]$project_catch, last_assess_yr = Data@Misc[[x]]$last_assess_yr,
                       next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr),
                       diagnostic = Data@Misc[[x]]$diagnostic)
    }
  }
  return(Rec)
  
}))
class(projection_MP) <- "MP"
environment(projection_MP) <- asNamespace("MSEtool")


make_projection_MP <- function(...) {
  fn <- projection_MP
  dots <- list(...)
  arg_ind <- pmatch(names(dots), names(formals(fn)))
  formals(fn)[arg_ind] <- dots
  class(fn) <- "MP"
  return(fn)
}

# Projection MPs
pMP_5 <- make_projection_MP(assessment_interval = 5)
pMP_10 <- make_projection_MP(assessment_interval = 10)

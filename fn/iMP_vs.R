
# Averaged Index and Buffered Index MPs with (numbers indicate assessment interval)
iMP_avg_5 <- make_iMP(assessment_interval = 5, I_smooth = "mean", smooth_par = 3,
                      SCA_arg = list(I_type = "VB", vulnerability = "dome", CAA_multiplier = 20))
iMP_avg_10 <- make_iMP(assessment_interval = 10, I_smooth = "mean", smooth_par = 3,
                       SCA_arg = list(I_type = "VB", vulnerability = "dome", CAA_multiplier = 20))

iMP_buffer_5 <- make_iMP(assessment_interval = 5, I_smooth = "buffer", smooth_par = 1,
                         SCA_arg = list(I_type = "VB", vulnerability = "dome", CAA_multiplier = 20))
iMP_buffer_10 <- make_iMP(assessment_interval = 10, I_smooth = "buffer", smooth_par = 1,
                          SCA_arg = list(I_type = "VB", vulnerability = "dome", CAA_multiplier = 20))


# Fixed TAC MPs (numbers indicate assessment interval) and annual assessment MPs
SCA_5 <- SCA_10 <- SCA_1 <- make_MP(SCA_Pope, HCR_MSY, diagnostic = "min", I_type = "VB", CAA_multiplier = 20,
                                    vulnerability = "dome")

# Projection MPs
pMP_5 <- make_projection_MP(assessment_interval = 5, SCA_arg = list(I_type = "VB", vulnerability = "dome", CAA_multiplier = 20))
pMP_10 <- make_projection_MP(assessment_interval = 10, SCA_arg = list(I_type = "VB", vulnerability = "dome", CAA_multiplier = 20))

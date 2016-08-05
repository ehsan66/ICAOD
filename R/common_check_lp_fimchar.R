

check_lp_fimchar <- function(fimchar, lp){
  # 1-check if the fimchar is existed in the list.
  # 2- match the character with right information matrix
  # 3- check whether the lp has the right length.
  # outpu: th FIM as a function


  ##check the lax also!!!
  ## for michaelis-menten menten model the lax is [0, x0]


  ModelNames <- c("FIM_logistic",
                  "FIM_logistic_4par",
                  "FIM_power_logistic",
                  "FIM_logistic_1par",
                  "FIM_michaelis",
                  #"FIM_michaelis2", "FIM_michaelis3",
                  #  "FIM_sig_emax",
                  "FIM_emax_3par",
                  # "FIM_poi1", "FIM_poi2", "FIM_second_poi",
                  "FIM_exp_3par", "FIM_exp_2par",
                  # "FIM_nbin",
                  "FIM_loglin",
                  # "FIM_richards", "FIM_Weibull",
                  "FIM_comp_inhibition", "FIM_noncomp_inhibition", "FIM_uncomp_inhibition", "FIM_mixed_inhibition"
                  #"FIM_compartmental",
                  #"FIM_sig_emax1", "FIM_sig_emax2", "FIM_sig_emax3" , "FIM_sig_emax4",
                  #"FIM_sig_emax5" , "FIM_sig_emax5", "FIM_sig_emax6",
                  #"FIM_inverse_quadratic1", "FIM_inverse_quadratic2",
                  #"FIM_exp_2par_censor1", "FIM_exp_2par_censor2",
                  #"FIM_exp_3par_censor1", "FIM_exp_3par_censor2",
                  #"FIM_chemical_kinetic",
                  #"FIM_compartmental_4par",
                  #"FIM_klimpel"
  )

  if (!(fimchar %in% ModelNames))
    stop("'",fimchar, "' is not in the list of available models.")

  fimfunc1 <- switch(fimchar, "FIM_logistic" = FIM_logistic, "FIM_power_logistic" = FIM_power_logistic,
                     "FIM_logistic_1par" = FIM_logisitic_1par,
                     "FIM_michaelis" = FIM_michaelis,
                     #"FIM_michaelis2" = FIM_michaelis2, "FIM_michaelis3" = FIM_michaelis3,
                     "FIM_emax_3par" = FIM_emax_3par,
                     # "FIM_poi1" = FIM_Poissson1, "FIM_poi2" = FIM_Poissson2, "FIM_second_poi" = FIM_second_poi,
                     "FIM_exp_3par" = FIM_exp_3par,
                     #                      "FIM_nbin" = FIM_nbin,
                     "FIM_loglin" = FIM_loglin,
                     # "FIM_richards" = FIM_richards,  "FIM_weibull" =  FIM_weibull,
                     "FIM_exp_2par" = FIM_exp_2par,

                     # because the deign is (S, I) we wrote some function in multi_dimensional_designs.R
                     # that makes the fisher information matrix only depends on x = (S, I), not S and I
                     "FIM_comp_inhibition" = FIM_comp_inhibition_x,
                     "FIM_noncomp_inhibition" = FIM_noncomp_inhibition_x,
                     "FIM_uncomp_inhibition" = FIM_uncomp_inhibition_x,
                     "FIM_mixed_inhibition" = FIM_mixed_inhibition_x,
                     #"FIM_compartmental" = FIM_compartmental
                     #                      "FIM_sig_emax" = FIM_sig_emax,
                     #                      "FIM_sig_emax1" = FIM_sig_emax1,
                     #                      "FIM_sig_emax2" = FIM_sig_emax2,
                     #                      "FIM_sig_emax3" = FIM_sig_emax3,
                     #                      "FIM_sig_emax4" = FIM_sig_emax4,
                     #                      "FIM_sig_emax5" = FIM_sig_emax5,
                     #                      "FIM_sig_emax6" = FIM_sig_emax6,
                     #"FIM_inverse_quadratic1" =  FIM_inverse_quadratic1,
                     #"FIM_inverse_quadratic2" = FIM_inverse_quadratic2,
                     "FIM_logistic_4par" = FIM_logistic_4par
                     #                      "FIM_exp_2par_censor1" = FIM_exp_2par_censor1,
                     #                      "FIM_exp_2par_censor2" = FIM_exp_2par_censor2,
                     #                      "FIM_exp_3par_censor1" = FIM_exp_3par_censor1,
                     #                      "FIM_exp_3par_censor2" = FIM_exp_3par_censor2,
                     #"FIM_chemical_kinetic" = FIM_chemical_kinetic,
                     #"FIM_compartmental_4par" = FIM_compartmental_4par,
                     #"FIM_klimpel" = FIM_klimpel
  )

  ModelNames_1par <- c("FIM_logistic_1par")
  ##check the length of lp
  ModelNames_2par <- c("FIM_logistic", "FIM_power_logistic",
                       "FIM_michaelis", "FIM_michaelis2", "FIM_michaelis3",
                       "FIM_poi1", "FIM_poi2",  "FIM_nbin", "FIM_exp_2par", "FIM_compartmental",
                       "FIM_exponential-type1", "FIM_exponential-type2")
                       #"FIM_klimpel", "FIM_chemical_kinetic")


  ModelNames_3par <- c("FIM_sig_emax", "FIM_emax_3par",  "FIM_second_poi",
                       "FIM_exp_3par", "FIM_loglin",
                       "FIM_comp_inhibition", "FIM_noncomp_inhibition", "FIM_uncomp_inhibition",
                       "FIM_inverse_quadratic1", "FIM_inverse_quadratic2")
  ModelNames_4par <- c("FIM_richards", "FIM_weibull", "FIM_mixed_inhibition",
                       "FIM_sig_emax1", "FIM_sig_emax3", "FIM_sig_emax5",
                       "FIM_logistic_4par")
                       #"FIM_compartmental_4par")
  ModelNames_5par <- c("FIM_sig_emax2", "FIM_sig_emax4", "FIM_sig_emax6")

  ## two-dimensional
  if (fimchar %in% ModelNames_1par)
    if (length(lp) != 1)
      stop("For model ", fimchar, ", the length of 'lp' and 'up' must be 1.")

  ## two-dimensional
  if (fimchar %in% ModelNames_2par)
    if (length(lp) != 2)
      stop("For model ", fimchar, ", the length of 'lp' and 'up' must be 2.")

  if (fimchar %in% ModelNames_3par)
    if (length(lp) != 3)
      stop("For model ", fimchar, ", the length of 'lp' and 'up' must be 3.")

  if (fimchar %in% ModelNames_4par)
    if (length(lp) != 4)
      stop("For model ", fimchar, ", the length of 'lp' and 'up' must be 4.")

  if (fimchar %in% ModelNames_5par)
    if (length(lp) != 5)
      stop("For model ", fimchar, ", the length of 'lp' and 'up' must be 5.")


  return(fimfunc1)
}


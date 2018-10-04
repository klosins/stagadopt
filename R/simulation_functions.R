# one simulation
#' @export
one_simulation <- function(sim_num = 1){

  all_sim_options <- data.table::expand.grid(design_letter = c("A", "B","C"),
                                 design_roman = c("I", "II"),
                                 N = params$num_people)


  tau_estimates <- data.table::data.table(design_letter = all_sim_options$design_letter,
                              design_roman = all_sim_options$design_roman,
                              N = all_sim_options$N,
                              tau_base = 0,
                              tau_1 = 0,
                              tau_2 = 0,
                              true_tau = 0,
                              bootstrap = 0,
                              bootstrap_CI_90 = 0,
                              bootstrap_CI_95 = 0,
                              V_hat_did = 0,
                              V_hat_did_CI_90 = 0,
                              V_hat_did_CI_95 = 0,
                              V_exact = 0,
                              V_exact_CI_90 = 0,
                              V_exact_CI_95 = 0,
                              liang_zeger = 0,
                              liang_zeger_CI_90 = 0,
                              liang_zeger_CI_95 = 0)

  for (i in c(1:nrow(tau_estimates))){
    #================
    # Make the data
    #================
    current_data <- design_spec_data(design_letter = tau_estimates$design_letter[i],
                                     num_people = tau_estimates$N[i],
                                     addoption_dates = params$addoption_dates ,
                                     num_T = length(params$time_periods),
                                     sigma_mat = sigma_mat_maker(sigma = 1, rho = 0)) # the base never chages
    # here we randomly add what year each individual has their policy given given the pi_vec
    current_data_with_epsilon <- add_roman_data_to_base_data(current_data,
                                                             tau_estimates$design_roman[i])
    current_data_long <- make_design_spec_data_long(current_data_with_epsilon)
    current_data_realized_long <- make_design_spec_data_long(current_data_with_epsilon[addoption_date
                                                                                       == addoption_date_true,])

    #================
    # Find the tau estimates
    #================
    tau_estimates$tau_1[i] <- tau_1(current_data_realized_long,
                                    get(paste0("pi_vec_",tau_estimates$design_roman[i])))

    tau_estimates$tau_2[i] <- tau_2(current_data_realized_long,
                                    get(paste0("pi_vec_",tau_estimates$design_roman[i])))

    current_mod <- lm(y_value ~ year + as.factor(person_id) + W,
                      data = current_data_realized_long)

    tau_estimates$tau_base[i] <- current_mod$coefficients[["W"]]


    #================
    # Find the var estimates
    #================



    tau_estimates$V_hat_did[i] <- V_hat_did(current_data_realized_long,
                                            get(paste0("pi_vec_", tau_estimates$design_roman[i])))

    tau_estimates$V_exact[i] <- V_exact(current_data_long,
                                        get(paste0("pi_vec_", tau_estimates$design_roman[i])) )

    tau_estimates$bootstrap[i] <- clustered_bootstrap(current_data_realized_long, B = 100)

    tau_estimates$liang_zeger[i] <- liang_zeger(current_data_realized_long)

    #================
    # Find whether or not there is coverage
    #================

    #=====
    # Find the ~true~ tau based on full data for confidence intervals
    #=====
    tau_estimates$true_tau[i] <- tau_2(current_data_long,
                                       get(paste0("pi_vec_",tau_estimates$design_roman[i])))

    #=====
    # save the coverage
    #=====

    tau_estimates$V_hat_did_CI_95[i] <- CI_test(.05, current_mod$coefficients[["W"]],
                                                tau_estimates$true_tau[i], tau_estimates$V_hat_did[i])

    tau_estimates$V_exact_CI_95[i] <- CI_test(.05, current_mod$coefficients[["W"]],
                                              tau_estimates$true_tau[i], tau_estimates$V_exact[i])

    tau_estimates$bootstrap_CI_95[i] <- CI_test(.05, current_mod$coefficients[["W"]],
                                                tau_estimates$true_tau[i], tau_estimates$bootstrap[i])

    tau_estimates$liang_zeger_CI_95[i] <- CI_test(.05, current_mod$coefficients[["W"]],
                                                  tau_estimates$true_tau[i], tau_estimates$liang_zeger[i])

    tau_estimates$V_hat_did_CI_90[i] <- CI_test(.1, current_mod$coefficients[["W"]],
                                                tau_estimates$true_tau[i], tau_estimates$V_hat_did[i])

    tau_estimates$V_exact_CI_90[i] <- CI_test(.1, current_mod$coefficients[["W"]],
                                              tau_estimates$true_tau[i], tau_estimates$V_exact[i])

    tau_estimates$bootstrap_CI_90[i] <- CI_test(.1, current_mod$coefficients[["W"]],
                                                tau_estimates$true_tau[i], tau_estimates$bootstrap[i])

    tau_estimates$liang_zeger_CI_90[i] <- CI_test(.1, current_mod$coefficients[["W"]],
                                                  tau_estimates$true_tau[i], tau_estimates$liang_zeger[i])
  }

  tau_estimates$sim_num <- sim_num # keep track of simulation data
  tau_estimates
}


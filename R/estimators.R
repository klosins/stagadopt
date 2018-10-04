# function 3.1 in paper
#' @export
g <- function(a,t, pi_vec){
  # params
  T <- nrow(pi_vec) - 1 # bc there is always one more adoption data than time period value
  g_val <- as.numeric(a <= t) - sum(pi_vec$pi[1:t]) +
    (1/T)*(ifelse(a <= T, a, 0) - sum(c(1:T)*pi_vec$pi[1:T])) +
    ((T+1)/T)*( as.numeric(a == Inf) - pi_vec[ad == Inf,(pi)] )
  g_val
}

# function 4.2 in paper
#' @export
gamma <- function(a, t, pi_vec, ad, yr){
  g_fix_vec <- function(a,t){
    g_fix_vec <- g(a,t,pi_vec)
  }
  numerator <- pi_vec[ad == a,(pi)]*g(a,t, pi_vec)
  denominator <- sum(pi_vec$pi[as.numeric(as.factor(ad))]*(map2_dbl(ad,  yr, g_fix_vec))^2)
  numerator/denominator
}


# estimate of treatment effect
#' @export
tau_2 <- function(data, pi_vec){
  unique_vals <- unique(data[,.(addoption_date, year_numeric)])
  mean_y <- map2_dbl(unique_vals$addoption_date,
                     unique_vals$year_numeric,
                     ~ mean(data[addoption_date == .x & year_numeric == .y,(y_value)]))
  gammas <- map2_dbl(unique_vals$addoption_date,
                     unique_vals$year_numeric,
                     ~ gamma(.x, .y, pi_vec = pi_vec,
                             ad = rep(params$addoption_dates,
                                      length(params$time_periods)),
                             yr = rep(params$time_periods,
                                      each = length(params$addoption_dates))
                     ))
  tau_2_value <- sum(mean_y*gammas)
  tau_2_value
}




# estimate of treatment effect
#' @export
tau_1 <- function(data, pi_vec){
  W_dot <- purrr::map2_dbl(data$addoption_date_true, data$year_numeric,
                    ~ g(.x, .y, pi_vec = pi_vec))
  tau_1_value <- sum(W_dot*data$y_value)/sum(W_dot^2)
  tau_1_value
}


#  one bootstrap sample
#' @export
one_boot <- function(sim_num, current_data_realized_long){
  num_units <- data.table::uniqueN(current_data_realized_long$person_id)
  sample_units <- data.table::data.table((table(sample(1:num_units, replace =  TRUE))))
  sample_units[, person_id := as.numeric(V1)]
  sample_units[, N := as.numeric(N)]

  boot_DT <- merge(current_data_realized_long, sample_units,
                   by = "person_id", all.x = TRUE)
  boot_DT <- boot_DT[!is.na(N),]
  boot_DT[, ID := .I]

  boot_DT <- boot_DT[rep(boot_DT$ID, boot_DT$N)]

  boot_mod <- lm(y_value ~ year + as.factor(person_id) + W,
                 data = boot_DT)

  boot_mod$coefficients[["W"]]
}

# bootstap with cluster
#' @export
clustered_bootstrap <- function(current_data_realized_long, B = 100){
  clustered_bootstrap_var <- var(unlist(lapply(c(1:B), one_boot, current_data_realized_long)), na.rm = TRUE)
  clustered_bootstrap_var
}

# to help make the clustered bootstrap sample where the number of people with each
# addoption date is fixed
#' @export
one_boot_subset_maker <- function(addoption_date_fixed, current_data_realized_long){
  current_data_realized_long_subset <- current_data_realized_long[addoption_date == addoption_date_fixed,]
  id_units <- unique(current_data_realized_long_subset$person_id)
  sample_units <- data.table::data.table((table(sample(id_units, replace =  TRUE))))
  sample_units[, person_id := as.numeric(V1)]
  sample_units[, N := as.numeric(N)]
  boot_DT <- merge(current_data_realized_long_subset, sample_units,
                   by = "person_id", all.x = TRUE)
  boot_DT <- boot_DT[!is.na(N),]
  boot_DT[, ID := .I]
  boot_DT <- boot_DT[rep(boot_DT$ID, boot_DT$N)]
  boot_DT
}

# Keeping the number of people in each addoption group fixed
#' @export
one_boot_addoption_fixed <- function(sim_num, current_data_realized_long){
  num_units <- data.table::uniqueN(current_data_realized_long$person_id)
  list_boot_subsets <- lapply(unique(current_data_realized_long$addoption_date),
                              one_boot_subset_maker, current_data_realized_long =
                                current_data_realized_long)
  boot_DT <- rbindlist(list_boot_subsets)

  boot_mod <- lm(y_value ~ year + as.factor(person_id) + W,
                 data = boot_DT)

  boot_mod$coefficients[["W"]]
}

# clustered bootstap jusst keeping assoption data fixed
#' @export
clustered_bootstrap_addoption_fixed <- function(current_data_realized_long, B = 100){
  clustered_bootstrap_var <- var(unlist(lapply(c(1:B), one_boot_addoption_fixed, current_data_realized_long)), na.rm = TRUE)
  clustered_bootstrap_var
}



# V hat did no adj DOES NOT have the T-1/N correction
#' @export
V_hat_did_no_adj <- function(current_data_realized_long, pi_vec){
  unique_vals <- unique(current_data_realized_long[,.(addoption_date, year_numeric)])
  unique_vals$mean_y_a_t <- purrr::map2_dbl(unique_vals$addoption_date,
                                     unique_vals$year_numeric,
                                     ~ mean(current_data_realized_long[addoption_date == .x &
                                                                         year_numeric == .y,(y_value)]))


  unique_vals$gammas <- purrr::map2_dbl(unique_vals$addoption_date,
                                 unique_vals$year_numeric,
                                 ~ gamma(.x, .y, pi_vec = pi_vec,
                                         ad = rep(params$addoption_dates,
                                                  length(params$time_periods)),
                                         yr = rep(params$time_periods,
                                                  each = length(params$addoption_dates))
                                 ))
  unique_vals[, mean_y_a_t_gammas := mean_y_a_t*gammas]
  y_bar_gamma_a_DT <- unique_vals[,.(y_bar_gamma_a = sum(mean_y_a_t_gammas)),.(addoption_date)]


  current_data_realized_long <- merge(current_data_realized_long, unique_vals,
                                      by = c("addoption_date", "year_numeric"), all.x = TRUE)
  current_data_realized_long[, y_i_gammas := y_value*gammas]
  y_gamma_a_DT <- current_data_realized_long[,.(y_gamma_a = sum(y_i_gammas)),
                                             .(person_id, addoption_date)] # keep track of addoption date
  y_gamma_a_total_DT <- merge(y_gamma_a_DT, y_bar_gamma_a_DT, by = "addoption_date", all.x = TRUE)
  y_gamma_a_total_DT[, minus_2 := (y_gamma_a - y_bar_gamma_a)^2]

  s_2_gamma_a_DT <- y_gamma_a_total_DT[, .(s_2_gamma_a = sum(minus_2)/(.N - 1),
                                           N_a = .N), .(addoption_date)]
  s_2_gamma_a_DT[, s_2_gamma_a_normalized := s_2_gamma_a*(1/N_a) ]

  V_hat_did_no_adj <- sum(s_2_gamma_a_DT$s_2_gamma_a_normalized)
}


# V hat did DOES have the T-1/N correction
#' @export
V_hat_did <- function(current_data_realized_long, pi_vec){
  unique_vals <- unique(current_data_realized_long[,.(addoption_date, year_numeric)])
  unique_vals$mean_y_a_t <- purrr::map2_dbl(unique_vals$addoption_date,
                                     unique_vals$year_numeric,
                                     ~ mean(current_data_realized_long[addoption_date == .x &
                                                                         year_numeric == .y,(y_value)]))


  unique_vals$gammas <- purrr::map2_dbl(unique_vals$addoption_date,
                                 unique_vals$year_numeric,
                                 ~ gamma(.x, .y, pi_vec = pi_vec,
                                         ad = rep(params$addoption_dates,
                                                  length(params$time_periods)),
                                         yr = rep(params$time_periods,
                                                  each = length(params$addoption_dates))
                                 ))
  unique_vals[, mean_y_a_t_gammas := mean_y_a_t*gammas]
  y_bar_gamma_a_DT <- unique_vals[,.(y_bar_gamma_a = sum(mean_y_a_t_gammas)),.(addoption_date)]


  current_data_realized_long <- merge(current_data_realized_long, unique_vals,
                                      by = c("addoption_date", "year_numeric"), all.x = TRUE)
  current_data_realized_long[, y_i_gammas := y_value*gammas]
  y_gamma_a_DT <- current_data_realized_long[,.(y_gamma_a = sum(y_i_gammas)),
                                             .(person_id, addoption_date)] # keep track of addoption date
  y_gamma_a_total_DT <- merge(y_gamma_a_DT, y_bar_gamma_a_DT, by = "addoption_date", all.x = TRUE)
  y_gamma_a_total_DT[, minus_2 := (y_gamma_a - y_bar_gamma_a)^2]

  s_2_gamma_a_DT <- y_gamma_a_total_DT[, .(s_2_gamma_a = sum(minus_2)/(.N - 1),
                                           N_a = .N), .(addoption_date)]
  s_2_gamma_a_DT[, s_2_gamma_a_normalized := s_2_gamma_a*(1/N_a +
                                                            (uniqueN(current_data_realized_long$year) - 1)/uniqueN(current_data_realized_long$person_id)) ]

  V_hat_did <- sum(s_2_gamma_a_DT$s_2_gamma_a_normalized)
}

# V exact
#' @export
V_exact <- function(current_data_long, pi_vec){
  unique_vals <- unique(current_data_long[,.(addoption_date, year_numeric)])
  unique_vals$mean_y_a_t <- purrr::map2_dbl(unique_vals$addoption_date,
                                     unique_vals$year_numeric,
                                     ~ sum(current_data_long[addoption_date == .x &
                                                               year_numeric == .y, (y_value)]) /uniqueN(current_data_long$person_id))



  unique_vals$gammas <- purrr::map2_dbl(unique_vals$addoption_date,
                                 unique_vals$year_numeric,
                                 ~ gamma(.x, .y, pi_vec = pi_vec,
                                         ad = rep(params$addoption_dates,
                                                  length(params$time_periods)),
                                         yr = rep(params$time_periods,
                                                  each = length(params$addoption_dates))))
  unique_vals[, mean_y_a_t_gammas := mean_y_a_t*gammas]
  y_bar_gamma_a_DT <- unique_vals[,.(y_bar_gamma_a =
                                       sum(mean_y_a_t_gammas)),.(addoption_date)]


  current_data_long <- merge(current_data_long, unique_vals,
                             by = c("addoption_date", "year_numeric"), all.x = TRUE)
  current_data_long[, y_i_gammas := y_value*gammas]
  y_gamma_a_DT <- current_data_long[,.(y_gamma_a = sum(y_i_gammas)),
                                    .(person_id, addoption_date, addoption_date_true)] # keep track of addoption date
  y_gamma_a_total_DT <- merge(y_gamma_a_DT, y_bar_gamma_a_DT,
                              by = "addoption_date", all.x = TRUE)
  y_gamma_a_total_DT[, minus := (y_gamma_a - y_bar_gamma_a)] # by minus I mean de-meaned
  y_gamma_a_total_DT[, minus_2 := (minus)^2]

  s_2_gamma_a_DT <- y_gamma_a_total_DT[, .(s_2_gamma_a = sum(minus_2)/(.N - 1),
                                           sum = sum(minus_2)),
                                       .(addoption_date)]
  # the N_a above will just bc the number of people bc every person gets one
  # of the options. If we want N_a to be real assignment we do the following
  real_N_a <- unique(y_gamma_a_total_DT[,.(person_id,addoption_date_true)])[, .(N_a = .N),
                                                                            .(addoption_date_true)]
  s_2_gamma_a_DT <- merge( s_2_gamma_a_DT, real_N_a, by.x = c("addoption_date"),
                           by.y = c("addoption_date_true"), all.x = TRUE)

  base_var <- sum(s_2_gamma_a_DT$s_2_gamma_a*(1/s_2_gamma_a_DT$N_a +
                                                (uniqueN(current_data_long$year) - 1)/uniqueN(current_data_long$person_id)),
                  na.rm = T)

  # finding covaraince
  s_2_gamma_a_cov_DT <- y_gamma_a_total_DT[, .(addoption_date, addoption_date_true,person_id, minus)]
  s_2_gamma_a_cov_DT <- data.table::dcast(s_2_gamma_a_cov_DT[addoption_date %in% params$addoption_date, ],
                              person_id ~ addoption_date, value.var = "minus")

  cov_combinations <- data.table::data.table(expand.grid( params$addoption_dates,
                                              params$addoption_dates))[Var1 < Var2,]

  cov_value <- purrr::map2_dbl(cov_combinations$Var1, cov_combinations$Var2,
                        ~ sum((s_2_gamma_a_cov_DT[[as.character(.x)]] +
                                 s_2_gamma_a_cov_DT[[as.character(.y)]] )^2)/
                          (uniqueN(current_data_long$person_id) - 1))

  cov_sum_val <- sum(cov_value)/uniqueN(current_data_long$person_id)

  V_exact_val <- base_var - cov_sum_val
  V_exact_val
}



# helper functions for liang_zeger
# g num is group number
#' @export
make_middle_mat <- function(g_num, X_mat, s_indicator_group_dummies, epsilon_vec, nrow_data){
  X_mat_g <- X_mat[s_indicator_group_dummies[,g_num]*c(1:nrow_data),]
  epsilon_vec_g <- epsilon_vec[s_indicator_group_dummies[,g_num]*c(1:nrow_data)]
  middle_mat <- (t(X_mat_g)%*%epsilon_vec_g%*%t(epsilon_vec_g)%*%X_mat_g)
  middle_mat
}


# https://www.princeton.edu/~mkolesar/papers/small-robust.pdf
#' @export
liang_zeger <- function(current_data_realized_long){

  current_data_realized_long <- current_data_realized_long[order(person_id)]
  nrow_data <- nrow(current_data_realized_long)
  Y_vec <- as.matrix(current_data_realized_long$y_value)

  X_mat <- fastDummies::dummy_cols(current_data_realized_long, select_columns = c("year", "person_id"),
                      remove_first_dummy = TRUE)

  X_mat <- X_mat[,-1:-6] # only keep dummies

  X_mat <- as.matrix(X_mat)

  beta_vec <- solve(t(X_mat)%*%X_mat)%*%t(X_mat)%*%Y_vec

  Y_hat_vec <- X_mat%*%beta_vec

  epsilon_vec <- Y_vec - Y_hat_vec

  s_indicator_group_dummies <- dummy_cols(current_data_realized_long$person_id)[,-1]

  base <-  make_middle_mat(1, X_mat = X_mat, s_indicator_group_dummies = s_indicator_group_dummies,
                           epsilon_vec = epsilon_vec, nrow_data = nrow_data )
  for (i in 2:uniqueN(current_data_realized_long$person_id)){
    new <-  make_middle_mat(i, X_mat = X_mat, s_indicator_group_dummies = s_indicator_group_dummies,
                            epsilon_vec = epsilon_vec, nrow_data = nrow_data)
    base <- base + new
  } # in the end base will be our final middle mat


  V_robust <- solve(t(X_mat)%*%X_mat)%*%(base)%*%solve(t(X_mat)%*%X_mat)
  as.numeric(V_robust[1,1]) # gives the W variance
}


# V adjusting
#' @export
V_adj <- function(current_data_realized_long, pi_vec){
  unique_vals <- unique(current_data_realized_long[,.(addoption_date, year_numeric)])
  unique_vals$mean_y_a_t <- purrr::map2_dbl(unique_vals$addoption_date,
                                     unique_vals$year_numeric,
                                     ~ mean(current_data_realized_long[addoption_date == .x &
                                                                         year_numeric == .y,(y_value)]))


  unique_vals$gammas <- purrr::map2_dbl(unique_vals$addoption_date,
                                 unique_vals$year_numeric,
                                 ~ gamma(.x, .y, pi_vec = pi_vec,
                                         ad = rep(params$addoption_dates,
                                                  length(params$time_periods)),
                                         yr = rep(params$time_periods,
                                                  each = length(params$addoption_dates))
                                 ))
  unique_vals[, mean_y_a_t_gammas := mean_y_a_t*gammas]
  y_bar_gamma_a_DT <- unique_vals[,.(y_bar_gamma_a = sum(mean_y_a_t_gammas)),.(addoption_date)]


  current_data_realized_long <- merge(current_data_realized_long, unique_vals,
                                      by = c("addoption_date", "year_numeric"), all.x = TRUE)
  current_data_realized_long[, y_i_gammas := y_value*gammas]
  y_gamma_a_DT <- current_data_realized_long[,.(y_gamma_a = sum(y_i_gammas)),
                                             .(person_id, addoption_date)] # keep track of addoption date
  y_gamma_a_total_DT <- merge(y_gamma_a_DT, y_bar_gamma_a_DT, by = "addoption_date", all.x = TRUE)
  y_gamma_a_total_DT[, minus_2 := (y_gamma_a - y_bar_gamma_a)^2]

  s_2_gamma_a_DT <- y_gamma_a_total_DT[, .(s_2_gamma_a = sum(minus_2)/(.N - 1),
                                           N_a = .N), .(addoption_date)]
  s_2_gamma_a_DT[, s_2_gamma_a_normalized := s_2_gamma_a*(1/N_a +
                                                            (uniqueN(current_data_realized_long$year) - 1)/uniqueN(current_data_realized_long$person_id)    ) ]

  base_var <- sum(s_2_gamma_a_DT$s_2_gamma_a_normalized)

  # covaraince

  s_2_gamma_a_cov_DT <- s_2_gamma_a_DT[, .(addoption_date, s_2_gamma_a)]
  s_2_gamma_a_cov_DT <- setNames(data.table::data.table(t(s_2_gamma_a_cov_DT[,-"addoption_date"])),
                                 as.character(s_2_gamma_a_cov_DT[["addoption_date"]]))

  cov_combinations <- data.table::data.table(expand.grid(as.numeric(colnames(s_2_gamma_a_cov_DT)),
                                             as.numeric(colnames(s_2_gamma_a_cov_DT))))[Var1 <
                                                                                          Var2,]

  cov_value <- purrr::map2_dbl(cov_combinations$Var1, cov_combinations$Var2,
                        ~ sum((s_2_gamma_a_cov_DT[[as.character(.x)]] +
                                 s_2_gamma_a_cov_DT[[as.character(.y)]] -
                                 2*sqrt(s_2_gamma_a_cov_DT[[as.character(.x)]]*s_2_gamma_a_cov_DT[[as.character(.y)]]) )/
                                uniqueN(current_data_realized_long$person_id)))

  cov_sum_val <- sum(cov_value)

  V_adj <- base_var - cov_sum_val
  V_adj
}

# V adjusted poplation takes the adjusted varaince originaly found on pg 20, but in the covaraince term uses population
# moments and not sample moments (like the classic var_adj) does
#' @export
V_adj_pop <- function(current_data_realized_long, current_data_long, pi_vec){
  unique_vals <- unique(current_data_realized_long[,.(addoption_date, year_numeric)])
  unique_vals$mean_y_a_t <- purrr::map2_dbl(unique_vals$addoption_date,
                                     unique_vals$year_numeric,
                                     ~ mean(current_data_realized_long[addoption_date == .x &
                                                                         year_numeric == .y,(y_value)]))


  unique_vals$gammas <- purrr::map2_dbl(unique_vals$addoption_date,
                                 unique_vals$year_numeric,
                                 ~ gamma(.x, .y, pi_vec = pi_vec,
                                         ad = rep(params$addoption_dates,
                                                  length(params$time_periods)),
                                         yr = rep(params$time_periods,
                                                  each = length(params$addoption_dates))
                                 ))
  unique_vals[, mean_y_a_t_gammas := mean_y_a_t*gammas]
  y_bar_gamma_a_DT <- unique_vals[,.(y_bar_gamma_a = sum(mean_y_a_t_gammas)),.(addoption_date)]


  current_data_realized_long <- merge(current_data_realized_long, unique_vals,
                                      by = c("addoption_date", "year_numeric"), all.x = TRUE)
  current_data_realized_long[, y_i_gammas := y_value*gammas]
  y_gamma_a_DT <- current_data_realized_long[,.(y_gamma_a = sum(y_i_gammas)),
                                             .(person_id, addoption_date)] # keep track of addoption date
  y_gamma_a_total_DT <- merge(y_gamma_a_DT, y_bar_gamma_a_DT, by = "addoption_date", all.x = TRUE)
  y_gamma_a_total_DT[, minus_2 := (y_gamma_a - y_bar_gamma_a)^2]

  s_2_gamma_a_DT <- y_gamma_a_total_DT[, .(s_2_gamma_a = sum(minus_2)/(.N - 1),
                                           N_a = .N), .(addoption_date)]
  s_2_gamma_a_DT[, s_2_gamma_a_normalized := s_2_gamma_a*(1/N_a +
                                                            (uniqueN(current_data_realized_long$year) - 1)/uniqueN(current_data_realized_long$person_id)    ) ]

  base_var <- sum(s_2_gamma_a_DT$s_2_gamma_a_normalized)
  #=======================================
  # covaraince NOW WITH POPULATION MOMENTS
  #=======================================
  unique_vals <- unique(current_data_long[,.(addoption_date, year_numeric)])
  unique_vals$mean_y_a_t <- purrr::map2_dbl(unique_vals$addoption_date,
                                     unique_vals$year_numeric,
                                     ~ sum(current_data_long[addoption_date == .x &
                                                               year_numeric == .y, (y_value)]) /uniqueN(current_data_long$person_id))

  unique_vals$gammas <- purrr::map2_dbl(unique_vals$addoption_date,
                                 unique_vals$year_numeric,
                                 ~ gamma(.x, .y, pi_vec = pi_vec,
                                         ad = rep(params$addoption_dates,
                                                  length(params$time_periods)),
                                         yr = rep(params$time_periods,
                                                  each = length(params$addoption_dates))))
  unique_vals[, mean_y_a_t_gammas := mean_y_a_t*gammas]
  y_bar_gamma_a_DT <- unique_vals[,.(y_bar_gamma_a =
                                       sum(mean_y_a_t_gammas)),.(addoption_date)]

  current_data_long <- merge(current_data_long, unique_vals,
                             by = c("addoption_date", "year_numeric"), all.x = TRUE)
  current_data_long[, y_i_gammas := y_value*gammas]
  y_gamma_a_DT <- current_data_long[,.(y_gamma_a = sum(y_i_gammas)),
                                    .(person_id, addoption_date, addoption_date_true)] # keep track of addoption date
  y_gamma_a_total_DT <- merge(y_gamma_a_DT, y_bar_gamma_a_DT,
                              by = "addoption_date", all.x = TRUE)
  y_gamma_a_total_DT[, minus := (y_gamma_a - y_bar_gamma_a)]
  y_gamma_a_total_DT[, minus_2 := (minus)^2]

  s_2_gamma_a_DT <- y_gamma_a_total_DT[, .(s_2_gamma_a = sum(minus_2)/(.N - 1),
                                           sum = sum(minus_2)),
                                       #N_a = .N),
                                       .(addoption_date)]
  # the N_a above will just bc the number of people bc every person gets one
  # of the options. If we want N_a to be real assignment we do the following
  real_N_a <- unique(y_gamma_a_total_DT[,.(person_id,addoption_date_true)])[, .(N_a = .N),
                                                                            .(addoption_date_true)]
  s_2_gamma_a_DT <- merge( s_2_gamma_a_DT, real_N_a, by.x = c("addoption_date"),
                           by.y = c("addoption_date_true"), all.x = TRUE)

  # finding covaraince

  s_2_gamma_a_cov_DT <- s_2_gamma_a_DT[, .(addoption_date, s_2_gamma_a)]
  s_2_gamma_a_cov_DT <- setNames(data.table::data.table(t(s_2_gamma_a_cov_DT[,-"addoption_date"])),
                                 as.character(s_2_gamma_a_cov_DT[["addoption_date"]]))

  cov_combinations <- data.table::data.table(expand.grid(as.numeric(colnames(s_2_gamma_a_cov_DT)),
                                             as.numeric(colnames(s_2_gamma_a_cov_DT))))[Var1 <
                                                                                          Var2,]

  cov_value <- purrr::map2_dbl(cov_combinations$Var1, cov_combinations$Var2,
                        ~ sum((s_2_gamma_a_cov_DT[[as.character(.x)]] +
                                 s_2_gamma_a_cov_DT[[as.character(.y)]] -
                                 2*sqrt(s_2_gamma_a_cov_DT[[as.character(.x)]]*s_2_gamma_a_cov_DT[[as.character(.y)]]) )/
                                uniqueN(current_data_realized_long$person_id)))

  cov_sum_val <- sum(cov_value)

  V_adj <- base_var - cov_sum_val
  V_adj
}



# confidence interval
#' @export
CI_test <- function(alpha, est_coefficent, real_coefficent, est_var){
  as.numeric(est_coefficent - qnorm(1 - alpha/2)*sqrt(est_var) <= real_coefficent &
               est_coefficent + qnorm(1 - .05/2)*sqrt(est_var) >= real_coefficent )
}







#' @importFrom rpart "rpart.control"
#' @importFrom import "from"


# creates base data table that all other columns will be merged to
#' @export
base_data_maker <- function(num_people, num_a, num_T){ # creates template for filling in data
  DT <- data.table::data.table(person_id = rep(1:num_people, each = num_a),
                   addoption_date = rep(c(1:(num_a-1), Inf), num_people)) # assumes last "a" is never
  DT
}


# returns mean potential outcomes given a design letter and addoption date
#' @export
design_letter_vec_maker <- function(design_letter, a){
  if (design_letter == "A"){
    mean_vec <- c(as.numeric(a <= 1),
                  as.numeric(a <= 2) + 3,
                  as.numeric(a <= 3) + 1
    )
  }
  if(design_letter == "B"){
    mean_vec <- c(0,
                  1 + as.numeric(a == 2),
                  1 + as.numeric(a == 2) + 10*as.numeric(a == 3)
    )
  }
  if(design_letter == "C"){
    mean_vec <- c(0,
                  1 + as.numeric(a == 2),
                  1 + 2*as.numeric(a == 2) + 2*as.numeric(a == 3)
    )
  }
  if(design_letter == "Z"){ # zero, the mean is zero for all potential outcome variable
    mean_vec <- c(0,0,0)
  }
  mean_vec
}

# creates data given design letterm addoption dates, the number of time periods
#' @export
design_letter_data_maker <- function(design_letter, addoption_dates, num_T){
  num_a <- length(addoption_dates)
  design_letter_data <- do.call(rbind,
                                (lapply(c(1:num_a), function(i) design_letter_vec_maker(design_letter = design_letter,
                                                                                        a = addoption_dates[i]))))
  design_letter_data <- data.table(design_letter_data)
  colnames(design_letter_data) <- paste0("Y_", 1:ncol(design_letter_data))
  design_letter_data[, addoption_date := addoption_dates]
}

design_roman_data_maker <- function(design_roman, num_people){
  design_roman_data <- data.table(person_id = (1:num_people))
  if (design_roman == "I"){
    design_roman_data[, addoption_date_true := sample(x = rep(params$addoption_dates,
                                                              num_people*pi_vec_I$pi),
                                                      size = num_people,
                                                      replace = FALSE)]
  }
  if (design_roman == "II"){
    design_roman_data[, addoption_date_true := sample(x = rep(params$addoption_dates,
                                                              num_people*pi_vec_II$pi),
                                                      size = num_people,
                                                      replace = FALSE)]
  }
  if (design_roman == "III"){
    design_roman_data[, addoption_date_true := sample(x = rep(params$addoption_dates,
                                                              num_people*pi_vec_III$pi),
                                                      size = num_people,
                                                      replace = FALSE)]
  }
  design_roman_data
}

# generates variance matrix for diagonal var
#' @export
sigma_mat_maker <- function(sigma = 1, rho = 0){
  sigma_mat <- matrix( rep( rho, len=9), nrow = 3)
  diag(sigma_mat) <- 1
  sigma_mat <- sigma^2 * sigma_mat
  sigma_mat
}


# generates noise in data
#' @export
epsilon_data_maker <- function(sigma_mat = sigma_mat_maker(sigma = 1, rho = 0),
                               num_people = 5){
  epsilon_data <- mvrnorm(n = num_people, c(0,0,0), sigma_mat)
  epsilon_data <- as.data.table(epsilon_data)
  colnames(epsilon_data) <- paste0("Y_", 1:ncol(epsilon_data))
  epsilon_data[, person_id := 1:num_people]
  epsilon_data
}


# uses functions above to make whole dataset
#' @export
# number people same as design number
design_spec_data <- function(design_letter, num_people, addoption_dates,
                             num_T, sigma_mat = sigma_mat_maker(sigma = 1, rho = 0)){
  #==============
  # params
  #==============
  num_a <- length(addoption_dates)
  #==============
  # creating raw datasets
  #==============
  base_data <- base_data_maker(num_people, num_a, num_T)
  epsilon_data <- epsilon_data_maker(sigma_mat, num_people)
  design_letter_data <- design_letter_data_maker(design_letter, addoption_dates, num_T)
  #==============
  # mergeing data together
  #==============
  # bring design letter values into the base
  base_data <- merge(base_data, design_letter_data, by = "addoption_date", all.x = TRUE)
  # add in epslion shocks
  base_data <- merge(base_data, epsilon_data, by = "person_id", all.x = TRUE)

  base_data[, Y_1 := Y_1.x + Y_1.y]
  base_data[, Y_2 :=  Y_2.x + Y_2.y]
  base_data[, Y_3 := Y_3.x + Y_3.y]
  base_data[, c("Y_1.x", "Y_1.y", "Y_2.x", "Y_2.y", "Y_3.x", "Y_3.y") := NULL]

  #==============
  # Return data
  #==============
  base_data
}

# I seperate this part of the data generating process because "we keep the set of potential outcomes fixed given a particular design, but in each simulation we change who adopts at that date"
add_roman_data_to_base_data <- function(base_data, design_roman ){
  num_people <- uniqueN(base_data$person_id)
  design_roman_data <- design_roman_data_maker(design_roman, num_people)
  # add in indicator that gives which option is relized
  base_data <- merge(base_data, design_roman_data, by = "person_id", all.x = TRUE)
  base_data
}


make_design_spec_data_long <- function(data){
  year_cols <- c(paste0("Y_",params$time_periods))
  data_long <- gather(data, year, y_value, year_cols, factor_key=TRUE) # make long for regressions
  setDT(data_long)
  data_long[, year_numeric := as.numeric(gsub("\\D", "", data_long$year)) ]
  data_long[, W := as.numeric(year_numeric >= addoption_date) ]
  data_long
}


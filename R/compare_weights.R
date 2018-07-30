#' Compare various methods to a simulated data set
#'
#' \code{compare_weights} will simulate a data set (according to \code{design})
#' and will calculate standardized z-statistics for a variety of methods
#' \param{weights_list_list} A list, where each element is a list describing a method. For example,
#' \code{list(method="fixed_c",delay=6)} or \code{list(method="landmark", time = 20)}.
#' \param{design} A list of design information, containing:
#' \code{med_c}, \code{rate_e_1}, \code{rate_e_2}, \code{rec_period}, \code{rec_power}, \code{delay}, \code{max_cal_t}
#' \return A vector of standardized test statistics, one for each method in \code{weights_list_list}.
#' @export

############################################

add_weights_list = function(weights_list, risk_table){

  if (weights_list$method == "fh"){
    weights_list$delay = NULL
  }
  else {
    weights_list$rho = NULL
    weights_list$gamma = NULL
  }

  add_weights(risk_table,
              method = weights_list$method,
              delay = weights_list$delay,
              rho = weights_list$rho,
              gamma = weights_list$gamma)

}

###########################################

compare_weights = function(weights_list_list, design){

  which_weights = unlist(lapply(weights_list_list,
                                function(x) x$method != "landmark"))

  real_weights = weights_list_list[which_weights]
  landmarks = weights_list_list[!which_weights]

  df = modestWLRT::get_df(med_c = design$med_c,
                          rate_e_1 = design$rate_e_1,
                          rate_e_2 = design$rate_e_2,
                          rec_period = design$rec_period,
                          rec_power = design$rec_power,
                          delay = design$delay,
                          max_cal_t = design$max_cal_t)


  risk_table = get_risk_table(df)
  
  if (length(real_weights) > 0){
  
    w_risk_table_list = lapply(real_weights,
                               add_weights_list,
                               risk_table = risk_table)
  
    w_z = lapply(w_risk_table_list, get_zs) %>%
            lapply(function(x)x[1]) %>%
            unlist()

  }
  else w_z = numeric(0)
  
  if (length(landmarks) > 0){
  
    land_times = unlist(lapply(landmarks, function(x)x$time))
    
    landmark_z = landmark(df = df,
                          time = land_times)
    
    names(landmark_z) = paste("landmark", land_times)
    
  }
  else landmark_z = numeric(0)
  
  c(w_z, landmark_z)

}
#' Get (standardized) score and weighted log-rank statistics.
#' 
#' \code{get_zs} Calculate (standardized) score and weighted log-rank statistics from risk table (with weights).
#' @param{risk_table} A risk table with weights produced by the functions \code{get_risk_table} and \code{add_weights}.
#' @param{method} If \code{method = "u"} then the logrank formula is used. If \code{method = "s"} then the score test formula is used.
#' @return The standardized weighted log-rank statistic.
#' @export
#' 
#' 
get_zs = function(risk_table, method = "u"){
  
  if (class(risk_table) == "list") risk_table = risk_table$risk_table
  
  n_e = max(risk_table$n_e)
  n_c = max(risk_table$n_c)
  
  # formulas for S and V[S] in Leton and Zuluaga pg. 595.
  
  s = with(risk_table, sum(d_c * c + l_c * C))
  v_s = with(risk_table, sum(d * c ^ 2 + l * C ^ 2)) * n_c * n_e / (n_c + n_e) / (n_c + n_e - 1)
  
  z_s = s / sqrt(v_s)
  
  # formulas for U and V[U] in Leton and Zuluaga pg. 596.
  
  u = with(risk_table, sum(w * (d_c - d * n_c / n)))
  v_u = with(risk_table, sum(w^2 * n_c * n_e * d * (n - d) / n / n / (n - 1), na.rm = TRUE))
  
  z_u = u / sqrt(v_u)

  if (method == "u") return(z_u)
  if (method == "s") return(z_s)
  c(z_u = z_u,
    z_s = z_s)
  
  
}

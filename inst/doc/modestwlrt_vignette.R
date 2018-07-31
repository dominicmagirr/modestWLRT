## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("dominicmagirr/modestWLRT")

## ----include=FALSE-------------------------------------------------------
library(dplyr)
library(ggplot2)
devtools::load_all("..")

## ------------------------------------------------------------------------
example_data = delayed_effect_sim(n_c = 10,
                                  n_e = 10,
                                  rec_period = 12,
                                  rec_power = 1,
                                  med_c = 15,
                                  rate_e_1 = log(2) / 15,
                                  rate_e_2 = 0.03,
                                  delay = 6,
                                  max_cal_t  = 36)

example_data

## ------------------------------------------------------------------------
example_risk_table = get_risk_table(example_data)

example_risk_table

## ------------------------------------------------------------------------
modest_weights = add_weights(example_risk_table, 
                             method = "fixed_c", 
                             delay = 12, 
                             plot_weights = TRUE)

## ------------------------------------------------------------------------
modest_weights$risk_table

## ----fig.cap = "Scores/weights from a modestWLRT."-----------------------
modest_weights$p

## ------------------------------------------------------------------------
get_zs(modest_weights)


#' Groundwater Lagged Return Flow Model 
#' 
#' Analytical lagged return flow model that calculates weekly groundwater accretion back to stream from 
#' weekly volumetric groundwater inputs. Model implements version of the Glover equation presented by Knight et al. (2005) 
#' where an infinite series of paired image wells are used to approximate the effect of an impermeable alluvial
#' valley boundary. Generally,  the first 20 pairs in the infinite series are probably sufficient and 
#' the user can assume all other effects are negligible.
#' 
#'
#'@author Seth Mason and Alex Brooks
#'@organization Lotic Hydrological
#' 
#'@param gw_input groundwater input time series [m/week]
#'@param inundation_area area inundated [m^2]
#'@param sim_years = number of simulation years
#'@param a  perpendicular distance between the augmentation structure and the stream [m]
#'@param b aquifer thickness [m]
#'@param aquifer_ksat saturated hydraulic conductivity of aquifer [m/wk]
#'@param Sy aquifer specific yield [-]
#'@param w alluvial valley width [m]
#'@param images the number of paired image wells to use from the infinite series [#]

#' @return A dataframe with timeseries of groundwater contributions [acre-feet] and weekly returns flow volumes [acre-feet] and rates [cfs]
require(tidyverse)
require(dplyr)

# function for iterating over parameters values
gw_lagged_return_flow <- function(gw_input, inundation_area, sim_years, a, b, aquifier_ksat, Sy, Tr, w, images=20){
  
  #calculate aquifer transmisivity 
  Tr = Ksat_s*b_s #[m^2/wk]
  
  #calculate the unit response of a one week irrigation impulse over a 10 year period

  #accretions to alluvial aquifer [m^3/week]
  gw_contribution = gw_input *inundation_area
  gw_contribution[gw_contribution <0] <- 0
  
  #calculate the unit response of a one week irrigation impulse over a 100 year period
  unit_response_function = unit.response(n = 7, ts = 365*100, a, Sy, Tr, w, images = 20)
  unit_response_function[unit_response_function<0] <- 0
  model_sim_year=10 # number of years over which to return outputs 
  
  #model return flow response using the full 10-year ramp up period
  return_flow = data.frame(rep(x = gw_contribution, times = sim_years))
  for (i in 1:nrow(return_flow)){
    column = i+1
    return_flow[,column] = lag(return_flow[i,1]*unit_response_function[1:nrow(return_flow)], n = i, default = 0)
  }
  
  return_flow_output= tibble(
    return_flow_af_wk = return_flow[,-1] %>%
      rowSums()*0.000810714 #sum impulse responses and convert from m3 to acre feet
  ) %>% 
    mutate(
      return_flow_cfs=return_flow_af_wk*0.0720237,
      return_flow_cume= cumsum(return_flow_af_wk)
    ) %>% 
    rowid_to_column(var='week')  %>% 
    mutate(
      gw_contribution_acre = sum(gw_contribution)*0.000247105
    ) 
  
}

#### Helper Functions ####

## internal function to compute the complementary error function using base R
erfc    <- function(x) {  # 1 - erf(x)
  2 * pnorm(-sqrt(2) * x)
}

glover <- function(t, a, Sy, Tr, w, images) {
  #' Streamflow return fractions for a fully-penetrating stream and no streambed resistance to flow.
  #' Function uses the version of the Glover equation presented by Knight et al. (2005) where an
  #' infinite series of paired image wells are used to approximate the effect of an impermeable alluvial
  #' valley boundary. Generally,  the first 10 pairs in the infinite series are probably sufficient and 
  #' the user can assume all other effects are negligible.
  #' 
  #' @param t solution time step of interest [T]
  #' @param a flowpath distance between the irrigation centroid and the stream [L]
  #' @param Sy aquifer specific yield [-]
  #' @param Tr aquifer transmissivity [L2/T]
  #' @param w alluvial valley width [L]
  #' @param images the number of paired image wells to use from the infinite series
  #' @return A numeric vector of return flow fractions \code{Qf} for a given groundwater recharge rate [-].
  #' If the pumping rate of the well (\code{Qw}; [L3/T]) is known, you can calculate volumetric streamflow depletion [L3/T] as \code{Qf*Qw}
  
  D = Tr/Sy  
  recharge = erfc(a / (2*(D*t)^(1/2)) ) #Glover Equation
  
  image_wells = 0
  for (n in 1:images) {
    image_wells = image_wells +
      (-1)^(n+1) * ( erfc((2*n*w - a) / (2*(D*t)^(1/2))  ) - erfc((2*n*w + a) / (2*(D*t)^(1/2)) ) )
  }
  
  out = recharge + image_wells
  
}

unit.response <- function(n, ts, a, Sy, Tr, w, images) {
  #' Calculate the return flow response pattern for a given impulse of infiltrated water. 
  #' 
  #' @param n  duration of the groundwater percolation impulse [T]
  #' @param ts number of time steps to compute groundwater returns over [T]
  #' @param a  perpendicular distance between the augmentation structure and the stream [L]
  #' @param Sy aquifer specific yield [-]
  #' @param Tr aquifer transmissivity [L2/T]
  #' @param w alluvial valley width [L]
  #' @param images the number of paired image wells to use from the infinite series
  #' @return A numeric vector of return flow fractions \code{Qf} for the time steps of interest [-].
  
  impulse_start = 1
  impulse_end = n
  impulse_rate = 1
  t = seq(from = 1, to = ts, by = 1)
  starts.all <- base::matrix(impulse_start, nrow = length(t), ncol = length(impulse_start), byrow = T)
  stops.all <- base::matrix(impulse_end, nrow = length(t), ncol = length(impulse_start), byrow = T)
  rates.all <- base::matrix(impulse_rate, nrow = length(t), ncol = length(impulse_start), byrow = T)
  t.all <- base::matrix(t, nrow = length(t), ncol = length(impulse_start))
  
  # calculate time since each pumping interval starts/stops, bounded at 0
  t.starts <- t.all+1 - starts.all
  t.starts[t.starts < 0] <- 0
  
  t.stops <- t.all - stops.all
  t.stops[t.stops < 0] <- 0
  
  # vectorize for calculations
  t.starts.vec <- c(t.starts)
  t.stops.vec <- c(t.stops)
  rates.all.vec <- c(rates.all)
  Qs.all.vec <- rep(0, length(t.starts.vec))
  
  # calculate depletion
  Qs.all.vec[t.starts.vec > 0] <-
    rates.all.vec[t.starts.vec > 0] *
    (glover(t = t.starts.vec[t.starts.vec > 0], a, Sy, Tr, w, images) -
       glover(t = t.stops.vec[t.starts.vec > 0], a, Sy, Tr, w, images))
  
  weeks = rep(1:(ts/n), each=n)
  out = data.frame(fraction = Qs.all.vec[1:length(weeks)], week = weeks)  %>%
    group_by(week) %>%
    dplyr::summarize(
      fraction = sum(fraction),
    ) %>%
    mutate(
      # fraction = ifelse(fraction<0 ,0,fraction)#,
      fraction = fraction/n
    ) %>% #normalize fractional response
    select(fraction)
  
  return(out$fraction)
  
}


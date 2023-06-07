#' Soil Water Balance Model
#' 
#' Hourly 1D soil column water balance model using input time series for irrigation, precipitation and evapotranspiration. 
#' The model partitions soil into five layers with ET partitioned between the top four layers based on rooting density.
#' The model assumes no groundwater is present within the rooting zone and that capillary action and lateral flows are negligible.
#' Model outputs include soil water content and fluxes including surface runoff, infiltration, deep percolation, and evapotranspiration.    
#'
#' @author Seth Mason and Alex Brooks
#' @organization Lotic Hydrological
#' 
#' @param model_inputs Dataframe with model time series inputs  [ft/hour]. Dataframe should contain columns labeled irrigation, precipitation, & ET. All values should be positive or zero.   
#' @param root_zone_depth depth of rooting zone (zone where ET is abstracted) [ft]
#' @param total_profile_depth total depth of vadose zone [ft]
#' @param HZ1_soil_porosity porosity of top soil layer [%]. Value should be between 0-1.
#' @param HZ2_soil_porosity porosity of 2nd soil layer [%]. Value should be between 0-1.
#' @param HZ3_soil_porosity porosity of 3rd soil layer [%]. Value should be between 0-1.
#' @param HZ4_soil_porosity porosity of 4th soil layer [%]. Value should be between 0-1.
#' @param HZ5_soil_porosity porosity of 5th soil layer [%]. Value should be between 0-1.
#' @param HZ1_Field_Capacity field capacity of top soil layer [%]. Value should be between 0-1.
#' @param HZ2_Field_Capacity field capacity of 2nd soil layer [%]. Value should be between 0-1.
#' @param HZ3_Field_Capacity field capacity of 3rd soil layer [%]. Value should be between 0-1.
#' @param HZ4_Field_Capacity field capacity of 4th soil layer [%]. Value should be between 0-1.
#' @param HZ5_Field_Capacity field capacity of 5th soil layer [%]. Value should be between 0-1.
#' @param HZ1_Wilting_Point wilting point of top soil layer [%]. Value should be between 0-1.
#' @param HZ2_Wilting_Point wilting point of 2nd soil layer [%]. Value should be between 0-1.
#' @param HZ3_Wilting_Point wilting point of 3rd soil layer [%]. Value should be between 0-1.
#' @param HZ4_Wilting_Point wilting point of 4th soil layer [%]. Value should be between 0-1.
#' @param Infiltration_Rate infiltration rate [ft/hour]
#' @param soil_KSAT saturated hydraulic conductivity  in soil [ft/hour]
#' @param HZ1_root_distribution proportion of roots in top soil layer [%]. Value should be between 0-1.
#' @param HZ2_root_distribution proportion of roots in 2nd soil layer [%]. Value should be between 0-1.
#' @param HZ3_root_distribution proportion of total roots in 3rd soil layer [%]. Value should be between 0-1.
#' @param HZ4_root_distribution proportion of total roots in 4th soil layer [%]. Value should be between 0-1.
#' @param HZ1_water_content initial volumetric water content in top soil layer [%]. Value should be between 0-1.
#' @param HZ2_water_content initial volumetric water content in 2nd soil layer [%]. Value should be between 0-1.
#' @param HZ3_water_content initial volumetric water content in 3rd soil layer [%]. Value should be between 0-1.
#' @param HZ4_water_content initial volumetric water content in 4th soil layer [%]. Value should be between 0-1.

#' @return A dataframe with model outputs including conditions (instantaneous values) and state variables (cumulative values). 

require(tidyverse)
require(deSolve)

soil_wb_model <- function(model_inputs, 
                         root_zone_depth =3,total_profile_depth = 6, 
                         HZ1_soil_porosity = 0.5,HZ2_soil_porosity = 0.5,
                         HZ3_soil_porosity = 0.5,HZ4_soil_porosity = 0.5,
                         HZ5_soil_porosity = 0.5,HZ1_Field_Capacity = 0.4,
                         HZ2_Field_Capacity = 0.4,HZ3_Field_Capacity = 0.4,
                         HZ4_Field_Capacity = 0.4,HZ5_Field_Capacity = 0.4,
                         HZ1_Wilting_Point = 0.15,HZ2_Wilting_Point = 0.15,
                         HZ3_Wilting_Point = 0.15,HZ4_Wilting_Point = 0.15,
                         Infiltration_Rate= 0.007, soil_KSAT = 0.01,
                         HZ1_root_distribution=0.5,
                         HZ2_root_distribution= 0.25, HZ3_root_distribution=0.15,
                         HZ4_root_distribution=0.1,HZ1_water_content =0.4,
                         HZ2_water_content=0.4,HZ3_water_content=0.4,
                         HZ4_water_content=0.4){  
  
  #Inputs
  Irr <- model_inputs$irrigation #ft/hour
  Ppt <- model_inputs$precipitation #ft/hour
  ET <- model_inputs$ET #ft/hour
  
  #set the model simulation time
  simtime <- seq(1,365*24) # hours
  
  parameters = c(
    Surface_Infiltration_Rate = Infiltration_Rate, # ft/hour
    #Average_Percolation_Rate = P_avg/t_step, # ft/day
    Total_Profile_Depth = total_profile_depth, # ft
    Root_Zone_Depth = root_zone_depth, # ft
    HZ1_Root_Distribution = HZ1_root_distribution, # %
    HZ2_Root_Distribution = HZ2_root_distribution, # %
    HZ3_Root_Distribution = HZ3_root_distribution, # %
    HZ4_Root_Distribution = HZ4_root_distribution # %
  )
  
 state_vars = c( 
    Ponded_Water = 0.0, # ft
    HZ1_Water_Depth = (root_zone_depth/4) * HZ1_water_content, # ft
    HZ2_Water_Depth = (root_zone_depth/4) * HZ2_water_content, # ft
    HZ3_Water_Depth = (root_zone_depth/4) * HZ3_water_content, # ft
    HZ4_Water_Depth = (root_zone_depth/4) * HZ4_water_content, # ft
    HZ5_Water_Depth = (total_profile_depth - root_zone_depth) * HZ5_water_content, # ft
    Consumptive_Use_Depth = 0.0, # ft
    GW_Depth = 0.0, # ft
    Infiltration_Depth = 0.0 #ft
  )
  
  #create the system of equations
  wb = function (time, state_vars, parameters) {
    Irrigation <- Irr[time] #time series of irrigation in ft/hour
    Precipitation <- Ppt[time] #time series of precipitation in ft/hour
    ET_Demand <- ET[time] #time series of ET demands in ft/hour
    with(as.list(c(state_vars, parameters)), {
      HZ1_Depth <- Root_Zone_Depth/4 # ft
      HZ1_WC <- pmin(HZ1_Water_Depth/HZ1_Depth, HZ1_soil_porosity) # %
      HZ2_Depth <- Root_Zone_Depth/4 # ft
      HZ2_WC <- pmin(HZ2_Water_Depth/HZ2_Depth, HZ2_soil_porosity) # %
      HZ3_Depth <- Root_Zone_Depth/4 # ft
      HZ3_WC <- pmin(HZ3_Water_Depth/HZ3_Depth, HZ3_soil_porosity) # %
      HZ4_Depth <- Root_Zone_Depth/4 # ft
      HZ4_WC <- pmin(HZ4_Water_Depth/HZ4_Depth, HZ4_soil_porosity) # %
      HZ5_Depth <- Total_Profile_Depth - Root_Zone_Depth # ft
      HZ5_WC <- pmin(HZ5_Water_Depth/HZ5_Depth, HZ5_soil_porosity) # %
      HZ1_TT <- ((HZ1_Depth*HZ1_soil_porosity)-(HZ1_Depth*HZ1_Field_Capacity))/soil_KSAT # travel time (days)
      HZ2_TT <- ((HZ2_Depth*HZ2_soil_porosity)-(HZ2_Depth*HZ2_Field_Capacity))/soil_KSAT # travel time (days)
      HZ3_TT <- ((HZ3_Depth*HZ3_soil_porosity)-(HZ3_Depth*HZ3_Field_Capacity))/soil_KSAT # travel time (days)
      HZ4_TT <- ((HZ4_Depth*HZ4_soil_porosity)-(HZ4_Depth*HZ4_Field_Capacity))/soil_KSAT # travel time (days)
      HZ5_TT <- ((HZ5_Depth*HZ5_soil_porosity)-(HZ5_Depth*HZ5_Field_Capacity))/soil_KSAT # travel time (days)
      Inflow <- Irrigation + (Precipitation) # ft/day
      HZ1_FC <- ifelse(HZ1_WC >= HZ1_Field_Capacity, 1, 0) # field capacity flag
      HZ1_PWP <- ifelse(HZ1_WC <= HZ1_Wilting_Point, 0, 1) # wilting point flag
      HZ1_Sat <- ifelse(HZ1_WC >= HZ1_soil_porosity, 0, 1) # saturation flag
      HZ2_FC <- ifelse(HZ2_WC >= HZ2_Field_Capacity, 1, 0) # field capacity flag
      HZ2_PWP <- ifelse(HZ2_WC <= HZ2_Wilting_Point, 0, 1) # wilting point flag
      HZ2_Sat <- ifelse(HZ2_WC >= HZ2_soil_porosity, 0, 1) # saturation flag
      HZ3_FC <- ifelse(HZ3_WC >= HZ3_Field_Capacity, 1, 0) # field capacity flag
      HZ3_PWP <- ifelse(HZ3_WC <= HZ3_Wilting_Point, 0, 1) # wilting point flag
      HZ3_Sat <- ifelse(HZ3_WC >= HZ3_soil_porosity, 0, 1) # saturation flag
      HZ4_FC <- ifelse(HZ4_WC >= HZ4_Field_Capacity, 1, 0) # field capacity flag
      HZ4_PWP <- ifelse(HZ4_WC <= HZ4_Wilting_Point, 0, 1) # wilting point flag
      HZ4_Sat <- ifelse(HZ4_WC >= HZ4_soil_porosity, 0, 1) # saturation flag
      HZ5_FC <- ifelse(HZ5_WC >= HZ5_Field_Capacity, 1, 0) # field capacity flag
      HZ5_Sat <- ifelse(HZ5_WC >= HZ5_soil_porosity, 0, 1) # saturation flag
      HZ1_Void <- (HZ1_Depth*HZ1_soil_porosity)-HZ1_Water_Depth
  
      Ponding_Exists <- ifelse(Ponded_Water > 0, 1, 0) #ponding flag
      ET1 <- ET_Demand * HZ1_PWP * HZ1_Root_Distribution # ft/hour
      ET2 <- ET_Demand * HZ2_PWP * HZ2_Root_Distribution # ft/hour
      ET3 <- ET_Demand * HZ3_PWP * HZ3_Root_Distribution # ft/hour
      ET4 <- ET_Demand * HZ4_PWP * HZ4_Root_Distribution # ft/hour
      Infiltration <- ifelse(Ponded_Water<Surface_Infiltration_Rate,pmin(Ponded_Water,HZ1_Void),
                             pmin(Surface_Infiltration_Rate,HZ1_Void)) * HZ1_Sat * Ponding_Exists # ft/hour
      Surface_Runoff <- pmax(1*(Ponded_Water)-Infiltration,0) # ft/hour
      Irrigation_Efficiency <- Infiltration/Irrigation # %
      HZ1_Excess <-  pmax(HZ1_Water_Depth - (HZ1_Depth*HZ1_Field_Capacity),0) # ft
      HZ2_Excess <-  pmax(HZ2_Water_Depth - (HZ2_Depth*HZ2_Field_Capacity),0) # ft
      HZ3_Excess <-  pmax(HZ3_Water_Depth - (HZ3_Depth*HZ3_Field_Capacity),0) # ft
      HZ4_Excess <-  pmax(HZ4_Water_Depth - (HZ4_Depth*HZ4_Field_Capacity),0) # ft
      HZ5_Excess <-  pmax(HZ5_Water_Depth - (HZ5_Depth*HZ5_Field_Capacity),0) # ft

      HZ1_Percolation <- HZ1_Excess * (1 - exp(-1/(HZ1_TT*24))) # ft/hour
      HZ2_Percolation <- HZ2_Excess * (1 - exp(-1/(HZ2_TT*24))) # ft/hour
      HZ3_Percolation <- HZ3_Excess * (1 - exp(-1/(HZ3_TT*24))) # ft/hour
      HZ4_Percolation <- HZ4_Excess * (1 - exp(-1/(HZ4_TT*24))) # ft/hour
      Deep_Percolation <- HZ5_Excess * (1 - exp(-1/(HZ5_TT*24))) # ft/hour

      #Deep_Percolation <- Average_Percolation_Rate * HZ5_FC # ft/day
      d_Consumptive_Use_Depth_dt <- ET1 + ET2 + ET3 + ET4 # d(ft)/d(time)
      d_GW_Depth_dt <- Deep_Percolation # d(ft)/d(time)
      d_Infiltration_Depth_dt <- Infiltration # d(ft)/d(time)
      d_HZ1_Water_Depth_dt <- Infiltration - ET1 - HZ1_Percolation # d(ft)/d(time)
      d_HZ2_Water_Depth_dt <- HZ1_Percolation - ET2 - HZ2_Percolation # d(ft)/d(time)
      d_HZ3_Water_Depth_dt <- HZ2_Percolation - ET3 - HZ3_Percolation# d(ft)/d(time)
      d_HZ4_Water_Depth_dt <- HZ3_Percolation - ET4 - HZ4_Percolation # d(ft)/d(time)
      d_HZ5_Water_Depth_dt <- HZ4_Percolation - Deep_Percolation # d(ft)/d(time)
      d_Ponded_Water_dt <- Inflow - Surface_Runoff - Infiltration  # d(ft)/d(time)
      
      return(list(
        c(
          d_Ponded_Water_dt, # d(ft)/d(time)
          d_HZ1_Water_Depth_dt,  # d(ft)/d(time)
          d_HZ2_Water_Depth_dt,  # d(ft)/d(time)
          d_HZ3_Water_Depth_dt,  # d(ft)/d(time)
          d_HZ4_Water_Depth_dt,  # d(ft)/d(time)
          d_HZ5_Water_Depth_dt,  # d(ft)/d(time)
          d_Consumptive_Use_Depth_dt,  # d(ft)/d(time)
          d_GW_Depth_dt, # d(ft)/d(time)
          d_Infiltration_Depth_dt# d(ft)/d(time)
          ),
        HZ1_WC = HZ1_WC, # water content of the surface layer expressed as %
        HZ2_WC = HZ2_WC, # water content of the second layer expressed as %
        HZ3_WC = HZ3_WC, # water content of the third layer expressed as %
        HZ4_WC = HZ4_WC, # water content of the fourth layer expressed as %
        HZ5_WC = HZ5_WC, # water content of the fourth layer expressed as %
        ET_Demand = ET_Demand,
        Precipitation = Precipitation,
        Irrigation = Irrigation,
        Infiltration = Infiltration,
        Surface_Runoff = Surface_Runoff,
        Deep_Percolation = Deep_Percolation,
        ET1 =ET1,
        ET2 =ET2,
        ET3 =ET3,
        ET4 =ET4,
        HZ1_Percolation=HZ1_Percolation,
        HZ2_Percolation=HZ2_Percolation,
        HZ3_Percolation=HZ3_Percolation,
        HZ4_Percolation=HZ4_Percolation,
        Ponded_Water=Ponded_Water,
        HZ1_Void =HZ1_Void
        ))
    })
  }
  
  # event function for controlling negative stock values
  positive.stocks <- function(time, state_vars, parameters){
      with(as.list(state_vars), {
        state_vars[which(state_vars < 0)] <- 0  
        return(state_vars)
      })
    }
  
  output_deSolve <- ode(y = state_vars,
                        times  = simtime,
                        func   = wb,
                        parms  = parameters, 
                        method = "rk4",
                        events=list(func = positive.stocks, time = simtime))
  
  result_df <- as_tibble(output_deSolve)
  return(result_df)
}

 
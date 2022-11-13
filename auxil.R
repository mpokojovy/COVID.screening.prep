####################################################################################################
# Copyright (C) 2022 Michael Pokojovy                                                              #
#                                                                                                  #
# Data preparation for COVID-19 patient screening using machine learning                           #                                                                                              #
####################################################################################################

str2days <- function(s) {
  if (s == "")
    return(NA)
  else {
    pos.days = as.integer(regexpr("days", s))
    pos.col  = as.integer(regexpr(":", s))
    
    days = as.numeric(substr(s, 1L, pos.days - 2L))
    hrs  = as.numeric(substr(s, pos.days + 5L, pos.col - 1L))
    mins = as.numeric(substr(s, pos.col + 1L, pos.col + 2L))
    
    return(days + hrs/24.0 + mins/(24.0*60.0))
  }
}
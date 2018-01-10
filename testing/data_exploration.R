rm(list = ls())
gc()

library(ncdf4)
library(dplyr)
library(plotly)
library(extdepth)

# open connection to TAS file
nc = nc_open('/Users/trevh/Research/extdepth/data/tas_ens_da_hydro_r.1000-2000_d.30-Nov-2017.nc')

# import and convert data
tas = ncvar_get(nc, attributes(nc$var)$names[1], count = c(-1, -1, 1, 100))
tas = lapply(seq(dim(tas)[3]), function(x) as.matrix(tas[ , , x]))

plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~tas[[1]])


# smooth surfaces

# find basis coefficients

# apply extdepth to basis coeff


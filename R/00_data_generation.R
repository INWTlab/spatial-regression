library(Kernelheaping)
library(gstat)
library(sp)

# Auswahl M1, M2, M3:
# M1: areas= FALSE, addSmooth = FALSE; 
# M2: areas= TRUE, addSmooth = FALSE;
# M3: areas= FALSE, addSmooth = TRUE.
generate_sim_data <- function(n = 500, areas = TRUE,
                              areaGridSize = 2.5,
                              mapBorder = c(-10, 10),
                              addSmooth = TRUE, 
                              y_sigma = 0.05){
  
  # Create true z coordinate values ------------------------------------------
  z_mu1 <- c(0, 0)
  z_mu2 <- c(5, 3)
  z_mu3 <- c(-4, 1)
  z_Sigma1 <- matrix(c(2, 0, 0, 2), 2, 2)
  z_Sigma2 <- matrix(c(1, 0, 0, 1), 2, 2)
  z_Sigma3 <- matrix(c(1, 0, 0, 3), 2, 2)
  
  z_mus <- rbind(z_mu1, z_mu2, z_mu3)
  z_Sigmas <- rbind(z_Sigma1, z_Sigma2, z_Sigma3)
  z_props <- c(1/3, 1/3, 1/3)
  ztrue=rmvnorm.mixt(n=n, mus=z_mus, Sigmas=z_Sigmas, props=z_props)
  
  # Create rounded z coordinate values ------------------------------------------
  roundvalue = 1.5 # choose from 0, 0.75, 1.5, 2.25
  zrounded = plyr::round_any(ztrue, roundvalue)
  
  # Create true x covariate values ------------------------------------------
  x_mu1 <- c(0, 0)
  x_mu2 <- c(5, 3)
  x_mu3 <- c(-4, 1)
  x_Sigma1 <- matrix(c(2, 0, 0, 2), 2, 2)
  x_Sigma2 <- matrix(c(1, 0, 0, 1), 2, 2)
  x_Sigma3 <- matrix(c(1, 0, 0, 3), 2, 2)
  
  x_mus <- rbind(x_mu1, x_mu2, x_mu3)
  x_Sigmas <- rbind(x_Sigma1, x_Sigma2, x_Sigma3)
  x_props <- c(1/3, 1/3, 1/3)
  
  sigma_x <- 0.1
  
  #Prof. Rendtels approach
  #xvalues <-  -5 * dmvnorm.mixt(ztrue, mus=x_mus, Sigmas=x_Sigmas, props=c(1,0,0)) +
  #  1 * dmvnorm.mixt(ztrue, mus=x_mus, Sigmas=x_Sigmas, props=c(0,1,0))+
  #  8 * dmvnorm.mixt(ztrue, mus=x_mus, Sigmas=x_Sigmas, props=c(0,0,1)) + 
  #  rnorm(n, 0, sigma_x)
  
  # generate regression data by inverse geo-weighting
  
  # Define 3 known points and z-values at the centers of the normal
  # distributions from the mixture above
  coords <- data.frame(x = c(0, 5, -4), y = c(0, 3, 1))
  z <- c(-5, 1, 8)
  data <- data.frame(coords, z = z)
  
  # Convert to spatial object
  coordinates(data) <- ~x + y
  
  # Need xvalues at the true z locations
  newdata <- data.frame(ztrue)
  colnames(newdata) <- c("x", "y")
  coordinates(newdata) <- ~x + y
  
  # Perform IDW interpolation
  # interpolate xvalues at ztrue locations with the only
  # known points specified in data
  idw_result <- idw(z ~ 1, data, newdata = newdata, idp = 2)
  
  # Convert to data frame for plotting
  xvalues <- as.data.frame(idw_result)[,3]
  
  if(addSmooth){
    # Use the same approach (IDW interpolation) 
    # as for the xvalues to generate a smooth spatial effect
    coords <- data.frame(x = c(-2, 8, -5, 9, -5), y = c(1, 4, 3, -1, -3))
    z <- c(3, -7, 11, 0.9, -3)
    data <- data.frame(coords, z = z)
    
    # Convert to spatial object
    coordinates(data) <- ~x + y
    
    newdata <- data.frame(ztrue)
    colnames(newdata) <- c("x", "y")
    coordinates(newdata) <- ~x + y
    
    # Perform IDW interpolation
    idw_result <- idw(z ~ 1, data, newdata = newdata, idp = 2)
    
    # Convert to data frame for plotting
    smoothvalues <- as.data.frame(idw_result)[,3]
  } else {
    # if no smooth effect, set to 0
    smoothvalues <- rep(0, n)
  }
  
  # Generate area values
  nAreas <- ceiling((mapBorder[2] -  mapBorder[1]) / areaGridSize) * ceiling((mapBorder[2] -  mapBorder[1]) / areaGridSize)
  if(areas == TRUE){
    u <- rnorm(nAreas, sd = 5)
  } else {
    # if no area effect, set u to 0
    u <- rep(0, nAreas)
  }
  # Assign area ids and u values to coordinates
  areaAssignment <- map_to_grid_with_u(ztrue[,1], ztrue[,2], u_values = u, 
                                       grid_size = areaGridSize,
                                       min_val = mapBorder[1],
                                       max_val = mapBorder[2])
  uvalues <- areaAssignment$u
  area_id <- areaAssignment$area_id
  
  #match areas with coordinates
  
  # Create y response values ------------------------------------------
  # uvalues and smoothvalues is zero if area/smooth effect is set to False
  # y_sigma <- 0.05
  yvalues <- 0 + 1*xvalues + uvalues + smoothvalues + rnorm(n, 0, y_sigma)
  
  return(list(zrounded = zrounded, roundvalue = roundvalue,
              yvalues = yvalues, xvalues = xvalues, areas = areas,
              areaGridSize = areaGridSize, mapBorder = mapBorder,
              addSmooth = addSmooth,u = u,
              z_mus=z_mus, z_Sigmas=z_Sigmas, z_props=z_props,
              uvalues = uvalues, smoothvalues = smoothvalues, ztrue = ztrue))
}

map_to_grid_with_u <- function(x, y, u_values, grid_size = 2.5, min_val = -10, max_val = 10) {
  # helper function to assign area id and u-values to coordinates
  # Number of cells per axis
  n_cells <- ((max_val - min_val) / grid_size)
  
  # Compute indices
  col_index <- floor((x - min_val) / grid_size) + 1
  row_index <- floor((y - min_val) / grid_size) + 1
  
  # Compute rectangle boundaries
  x_min <- min_val + (col_index - 1) * grid_size
  x_max <- x_min + grid_size
  y_min <- min_val + (row_index - 1) * grid_size
  y_max <- y_min + grid_size
  
  # Compute cell ID (row-major order)
  area_id <- (row_index - 1) * n_cells + col_index
  
  # Assign z-values based on cell_id
  if(!is.null(u_values)){
    u_assigned <- u_values[area_id]
  } else {
    # if no u-values are provided, set to 0
    u_assigned <- rep(0, length(area_id))
  }
  
  # Return as data frame
  data.frame(
    x = x, y = y,
    row = row_index, col = col_index,
    area_id = factor(area_id, levels = as.character(1:n_cells**2)),
    u = u_assigned,
    x_min = x_min, x_max = x_max,
    y_min = y_min, y_max = y_max
  )
}

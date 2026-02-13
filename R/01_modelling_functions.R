library(mgcv)
# Achtung: Variablenbezeichnung folgt den Bezeichnungen aus KernelHeaping Paket.
# Die Funktion dbivr_X erweitert die dbivr Funktion aus dem KernelHeaping Paket.
# Diese findet sich auch weiter unten in diesem Skript, nach der dbivr_X Funktion.

# Auswahl M1, M2, M3:
# M1: areas= FALSE, addSmooth = FALSE; 
# M2: areas= TRUE, addSmooth = FALSE;
# M3: areas= FALSE, addSmooth = TRUE.
dbivr_X <- function (xrounded, roundvalue, xvalues, yvalues, burnin = 2, samples = 5, 
          gridsize = 200, areas = TRUE, areaGridSize = 2.5, mapBorder = c(-10, 10), addSmooth = TRUE)
{
  gridx = seq(min(xrounded[, 1]) - 0.5 * roundvalue, max(xrounded[, 1]) + 0.5 * roundvalue, length = gridsize)
  gridy = seq(min(xrounded[, 2]) - 0.5 * roundvalue, max(xrounded[, 2]) + 0.5 * roundvalue, length = gridsize)

  # Initial density estimate based on rounded coordinates -----
  Mestimates <- ks::kde(x = xrounded, H = diag(c(roundvalue, 
                                                 roundvalue))^2, gridsize = c(length(gridx), length(gridy)), 
                        xmin = c(min(gridx), min(gridy)), xmax = c(max(gridx), 
                                                                   max(gridy)))
  resultDensity = array(dim = c(burnin + samples, length(gridx), 
                                length(gridy)))
  resultX = array(dim = c(samples + burnin, nrow(xrounded), 
                          2))
  reg_models = vector("list", samples + burnin)
  rvalues = unique(xrounded)
  result_m_reg <- list()

  # Find all density eval points around a unique rounded coordinate value
  selectionGrid <- lapply(1:nrow(rvalues), function(k) {
    selectionX = which(Mestimates$eval.points[[1]] >= rvalues[k, 1] - roundvalue * 0.5 & Mestimates$eval.points[[1]] < 
                         rvalues[k, 1] + roundvalue * 0.5)
    selectionY = which(Mestimates$eval.points[[2]] >= rvalues[k, 2] - roundvalue * 0.5 & Mestimates$eval.points[[2]] < 
                         rvalues[k, 2] + roundvalue * 0.5)
    list(selectionX, selectionY)
  })
  
  # Initial estimate of the regression model m_reg and m_y_smooth -----
  #Step 1+2: Get values of y on grid by some smoothing routine:
  if(areas){
    area_id <- map_to_grid_with_u(xrounded[,1], xrounded[,2], u_values = NULL, 
                                  grid_size = areaGridSize,
                                  min_val = mapBorder[1],
                                  max_val = mapBorder[2])$area_id
    dat_reg = data.frame(xvalues = xvalues, yvalues = yvalues, area_id = area_id, x = xrounded[,1], y = xrounded[,2])
    #an area random effect is not included here yet as it severely worsens results (likely due to only a few areas attributable to the rounded data)
    if(addSmooth){
      m_reg <- mgcv::gam(yvalues ~ xvalues  + s(x,y), data = dat_reg)
    } else {
      m_reg <- mgcv::gam(yvalues ~ xvalues , data = dat_reg)
    }
    dat_tmp <- data.frame(x = xrounded[,1], y = xrounded[,2], z = yvalues, area_id = area_id)
    m_y_smooth <- mgcv::gam(z~s(x,y) , data = dat_tmp)
  } else {
    dat_reg = data.frame(xvalues = xvalues, yvalues = yvalues, x = xrounded[,1], y = xrounded[,2])
    if(addSmooth){
      m_reg <- mgcv::gam(yvalues ~ xvalues + s(x,y), data = dat_reg)
    } else {
      m_reg <- mgcv::gam(yvalues ~ xvalues, data = dat_reg)
    }
    dat_tmp <- data.frame(x = xrounded[,1], y = xrounded[,2], z = yvalues)
    m_y_smooth <- mgcv::gam(z~s(x,y), data = dat_tmp)
  }
  # Start sampling loop ---------
  for (j in 1:(burnin + samples)) {
    new = c()
    indicators_all = c()
    for (i in 1:nrow(rvalues)) {
      ## get the density values around the rounded coordinate --------
      probs = as.vector(Mestimates$estimate[selectionGrid[[i]][[1]], 
                                            selectionGrid[[i]][[2]]])
      # get the candidate coordinate points around the rounded coordinate ------
      points = cbind(rep(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]], 
                         times = length(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]])), 
                     rep(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]], 
                         each = length(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]])))
      
      # indicators determines observations that belong to a specific rounded value -----
      indicators = which(xrounded[, 1] == rvalues[i, 1] & xrounded[, 2] == rvalues[i, 2])
      
      if(areas){
        # assign area ids to points ----
        area_id <- map_to_grid_with_u(points[,1], points[,2], u_values = NULL, 
                                      grid_size = areaGridSize,
                                      min_val = mapBorder[1],
                                      max_val = mapBorder[2])$area_id
        data_pred <- data.frame(x = points[,1], y = points[,2], area_id = area_id)
      } else {
        data_pred <- data.frame(x = points[,1], y = points[,2])
      }
      ## predict y values on candidate points -----------
      smoothed_y <- predict(m_y_smooth, newdata = data_pred, newdata.guaranteed=TRUE)
      
      for (k in 1:length(indicators)) {
        # weight the sampling probabilities based on m_y_smooth predictions and true y -----
        probs_combined <- dnorm((yvalues[indicators[k]] - smoothed_y), 0, sd = sqrt(m_y_smooth$sig2)) * probs
        # probs_combined <- (abs(yvalues[indicators[k]] - smoothed_y) <= (1 * sqrt(m_reg$sig2))) * probs
        new = rbind(new, points[sample(1:nrow(points), size = 1, replace = T, prob = probs_combined), ])
        indicators_all <- rbind(indicators_all, indicators[k])
      }
    }
    # Estimation of m_reg and m_y_smooth with newly generated coordinates -----
    if(areas){
      area_id <- map_to_grid_with_u(new[,1], new[,2], u_values = NULL, 
                                    grid_size = areaGridSize,
                                    min_val = mapBorder[1],
                                    max_val = mapBorder[2])$area_id
      dat_reg = data.frame(xvalues = xvalues[indicators_all], yvalues = yvalues[indicators_all], area_id = area_id, x = new[,1], y = new[,2])
      if(addSmooth){
        m_reg <- mgcv::gam(yvalues ~ xvalues + s(area_id, bs = "re") + s(x,y), data = dat_reg)
      } else {
        m_reg <- mgcv::gam(yvalues ~ xvalues + s(area_id, bs = "re"), data = dat_reg)
      }
      dat_tmp <- data.frame(x = new[,1], y = new[,2], z = yvalues[indicators_all], area_id = area_id)
      m_y_smooth <- mgcv::gam(z~s(x,y) + s(area_id, bs = "re"), data = dat_tmp)
    } else{
      dat_reg = data.frame(xvalues = xvalues[indicators_all], yvalues = yvalues[indicators_all], x = new[,1], y = new[,2])
      if(addSmooth){
        m_reg <- mgcv::gam(yvalues ~ xvalues + s(x,y), data = dat_reg)
      } else {
        m_reg <- mgcv::gam(yvalues ~ xvalues, data = dat_reg)
      }
      dat_tmp <- data.frame(x = new[,1], y = new[,2], z = yvalues[indicators_all])
      m_y_smooth <- mgcv::gam(z~s(x,y), data = dat_tmp)
    }
    
    H <- ks::Hpi(x = new, binned = TRUE)
    
    # new sind die neu generierten z-Koordinaten
    # Density estimate based on newly sampled coordinates --------
    Mestimates <- ks::kde(x = new, H = H, gridsize = c(length(gridx), 
                                                         length(gridy)), bgridsize = c(length(gridx), 
                                                                                       length(gridy)), xmin = c(min(gridx), min(gridy)), 
                            xmax = c(max(gridx), max(gridy)), binned = TRUE)
    Mestimates$estimate[is.na(Mestimates$estimate)] = 1e-96
    resultDensity[j, , ] = Mestimates$estimate
    reg_models[[j]] = m_reg
    resultX[j, , ] = new
    result_m_reg[[j]] <- m_reg[c("coefficients", "sig2")]
    print(paste("Iteration:", j, "of", burnin + samples))
    # End of Sampling Loop -----------
  }
  # Collect and return results ------
  Mestimates$estimate = apply(resultDensity[-c(1:burnin), , 
  ], c(2, 3), mean)
  
  result_m_reg_xvalues_samples <- sapply(result_m_reg, function(x){x$coefficients["xvalues"]})
  result_m_reg_xvalues <- mean(result_m_reg_xvalues_samples[-c(1:burnin)])
  
  est <- list(Mestimates = Mestimates, resultDensity = resultDensity, 
              resultX = resultX, xrounded = xrounded, gridx = gridx, 
              gridy = gridy, roundvalue = roundvalue, burnin = burnin, 
              samples = samples, adaptive = TRUE, reg_models = reg_models)
  class(est) <- "bivrounding"
  return(est)
  # End of function -----------
}


# Old method without regression model, only density estimation based on rounded coordinates ----

dbivr <- function (xrounded, roundvalue, burnin = 2, samples = 5, 
                            gridsize = 200) 
{
  gridx = seq(min(xrounded[, 1]) - 0.5 * roundvalue, max(xrounded[, 1]) + 0.5 * roundvalue, length = gridsize)
  gridy = seq(min(xrounded[, 2]) - 0.5 * roundvalue, max(xrounded[, 2]) + 0.5 * roundvalue, length = gridsize)
  Mestimates <- ks::kde(x = xrounded, H = diag(c(roundvalue, 
                                                 roundvalue))^2, gridsize = c(length(gridx), length(gridy)), 
                        xmin = c(min(gridx), min(gridy)), xmax = c(max(gridx), 
                                                                   max(gridy)))
  resultDensity = array(dim = c(burnin + samples, length(gridx), 
                                length(gridy)))
  resultX = array(dim = c(samples + burnin, nrow(xrounded), 
                          2))
  rvalues = unique(xrounded)
  selectionGrid <- lapply(1:nrow(rvalues), function(k) {
    selectionX = which(Mestimates$eval.points[[1]] >= rvalues[k, 1] - roundvalue * 0.5 & Mestimates$eval.points[[1]] < 
                         rvalues[k, 1] + roundvalue * 0.5)
    selectionY = which(Mestimates$eval.points[[2]] >= rvalues[k, 2] - roundvalue * 0.5 & Mestimates$eval.points[[2]] < 
                         rvalues[k, 2] + roundvalue * 0.5)
    list(selectionX, selectionY)
  })
  
  for (j in 1:(burnin + samples)) {
    new = c()
    for (i in 1:nrow(rvalues)) {
      probs = as.vector(Mestimates$estimate[selectionGrid[[i]][[1]], 
                                            selectionGrid[[i]][[2]]])
      points = cbind(rep(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]], 
                         times = length(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]])), 
                     rep(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]], 
                         each = length(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]])))
      npoints = length(which(xrounded[, 1] == rvalues[i, 
                                                      1] & xrounded[, 2] == rvalues[i, 2]))
      new = rbind(new, points[sample(1:nrow(points), size = npoints, 
                                     replace = T, prob = probs), ])
    }
    H <- ks::Hpi(x = new, binned = TRUE)
    
    Mestimates <- ks::kde(x = new, H = H, gridsize = c(length(gridx), 
                                                       length(gridy)), bgridsize = c(length(gridx), 
                                                                                     length(gridy)), xmin = c(min(gridx), min(gridy)), 
                          xmax = c(max(gridx), max(gridy)), binned = TRUE)
    Mestimates$estimate[is.na(Mestimates$estimate)] = 1e-96
    resultDensity[j, , ] = Mestimates$estimate
    resultX[j, , ] = new
    print(paste("Iteration:", j, "of", burnin + samples))
  }
  Mestimates$estimate = apply(resultDensity[-c(1:burnin), , 
  ], c(2, 3), mean)
  est <- list(Mestimates = Mestimates, resultDensity = resultDensity, 
              resultX = resultX, xrounded = xrounded, gridx = gridx, 
              gridy = gridy, roundvalue = roundvalue, burnin = burnin, 
              samples = samples, adaptive = TRUE)
  class(est) <- "bivrounding"
  return(est)
}
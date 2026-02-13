run_sim <- function(NSIM, n, y_sigma = 0.05, areas = F, addSmooth = F){
  new_method <- rep(0,NSIM)
  old_method <- rep(0,NSIM)
  
  new_method_RE <- rep(0,NSIM)
  old_method_RE <- rep(0,NSIM)
  new_method_FX <- rep(0,NSIM)
  old_method_FX <- rep(0,NSIM)
  
  for (l in 1:NSIM){
    # Daten generieren
    
    # Auswahl M1, M2, M3:
    # M1: areas= FALSE, addSmooth = FALSE; 
    # M2: areas= TRUE, addSmooth = FALSE;
    # M3: areas= FALSE, addSmooth = TRUE.
    sim_data <- generate_sim_data(n, areas = areas, addSmooth = addSmooth, y_sigma = y_sigma)
    
    # Neues Modell fitten
    burnin <- 5
    # Auswahl von M1, M2, M3 wird weitergegeben über sim_data Objekt
    res <- dbivr_X(xrounded = sim_data$zrounded, roundvalue = sim_data$roundvalue, yvalues = sim_data$yvalues,
                   xvalues = sim_data$xvalues,burnin = burnin, samples = 10, gridsize = 200,
                   areas = sim_data$areas, areaGridSize = sim_data$areaGridSize,
                   mapBorder = sim_data$mapBorder, addSmooth = sim_data$addSmooth)
    
    # Altes Modell fitten
    res2 <- dbivr(sim_data$zrounded, sim_data$roundvalue, burnin = burnin, samples = 10, gridsize = 200)
    #evaluation:
    # Wahre Dichte erstellen
    dens=dmvnorm.mixt(x=expand.grid(res$Mestimates$eval.points[[1]],res$Mestimates$eval.points[[2]]),
                      mus=sim_data$z_mus, Sigmas=sim_data$z_Sigmas, props=sim_data$z_props)
    # In gleiches Format umwandeln wie aus dbivr_X
    dens=matrix(dens,nrow=length(res$gridx),ncol=length(res$gridy))
    
    # RMSE berechnen für jeden Dichte-Evaluationspunkt
    new_method[l] <- mean((c((dens)-(res$Mestimates$estimate)))**2)**0.5
    old_method[l] <- mean((c((dens)-(res2$Mestimates$estimate)))**2)**0.5
    
    #compare RE if M2
    models <- res$reg_models[-c(1:burnin)]
    fittedRE <- rowMeans(sapply(1:length(models), function(x) predict(models[[x]], data.frame(xvalues = 0, area_id = 1:64))-coefficients(models[[x]])[1]))
    area_id <- map_to_grid_with_u(sim_data$zrounded[,1], sim_data$zrounded[,2], u_values = NULL, 
                                  grid_size = sim_data$areaGridSize,
                                  min_val = sim_data$mapBorder[1],
                                  max_val = sim_data$mapBorder[2])$area_id
    dat_reg = data.frame(xvalues = sim_data$xvalues, yvalues = sim_data$yvalues, area_id = area_id, x = sim_data$zrounded[,1], y = sim_data$zrounded[,2])
    m_reg <- mgcv::gam(yvalues ~ xvalues + s(area_id, bs = "re"), data = dat_reg)
    fitted_RE_OLD <- predict(m_reg, data.frame(xvalues = 0, area_id = 1:64))-coefficients(m_reg)[1]
    new_method_RE[l] <- sqrt(mean((sim_data$u[fittedRE!=0] - fittedRE[fittedRE!=0] )^2))
    old_method_RE[l] <- sqrt(mean((sim_data$u[fittedRE!=0] - fitted_RE_OLD[fittedRE!=0] )^2))
    new_method_FX[l] <- mean(sapply(1:length(models), function(y) coefficients(models[[y]])[2]))
    old_method_FX[l] <- coefficients(m_reg)[2]
  }
  
  result = cbind(new_method, old_method)
  result_FX = cbind(new_method_FX, old_method_FX)
  result_RE = cbind(new_method_RE, old_method_RE)
  
  return(list(result = result, result_FX = result_FX, result_RE = result_RE))
}
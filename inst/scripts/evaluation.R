# Setze Arbeitsverzeichnis zum Projektordner
setwd("~/Github/spatial-regression")

# Lade Funktionen
source("R/00_data_generation.R")
source("R/01_modelling_functions.R")
source("R/02_result_tables.R")
source("R/03_run_simulation.R")

output_dir <- "sim_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Simulation mit einem n und y_sigma ----
set.seed(98)

# Auswahl M1, M2, M3 für das Regressionsmodell in der run_sim Funktion:
# M1: areas= FALSE, addSmooth = FALSE; 
# M2: areas= TRUE, addSmooth = FALSE;
# M3: areas= FALSE, addSmooth = TRUE.
# Per default ist M1 eingestellt

results <- run_sim(NSIM = 10, n = 200, y_sigma = 0.05, areas = F, addSmooth = F)
res <- results$result
summary(res)

# Schätzung der Regressionkoeffizienten ----
res_FX <- results$result_FX
summary(res_FX)
# berechne RMSE wenn wahrer Wert immer 1 ist
sqrt(mean((res_FX[,1] - 1)^2)) # neu
sqrt(mean((res_FX[,2] - 1)^2)) # alt

res_RE <- results$result_RE
summary(res_RE)

# Verschiedene n ausprobieren -----
# Nutze diesen Abschnitt, um die Simulation mit verschiedenen n laufen zu lassen und die Ergebnisse zu speichern.
NSIM <- 20
n_list = c(500, 600, 700, 800, 900, 1000)

results_list <- list()

for (n in n_list){
  res <- run_sim(NSIM, n, areas = F, addSmooth = F)
  print(paste0("n=", n))
  results_list[[as.character(n)]] <- res
}

save(results_list, file = "sim_results/evaluation_results.RData")

print(mean_results(results_list, n_list, "n"))
print(median_results(results_list, n_list, "n"))

load("sim_results/evaluation_results.RData")


# Residuenvarianz variieren ----
# Nutze diesen Abschnitt, um die Simulation mit verschiedenen y_sigma laufen zu lassen und die Ergebnisse zu speichern.
y_sigma_list = c(0.01, 0.05, 0.5)

results_var_list <- list()
for (y_sigma in y_sigma_list){
  res <- run_sim(NSIM = 20, n = 500, y_sigma = y_sigma, areas = F, addSmooth = F)
  print(paste0("y_sigma=", y_sigma))
  results_var_list[[as.character(y_sigma)]] <- res
}

save(results_var_list, file = "sim_results/evaluation_var_results.RData")

print(mean_results(results_var_list, n_list, "y_sigma"))
print(median_results(results_var_list, n_list, "y_sigma"))

load("sim_results/evaluation_var_results.RData")
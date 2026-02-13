# Notizen: ggplot kann nur mit den zrounded umgehen nicht mit den ztrue
# Dadurch werden manche Werte gedropped, weil die Koordinatenwerte doppelt sind

source("R/00_data_generation.R")
source("R/01_modelling_functions.R")
set.seed(98)
n=500
y_sigma=0.05
# Daten generieren

# Auswahl M1, M2, M3:
# M1: areas= FALSE, addSmooth = FALSE; 
# M2: areas= TRUE, addSmooth = FALSE;
# M3: areas= FALSE, addSmooth = TRUE.
sim_data <- generate_sim_data(n, areas = T, addSmooth = T, y_sigma = y_sigma)

z1 <- sim_data$zrounded[,1]
z2 <- sim_data$zrounded[,2]

df <- data.frame(q = sim_data$xvalues, z1 = z1, z2 = z2, smooth = sim_data$smoothvalues)

# Load the library
library(ggplot2)

# Create the plot
ggplot(df, aes(x = z1, y = z2, z = q)) +
  # 1. Fill the background based on value q
  geom_contour_filled() +
  # 2. Add the specific contour lines (Höhenlinien)
  geom_contour(color = "blue", linewidth = 0.3) +
  # 3. Use a nice color palette
  scale_fill_viridis_d() + 
  theme_minimal() +
  labs(fill = "Value (xvalues)")

ggplot(df, aes(x = z1, y = z2, z = smooth)) +
  # 1. Fill the background based on value q
  geom_contour_filled() +
  # 2. Add the specific contour lines (Höhenlinien)
  geom_contour(color = "blue", linewidth = 0.3) +
  # 3. Use a nice color palette
  scale_fill_viridis_d() + 
  theme_minimal() +
  labs(fill = "Value (smoothvalues)")


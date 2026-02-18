# Extrahiere nur den Mittelwert der RMSEs und speichere in einer Tabelle
# Und runden auf 5 Stellen nach dem Komma
mean_results <- function(results_list, n_list, name, which_table = "result"){dat <- data.frame(
  New_Method_Mean_RMSE = sapply(n_list, function(n) round(mean(results_list[[as.character(n)]][[which_table]][,1]), 5)),
  Old_Method_Mean_RMSE = sapply(n_list, function(n) round(mean(results_list[[as.character(n)]][[which_table]][,2]), 5)))
dat["differenz"] <- dat$Old_Method_Mean_RMSE - dat$New_Method_Mean_RMSE
dat["pct_diff"] <- round(100 * dat$differenz / dat$Old_Method_Mean_RMSE, 2)
dat[name] <- n_list
return(dat)
}

# Das gleiche für den Median
median_results <- function(results_list, n_list, name, which_table = "result"){dat <- data.frame(
  New_Method_Median_RMSE = sapply(n_list, function(n) round(median(results_list[[as.character(n)]][[which_table]][,1]), 5)),
  Old_Method_Median_RMSE = sapply(n_list, function(n) round(median(results_list[[as.character(n)]][[which_table]][,2]), 5)))
dat["differenz"] <- dat$Old_Method_Median_RMSE - dat$New_Method_Median_RMSE
dat["pct_diff"] <- round(100 * dat$differenz / dat$Old_Method_Median_RMSE, 2)
dat[name] <- n_list
return(dat)
}
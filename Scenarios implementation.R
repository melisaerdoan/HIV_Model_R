# Here dynamic number of scenarios created and plotted.
# IMPORTANT :lines between 113 and 116 should be manually determined. 
library(deSolve)
library(MASS)
library(mvtnorm)

source("run_model_HIV2.R")

beginning = 0
intermediate = 4
ending = 15

run_model <- function(scenario_number, scenario_change_in_VCT, scenario_change_in_others) {   

  initial_values1 <-  c(SH=0.02643, SL=0.97295, SVCT= 0, SO=0,  U=0.0003225, IVCT=1.67*10^-6, IO= 7.318*10^-5, ART=0.00022)
  
  beta_1 = 0.5247                    # peak of the posterior N=500
  beta_2 = 0.2339                    # peak of the posterior N=500
  
  # Base case Part 1
  
  params1 <- c( beta1 = beta_1, beta2 = beta_2, sigma1 = 0.00206, sigma2 = 0.13241, delta1 = 0.00434, delta2 = 0.19086, epsilon1 = 0.53
               , epsilon2 = 0.96,  omega = 2, psi = 0.219,  alpha = 0.747)
  times1 <- seq(from = beginning, to = intermediate, by=1)
  
  total_numberof_people = 62819553
  
  D1 <- ode(
    y = initial_values1,
    times = times1,
    func = closed.sir.model,
    parms = params1
  )  
  
  colnames(D1)[10] <- "incidence"
  
  prevalence <- matrix(0, nrow=nrow(D1))
  for (i in 1:nrow(D1)){
    prevalence[i] = sum(D1[i,6:9])
  }
  D1 <- cbind(D1, prevalence)
  colnames(D1)[11] <- "prevalence"

  
  diagnosed_cases <- c() 
  diagnosed_cases <- (params1["delta1"] + params1["delta2"]) * D1[,6] * total_numberof_people   # Do not need to consider the change for the parameter because it is in the formula
  D1 <- cbind(D1, diagnosed_cases)
  colnames(D1)[12] <- "diagnosed_cases"
  
  
  applied_test_rates_VCT <- c() 
  applied_test_rates_VCT <- (params1["sigma1"] * D1[,2]) + (params1["delta1"] * D1[,6])   # S_High and U compartments
  D1 <- cbind(D1, applied_test_rates_VCT)
  colnames(D1)[13] <- "applied_tests_VCT"
  
  applied_test_rates_others <- c() 
  applied_test_rates_others <- (params1["sigma2"] * D1[,3]) + (params1["delta2"] * D1[,6])   # S_Low and U compartments
  D1 <- cbind(D1, applied_test_rates_others)
  colnames(D1)[14] <- "applied_tests_others"
  
  
  # Change in the testing & diagnosis - Part 2
  
  change_in_VCT = scenario_change_in_VCT      
  change_in_others = scenario_change_in_others   
  
  initial_values2 <- D1[(intermediate+1), 2:9]  # Number of columns should be compatible with inputs/outputs of ode
  
  params2 <- c( beta1 = beta_1, beta2 = beta_2, sigma1 = 0.00206*change_in_VCT, sigma2 = 0.13241*change_in_others, delta1 = 0.00434*change_in_VCT, delta2 = 0.19086*change_in_others, epsilon1 = 0.53
               , epsilon2 = 0.96,  omega = 2, psi = 0.219,   alpha = 0.747)
  times2 <- seq(from = intermediate, to = ending, by=1)
  
  total_numberof_people = 62819553
  
  D2 <- ode(
    y = initial_values2,
    times = times2,
    func = closed.sir.model,
    parms = params2
  )  
  
  colnames(D2)[10] <- "incidence"
  
  prevalence <- matrix(0,nrow=nrow(D2))
  for (i in 1:nrow(D2)){
    prevalence[i] = sum(D2[i,6:9])
  }
  D2 <- cbind(D2, prevalence)
  colnames(D2)[11] <- "prevalence"

  
  diagnosed_cases <- c() 
  diagnosed_cases <- (params2["delta1"] + params2["delta2"]) * D2[,6] * total_numberof_people   # Do not need to consider the change in the parameter because it is in the formula
  D2 <- cbind(D2, diagnosed_cases)
  colnames(D2)[12] <- "diagnosed_cases"
  
  
  applied_test_rates_VCT <- c() 
  applied_test_rates_VCT <- (params2["sigma1"] * D2[,2]) + (params2["delta1"] * D2[,6])   # S_High and U compartments
  D2 <- cbind(D2, applied_test_rates_VCT)
  colnames(D2)[13] <- "applied_tests_VCT"
  
  applied_test_rates_others <- c() 
  applied_test_rates_others <- (params2["sigma2"] * D2[,3]) + (params2["delta2"] * D2[,6])   # S_Low and U compartments
  D2 <- cbind(D2, applied_test_rates_others)
  colnames(D2)[14] <- "applied_tests_others"
  
  
  D <- rbind(D1,D2[(2:(ending-intermediate+1)),])
  return (D)

}

scenario_count <- 5 
main_scenario <- 2
scenario_change_in_VCT <- c(1, 2, 4, 2.6, 3.8)       # 5-2 c(1, 2, 4, 2.6, 3.8)  # 5-3 c(1, 2, 4, 1, 1)        #3-4 c(1,2,4)
scenario_change_in_others <- c(1, 1, 1, 1, 1)                # 5-2 c(1, 1, 1, 1, 1)        # 5-3 c(1, 1, 1, 2, 4)  # 3-4 c(1,2,4)
scenario_results <- list()
prevented_incidence_per_increased_test_VCT <- c()
prevented_incidence_per_increased_test_others <- c()
all_results <- c()
col <- c( "red", "blue", "orange","limegreen", "gray" , "brown", "navyblue")
color <- col[1:scenario_count]
name_legend = c("Base")
par(mar = c(5,4,4,5)+1 , mai= c(0.6, 0.6, 0.5, 0.3))
nf <- layout( matrix(rbind(c(1, 2, 3),  c(4, 5, 6), c(7, 8, 9), c(10, 11, 12)), nrow = 4, ncol=3) )
#layout.show(nf)

#for (i in 1:(scenario_count-1)) {
#  name_legend <- c(name_legend, paste0("Scenario ", main_scenario, ".", i))
#}

if (main_scenario == 2) {
  name_legend = c("No change", "2-fold increase (VCT)", "4-fold increase (VCT)", "2.6-fold increase (VCT)", "3.8-fold increase (VCT)")
} 

if (main_scenario == 3) {
  name_legend = c("No change", "2-fold increase (VCT)", "4-fold increase (VCT)", "2-fold increase (others)", "4-fold increase (others)" )
} 


if (main_scenario == 4) {
  name_legend = c("No change", "2-fold increase (Both)", "4-fold increase (Both)")
} 


# Model run is applied here
for (i in 1:scenario_count) {
  scenario_results[[i]] <- run_model(i, scenario_change_in_VCT[i], scenario_change_in_others[i])
  all_results <- rbind(all_results, scenario_results[[i]])
}

# This list includes only the ones that you want to graph
result_name_list <- list("time", "SH", "SL", "SVCT", "SO", "U", "IVCT", "IO", "ART", "incidence", "prevalence", "diagnosed_cases")
fold_list <- list(1,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1)
names(fold_list) <- result_name_list 
# data without fold_list to calculate prevented per increased test
data <- data.frame(time = c(scenario_results[[1]][, 1]), SH = c(scenario_results[[1]][, 2]), 
                   SL = c(scenario_results[[1]][, 3]), SVCT = c(scenario_results[[1]][, 4]), 
                   SO = c(scenario_results[[1]][, 5]), U = c(scenario_results[[1]][, 6]), 
                   IVCT = c(scenario_results[[1]][, 7]), IO = c(scenario_results[[1]][, 8]), 
                   ART = c(scenario_results[[1]][, 9]), incidence = c(scenario_results[[1]][, 10]), 
                   prevalence = c(scenario_results[[1]][, 11]), diagnosed_cases = c(scenario_results[[1]][, 12]),
                   applied_tests_VCT = c(scenario_results[[1]][, 13]), applied_tests_others = c(scenario_results[[1]][, 14])
                   )
#data_for_plot is with fold_list to use only for  plots
data_for_base_plot <- data.frame(time = c(scenario_results[[1]][, 1])*fold_list[[1]], SH = c(scenario_results[[1]][, 2])*fold_list[[2]], 
                   SL = c(scenario_results[[1]][, 3])*fold_list[[3]], SVCT = c(scenario_results[[1]][, 4])*fold_list[[4]], 
                   SO = c(scenario_results[[1]][, 5])*fold_list[[5]], U = c(scenario_results[[1]][, 6])*fold_list[[6]], 
                   IVCT = c(scenario_results[[1]][, 7])*fold_list[[7]], IO = c(scenario_results[[1]][, 8])*fold_list[[8]], 
                   ART = c(scenario_results[[1]][, 9])*fold_list[[9]], incidence = c(scenario_results[[1]][, 10])*fold_list[[10]], 
                   prevalence = c(scenario_results[[1]][, 11])*fold_list[[11]], diagnosed_cases = c(scenario_results[[1]][, 12])*fold_list[[12]],
                   applied_tests_VCT = c(scenario_results[[1]][, 13]), applied_tests_others = c(scenario_results[[1]][, 14])
)
for (i in 1:scenario_count) {
# Prevented incidence & increased tests calculation
  column_of_incidence <- 10
  prevented <- scenario_results[[i]][, column_of_incidence] - data[["incidence"]]   # Related scenario results - Base results
  scenario_results[[i]] <- cbind(scenario_results[[i]], prevented)
  colnames(scenario_results[[i]])[15] <- "prevented_incidence"
  
  column_of_applied_test_rates_VCT <- 13
  increased_test_rates_VCT <- scenario_results[[i]][, column_of_applied_test_rates_VCT] - data[["applied_tests_VCT"]] # Related scenario results - Base results
  scenario_results[[i]] <- cbind(scenario_results[[i]], increased_test_rates_VCT)
  colnames(scenario_results[[i]])[16] <- "increased_tests_VCT"
  
  column_of_applied_test_rates_others <- 14
  increased_test_rates_others <- scenario_results[[i]][, column_of_applied_test_rates_others] - data[["applied_tests_others"]]
  scenario_results[[i]] <- cbind(scenario_results[[i]], increased_test_rates_others)
  colnames(scenario_results[[i]])[17] <- "increased_tests_others"
  
  if (main_scenario == 3) {
    prevented_incidence_per_increased_test_VCT[i] <- sum(scenario_results[[i]][,15]) / sum(scenario_results[[i]][6:16,16])
    prevented_incidence_per_increased_test_others[i] <- sum(scenario_results[[i]][,15]) / sum(scenario_results[[i]][6:16,17])
    prevented_incidence_per_increased_test <- cbind(prevented_incidence_per_increased_test_VCT, prevented_incidence_per_increased_test_others)
    write.csv(x = prevented_incidence_per_increased_test, file = paste0(main_scenario,".", "prevented_per_increased_tests", ".csv"))
    }  
}

# Plots the results of the base case. It will take the values for first and second part. 
# When you want to see all the results in one page run the following:
#
pdf(file = paste0("My_Plot_", main_scenario, ".pdf" ),  width = 4, height = 4) # Exporting the plot directly to a pdf but do not forget dev.off !
par(mfrow = c(1,1))

for (s in 2:length(result_name_list)) {     # s starts from 2 because there is time in the first column and SH starts from the second column
  # plot for base is created (with fold_list) 
  if (s == 12) {   # The y axis label for diagnosed cases is different
    plot(c((1+2017):(ending+1+2017)), data_for_base_plot[[result_name_list[[s]]]], 
         xlim = c((1+2017),(ending+1+2017)), 
         ylim = c(min(all_results[,s]*fold_list[[s]]), (max(all_results[,s])+sd(all_results[,s]))*fold_list[[s]]), 
         type = "l", main = paste("Results -", result_name_list[[s]]) , 
         xlab = "Time", ylab = paste(""),  
         xaxt='n', col= col[1], las=2 )
    title(ylab="Cases", line=3.3, cex.lab=1)
  } else {
    plot(c((1+2017):(ending+1+2017)), data_for_base_plot[[result_name_list[[s]]]], 
          xlim = c((1+2017),(ending+1+2017)), 
          ylim = c(min(all_results[,s]*fold_list[[s]]), (max(all_results[,s])+sd(all_results[,s]))*fold_list[[s]]), 
          type = "l", main = paste("Results -", result_name_list[s]) , 
          xlab = "Time", ylab = paste("Rates per", fold_list[[s]]),  
          xaxt='n', col= col[1], las=2 )
  }
  # plot for scenarios other than base is created (with fold_list) 
      for (i in (2:scenario_count)){
        lines(c((1+intermediate+2017):(ending+1+2017)), scenario_results[[i]][(intermediate+1):(ending+1), s]*fold_list[[s]], 
              xlim = c((1+intermediate+2017),(ending+1+2017)), 
              ylim = c(min(all_results[,s])*fold_list[[s]], (max(all_results[,s])+sd(all_results[,s]))*fold_list[[s]]), 
              type = "l", main = "Results", xlab = "Time", ylab = "Rates",  
              xaxt='n',col= col[i]  )
        axis(1, at= c((1+2017):(ending+1+2017)), labels=c((1+2017):(ending+1+2017)), las=2)
        #legend("topright", legend= name_legend, cex= 0.6 , lty = 1, col= color, bty = "n", inset = c(-0.30, 0.4), xpd = TRUE)
        }
}

plot.new()
legend("center", legend= name_legend, cex= 1 , lty = 1, col= color, bty = "n", inset = c(0.7, 0.2), xpd = TRUE)

dev.off()  # Exporting the plot directly to a pdf but do not forget dev.off !


for (i in 1:scenario_count) { 
  write.csv(x = scenario_results[[i]], file = paste0(main_scenario,".", name_legend[i],".csv"))
}


# Change in the results depends on the followings:
# Scenario_count 
# scenario_change_in_VCT and scenario_change_in_others
# main_scenario


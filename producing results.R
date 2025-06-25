# This file is for posterior predictive check 
# Line 40 is decision to which dataset to be used N=500 or N=5000. 
# IMPORTANT!! First find D_all then draw the plots. Because D_all is random everytime you run.
# Line 98 shows the compartment to be plotted for predictive check
library(deSolve)
library(MASS)
library(mvtnorm)

source("run_model_HIV2.R")

sample_number <- 1000+1  # +1 is base case. We will show the base case
number_of_years = 6
diagnosed_results <- matrix(ncol = number_of_years, nrow = sample_number) # includes base case that's why incresed by 1
D_all <- c()
random_beta1 <- 0.5247  # The first element of D_all will be the base case. That's why the order of codes is changed.
random_beta2 <- 0.2339
begin = 0
middle = 4
end = 16

# For N = 5000 Run these N = 500 and N=5000 separately bc they use some common variables that differ such as N.
data_for_check_N500 <- read.csv("read_posterior_N_500_first10.csv", header = FALSE)

posterior_density_beta1 <- density(data_for_check_N500[,1]) 
peak_beta1 <- posterior_density_beta1$x[which.max(posterior_density_beta1$y)]

posterior_density_beta2 <- density(data_for_check_N500[,2]) 
peak_beta2 <- posterior_density_beta2$x[which.max(posterior_density_beta2$y)]

mean(data_for_check_N500[,1])
mean(data_for_check_N500[,2])
peak_beta1
peak_beta2

plot(posterior_density_beta1, ylim=c(min(posterior_density_beta1$y),max(posterior_density_beta1$y)+0.3),las=1, main="HIV model" , xlab=expression(beta[1]), ylab="Density") 
plot(posterior_density_beta2, ylim= c(min(posterior_density_beta2$y),max(posterior_density_beta2$y)+5),las=1, main="HIV model" , xlab=expression(beta[2]), ylab="Density")
 # For N = 5000 Plots
# data_for_check_N5000 <- read.csv("read_posterior_N_5000.csv", header = FALSE)
# posterior_density_beta1 <- density(data_for_check_N5000[,1]) 
# peak_beta1 <- posterior_density_beta1$x[which.max(posterior_density_beta1$y)]
# 
# posterior_density_beta2 <- density(data_for_check_N5000[,2]) 
# peak_beta2 <- posterior_density_beta2$x[which.max(posterior_density_beta2$y)]
# plot(posterior_density_beta1, ylim=c(min(posterior_density_beta1$y),max(posterior_density_beta1$y)+0.3),las=1, main="HIV model" , xlab=expression(beta[1]), ylab="Density") 
# plot(posterior_density_beta2, ylim= c(min(posterior_density_beta2$y),max(posterior_density_beta2$y)+5),las=1, main="HIV model" , xlab=expression(beta[2]), ylab="Density")
# 
# mean(data_for_check_N5000[,1])
# mean(data_for_check_N5000[,2])
# peak_beta1
# peak_beta2

# These changes apply the second part (D2)
VCT_fold <- 4      #1 if no change 
others_fold <- 1   #1 if no change
  
for (i in 1:sample_number) { 
  initial_values <-  c(SH=0.02643, SL=0.97295, SVCT= 0, SO=0,  U=0.0003225, IVCT=1.67*10^-6, IO= 7.318*10^-5, ART=0.00022)
  params <- c( beta1 = random_beta1, beta2 = random_beta2, sigma1 = 0.00206, sigma2 = 0.13241, delta1 = 0.00434, delta2 = 0.19086, epsilon1 = 0.53
               , epsilon2 = 0.96,  omega = 2, psi = 0.219,   alpha = 0.747)
  times <- seq(from=begin, to=middle , by=1)
  
  total_numberof_people = 62819553
  
 
  D1 <- ode(
    y = initial_values,
    times = times,
    func = closed.sir.model,
    parms = params
  )      
  
  #since we changed the test rates between 2023-2034
  initial_values2 <-  D1[nrow(D1), 2:(ncol(D1)-1)]
  times2 <-  seq(from=middle, to=end , by=1)
  params2 <- c( beta1 = random_beta1, beta2 = random_beta2, sigma1 = 0.00206*VCT_fold, sigma2 = 0.13241*others_fold , delta1 = 0.00434*VCT_fold, delta2 = 0.19086*others_fold, epsilon1 = 0.53
               , epsilon2 = 0.96,  omega = 2, psi = 0.219,   alpha = 0.747)
  
  D2 <- ode(
    y = initial_values2,
    times = times2,
    func = closed.sir.model,
    parms = params2
  )     
    
  D = rbind(D1,D2[-1,]) # to prevent the initial2 to be repeated twice, we used -1.
  
  # Then prevalence and diagnosed cases are added to D matrix
  prevalence <- matrix(0, nrow=nrow(D))
  for (r in 1:nrow(D)){
    prevalence[r] = sum(D[r,6:9])
  }
  D <- cbind(D, prevalence)
  colnames(D)[11] <- "prevalence"
  

  before_change_nod <- (params["delta1"] + params["delta2"]) * D[,6] * total_numberof_people
  after_change_nod <-  (params["delta1"]*VCT_fold + params["delta2"]*others_fold) * D[,6]* total_numberof_people
  numberof_diagnosed <- before_change_nod[1:(middle+1)]
  numberof_diagnosed[(middle+2):(end+1)] <- after_change_nod[(middle+2):(end+1)]

  
  D <- cbind(D, numberof_diagnosed)
  colnames(D)[12] <- "diagnosed_cases"
  
  D_all <- rbind(D_all, D)  # for predictive check plots
  
  # IMPORTANT NOTE !!!! Which data is used N=500 or N = 5000. data_for_check changes according to this decision!!
  # Decision to be made N=500 or N=5000 (For now N=500)
  set.seed(123+i) 
  random_sample <- data_for_check_N500[sample(nrow(data_for_check_N500), 1),] 
  #random_sample <- data_for_check[sample(nrow(data_for_check_N5000), 1),]
  random_beta1  <-as.numeric(random_sample[1])
  random_beta2  <- as.numeric(random_sample[2])
}

for (j in 1:6) {   # To get each year's prediction from sampled results
  for (i in 1:sample_number) {
    diagnosed_results[i,j] = D_all[(j+(i-1)*17),12]
  }
}

write.csv(x = diagnosed_results, file = "predictive_check.csv" )


#Plotting the prediction intervals from 2018 till end
column <- 12 #10   #which compartment you want to see, or which output
fold  <- 100000 #10^7   # if you want to calculate the numbers or 
result_name_list <- list("time", "SH", "SL", "SVCT", "SO", "U", "IVCT", "IO", "ART", "incidence", "prevalence", "diagnosed_cases")
result_name <- paste0(result_name_list[column])
pdf(file = paste0("Preditive_check_plot", ".pdf" ),  width = 4, height = 4) # Exporting the plot directly to a pdf
par(mfrow = c(1,1))
if (column == 12) 
{ fold <- 1 }

mean_D_all <- c()
for (j in 1:17) {
  sum = 0
  for (i in 1:sample_number) {
    sum = sum + D_all[(j+(i-1)*17),column]
  }
  mean_D_all <- rbind(mean_D_all, sum/sample_number)
}
years <- c(2018:(end+2018))
par(mar = c(3, 4, 2, 3), xpd = TRUE) 
y_limit <- c(min(D_all[,column])*fold, (max(D_all[,column])+sd(D_all[,column]))*fold)

#Plot the Base case first:
plot( years , D_all[1:(end+1), column]*fold, xlim = c(2018,(end+2018)), 
      ylim = y_limit,
      cex.axis = 0.9,
      type = "l", main = result_name , xlab = "", ylab = "", 
      las=2, col = "Red", xaxt = "n") 
title(xlab="Time", line = 3, cex.lab=1)
title(ylab="Cases", line = 3, cex.lab=1) 
axis(1, at = years, labels = years, las = 2, cex.axis = 0.9)

#Plot the other random parameter sets' results
for (i in 2:(sample_number)) {
  lines(years, 
        D_all[(1+((i-1)*(end+1))):(i*(end+1)), column]*fold, 
        xlim = c(1,(end+1)), 
        ylim = y_limit,
        type = "l",
        main = "SH", xlab = "Time",  xaxt='n',col= "Gray"  )
}
lines( years , mean_D_all*fold, xlim = c(2018,(end+2018)), 
      ylim = y_limit, 
      type = "l", main = result_name , xlab = "Time",
      las=2, col = "Blue", xaxt = "n")
# Again the base case is repeated because after adding lines they disappeared
lines( years , D[1:(end+1), column]*fold, xlim = c(2018,(end+2018)), 
      ylim = y_limit, 
      type = "l", main = result_name , xlab = "Time",
      las=2, col = "Red", xaxt = "n")


# obtain the data for 95% Conficence Intervals
all_years = 17
per_pop = 1000
incidence_all_randoms <- matrix(ncol = all_years, nrow = sample_number)
prevalence_all_randoms <- matrix(ncol = all_years, nrow = sample_number)
diagnosed_all_randoms <- matrix(ncol = all_years, nrow = sample_number) # includes base case that's why incresed by 1

for (j in 1:all_years) {   # To get each year's prediction from sampled results
  for (i in 1:sample_number) {
    incidence_all_randoms[i,j] = D_all[(j+(i-1)*17),10]*per_pop
  }
}
colnames(incidence_all_randoms) <- 2018:(2018 + all_years - 1)

for (j in 1:all_years) {   # To get each year's prediction from sampled results
  for (i in 1:sample_number) {
    prevalence_all_randoms[i,j] = D_all[(j+(i-1)*17),11]*per_pop
  }
}
colnames(prevalence_all_randoms) <- 2018:(2018 + all_years - 1)
for (j in 1:all_years) {   # To get each year's prediction from sampled results
  for (i in 1:sample_number) {
    diagnosed_all_randoms[i,j] = D_all[(j+(i-1)*17),12]
  }
}
colnames(diagnosed_all_randoms) <- 2018:(2018 + all_years - 1)

# For Credible Intervals (CI) 
# Applied for each year seperately 
CI_lower_incidence <- apply(incidence_all_randoms, 2, quantile, probs = 0.025)
CI_upper_incidence <- apply(incidence_all_randoms, 2, quantile, probs = 0.975)

CI_lower_prevalence <- apply(prevalence_all_randoms, 2, quantile, probs = 0.025)
CI_upper_prevalence <- apply(prevalence_all_randoms, 2, quantile, probs = 0.975)

CI_lower_diagnosed <- apply(diagnosed_all_randoms, 2, quantile, probs = 0.025)
CI_upper_diagnosed <- apply(diagnosed_all_randoms, 2, quantile, probs = 0.975)


plot.new()
legend("center", 
       #inset = c(-0.7, -0.03), 
       legend = c("Base","Mean"), 
       x.intersp = 0.4,  #distance between legend text and legend line
       y.intersp = 0.8,  #distance between legend elements
       lty = c(1,1), col = c("Red", "Blue"), 
       cex = 1, 
       bty = "n") # nobox

dev.off()  # Exporting the plot directly to a pdf

#Plotting validation part 
#Example data: counts of some discrete events
data <- data.frame(
  category = c("2019", "2020", "2021", "2022", "2023"),
  count = c(mean(diagnosed_results[,2]), mean(diagnosed_results[,3]), mean(diagnosed_results[,4]), mean(diagnosed_results[,5]), mean(diagnosed_results[,6])  ) )
# Function to calculate Poisson confidence intervals
poisson_ci <- function(count, conf_level = 0.95) {
  alpha <- 1 - conf_level
  lower <- qpois(alpha / 2, count)
  upper <- qpois(1 - alpha / 2, count)
  return(c(lower, upper))
}
# Calculate confidence intervals for each count
c1_l <- mean(diagnosed_results[,2]) - sd(diagnosed_results[,2]) * (1.96 / sample_number)
c2_l <- mean(diagnosed_results[,3]) - sd(diagnosed_results[,3]) * (1.96 / sample_number)
c3_l <- mean(diagnosed_results[,4]) - sd(diagnosed_results[,4]) * (1.96 / sample_number)
c4_l <- mean(diagnosed_results[,5]) - sd(diagnosed_results[,5]) * (1.96 / sample_number)

c1_u <- mean(diagnosed_results[,2]) + sd(diagnosed_results[,2]) * (1.96 / sample_number)
c2_u <- mean(diagnosed_results[,3]) + sd(diagnosed_results[,3]) * (1.96 / sample_number)
c3_u <- mean(diagnosed_results[,4]) + sd(diagnosed_results[,4]) * (1.96 / sample_number)
c4_u <- mean(diagnosed_results[,5]) + sd(diagnosed_results[,5]) * (1.96 / sample_number)


data$ci_lower <- c(c1_l, c2_l, c3_l, c4_l)
data$ci_upper <- c(c1_u, c2_u, c3_u, c4_u)

# Print the data with confidence intervals
print(data)

# Plot points with error bars
par(mar = c(3,5,2,3) )
#nf <- layout(matrix(c(1, 2), ncol=2))
#layout.show(nf)

plot(1:length(data$category), data$count, ylim = c(3000, max(data$count) + 500),
     main = "Diagnosed cases vs. Predictions", xlab = "Category", ylab = "Cases",
     pch = 16, col = "blue", xaxt = "n", xlim = c(0.5, length(data$category) + 0.5),
     cex = 1.5)
# Add custom x-axis labels
axis(1, at = 1:length(data$category), labels = data$category)
# Add error bars
#arrows(x0 = 1:length(data$category), y0 = data$ci_lower, x1 = 1:length(data$category),
#       y1 = data$ci_upper, angle = 90, code = 3, length = 0.05, col = "darkblue")
target <- c(3955, 4298, 3203, 4285, 5710, 6029)
points(1,target[2], col="red", pch = 17, cex = 1.5)
points(2,target[3], col="red", pch = 17, cex = 1.5)
points(3,target[4], col="red", pch = 17, cex = 1.5)
points(4,target[5], col="red", pch = 17, cex = 1.5)
points(5,target[6], col="red", pch = 17, cex = 1.5)

legend("bottomright" , inset = c(-0.28, 0), legend = c("Cases", "Predicted"),
       pch = c(17, 16), col = c("red", "blue"),
       bty = "n", pt.cex = 1, cex = 0.9)

# old legend with "=" sign: legend("bottomright" , inset = c(-0.3, 0), legend = c("Cases", "Predicted","95% CI"),
#       pch = c(17, 16, (charToRaw("=="))), col = c("red", "blue", "blue"),
#       bty = "n", pt.cex = 1, cex = 0.7)


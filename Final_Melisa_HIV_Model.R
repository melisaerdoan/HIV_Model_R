#Version: Copy of General version of version 17
# initial S is converted into uniform dist. # rpo,s deleted and rounds and dpois deleted.

#when you use truncated productivities: dtmvnorm & rtmvnorm
#install.packages("tmvtnorm",dependencies=TRUE, repos='http://cran.rstudio.com/')
#library(tmvtnorm)

#final_my_func17 <- function() {
start_time <- Sys.time()

# This is ABC_SMC algorithm
library(deSolve)
library(MASS)
library(mvtnorm)
library(matlib)
#install.packages('mvtnorm')
modeltype <- 'HIV_model'  #'model2' , 'model3' #other way modeltype <- readline(prompt='Enter model_type:') 
epsilon_type <- 'epsilon3'  #epsilon1,2 or 3
calib_initial <- c('S')  # calib_initial HIV_modelde yok olmayacak!

# Based on the modeltype we import the closed_sir_model function to use the differential equations
if (modeltype == 'model1') {source("run_model1.R")}  
if (modeltype == 'model2') {source("run_model2.R")}  
if (modeltype == 'model3') {source("run_model3.R")}
#if (modeltype == 'HIV_model') {source("run_model_HIV.R")}  ## progression-6'da sundu�um hali v1
if (modeltype == 'HIV_model') {source("run_model_HIV.R")}   ## run_model_HIV2 ise (1-psi)*VCT'lisi v2

#Calculating the actual results
if (modeltype=='model1') {
calib_params <- c('gen','no', 'beta', 'gamma', 'S', 'distance', 'weight')
params <- c(beta=1.5,gamma=0.5)   # define the model parameters  
times <- seq(from=0, to=17, by=1)
initial_values <- c(S=99, I=1, R=0)
true_initial <- c(S=99, I=1, R=0)
n_par <- length(calib_params)-4      # Number of parameters to be estimated
} else if (modeltype=='model2') {
calib_params <- c('gen','no', 'beta', 'sigma', 'gamma', 'S', 'distance', 'weight')
params <- c(beta=1.5,sigma=0.75, gamma=0.5)   # define the model parameters  
times <- seq(from=0, to=17, by=1)
initial_values <- c(S=99, E=2,  I=1, R=0)
true_initial <- c(S=99, E=2,  I=1, R=0)
n_par <- length(calib_params)-4      
} else if (modeltype=='model3') {
calib_params <- c('gen','no', 'delta', 'beta', 'sigma', 'gamma', 'S', 'distance', 'weight')
params <- c(delta= 2, beta=1.5,sigma=0.75, gamma=0.5, mu=0.003, ni= 0.0015)   # define the model parameters  
times <- seq(from=0, to=17, by=1)
initial_values <- c(M= 20, S=99, E=2,  I=1, R=0) 
true_initial <- c(M= 20, S=99, E=2,  I=1, R=0) 
n_par <- length(calib_params)-4      
} else {
  calib_params <- c('gen','no','beta1','beta2', 'distance', 'weight')
  params <- c( beta1 = 0, beta2 = 0, sigma1 = 0.00206, sigma2 = 0.13241, delta1 = 0.00434, delta2 = 0.19086, epsilon1 = 0.53
                  , epsilon2 = 0.96,  omega = 2, psi = 0.219,   alpha = 0.747)  # define the model parameters  # put 0 for beta1 and beta2 to understand we do not use them from here
  times <- seq(from=0, to=17, by=1)
  initial_values <-c(SH=0.02643, SL=0.97295, SVCT= 0, SO=0,  U=0.0003225, IVCT=1.67*10^-6, IO= 7.318*10^-5, ART=0.00022)
  true_initial <-  c(SH=0.02643, SL=0.97295, SVCT= 0, SO=0,  U=0.0003225, IVCT=1.67*10^-6, IO= 7.318*10^-5, ART=0.00022)

  n_par <- length(calib_params)-4 # number of parameters that will be calibrated.      
}

print(paste('# Calibrated:', length(calib_params)-4, '  # All parameters (except initial):', length(params)))
print(modeltype)
print(epsilon_type)
  
# D  <- ode(
#   y = initial_values,
#   times = times,
#   func = closed.sir.model,
#   parms = params
# )

#ABC_SMC set-up
#G <- 6
N <- 500      # Number of accepted particles
#epsilon <- seq(400,400-(G-1)*10,-10)   # Tolerance value 
#epsilon <- c(400,300,200,100,80,30,20,10,9,8)   # Tolerance value 
epsilon <- c(20000)  #Can be adapted for  

accepted <- matrix(ncol = n_par+4, nrow = N)    # This is an empty matrix to store results
accepted2 <- matrix(ncol = n_par+4, nrow = N)    # This is an empty matrix to store results
accepteds <- c()    #store all accepted matrices (generations) in accepteds 
rejecteds<- c()       #store all rejecteds
acc_rate_pop <- c()   # acc_rate per population
colnames(accepted) <- calib_params
g <- 1   # Initializing the counter of generation
i <- 1   # Initializing the counter of accepted particles
j <- 1   # Initializing counter of proposed particles 
i2<- 1   # the same as i and j but this time used for after first generation        
j2 <-1
priorIsZero <- 0

totalNumberOfPeople = 62819553   # old one : 62858367
target <- c(3955, 4298, 3203, 4285, 5710)   #2018,2019,2020,2021,2022
#function to calculate distance
calc_distance <- function(target, estimates){    #before: function(D, D_star)
  #dist = sum((D-D_star)^2)   #Euclidean distance
  dist = sum(abs(target-estimates))   #Absolute distance
  dist
}
  

if (g==1) {               # this block is for the 1st generation
  while(i <= N)  {        # While accepted parameters are less then or equal to N. That is, run until accepting N particles
    
    theta_star <- c()                   # defined as empty and added based on the ORDER of the names in model schematic
    if (modeltype=='model3') {
      theta_star['delta'] <- runif(1,1,4)
    }
    #theta_star['beta'] <- runif(1,0,3)  # Sample from prior distribution
    if (modeltype=='model2' | modeltype=='model3') {
      theta_star['sigma'] <- runif(1,0,1.25)
    }
    #theta_star['gamma'] <- runif(1,0,1)  
   # initial_values['S'] <- runif(1,90,110)   #random numbers from poisson distribution 
    
    if (modeltype=='HIV_model') {   #Calibrated parameters are beta1 and beta2 for rate of movement from S_high and S_low to HIV+ Unaware compartment ,respectively
      theta_star['beta1'] <- runif(1,0,1)
      theta_star['beta2'] <- runif(1,0,1)
    }

    # Simulate data set w new random parameter sets ( using new random beta1 and)
    
 if (modeltype=='HIV_model') {
       sir_values <- ode(
         y = initial_values,
         times = times,
         func = closed.sir.model,
         parms = c(theta_star, params['sigma1'], params['sigma2'], params['delta1'], params['delta2'], params['epsilon1'],
         params['epsilon2'],  params['omega'], params['psi'],  params['alpha']) )
D_star <- as.data.frame(sir_values)  } # record as data-frame
    
    
    estimated_cases1 <- (params["delta1"] + params["delta2"]) * D_star[,6] * totalNumberOfPeople
    estimated_cases1d <- estimated_cases1[1:5]
    estimated_cases2 <- c()
    estimated_cases3 <- c()
    
    for (e in 2:nrow(D_star)) {
      estimated_cases2[e-1] <- D_star[e,7] - D_star[e-1,7] + D_star[e,8] - D_star[e-1,8]
      estimated_cases2[e-1] <- estimated_cases2[e-1] * totalNumberOfPeople
    } 
    
    for (e in 2:nrow(D_star)) {
      estimated_cases3[e-1] <- D_star[e,7] - D_star[e-1,7] + D_star[e,8] - D_star[e-1,8] + D_star[e,9] - D_star[e-1,9]
      estimated_cases3[e-1] <- estimated_cases3[e-1] * totalNumberOfPeople
    }
    #write.csv(cbind(estimated_cases1,c(estimated_cases2,"end"),c(estimated_cases3,"end")), file="estimates.csv")
    #estimated_cases2d <-  estimated_cases2[1:5]
    
    
    distance <- calc_distance(target, estimated_cases1d) #OR
    #distance <- calc_distance(target, estimated_cases2d)
  
    if (distance <= epsilon[g]) {           #check if the distance is less than the tolerance or not
        accepted[i,] = c(g,i,theta_star, distance, NA)  #store accepted parameters
        #accepted <- accepted[order(accepted[,6], decreasing=FALSE),]     #sorting the accepted particles based on distance
        accepted[,n_par+4] <- 1/N                           #assigning weights to the particles
        ##note that [i,n_par+4] runs wrongly!!
        #accepted[,n_par+4]<- cumsum(accepted[,n_par+4])     #assigning cumulative distribution based on weights
        i <- i+1 
        # Update the counter
    } else {
      rejecteds <- rbind(rejecteds,c(g,i,theta_star,distance))
    } 
    j <- j +1    # update counter
    #print(j)
  } 
accepteds <- rbind(accepteds,accepted)

threshold_mean <- mean(accepted[,n_par+3]) #creating new threshold
threshold_sd <- sd(accepted[,n_par+3])

acc_rate1 <- (i-1)/(j-1)  #acceptence rate of first population
acc_rate_pop <- rbind(acc_rate_pop,data.frame(gen=g,acc_rate=acc_rate1))
}

g <- g+1       #going to next generation
#print(acc_rate2)
if (epsilon_type == 'epsilon1') {
  #print('epsilon1')
  epsilon[g] <- threshold_mean - (threshold_sd)
  } else if (epsilon_type == 'epsilon2'){
  #print('epsilon2')
  epsilon[g] <- threshold_mean - ((threshold_sd)/2)
  } else {
  #print('epsilon3')
  epsilon[g] <- median(accepted[,n_par+3])  
} 
#epsilon[g] <- (threshold_mean*(1-acc_rate2))

#From now on, variable names are defined by adding them 2 to separate them from previous ones.
 
lowcov <- 0
highcov <- 0

while (tail(epsilon,1)>=2450){
  i2 <- 1
  print(g)
  if (g==2) {
    prob1=0.1
    prob2=0.9
  }else{
    prob1 <- lowcov/(lowcov+highcov)
    prob2 <- highcov/(lowcov+highcov)
  } 
  while(i2 <= N){
    #sampling from the previous accepted based on weights
    theta_star2 <- c() 
    my_sample <- sample(x=accepted[,2],size=1,replace = TRUE, prob= accepted[,n_par+4])
    #print(accepted[my_sample,])   #previous particle 

    #### mean and covariance determined by Covariance intersection 
    
    intermediate <- accepteds[(((g-2)*N)+1):((g-1)*N),]   # just previous pop
    intermediate_a <- intermediate[intermediate[,n_par+3]<median(intermediate[,n_par+3]),3:(n_par+3)] #just points with calibrated parameters
    intermediate_b <- intermediate[intermediate[,n_par+3]<quantile(intermediate[,n_par+3])[2],3:(n_par+3)]
    randomsample_c <- intermediate[,3:(n_par+3)]

    columns <- colnames(intermediate_a)
    xA <- c()
    xB <- c()
    for (co in 1:n_par) {
      xA <- rbind(xA,mean(intermediate_a[,co]))
      xB <- rbind(xB,mean(intermediate_b[,co]))
    }   
    rownames(xA) <- columns[-n_par-1]
    rownames(xB) <- columns[-n_par-1]
    xC <- accepted[my_sample,3:(n_par+2)]
       
    CA <- cov(intermediate_a[,1:n_par])
    CB <- cov(intermediate_b[,1:n_par])
    CC <- cov(randomsample_c[,1:n_par])

    omega1 <- 0.333
    omega2 <- 0.333
    omega3 <- 0.334
    
    xA <- t(t(xA))
    xB <- t(t(xB))
    xC <- t(t(xC))
    C  <-  Ginv((omega1*Ginv(CA))+(omega2*Ginv(CB))+(omega3*Ginv(CC))) 
    xCov <- C%*%((omega1*Ginv(CA)%*%xA) + (omega2*Ginv(CB)%*%xB) + (omega3*Ginv(CC)%*%xC))
     
    
    bias <- cov(rejecteds[,3:(n_par+2)])
    rcoef <- sample(x=c(0,1),size= 1,replace = TRUE, prob=c(prob1,prob2))
    
    sigma <- C +  rcoef * bias  
    mean <- xCov # Covariance intersection mean

    theta_double_star <- mvrnorm(1, mu = mean, Sigma = sigma) #  this is our new particle
    #theta_double_star <- rtmvnorm(1, mean=mean, sigma =sigma, lower=c(0,0,85), upper=c(3,1,105))  #if you use truncated mvnormal Kernel
    #check: hist(mvrnorm(1000, mu = mean, Sigma = sigma)[,1])
    theta_double_star <- t(as.matrix(theta_double_star))
    colnames(theta_double_star) <- columns[-n_par-1]  #names of calibrated parameters
    
    #ORDER  MATTERS while constructing theta_star_2
    if (modeltype=='model3') {
      theta_star2['delta'] <- theta_double_star[,'delta'] }
    
    #theta_star2['beta'] <-  theta_double_star[,'beta']  #BIG CHANGE::: ['beta'] is added instead of [1] for theta_double_star part
    if (modeltype=='model2' | modeltype=='model3') {
      theta_star2['sigma'] <- theta_double_star[,'sigma'] }
    
    if (modeltype=='HIV_model') {
      theta_star2['beta1'] <- theta_double_star[,'beta1'] 
      theta_star2['beta2'] <- theta_double_star[,'beta2'] }

    #theta_star2['gamma'] <- theta_double_star[,'gamma']
    #initial_values['S'] <-  theta_double_star[,'S']
    #prior calculation & control whether prior is non-zero or not  # it should be revised.
    
    if (modeltype=='model1'){
      prior <- dunif(theta_star2['beta'],0,3)*dunif(theta_star2['gamma'],0,1)*dunif(initial_values['S'],90,110)
    }
    if (modeltype=='model2'){
      prior <- dunif(theta_star2['beta'],0,3)*dunif(theta_star2['sigma'],0,1.25)*dunif(theta_star2['gamma'],0,1)*dunif(initial_values['S'],90,110)
    }
    if (modeltype=='model3'){
      prior <- dunif(theta_star2['delta'],1,4)*dunif(theta_star2['beta'],0,3)*dunif(theta_star2['sigma'],0,1.25)*dunif(theta_star2['gamma'],0,1)*dunif(initial_values['S'],90,110)
    }
    if (modeltype=='HIV_model'){
      prior <- dunif(theta_star2['beta1'],0,1)*dunif(theta_star2['beta2'],0,1)
    }
    
    if (prior != 0) {         
      
      # if (modeltype == 'model3') {
      #   sir_values2 <- ode(
      #     y = initial_values,
      #     times = times,
      #     func = closed.sir.model,
      #     parms = c(theta_star2,params['mu'], params['ni']) )
      # } else if  (modeltype=='HIV_model') {
      #   sir_values2 <- ode(
      #     y = initial_values,
      #     times = times,
      #     func = closed.sir.model,
      #     parms = c(theta_star2, params['sigma1'], params['sigma2'], params['delta1'], params['delta2'], params['epsilon1'],
      #               params['epsilon2'],  params['omega'], params['psi'],  params['alpha']) )
      # } else {
      #   sir_values2 <- ode(
      #    y = initial_values,
      #    times = times,
      #    func = closed.sir.model,
      #    parms = theta_star2)
      # }
      
      if (modeltype=='HIV_model') {
        sir_values2 <- ode(
          y = initial_values,
          times = times,
          func = closed.sir.model,
          parms = c(theta_star2, params['sigma1'], params['sigma2'], params['delta1'], params['delta2'], params['epsilon1'],
                    params['epsilon2'],  params['omega'], params['psi'],  params['alpha']) )
          } # record as data-frame
      D_star <- as.data.frame(sir_values2)
      
      estimated_cases1 <- (params["delta1"] + params["delta2"]) * D_star[,6] * totalNumberOfPeople
      estimated_cases1d <- estimated_cases1[1:5]
      estimated_cases2 <- c()
      for (e in 2:nrow(D_star)) {
        estimated_cases2[e-1] <- D_star[e,7] - D_star[e-1,7] + D_star[e,8] - D_star[e-1,8] 
        estimated_cases2[e-1] <- estimated_cases2[e-1] * totalNumberOfPeople
      } 
      #write.csv(cbind(estimated_cases1,c(estimated_cases2,"end")), file="estimates.csv")
      estimated_cases2d <-  estimated_cases2[1:5]
      
      
      distance2 <- calc_distance(target, estimated_cases1d) #OR
      #distance2 <- calc_distance(target, estimated_cases2d)
      
      
      if (distance2 <= epsilon[g]) {           #checks the distance is less than the tolerance or not
        accepted2[i2,] = c(g,i2,theta_star2, distance2, NA)  #store accepted parameters
        lowcov <- lowcov + 1 * (1-rcoef)  # If rcoef=1 then it is highcov so lowcov stays the same
        highcov <- highcov + 1 * rcoef # If rcoef=0 then it is lowcov so highcov stays the same
        
        #assigning weights to the particles       
        weight_num <- prior                
        
        weight_den <- 0
        for (a in 1:N) {     #calculation of denominator of weight calculation
          weight_den <-weight_den + accepted[a,n_par+4]*dmvnorm(accepted2[i2,3:(n_par+2)], mean=accepted[a,3:(n_par+2)], sigma=sigma) #previous pop. weight * (current point, mean=previous particle, sigma=from previous particle
          #weight_den <-weight_den + accepted[a,7]*dtmvnorm(accepted2[i2,3:5], mean=accepted[a,3:5], sigma=sigma,lower=c(0,0,85),upper=c(3,1,105)) #if you use truncated mvnorm
          }
        accepted2[i2,n_par+4] <- weight_num/weight_den   #weight calculation PART I
        i2 <- i2+1  
        # Update the counter
      } else {
        rejecteds <- rbind(rejecteds,c(g,i2,theta_star2,distance2))
      } 
        
      } else { 
        priorIsZero=priorIsZero+1
      }
    
    j2 <- j2 +1    # update counter #when prior=0 or not 0 it increases by 1
    #print(j2)
    #acc_rate <- (N*G)/(j-1+j2-1)   #Calculate acceptance rate  (# of accepted particle/ all runs)
    acc_rate <- (N*g)/(j-1+j2-1) #with prior=0
    #acc_rate <- (N*g)/(length(accepteds[,1])+ length(rejecteds[,1])) #without prior=0
    #print(acc_rate)
    
  }
  #this part comes after accepting N particles in one generation
  #accepted2 <- accepted2[order(accepted2[,6], decreasing=FALSE),]     #sorting the accepted particles
  #accepted2[1:N,2] <- as.matrix(seq(1,N,1))  #giving number to each sorted particle
  numberOfParticle <- N+sum(rejecteds[,1]==g)  #because here there is no current g accepteds
  acc_rate_pop <- rbind(acc_rate_pop,data.frame(gen=g,acc_rate=N/numberOfParticle))
  #weight calculation PART II (normalization)
  sum_weights <- sum(accepted2[,n_par+4])                 # this part will be done after all weights are calculated
  accepted2[,n_par+4] <- accepted2[,n_par+4]/sum_weights  #  new normalized weights   
  #accepted2[,n_par+4] <- cumsum(accepted2[,n_par+4])                           
  
  accepteds <- rbind(accepteds,accepted2)
  accepted <- accepted2    #update accepted to sample from last accepted particles
  accepted2 <- matrix(nrow = N, ncol = n_par+4)   #reset accepted2
  g <- g+1 
  
  threshold_mean <- mean(accepted[,n_par+3]) #creating new threshold
  threshold_sd <- sd(accepted[,n_par+3])
  
  #epsilon[g] <- threshold_mean - (threshold_sd)
  if (epsilon_type == 'epsilon1') {
    #print('epsilon1')
    epsilon[g] <- threshold_mean - (threshold_sd)
  } else if (epsilon_type == 'epsilon2'){
    #print('epsilon2')
    epsilon[g] <- threshold_mean - ((threshold_sd)/2)
  } else {
    #print('epsilon3')
    epsilon[g] <- median(accepted[,n_par+3])  
  }
  
}

end_time <- Sys.time()
time_spent <- end_time-start_time

posterior <- tail(accepteds[,3:(n_par+2)],N)

#Plotting the posterior density
plot(density(posterior[1:N,1]), ylim=c(0,2.5),las=1, main="HIV model" , xlab=expression(beta1), ylab="Density") 
plot(density(posterior[1:N,2]), ylim=c(0,100),las=1, main="HIV model" , xlab=expression(beta2), ylab="Density")
#plot(hist(posterior[1:50,2]))

#extras:
#lines(density(postfigs[1:7500,6]), col="red")
#lines(density(postfigs[1:7500,31]), col="blue", lty = 3)
#abline(v=99, col="gray", lty=3)

write.csv(x=posterior, file = "Dist.csv")
print(acc_rate)
print(time_spent)
print(acc_rate_pop)
print(priorIsZero)
print(j)
print(j2)


#Calculation of Mean Absolute Error:
# MAE_use <- tail(accepteds, N)
# 
# MAE_sum <- cbind(c(abs(MAE_use[,3]-params[1])),c(abs(MAE_use[,4]-params[2])), c(abs(MAE_use[,5]-99)))
# MAE <- (colSums(MAE_sum))/N
# MAE <- t(MAE)

          #Calculation of Mean Absolute Error:
          ####   
# true_values <- c()
# MAE_sum <- c()
# MAE_use <- tail(accepteds, N)
# n_par_initial <- length(calib_initial)
# for (t in 1:(n_par-n_par_initial)) {
#   true_values[t] <- c(params[t])
# }
# for (t in calib_initial) {
#   true_values <- append(true_values,true_initial[t])
# }
# 
# for(m in 1:n_par) {
#   MAE_sum <- cbind(MAE_sum,(abs(MAE_use[,(2+m)]-true_values[m])))
# }
# MAE <- (colSums(MAE_sum))/N
# MAE <- t(MAE)

# return(list('Acceptence'= format(acc_rate, digits=3),'Acc_Pop'= acc_rate_pop, 'Time'= format(time_spent, digits=3),'Accuracy'=MAE, 'Posterior'=posterior))
# } 


# numrow1 <- (G-1)*N+1
# numrow2 <- G*N
# hist(accepteds[numrow1:numrow2,3])
# hist(accepteds[numrow1:numrow2,4])
# hist(accepteds[numrow1:numrow2,5])
# 
# 
# # #Visualizing Tool
# pdf(file='Faster9.pdf')
# #or visualize <- read.csv('accepteds_visual.csv')
# ############## G <- length(epsilon)-1  # if it has an adaptive threshold
# visualize <- accepteds[,n_par+3]
# maximum <- max(visualize)   #or max(visualize[,7])
# plot.new()
# plot.window(ylim=c(0,maximum+120), xlim=c(0,G*N))
# axis(side=1, at= seq(0,G*N,10), labels=seq(0,G*N,10))
# axis(side=2, at=seq(0,maximum+120,50), labels= seq(0,maximum+120,50), las=2)
# my_data <- visualize
# points(x=my_data, col='red')
# clip(0,120,0,maximum+100)
# abline(v=seq(N,G*N,N), col='grey', lty=3)
# mtext('Accepted particle number', side=1, line=2.3)
# mtext('Distance', side=2, line=2.3)
# 
# for (pop in 1:G) {
#   #print(N/2+(pop-1)*N)
#   text(N/2+(pop-1)*N, maximum+100, paste0('Pop #',pop))
# }
# 
# dev.off()


# 
# #calculating distance:
# melisa  <- c()
# for (im in 21:40)
# {# print(im, im-20)
#   melisa <- rbind(melisa,sum(abs(accepteds[im,3:5]-accepteds[im-20,3:5])))
# }
# 
# part1 <- accepteds[21:40,]
# part1 <- cbind(part1,melisa)

#Adding plots for each population with different colours

#plot(accepteds[1:500,3],accepteds[1:500,4], col=c('red'))
#points( accepteds[501:1000,3],accepteds[501:1000,4], col=c('green'))
#points( accepteds[1001:1500,3],accepteds[1001:1500,4], col=c('blue'))
#points( accepteds[1501:2000,3],accepteds[1501:2000,4], col=c('gray'))
#points( accepteds[2001:2500,3],accepteds[2001:2500,4], col=c('pink'))

#sum(accepteds[which(accepteds[,1]==1),6]<147)

#high_weights <- accepted[accepted[,7]>median(accepted[,7]),]
#min_dist_pool <- accepteds[accepteds[, n_par+3]<(1-acc_rate_pop[g-1,2])*epsilon[g]]


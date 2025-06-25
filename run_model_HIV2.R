#This is the model
closed.sir.model <- function(t,x,params) {
  with(as.list(c(x,params)), {
    
    SH   <- x[1] 
    SL   <- x[2]
    SVCT <- x[3]
    SO   <- x[4]
    U    <- x[5]
    IVCT <- x[6]
    IO   <- x[7]
    ART  <- x[8]

    #parameters
    beta1 <- params['beta1']           
    beta2 <- params['beta2']
    sigma1 <- params['sigma1']
    sigma2 <- params['sigma2']      
    delta1 <- params['delta1']
    delta2 <- params['delta2'] 
    epsilon1 <- params['epsilon1']  
    epsilon2 <- params['epsilon2']  
    omega <- params['omega']        
    psi <- params['psi']
    alpha <- params['alpha']

    N  <- SH + SL + SVCT + SO + U + IVCT + IO + ART 
    N <- 1
    
    #model equations
    dSHdt  <-  (1/omega)*SVCT - sigma1*SH - beta1*SH*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART)       
    dSLdt  <-  (1/omega)*SO - sigma2*SL - beta2*SL*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART)
    dSVCTdt <-  sigma1*SH - (1/omega)*SVCT - (1-psi)*beta1*SVCT*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART)
    dSOdt <-    sigma2*SL - (1/omega)*SO - beta2*SO*( U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART)
    dUdt  <-  beta1*SH*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART) +  beta2*SL*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART) + (1-psi)*beta1*SVCT*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART) + beta2*SO*( U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART)-(delta1+delta2)*U 
 
       
    dIVCTdt  <-  delta1*U - alpha*IVCT
    dIOdt  <-   delta2*U - alpha*IO
    dARTdt  <-  alpha*(IVCT + IO)

    incidence <- beta1*SH*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART) + beta2*SL*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART) + (1-psi)*beta1*SVCT*(U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART) + beta2*SO*( U + (1-epsilon1)*(1-psi)*IVCT + (1-epsilon1)*IO + (1-epsilon2)*ART)  
    dxdt <- c(dSHdt, dSLdt, dSVCTdt, dSOdt, dUdt, dIVCTdt, dIOdt, dARTdt)   #put the results into a single vector
    list(dxdt,incidence)    #return should be a list, so result is a list ! ,incidence
  }) 
}
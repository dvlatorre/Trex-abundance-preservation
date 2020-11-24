# Logistic function for growth
logistic <- function(M_max,b = 0.55,c = 16.2,d = 5,t){
  (M_max/(1+exp(-b*(t-c))))+d
}

# Survivorship - from Erickson et al 2006 (fig 2)
surv <- function(t, a, g, ni = 1){
  ni*exp((a/g)*(1-exp(g*t)))
}
# weighted survivorship
w_surv <- function(t, a, g, ni = 1){surv(t, a, g)/sum(surv(t, a, g))}

bootstrap_survivorship <- function(n, obsAges = c(2,6,8,9,11,14,14,15,16,17,18,18,18,18,18,19,21,21,21,22,22,22,22,22,22,23,23,24,24,28)){
  aEstimates <- vector()
  gEstimates <- vector()
  for(i in 1:n){
    while(TRUE){
      sampledAges <- sample(obsAges, size = length(obsAges), replace = T)
      propSurvSamp <- vector()
      for(j in 1:length(sampledAges)){
        propSurvSamp[j] <- 1-sum(sampledAges<sampledAges[j])/length(sampledAges)
      }
      
      gompertzFit <- try(nls(propSurvSamp ~ exp( (a/g)*(1-exp(g*sampledAges))), start = list(a = 0.002, g = 0.2214)), silent = TRUE)
      if(class(gompertzFit) != "try-error"){
        aEstimates[i] <- summary(gompertzFit)$parameters[1,1]
        gEstimates[i] <- summary(gompertzFit)$parameters[2,1]
        break()
      }
    }
  }
  return(cbind(aEstimates, gEstimates))
}

#Simulations to calculate expected number of T. rex to ever live
trex_ever_lived2 <- function(n, AsymBMmean, AsymBMsdev, damuth_slope = -0.75, logA_mid, logA_sd, SMmean, SMsdev, TRmin, TRmax, GRmean, GRsdev,  agesResolution = 0.1, maxAge = 28){

  #--------------------------------------
  #--- Age of sexual maturity (years) ---
  #--------------------------------------
  
  sex_mat <- round(rnorm(n, SMmean, SMsdev), digits = 1)
  
  #--------------------------------------
  #--- Ages from sexMat to 28 (years) ---
  #--- generates a list of vectors ------
  #--- to be used for growth and lx -----
  #--------------------------------------
  
  max_age <- maxAge
  ages <- lapply(sex_mat, seq, to = max_age, by = agesResolution)
  
  #-------------------------------------
  #--- Asymptotic_mass body mass (kg) ---
  #-------------------------------------
  
  asymp_mass <- rnorm(n, AsymBMmean, AsymBMsdev)
  
  #--------------------------------------------
  #--- Bootstrap gompertz survivorship fits ---
  #--------------------------------------------
  
  AxG <- bootstrap_survivorship(n)
  
  #-----------------------
  #---  Mean Body Mass ---
  #--- and lx ------------
  #-----------------------
  
  lx <- list()
  mean_body_mass <- vector()
  
  for(i in 1:n){
    
    foo <- w_surv(t = ages[[i]], a = AxG[i,1], g = AxG[i,2])
    mean_body_mass[i] <- sum(logistic(asymp_mass[i], t = ages[[i]])*foo)*1000
    
    lx[[i]] <- surv(t = ages[[i]], a = AxG[i,1], g = AxG[i,2])
    
  }

#--- Generation time (yrs) ---
#-----------------------------
  
  b_i <- sapply(lx, function(x){1/sum(x)})
  
  # ages * lx * b_i - but there are lists, so using Map
  i_bi_li <- Map('*', Map('*', ages, lx), b_i)
  
  gen_time <- sapply(i_bi_li, sum)

  # gen_time <- vector()
  # i_bi_li <- vector()
  # b_i <- vector()
  
#  for (i in 1:n){
    
    
#    b_i <- 1/sum(lx[sex_mat[i]:length(lx)])
    
#    i_bi_li[1:(sex_mat[i] - 1)] <- 0
    
#    for (j in sex_mat[i]:length(lx)){
#      i_bi_li[j] <- j*lx[j]*b_i
#    }
    
#    gen_time[i] <- sum(i_bi_li)
#    i_bi_li <- vector()
#  }
  
  #------------------------------------------
  #--- Mean population density (ind/km^2) ---
  #------------------------------------------
  
  pop_density <- (10**(rnorm(n, logA_mid, logA_sd) + damuth_slope * log(mean_body_mass, base = 10)))
  
  #---------------------------
  #--- Temporal range (ma) ---
  #---------------------------
  
  temporal_range <- runif(n, TRmin, TRmax)
  
  #-----------------------------
  #--- Number of generations ---
  #-----------------------------
  
  num_generations <- (temporal_range * 1e+06)/gen_time
  
  #----------------------------------
  #--- Geographic range (M. km^2) ---
  #----------------------------------
  
  geographic_range <- rnorm(n, GRmean, GRsdev)
  
  #--------------------------------------
  #--- Standing number of individuals ---
  #--------------------------------------
  
  standing_number <- geographic_range * 1e+06 * pop_density
  
  #-----------------------------------
  #--- Total number of individuals ---
  #-----------------------------------
  
  total_number <- standing_number * num_generations
  
  return(list(total_number = total_number, standing_number = standing_number, geographic_range = geographic_range, num_generations = num_generations, temporal_range = temporal_range, gen_time = gen_time, pop_density = pop_density, mean_body_mass = mean_body_mass))
}

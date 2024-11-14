## Simulation study to see if WAIC can correctly detect Winner-take-all scenarios

library(NeuralComp)
library(MASS)
library(truncnorm)
library(ggplot2)
library(latex2exp)
library(truncnorm)
library(future.apply)
library(ggpubr)
library(scales)
library(transport)
library(ggallin)
library(statmod)

save_dir <- "/Users/ndm34/Documents/WAIC_Sim4"

#### Generate WTA Scenarios

### Function to generate data from competition model
generate_data_TI <- function(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, N_A, N_B, N_AB, seed, time, basis_degree,
                             boundary_knots, internal_knots){
  set.seed(seed)
  n_A <- rep(NA, N_A)
  n_B <- rep(NA, N_B)
  n_AB <- rep(NA, N_AB)
  L_AB <- list()
  X_A <- list()
  X_B <- list()
  X_AB <- list()
  
  
  ## Generate A data
  for(i in 1:N_A){
    total_time <- 0
    spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
    spike <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
    total_time <- spike
    while(total_time < time){
      spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
      spike_i <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
      total_time <- total_time + spike_i
      spike <- c(spike, spike_i)
    }
    total_time <- total_time - spike[length(spike)]
    spike <- spike[-length(spike)]
    X_A[[i]] <- spike
    n_A[i] <- length(spike)
  }
  
  ## Generate B data
  for(i in 1:N_B){
    total_time <- 0
    spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
    spike <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
    total_time <- spike
    while(total_time < time){
      spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
      spike_i <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
      total_time <- total_time + spike_i
      spike <- c(spike, spike_i)
    }
    total_time <- total_time - spike[length(spike)]
    spike <- spike[-length(spike)]
    X_B[[i]] <- spike
    n_B[i] <- length(spike)
  }
  
  ## Generate AB data
  
  for(i in 1:N_AB){
    total_time <- 0
    spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
    spike_A <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
    spike_B <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
    if(spike_A < spike_B){
      spike <- spike_A
      L_i <- 0
      total_time <- spike
    }else{
      spike <- spike_B
      L_i <- 1
      total_time <- spike
    }
    
    iter <- 1
    while(total_time < time){
      spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
      spike_A <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
      spike_B <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
      if(L_i[iter] == 0){
        if(spike_A < spike_B + delta){
          spike <- c(spike, spike_A)
          L_i <- c(L_i, 0)
        }else{
          spike <- c(spike, spike_B + delta)
          L_i <- c(L_i, 1)
        }
      }
      if(L_i[iter] == 1){
        if(spike_A + delta < spike_B){
          spike <- c(spike, spike_A + delta)
          L_i <- c(L_i, 0)
        }else{
          spike <- c(spike, spike_B)
          L_i <- c(L_i, 1)
        }
      }
      iter <- iter + 1
      total_time <- sum(spike)
    }
    total_time <- total_time - spike[length(spike)]
    L_i <- L_i[-length(spike)]
    spike <- spike[-length(spike)]
    n_AB[i] <- length(spike)
    L_AB[[i]] <- L_i
    X_AB[[i]] <- spike
  }
  
  spline <- GetBSpline(seq(0, 1.0, 0.01), basis_degree, boundary_knots, internal_knots)
  return(list("X_A" = X_A, "X_B" = X_B, "X_AB" = X_AB, "L_AB" = L_AB, "n_A" = n_A,
              "n_B" = n_B, "n_AB" = n_AB, "A_CI" = (I_A * exp(spline %*% basis_coef_A)), "B_CI" = (I_B * exp(spline %*% basis_coef_B))))
}



run_sim5 <- function(iter){
  set.seed(iter)
  
  I_A <- rtruncnorm(1, a = 0, mean = 40, sd = 6)
  I_B <- rtruncnorm(1, a = 0, mean = 80, sd = 6)
  sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
  sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  delta <- rlnorm(1, -2.5, 0.5)
  
  dat <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 25, 25, 25, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  dat1 <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 25, 25, 25, iter + 1, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  WTA_truth <- "A"
  if(runif(1) > 0.5){
    dat$X_AB <- dat1$X_A
    dat$n_AB <- dat1$n_A
  }else{
    dat$X_AB <- dat1$X_B
    dat$n_AB <- dat1$n_B
    WTA_truth <- "B"
  }
  
  res <- Sampler_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  ## Run IIGPP Model
  res_A <- Sampler_IIGPP(dat$X_A, dat$n_A, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_B <- Sampler_IIGPP(dat$X_B, dat$n_B, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_AB <- Sampler_IIGPP(dat$X_AB, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  dat_A_AB <- c(dat$X_A, dat$X_AB)
  n_A_AB <- c(dat$n_A, dat$n_AB)
  res_A_AB <- Sampler_IIGPP(dat_A_AB, n_A_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  dat_B_AB <- c(dat$X_B, dat$X_AB)
  n_B_AB <- c(dat$n_B, dat$n_AB)
  res_B_AB <- Sampler_IIGPP(dat_B_AB, n_B_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  WAIC_WTA_A <- WAIC_Winner_Take_All(dat_A_AB, dat$X_B, n_A_AB, dat$n_B, res_A_AB, res_B, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  WAIC_WTA_B <- WAIC_Winner_Take_All(dat$X_A, dat_B_AB, dat$n_A, n_B_AB, res_A, res_B_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  WAIC_Comp <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  waic_iigpp <- WAIC_IIGPP(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res_A, res_B, res_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  
  params <- list("I_A" = I_A, "I_B" = I_B, "sigma_A" = sigma_A, "sigma_B" = sigma_B,
                 "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B)
  output <- list("WAIC_Comp" = WAIC_Comp, "WAIC_IIGPP" = waic_iigpp, 
                 "WAIC_WTA_A" = WAIC_WTA_A, "WAIC_WTA_B" = WAIC_WTA_B, "WTA_truth" = WTA_truth,
                 "res" = res, "res_A" = res_A, "res_B" = res_B, "res_AB" = res_AB, "res_A_AB" = res_A_AB, 
                 "res_B_AB" = res_B_AB, "data" = dat, "params" = params)
  saveRDS(output, paste0(save_dir, "/WTA/output", iter,".RDS"))
}

library(future.apply)

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
dir.create(paste0(save_dir, "/WTA"))
already_ran <- dir(paste0(save_dir, "/WTA"))
to_run <- which(!paste0("output", 1:100, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim5(this_seed))



### Post Processing
files <- dir(paste0(save_dir, "/WTA"))
WAIC <- matrix(0, nrow = length(files), 4)
wasserstein_dist <- rep(0, length(files))
truth <- rep(0, length(files))
for(i in 1:length(files)){
  output <- readRDS(paste0(save_dir, "/WTA/", files[i]))
  WAIC[i,1] <- output$WAIC_Comp$WAIC
  WAIC[i,2] <- output$WAIC_IIGPP$WAIC
  WAIC[i,3] <- output$WAIC_WTA_A$WAIC
  WAIC[i,4] <- output$WAIC_WTA_B$WAIC
  wasserstein_dist[i] <- wasserstein1d(output$data$n_A, output$data$n_B)
  truth[i] <- output$WTA_truth
  print(i)
}

min_index <- rep(0, 100)
for(i in 1:length(files)){
  min_index[i] <- which.min(WAIC[i,])
}

table(min_index, truth)







#####################################################
########## Generate from Competition Model ##########
#####################################################


### Function to generate data from competition model
generate_data_TI <- function(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, N_A, N_B, N_AB, seed, time, basis_degree,
                             boundary_knots, internal_knots){
  set.seed(seed)
  n_A <- rep(NA, N_A)
  n_B <- rep(NA, N_B)
  n_AB <- rep(NA, N_AB)
  L_AB <- list()
  X_A <- list()
  X_B <- list()
  X_AB <- list()
  
  
  ## Generate A data
  for(i in 1:N_A){
    total_time <- 0
    spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
    spike <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
    total_time <- spike
    while(total_time < time){
      spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
      spike_i <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
      total_time <- total_time + spike_i
      spike <- c(spike, spike_i)
    }
    total_time <- total_time - spike[length(spike)]
    spike <- spike[-length(spike)]
    X_A[[i]] <- spike
    n_A[i] <- length(spike)
  }
  
  ## Generate B data
  for(i in 1:N_B){
    total_time <- 0
    spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
    spike <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
    total_time <- spike
    while(total_time < time){
      spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
      spike_i <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
      total_time <- total_time + spike_i
      spike <- c(spike, spike_i)
    }
    total_time <- total_time - spike[length(spike)]
    spike <- spike[-length(spike)]
    X_B[[i]] <- spike
    n_B[i] <- length(spike)
  }
  
  ## Generate AB data
  
  for(i in 1:N_AB){
    total_time <- 0
    spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
    spike_A <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
    spike_B <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
    if(spike_A < spike_B){
      spike <- spike_A
      L_i <- 0
      total_time <- spike
    }else{
      spike <- spike_B
      L_i <- 1
      total_time <- spike
    }
    
    iter <- 1
    while(total_time < time){
      spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
      spike_A <- rinvgauss(1, mean = 1 / (I_A * exp(spline %*% basis_coef_A)), shape = (1 / sigma_A)^2)
      spike_B <- rinvgauss(1, mean = 1 / (I_B * exp(spline %*% basis_coef_B)), shape = (1 / sigma_B)^2)
      if(L_i[iter] == 0){
        if(spike_A < spike_B + delta){
          spike <- c(spike, spike_A)
          L_i <- c(L_i, 0)
        }else{
          spike <- c(spike, spike_B + delta)
          L_i <- c(L_i, 1)
        }
      }
      if(L_i[iter] == 1){
        if(spike_A + delta < spike_B){
          spike <- c(spike, spike_A + delta)
          L_i <- c(L_i, 0)
        }else{
          spike <- c(spike, spike_B)
          L_i <- c(L_i, 1)
        }
      }
      iter <- iter + 1
      total_time <- sum(spike)
    }
    total_time <- total_time - spike[length(spike)]
    L_i <- L_i[-length(spike)]
    spike <- spike[-length(spike)]
    n_AB[i] <- length(spike)
    L_AB[[i]] <- L_i
    X_AB[[i]] <- spike
  }
  
  spline <- GetBSpline(seq(0, 1.0, 0.01), basis_degree, boundary_knots, internal_knots)
  return(list("X_A" = X_A, "X_B" = X_B, "X_AB" = X_AB, "L_AB" = L_AB, "n_A" = n_A,
              "n_B" = n_B, "n_AB" = n_AB, "A_CI" = (I_A * exp(spline %*% basis_coef_A)), "B_CI" = (I_B * exp(spline %*% basis_coef_B))))
}



run_sim_comp <- function(iter){
  set.seed(iter)
  
  I_A <- rtruncnorm(1, a = 0, mean = 40, sd = 6)
  I_B <- rtruncnorm(1, a = 0, mean = 80, sd = 6)
  sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
  sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  delta <- rlnorm(1, -2.5, 0.5)
  
  dat <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 25, 25, 25, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  
  res <- Sampler_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  ## Run IIGPP Model
  res_A <- Sampler_IIGPP(dat$X_A, dat$n_A, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_B <- Sampler_IIGPP(dat$X_B, dat$n_B, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_AB <- Sampler_IIGPP(dat$X_AB, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  dat_A_AB <- c(dat$X_A, dat$X_AB)
  n_A_AB <- c(dat$n_A, dat$n_AB)
  res_A_AB <- Sampler_IIGPP(dat_A_AB, n_A_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  dat_B_AB <- c(dat$X_B, dat$X_AB)
  n_B_AB <- c(dat$n_B, dat$n_AB)
  res_B_AB <- Sampler_IIGPP(dat_B_AB, n_B_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  WAIC_WTA_A <- WAIC_Winner_Take_All(dat_A_AB, dat$X_B, n_A_AB, dat$n_B, res_A_AB, res_B, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  WAIC_WTA_B <- WAIC_Winner_Take_All(dat$X_A, dat_B_AB, dat$n_A, n_B_AB, res_A, res_B_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  WAIC_Comp <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  waic_iigpp <- WAIC_IIGPP(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res_A, res_B, res_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  
  params <- list("I_A" = I_A, "I_B" = I_B, "sigma_A" = sigma_A, "sigma_B" = sigma_B,
                 "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B, "delta" = delta)
  output <- list("WAIC_Comp" = WAIC_Comp, "WAIC_IIGPP" = waic_iigpp, 
                 "WAIC_WTA_A" = WAIC_WTA_A, "WAIC_WTA_B" = WAIC_WTA_B, 
                 "res" = res, "res_A" = res_A, "res_B" = res_B, "res_AB" = res_AB, "res_A_AB" = res_A_AB, 
                 "res_B_AB" = res_B_AB, "data" = dat, "params" = params)
  saveRDS(output, paste0(save_dir, "/Comp/output", iter,".RDS"))
}

library(future.apply)

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
dir.create(paste0(save_dir, "/Comp"))
already_ran <- dir(paste0(save_dir, "/Comp"))
to_run <- which(!paste0("output", 1:100, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim_comp(this_seed))


### Post Processing
files <- dir(paste0(save_dir, "/Comp"))
WAIC <- matrix(0, nrow = length(files), 4)
wasserstein_dist <- rep(0, length(files))
for(i in 1:length(files)){
  output <- readRDS(paste0(save_dir, "/Comp/", files[i]))
  WAIC[i,1] <- output$WAIC_Comp$WAIC
  WAIC[i,2] <- output$WAIC_IIGPP$WAIC
  WAIC[i,3] <- output$WAIC_WTA_A$WAIC
  WAIC[i,4] <- output$WAIC_WTA_B$WAIC
  wasserstein_dist[i] <- wasserstein1d(output$data$n_A, output$data$n_B)
  print(i)
}

min_index <- rep(0, 100)
for(i in 1:length(files)){
  min_index[i] <- which.min(WAIC[i,])
}

table(min_index)


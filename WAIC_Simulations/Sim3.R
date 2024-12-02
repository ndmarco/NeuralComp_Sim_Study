## This simulation is to check which model is chosen under an oscillatory process
## that does not align with the specified multiplexing switching. Specifically, we
## will simulate data that randomly switches between the A and B process. 

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

save_dir <- "."

################################################################################
################################################################################
###################### Sim from random switching Model #########################
################################################################################
################################################################################

### Function to generate data from a model that randomly switches between the two processes
generate_data_random_switching <- function(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, 
                              prob_stay, N_A, N_B, N_AB, seed, time, basis_degree,
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
    if(runif(1) > 0.5){
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
      if(runif(1) > prob_stay){
        L_i = c(L_i, abs(L_i[iter] - 1))
        if(abs(L_i[iter] - 1) == 0){
          spike <- c(spike, spike_A)
        }else{
          spike <- c(spike, spike_B)
        }
      }else{
        L_i = c(L_i, L_i[iter])
        if(L_i[iter] == 0){
          spike <- c(spike, spike_A)
        }else{
          spike <- c(spike, spike_B)
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

run_sim4 <- function(iter){
  set.seed(iter)
  
  I_A <- rtruncnorm(1, a = 0, mean = 40, sd = 4)
  I_B <- rtruncnorm(1, a = 0, mean = 80, sd = 4)
  sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
  sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
  prob_stay <- rbeta(1, 10, 2)

  
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  dat <- generate_data_random_switching(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, prob_stay, 25, 25, 25, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  
  ## Run Competition Model 
  res <- Sampler_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  ## Run IIGPP Model
  res_A <- Sampler_IIGPP(dat$X_A, dat$n_A, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_B <- Sampler_IIGPP(dat$X_B, dat$n_B, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_AB <- Sampler_IIGPP(dat$X_AB, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  WAIC_Comp <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  waic_iigpp <- WAIC_IIGPP(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res_A, res_B, res_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)

  params <- list("I_A" = I_A, "I_B" = I_B, "sigma_A" = sigma_A, "sigma_B" = sigma_B, "prob_stay" = prob_stay,
                 "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B)
  output <- list("WAIC_Comp" = WAIC_Comp, "WAIC_IIGPP" = waic_iigpp, 
                 "res" = res, "res_A" = res_A, "res_B" = res_B, "res_AB" = res_AB, "data" = dat,
                 "params" = params)
  saveRDS(output, paste0(save_dir, "/output", iter,".RDS"))
}



ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
already_ran <- dir(save_dir)
to_run <- which(!paste0("output", 1:100, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim4(this_seed))

files <- dir(save_dir)
WAIC <- matrix(0, nrow = length(files), 2)
lppd <- matrix(0, nrow = length(files), 2)
p_theta <- matrix(0, nrow = length(files), 2)
prob_stay <- rep(0, length(files))
switches <- rep(0, length(files))
dist_FR <- rep(0, length(files))
rel_A_B <- rep(0, length(files))
AB_dist_scale <- rep(0, length(files))
number_changes <- function(labels, n){
  output <- 0
  for(i in 1:n){
    for(j in 2:length(labels[[i]])){
      if(labels[[i]][j-1] != labels[[i]][j]){
        output <- output + 1
      }
    }
  }
  return(output)
}

for(i in 1:length(files)){
  output <- readRDS(paste0(save_dir, "/", files[i]))
  WAIC[i,1] <- output$WAIC_Comp$WAIC
  WAIC[i,2] <- output$WAIC_IIGPP$WAIC
  lppd[i,1] <- output$WAIC_Comp$lppd
  lppd[i,2] <- output$WAIC_IIGPP$lppd
  p_theta[i,1] <- output$WAIC_Comp$Effective_pars
  p_theta[i,2] <- output$WAIC_IIGPP$Effective_pars
  prob_stay[i] <- output$params$prob_stay
  switches[i] <- number_changes(output$data$L_AB, length(output$data$n_AB))
  rel_A_B[i] <- (mean(output$data$n_AB) - mean(output$data$n_A)) / (mean(output$data$n_B) - mean(output$data$n_A))
  AB_dist_scale[i] <- abs(output$params$sigma_A - output$params$sigma_B)
  print(i)
}

# for(i in 1:length(files)){
#   output <- readRDS(paste0(save_dir, "/", files[i]))
#   numbers <- gregexpr("[0-9]+", files[i])
#   result <- as.numeric(regmatches(files[i], numbers))
#   set.seed(result)
#   I_A <- rtruncnorm(1, a = 0, mean = 40, sd = 4)
#   I_B <- rtruncnorm(1, a = 0, mean = 80, sd = 4)
#   sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
#   sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
#   prob_stay <- rbeta(1, 10, 2)
#   output$prob_stay <- prob_stay
#   saveRDS(output, paste0(save_dir, "/", files[i]))
#   print(i)
# }

time_df <- matrix(0, length(files), 3)
time_df[1:length(files),1] <- WAIC[,1] - WAIC[,2]
time_df[1:length(files),2] <- "WAIC Comp Conditional"
time_df[1:length(files),3] <- prob_stay
time_df <- as.data.frame(time_df)
colnames(time_df) <- c("WAIC", "Method", "Prob_Stay")
time_df$WAIC <- as.numeric(time_df$WAIC)
time_df$Method <- as.factor(time_df$Method)
time_df$Prob_Stay <- as.numeric(time_df$Prob_Stay)
p1 <- ggplot(time_df, aes(x=Prob_Stay, y=WAIC)) + scale_y_continuous(breaks = c(-225,  0,  225, 900, 2500, 5000, 7500, 10000, 12500), trans = ssqrt_trans) + ylab(TeX("WAIC$_{comp}$ - WAIC$_{IIGPP}$"))+
  geom_point() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab(TeX("$p_s$")) +
  geom_hline(yintercept = 0,  color="red", linetype="dashed")

### Simulation to compare accuracy and speed of various methods to estimate Conditional WAIC
library(NeuralComp)
library(MASS)
library(truncnorm)
library(ggplot2)
library(latex2exp)
library(truncnorm)
library(future.apply)
library(scales)
library(transport)
library(ggallin)
library(statmod)


save_dir <- "."

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
    n_AB[i] <- length(spike)
    L_AB[[i]] <- L_i
    X_AB[[i]] <- spike
  }
  
  spline <- GetBSpline(seq(0, 1.0, 0.01), basis_degree, boundary_knots, internal_knots)
  return(list("X_A" = X_A, "X_B" = X_B, "X_AB" = X_AB, "L_AB" = L_AB, "n_A" = n_A,
              "n_B" = n_B, "n_AB" = n_AB, "A_CI" = (I_A * exp(spline %*% basis_coef_A)), "B_CI" = (I_B * exp(spline %*% basis_coef_B))))
}


run_sim1 <- function(iter){
  set.seed(iter)
  
  I_A <- rtruncnorm(1, a = 0, mean = 40, sd = 4)
  I_B <- rtruncnorm(1, a = 0, mean = 80, sd = 4)
  sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
  sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
  delta <- rlnorm(1, -2.5, 0.5)
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  dat <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 50, 50, 50, 1, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  
  ## Run Competition Model 
  res <- Sampler_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, 2000, 3, c(0,1), c(0.25, 0.5, 0.75))
  
  time_1_start <- Sys.time();
  WAIC_Comp_val <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), method = "numerical_approx")
  time_1_end <- Sys.time();
  
  time_2_start <- Sys.time();
  WAIC_Comp_val2 <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), method = "numerical_approx", n_spike_evals = 50)
  time_2_end <- Sys.time();
  
  time_5_start <- Sys.time();
  WAIC_Comp_val3 <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), method = "numerical_approx", n_spike_evals = 100)
  time_5_end <- Sys.time();
  
  time_3_start <- Sys.time();
  WAIC_Comp_dir <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75))
  time_3_end <- Sys.time();
  
  time_4_start <- Sys.time();
  WAIC_Comp_approx <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), method = "sampling", n_MCMC_approx2 = 40)
  time_4_end <- Sys.time();
  
  params <- list("I_A" = I_A, "I_B" = I_B, "sigma_A" = sigma_A, "sigma_B" = sigma_B, "delta" = delta,
                 "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B)
  output <- list("WAIC_Comp1" = WAIC_Comp_val, "WAIC_Comp2" = WAIC_Comp_val2, "WAIC_Comp3" = WAIC_Comp_val3,  "WAIC_Comp_dir" = WAIC_Comp_dir, "WAIC_Comp_approx" = WAIC_Comp_approx, "dat" = dat,
                 "res" = res, "Comp1_time" = time_1_end - time_1_start, "Comp2_time" = time_2_end - time_2_start, "Comp3_time" = time_5_end - time_5_start, "Comp_dir_time" = time_3_end - time_3_start,
                 "Comp_approx_time" = time_4_end - time_4_start, "params" = params)
  saveRDS(output, paste0(save_dir, "/output", iter,".RDS"))
}

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
already_ran <- dir(save_dir)
to_run <- which(!paste0("output", 1:50, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim1(this_seed))

## Post processing
files <- dir(save_dir)
WAIC <- matrix(0, nrow = length(files), 5)
llpd <- matrix(0, nrow = length(files), 5)
p_theta <- matrix(0, nrow = length(files), 5)
time <- matrix(0, nrow = length(files), 5)
for(i in 1:length(files)){
  output <- readRDS(paste0(save_dir, "/", files[i]))
  WAIC[i,1] <- output$WAIC_Comp1$WAIC
  WAIC[i,2] <- output$WAIC_Comp2$WAIC
  WAIC[i,3] <- output$WAIC_Comp3$WAIC
  WAIC[i,4] <- output$WAIC_Comp_dir$WAIC
  WAIC[i,5] <- output$WAIC_Comp_approx$WAIC
  llpd[i,1] <- output$WAIC_Comp1$llpd
  llpd[i,2] <- output$WAIC_Comp2$llpd
  llpd[i,3] <- output$WAIC_Comp3$llpd
  llpd[i,4] <- output$WAIC_Comp_dir$llpd
  llpd[i,5] <- output$WAIC_Comp_approx$llpd
  p_theta[i,1] <- output$WAIC_Comp1$Effective_pars
  p_theta[i,2] <- output$WAIC_Comp2$Effective_pars
  p_theta[i,3] <- output$WAIC_Comp3$Effective_pars
  p_theta[i,4] <- output$WAIC_Comp_dir$Effective_pars
  p_theta[i,5] <- output$WAIC_Comp_approx$Effective_pars
  time[i,1] <- output$Comp1_time
  time[i,2] <- output$Comp2_time
  time[i,3] <- output$Comp3_time
  time[i,4] <- output$Comp_dir_time
  time[i,5] <- output$Comp_approx_time
  print(i)
}


waic_df <- matrix(0, length(files)*4, 2)
waic_df[1:length(files),1] <- (WAIC[,1] - WAIC[,3])/ WAIC[,3]
waic_df[1:length(files),2] <- 1
waic_df[(length(files) + 1):(2*length(files)),1] <- (WAIC[,2] - WAIC[,3])/ WAIC[,3]
waic_df[(length(files) + 1):(2*length(files)),2] <- 2
waic_df[(2*length(files) + 1):(3*length(files)),1] <- (WAIC[,4] - WAIC[,3])/ WAIC[,3]
waic_df[(2*length(files) + 1):(3*length(files)),2] <- 4
waic_df[(3*length(files) + 1):(4*length(files)),1] <- (WAIC[,5] - WAIC[,3])/ WAIC[,3]
waic_df[(3*length(files) + 1):(4*length(files)),2] <- 5
waic_df <- as.data.frame(waic_df)
colnames(waic_df) <- c("Relative_Error", "Method")
waic_df$Method <- as.factor(waic_df$Method)
p1 <- ggplot(waic_df, aes(x=Method, y=`Relative_Error`)) + scale_y_continuous(labels = scales::percent) + ggtitle(TeX("WAIC")) + ylab("Relative Error") + 
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))


llpd_df <- matrix(0, length(files)*4, 2)
llpd_df[1:length(files),1] <- (llpd[,1] - llpd[,3])/ llpd[,3]
llpd_df[1:length(files),2] <- 1
llpd_df[(length(files) + 1):(2*length(files)),1] <- (llpd[,2] - llpd[,3])/ llpd[,3]
llpd_df[(length(files) + 1):(2*length(files)),2] <- 2
llpd_df[(2*length(files) + 1):(3*length(files)),1] <- (llpd[,4] - llpd[,3])/ llpd[,3]
llpd_df[(2*length(files) + 1):(3*length(files)),2] <- 4
llpd_df[(3*length(files) + 1):(4*length(files)),1] <- (llpd[,5] - llpd[,3])/ llpd[,3]
llpd_df[(3*length(files) + 1):(4*length(files)),2] <- 5
llpd_df <- as.data.frame(llpd_df)
colnames(llpd_df) <- c("Relative_Error", "Method")
llpd_df$Method <- as.factor(llpd_df$Method)
p2 <- ggplot(llpd_df, aes(x=Method, y=`Relative_Error`)) + scale_y_continuous(labels = scales::percent) + ggtitle(TeX("LPPD")) + ylab("Relative Error") + 
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))


p_theta_df <- matrix(0, length(files)*4, 2)
p_theta_df[1:length(files),1] <- (p_theta[,1] - p_theta[,3])/ p_theta[,3]
p_theta_df[1:length(files),2] <- 1
p_theta_df[(length(files) + 1):(2*length(files)),1] <- (p_theta[,2] - p_theta[,3])/ p_theta[,3]
p_theta_df[(length(files) + 1):(2*length(files)),2] <- 2
p_theta_df[(2*length(files) + 1):(3*length(files)),1] <- (p_theta[,4] - p_theta[,3])/ p_theta[,3]
p_theta_df[(2*length(files) + 1):(3*length(files)),2] <- 4
p_theta_df[(3*length(files) + 1):(4*length(files)),1] <- (p_theta[,5] - p_theta[,3])/ p_theta[,3]
p_theta_df[(3*length(files) + 1):(4*length(files)),2] <- 5
p_theta_df <- as.data.frame(p_theta_df)
colnames(p_theta_df) <- c("Relative_Error", "Method")
p_theta_df$Method <- as.factor(p_theta_df$Method)
p3 <- ggplot(p_theta_df, aes(x=Method, y=`Relative_Error`)) + scale_y_continuous(labels = scales::percent) + ggtitle(TeX("$p_{\\theta}$")) + ylab("Relative Error") + 
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

time_df <- matrix(0, length(files)*5, 2)
time_df[1:length(files),1] <- time[,1]
time_df[1:length(files),2] <- 1
time_df[(length(files) + 1):(2*length(files)),1] <- time[,2]
time_df[(length(files) + 1):(2*length(files)),2] <- 2
time_df[(2*length(files) + 1):(3*length(files)),1] <- time[,3]
time_df[(2*length(files) + 1):(3*length(files)),2] <- 3
time_df[(3*length(files) + 1):(4*length(files)),1] <- time[,4]
time_df[(3*length(files) + 1):(4*length(files)),2] <- 4
time_df[(4*length(files) + 1):(5*length(files)),1] <- time[,5]
time_df[(4*length(files) + 1):(5*length(files)),2] <- 5
time_df <- as.data.frame(time_df)
colnames(time_df) <- c("Time", "Method")
time_df$Method <- as.factor(time_df$Method)
p4 <- ggplot(time_df, aes(x=Method, y=`Time`)) + ggtitle(TeX("Time")) + ylab("Time (Minutes)") + 
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))
library(gridExtra)
grid.arrange(p1,p2,p3,p4, ncol = 2)
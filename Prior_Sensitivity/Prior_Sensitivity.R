### Simulation to examine how well we can recover various parameters of interest as sample size increases
library(NeuralComp)
library(MASS)
library(statmod)
library(truncnorm)
library(ggplot2)
library(latex2exp)
library(truncnorm)
library(future.apply)
library(ggpubr)
library(grid)
library(transport)
library(ggallin)
library(gridExtra)

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

run_sim <- function(iter, I_A_mean1, I_B_mean1, sigma_A_mean1, sigma_B_mean1, delta_shape1, delta_rate1, prior_set){
  set.seed(iter)
  I_A <- rtruncnorm(1, a = 10, mean = 40, sd = 20)
  I_B <- rtruncnorm(1, a = 10, mean = 80, sd = 20)
  sigma_A <- rtruncnorm(1, a = 3, mean = sqrt(40), sd = 5)
  sigma_B <- rtruncnorm(1, a = 3, mean = sqrt(80), sd = 5)
  delta <- rlnorm(1, -4, 1)
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  dat <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 25, 25, 25, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  
  ## Run Competition Model 
  res <- Sampler_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, 
                             10000, 3, c(0,1), c(0.25, 0.5, 0.75), 1, I_A_mean = I_A_mean1, 
                             I_B_mean = I_B_mean1, sigma_A_mean = sigma_A_mean1, 
                             sigma_B_mean = sigma_B_mean1, delta_shape = delta_shape1,
                             delta_rate = delta_rate1)
  
  CI <- FR_CI_Competition(seq(0,1,0.01), 3, c(0, 1), c(0.25, 0.5, 0.75), res, burnin_prop = 0.33)
  
  coverage_A <- sum((CI$A_FR_CI[,1] < dat$A_CI) & (dat$A_CI < CI$A_FR_CI[,2])) / length(dat$A_CI)
  coverage_B <- sum((CI$B_FR_CI[,1] < dat$B_CI) & (dat$B_CI < CI$B_FR_CI[,2])) / length(dat$B_CI)
  
  rel_error_A <- sum((dat$A_CI - CI$A_FR_median)^2) / sum(dat$A_CI^2)
  rel_error_B <- sum((dat$B_CI - CI$B_FR_median)^2) / sum(dat$B_CI^2)
  area_CI_A <- sum((CI$A_FR_CI[,1] - CI$A_FR_CI[,2])^2)
  area_CI_B <- sum((CI$B_FR_CI[,1] - CI$B_FR_CI[,2])^2)
  
  CI_sigma_A <- quantile(res$theta[4501:14500,3], c(0.025, 0.5, 0.975))
  CI_sigma_B <- quantile(res$theta[4501:14500,4], c(0.025, 0.5, 0.975))
  rel_error_sigma_A <- (CI_sigma_A[2] - sigma_A)^2 / (sigma_A^2)
  rel_error_sigma_B <- (CI_sigma_B[2] - sigma_B)^2 / (sigma_B^2)
  CI_width_sigma_A <- CI_sigma_A[3] - CI_sigma_A[1]
  CI_width_sigma_B <- CI_sigma_B[3] - CI_sigma_B[1]
  coverage_sigma_A <- (CI_sigma_A[1] < sigma_A) & (sigma_A < CI_sigma_A[3])
  coverage_sigma_B <- (CI_sigma_B[1] < sigma_B) & (sigma_B < CI_sigma_B[3])
  
  CI_delta <- quantile(res$theta[4501:14500,5], c(0.025, 0.5, 0.975))
  max_isi <- max(unlist(dat$X_AB))
  if(delta > max_isi){
    rel_error_delta <- NA
    CI_width_delta <- NA
    coverage_delta <- NA
  }else{
    rel_error_delta <- (CI_delta[2] - delta)^2 / (delta^2)
    CI_width_delta <- CI_delta[3] - CI_delta[1]
    coverage_delta <- (CI_delta[1] < delta) & (delta < CI_delta[3])
  }
  
  
  ### posterior_predictive
  post_pred_sample <- Competition_Posterior_Predictive(1, 3, c(0,1), c(0.25, 0.5, 0.75), res, burnin_prop =  0.33, n_samples = 100000)
  prob_switches_sample <- table(post_pred_sample$n_switches) / length(post_pred_sample$n_switches)
  density_spikes_AB_sample <- density(post_pred_sample$n_AB, from = min(post_pred_sample$n_AB) - 10, to = max(post_pred_sample$n_AB) + 10, n = 512)
  
  prop_time_A_sample <- rep(0, length(post_pred_sample$n_AB))
  for(i in 1:length(post_pred_sample$n_AB)){
    prop_time_A_sample[i] = sum(post_pred_sample$posterior_pred_samples_AB[[i]][post_pred_sample$posterior_pred_labels[[i]] == 0]) / sum(post_pred_sample$posterior_pred_samples_AB[[i]])
  }
  
  
  ## From true distribution
  res1 <- res
  for(i in 1:(length(res1$theta[,1]))){
    res1$theta[i,] <- c(I_A, I_B, sigma_A, sigma_B, delta)
    res1$basis_coef_A[i,] <- basis_coef_A
    res1$basis_coef_B[i,] <- basis_coef_B
  }
  
  post_pred_true <- Competition_Posterior_Predictive(1, 3, c(0,1), c(0.25, 0.5, 0.75), res1, burnin_prop =  0, n_samples = 100000)
  prob_switches_true <- table(post_pred_true$n_switches) / length(post_pred_true$n_switches)
  
  prop_time_A_true <- rep(0, length(post_pred_true$n_AB))
  for(i in 1:length(post_pred_sample$n_AB)){
    prop_time_A_true[i] = sum(post_pred_true$posterior_pred_samples_AB[[i]][post_pred_true$posterior_pred_labels[[i]] == 0]) / sum(post_pred_true$posterior_pred_samples_AB[[i]])
  }
  
  
  # Distance from posterior predictive distributions
  rel_dist_n_switches <- sum(abs(post_pred_sample$n_switches[order(post_pred_sample$n_switches)] - post_pred_true$n_switches[order(post_pred_true$n_switches)]))/ abs(sum(post_pred_true$n_switches))
  rel_dist_prop_A <- sum(abs(prop_time_A_sample[order(prop_time_A_sample)] - prop_time_A_true[order(prop_time_A_true)]))/ abs(sum(prop_time_A_true))
  rel_dist_n_AB <- sum(abs(post_pred_sample$n_AB[order(post_pred_sample$n_AB)] - post_pred_true$n_AB[order(post_pred_true$n_AB)]))/ abs(sum(post_pred_true$n_AB))
  
  params <- list("I_A" = I_A, "I_B" = I_B, "sigma_A" = sigma_A, "sigma_B" = sigma_B,
                 "delta" = delta, "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B)
  output <- list("res" = res, "data" = dat, "coverage_A_FR" = coverage_A , "coverage_B_FR" = coverage_B,
                 "rel_error_A_FR" = rel_error_A, "rel_error_B_FR" = rel_error_B, "area_CI_A_FR" = area_CI_A,
                 "area_CI_B_FR" = area_CI_B, "rel_error_sigma_A" = rel_error_sigma_A, "rel_error_sigma_B" = rel_error_sigma_B,
                 "coverage_sigma_A" = coverage_sigma_A, "coverage_sigma_B" = coverage_sigma_B,
                 "CI_width_sigma_A" = CI_width_sigma_A, "CI_width_sigma_B" = CI_width_sigma_B,
                 "rel_error_delta" = rel_error_delta, "CI_width_delta" = CI_width_delta, 
                 "coverage_delta" = coverage_delta, "rel_dist_n_switches" = rel_dist_n_switches, "rel_dist_prop_A" = rel_dist_prop_A,
                 "rel_dist_n_AB" = rel_dist_n_AB, "params" = params)
  saveRDS(output, paste0(save_dir, "/", prior_set, "/output", iter,".RDS"))
}

I_A_mean1 = c(20, 40, 80)
I_B_mean1 = c(20, 40, 80)
sigma_A_mean1 = c(sqrt(20), sqrt(40), sqrt(80))
sigma_B_mean1 = c(sqrt(20), sqrt(40), sqrt(80))
delta_shape1 = c(0.1, 0.01, 0.001)
delta_rate1 = c(0.1, 0.1, 0.01)
prior_set = c("A", "B", "C")
for(i in 1:3){
  ncpu <- min(5, availableCores())
  plan(multisession, workers = ncpu)
  dir.create(paste0(save_dir, "/", prior_set[i])) 
  already_ran <- dir(paste0(save_dir, "/", prior_set[i]))
  to_run <- which(!paste0("output", 1:25, ".RDS") %in% already_ran)
  future_lapply(to_run, function(this_seed) run_sim(this_seed, I_A_mean1[i], I_B_mean1[i], sigma_A_mean1[i], sigma_B_mean1[i], delta_shape1[i], delta_rate1[i], prior_set[i]))
}



prior_set = c("A", "B", "C")
files <- dir(paste0(save_dir, "/", prior_set[1]))
coverage_A_FR <- matrix(NA, nrow = length(files), ncol = length(prior_set))
coverage_B_FR <- matrix(NA, nrow = length(files), ncol = length(prior_set))
rel_error_A <- matrix(NA, nrow = length(files), ncol = length(prior_set))
rel_error_B <- matrix(NA, nrow = length(files), ncol = length(prior_set))
area_CI_A <- matrix(NA, nrow = length(files), ncol = length(prior_set))
area_CI_B <- matrix(NA, nrow = length(files), ncol = length(prior_set))
rel_error_sigma_A <- matrix(NA, nrow = length(files), ncol = length(prior_set))
rel_error_sigma_B <- matrix(NA, nrow = length(files), ncol = length(prior_set))
CI_width_sigma_A = matrix(NA, nrow = length(files), ncol = length(prior_set))
CI_width_sigma_B = matrix(NA, nrow = length(files), ncol = length(prior_set))
coverage_sigma_A <- matrix(NA, nrow = length(files), ncol = length(prior_set))
coverage_sigma_B <- matrix(NA, nrow = length(files), ncol = length(prior_set))
rel_error_delta <- matrix(NA, nrow = length(files), ncol = length(prior_set))
CI_width_delta <- matrix(NA, nrow = length(files), ncol = length(prior_set))
coverage_delta <- matrix(NA, nrow = length(files), ncol = length(prior_set))
dist_switches <- matrix(NA, nrow = length(files), ncol = length(prior_set))
dist_prop_A <- matrix(NA, nrow = length(files), ncol = length(prior_set))
dist_n_AB <- matrix(NA, nrow = length(files), ncol = length(prior_set))
for(j in 1:length(prior_set)){
  for(i in 1:length(files)){
    output <- readRDS(paste0(save_dir, "/", prior_set[j], "/", files[i]))
    coverage_A_FR[i,j] <- output$coverage_A_FR
    coverage_B_FR[i,j] <- output$coverage_B_FR
    rel_error_A[i,j] <- output$rel_error_A
    rel_error_B[i,j] <- output$rel_error_B
    area_CI_A[i,j] <- output$area_CI_A
    area_CI_B[i,j] <- output$area_CI_B
    rel_error_sigma_A[i,j] <- output$rel_error_sigma_A
    rel_error_sigma_B[i,j] <- output$rel_error_sigma_B
    coverage_sigma_A[i,j] <- output$coverage_sigma_A
    coverage_sigma_B[i,j] <- output$coverage_sigma_B
    CI_sigma_A <- quantile(output$res$theta[4501:14500,3], c(0.025, 0.5, 0.975))
    CI_sigma_B <- quantile(output$res$theta[4501:14500,4], c(0.025, 0.5, 0.975))
    CI_width_sigma_A[i,j] <- CI_sigma_A[3] - CI_sigma_A[1]
    CI_width_sigma_B[i,j] <- CI_sigma_B[3] - CI_sigma_B[1]
    rel_error_delta[i,j] <- output$rel_error_delta
    CI_width_delta[i,j] <- output$CI_width_delta
    coverage_delta[i,j] <- output$coverage_delta
    dist_switches[i,j] <- output$rel_dist_n_switches
    dist_prop_A[i,j] <- output$rel_dist_prop_A
    dist_n_AB[i,j] <- output$rel_dist_n_AB
    print(i)
  }
}


# Plot R-MISE for the
FR_RSE <- matrix(0, length(files) * length(prior_set), 2)
FR_RSE[1:length(files),1] <- (rel_error_A[,1] + rel_error_B[,1]) / 2
FR_RSE[1:length(files),2] <- "A"
FR_RSE[(length(files) + 1):(2 * length(files)),1] <- (rel_error_A[,2] + rel_error_B[,2]) / 2
FR_RSE[(length(files) + 1):(2 * length(files)),2] <- "B"
FR_RSE[(2*length(files) + 1):(3 * length(files)),1] <- (rel_error_A[,3] + rel_error_B[,3]) / 2
FR_RSE[(2*length(files) + 1):(3 * length(files)),2] <- "C"
FR_RSE <- as.data.frame(FR_RSE)
colnames(FR_RSE) <- c("RSE", "Prior")
FR_RSE$RSE <- as.numeric(FR_RSE$RSE)
FR_RSE$Prior <- as.factor(FR_RSE$Prior)
p1 <- ggplot(FR_RSE, aes(x=Prior, y=`RSE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.1, 0.01, 0.001, 0.0001))  + 
  ggtitle(TeX("$I\\exp\\{\\phi' b(t) \\}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))+ xlab("Prior Set")

sigma_RSE <- matrix(0, length(files) * length(prior_set), 2)
sigma_RSE[1:length(files),1] <- (rel_error_sigma_A[,1] + rel_error_sigma_B[,1]) / 2
sigma_RSE[1:length(files),2] <- "A"
sigma_RSE[(length(files) + 1):(2 * length(files)),1] <- (rel_error_sigma_A[,2] + rel_error_sigma_B[,2]) / 2
sigma_RSE[(length(files) + 1):(2 * length(files)),2] <- "B"
sigma_RSE[(2*length(files) + 1):(3 * length(files)),1] <- (rel_error_sigma_A[,3] + rel_error_sigma_B[,3]) / 2
sigma_RSE[(2*length(files) + 1):(3 * length(files)),2] <- "C"
sigma_RSE <- as.data.frame(sigma_RSE)
colnames(sigma_RSE) <- c("RSE", "Prior")
sigma_RSE$Prior <- as.factor(sigma_RSE$Prior)
sigma_RSE$RSE <- as.numeric(sigma_RSE$RSE)
p2 <- ggplot(sigma_RSE, aes(x=Prior, y=`RSE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.01, 0.001, 0.0001, 0.00001, 0.000001))  + 
  ggtitle(TeX("$\\sigma$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))+ xlab("Prior Set")

delta_RSE <- matrix(0, length(files) * length(prior_set), 2)
delta_RSE[1:length(files),1] <- rel_error_delta[,1]
delta_RSE[1:length(files),2] <- "A"
delta_RSE[(length(files) + 1):(2 * length(files)),1] <- rel_error_delta[,2]
delta_RSE[(length(files) + 1):(2 * length(files)),2] <- "B"
delta_RSE[(2*length(files) + 1):(3 * length(files)),1] <- rel_error_delta[,3]
delta_RSE[(2*length(files) + 1):(3 * length(files)),2] <- "C"
delta_RSE <- as.data.frame(delta_RSE)
colnames(delta_RSE) <- c("RSE", "Prior")
delta_RSE$Prior <- as.factor(delta_RSE$Prior)
delta_RSE$RSE <- as.numeric(delta_RSE$RSE)
p3 <- ggplot(delta_RSE, aes(x=Prior, y=`RSE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(1, 0.01, 0.0001, 0.000001, 0.00000001))  + 
  ggtitle(TeX("$\\delta$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))+ xlab("Prior Set")

options(scipen = 999)
df <- matrix(0, length(files) * length(prior_set), 2)
df[1:length(files),1] <- dist_switches[,1]
df[1:length(files),2] <- "A"
df[(length(files) + 1):(2 * length(files)),1] <- dist_switches[,2]
df[(length(files) + 1):(2 * length(files)),2] <- "B"
df[(2*length(files) + 1):(3 * length(files)),1] <- dist_switches[,3]
df[(2*length(files) + 1):(3 * length(files)),2] <- "C"
df <- as.data.frame(df)
colnames(df) <- c("Wasserstein Distance", "Prior")
df$Prior <- as.factor(df$Prior)
df$`Wasserstein Distance` <- as.numeric(df$`Wasserstein Distance`)
p4 <- ggplot(df, aes(x=Prior, y=`Wasserstein Distance`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(1, 0.1, 0.01, 0.001))  + 
  ggtitle(TeX("Distribution of the # of Switches")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))+ xlab("Prior Set")


df <- matrix(0, length(files) * length(prior_set), 2)
df[1:length(files),1] <- dist_prop_A[,1]
df[1:length(files),2] <- "A"
df[(length(files) + 1):(2 * length(files)),1] <- dist_prop_A[,2]
df[(length(files) + 1):(2 * length(files)),2] <- "B"
df[(2*length(files) + 1):(3 * length(files)),1] <- dist_prop_A[,3]
df[(2*length(files) + 1):(3 * length(files)),2] <- "C"
df <- as.data.frame(df)
colnames(df) <- c("Wasserstein Distance", "Prior")
df$Prior <- as.factor(df$Prior)
df$`Wasserstein Distance` <- as.numeric(df$`Wasserstein Distance`)
p5 <- ggplot(df, aes(x=Prior, y=`Wasserstein Distance`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(1, 0.1, 0.01, 0.001))  + 
  ggtitle(TeX("Distribution of Time Encoding for A")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))+ xlab("Prior Set")


df <- matrix(0, length(files) * length(prior_set), 2)
df[1:length(files),1] <- dist_n_AB[,1]
df[1:length(files),2] <- "A"
df[(length(files) + 1):(2 * length(files)),1] <- dist_n_AB[,2]
df[(length(files) + 1):(2 * length(files)),2] <- "B"
df[(2*length(files) + 1):(3 * length(files)),1] <- dist_n_AB[,3]
df[(2*length(files) + 1):(3 * length(files)),2] <- "C"
df <- as.data.frame(df)
colnames(df) <- c("Wasserstein Distance", "Prior")
df$Prior <- as.factor(df$Prior)
df$`Wasserstein Distance` <- as.numeric(df$`Wasserstein Distance`)
p6 <- ggplot(df, aes(x=Prior, y=`Wasserstein Distance`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.1, 0.01,0.001,0.0001))  + 
  ggtitle(TeX("Distribution of Spike Counts in AB Trials")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Prior Set")

### Note some points may be removed in p3 if max_ISI < delta (see simulation code above) -- This is expected
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
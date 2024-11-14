### Simulation to examine how WAIC performs when some of the trials come from the competition framework,
### and the rest come from the IIGPP model
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
library(gridExtra)

save_dir <- "/Users/ndm34/Documents/Post_Pval_Sim"



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

run_sim6 <- function(iter){
  set.seed(iter)
  
  I_A <- rtruncnorm(1, a = 0, mean = 40, sd = 4)
  I_B <- rtruncnorm(1, a = 0, mean = 80, sd = 4)
  sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
  sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
  if(iter <= 80){
    delta <- rlnorm(1, -2.5, 0.5)
  }else{
    delta <- 0
  }
  
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  dat <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 50, 50, 50, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))

  
  ## Run IIGPP Model
  res_A <- Sampler_IIGPP(dat$X_A, dat$n_A, 2000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_B <- Sampler_IIGPP(dat$X_B, dat$n_B, 2000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_AB <- Sampler_IIGPP(dat$X_AB, dat$n_AB, 2000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  IIGPP_Test_AB <- Test_IIGPP_Fit(dat$X_AB, dat$n_AB, res_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  IIGPP_Test_A <- Test_IIGPP_Fit(dat$X_A, dat$n_A, res_A, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  IIGPP_Test_B <- Test_IIGPP_Fit(dat$X_B, dat$n_B, res_B, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  
  params <- list("I_A" = I_A, "I_B" = I_B, "sigma_A" = sigma_A, "sigma_B" = sigma_B, "delta" = delta,
                 "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B)
  output <- list("IIGPP_Test_A" = IIGPP_Test_A, "IIGPP_Test_B" = IIGPP_Test_B, "IIGPP_Test_AB" = IIGPP_Test_AB, 
                 "res_A" = res_A, "res_B" = res_B, "res_AB" = res_AB, "data" = dat,
                 "params" = params)
  saveRDS(output, paste0(save_dir, "/output", iter,".RDS"))
}

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
already_ran <- dir(save_dir)
to_run <- which(!paste0("output", 1:100, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim6(this_seed))


### Post Processing
files <- dir(save_dir)
IIGPP_Test <- matrix(0, nrow = length(files), 3)
mean_test <- matrix(0, nrow = length(files), 3)
var_test <- matrix(0, nrow = length(files), 3)
wasserstein_dist <- rep(0, length(files))
delta <- rep(0, length(files))
for(i in 1:length(files)){
  output <- readRDS(paste0(save_dir, "/", files[i]))
  IIGPP_Test[i,1] <- output$IIGPP_Test_A$p_val
  IIGPP_Test[i,2] <- output$IIGPP_Test_B$p_val
  IIGPP_Test[i,3] <- output$IIGPP_Test_AB$p_val
  mean_test[i,1] <- output$IIGPP_Test_A$p_val_mean_SC
  mean_test[i,2] <- output$IIGPP_Test_B$p_val_mean_SC
  mean_test[i,3] <- output$IIGPP_Test_AB$p_val_mean_SC
  var_test[i,1] <- output$IIGPP_Test_A$p_val_var_SC
  var_test[i,2] <- output$IIGPP_Test_B$p_val_var_SC
  var_test[i,3] <- output$IIGPP_Test_AB$p_val_var_SC
  wasserstein_dist[i] <- wasserstein1d(output$data$n_A, output$data$n_B)
  delta[i] <- output$params$delta
  print(i)
}

df <- matrix(0, length(files)*3, 3)
df[1:length(files),1] <- IIGPP_Test[,1]
df[1:length(files),2] <-  "A"
df[1:length(files),3] <- delta
df[(length(files) + 1):(2*length(files)),1] <- IIGPP_Test[,2]
df[(length(files) + 1):(2*length(files)),2] <- "B"
df[(length(files) + 1):(2*length(files)),3] <- delta
df[(2*length(files) + 1):(3*length(files)),1] <- IIGPP_Test[,3]
df[(2*length(files) + 1):(3*length(files)),2] <- "AB"
df[(2*length(files) + 1):(3*length(files)),3] <- delta
df <- as.data.frame(df)
colnames(df) <- c("Pval", "Condition", "delta")
df$Pval <- as.numeric(df$Pval)
df$Condition <- as.factor(df$Condition)
df$delta <- as.numeric(df$delta)
IIGPP_plot <- ggplot(df, aes(x=delta, y=Pval, colour = Condition)) + ylab(TeX("P-value"))+
  geom_point(aes(colour = Condition)) +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab(TeX("$\\delta$")) + 
  scale_colour_manual(values = c("A" = "#516EB4", "B" = "#E76D8E", "AB" = "#7F3F98")) + ggtitle("Discrepency - Log Likelihood")


df <- matrix(0, length(files)*3, 3)
df[1:length(files),1] <- mean_test[,1]
df[1:length(files),2] <-  "A"
df[1:length(files),3] <- delta
df[(length(files) + 1):(2*length(files)),1] <- mean_test[,2]
df[(length(files) + 1):(2*length(files)),2] <- "B"
df[(length(files) + 1):(2*length(files)),3] <- delta
df[(2*length(files) + 1):(3*length(files)),1] <- mean_test[,3]
df[(2*length(files) + 1):(3*length(files)),2] <- "AB"
df[(2*length(files) + 1):(3*length(files)),3] <- delta
df <- as.data.frame(df)
colnames(df) <- c("Pval", "Condition", "delta")
df$Pval <- as.numeric(df$Pval)
df$Condition <- as.factor(df$Condition)
df$delta <- as.numeric(df$delta)
Mean_plot <- ggplot(df, aes(x=delta, y=Pval, colour = Condition)) + ylab(TeX("P-value"))+
  geom_point(aes(colour = Condition)) +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab(TeX("$\\delta$")) + 
  scale_colour_manual(values = c("A" = "#516EB4", "B" = "#E76D8E", "AB" = "#7F3F98")) + ggtitle("Discreprency - Mean Spike Count")


df <- matrix(0, length(files)*3, 3)
df[1:length(files),1] <- var_test[,1]
df[1:length(files),2] <-  "A"
df[1:length(files),3] <- delta
df[(length(files) + 1):(2*length(files)),1] <- var_test[,2]
df[(length(files) + 1):(2*length(files)),2] <- "B"
df[(length(files) + 1):(2*length(files)),3] <- delta
df[(2*length(files) + 1):(3*length(files)),1] <- var_test[,3]
df[(2*length(files) + 1):(3*length(files)),2] <- "AB"
df[(2*length(files) + 1):(3*length(files)),3] <- delta
df <- as.data.frame(df)
colnames(df) <- c("Pval", "Condition", "delta")
df$Pval <- as.numeric(df$Pval)
df$Condition <- as.factor(df$Condition)
df$delta <- as.numeric(df$delta)
Var_plot <- ggplot(df, aes(x=delta, y=Pval, colour = Condition)) + ylab(TeX("P-value"))+
  geom_point(aes(colour = Condition)) +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab(TeX("$\\delta$")) + 
  scale_colour_manual(values = c("A" = "#516EB4", "B" = "#E76D8E", "AB" = "#7F3F98")) + ggtitle("Discrepency - Var Spike Count")

grid.arrange(IIGPP_plot, Mean_plot, Var_plot, ncol = 2)

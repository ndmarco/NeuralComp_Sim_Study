### Simulation to examine WAIC performance for data generated from Competition model
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

save_dir <- "/Users/ndm34/Documents/WAIC_Sim1"

################################################################################
################################################################################
########################## Sim from Competition Model ##########################
################################################################################
################################################################################

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

run_sim2 <- function(iter){
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
  dat <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 25, 25, 25, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  
  ## Run Competition Model 
  res <- Sampler_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  ## Run IIGPP Model
  res_A <- Sampler_IIGPP(dat$X_A, dat$n_A, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_B <- Sampler_IIGPP(dat$X_B, dat$n_B, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_AB <- Sampler_IIGPP(dat$X_AB, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  WAIC_Comp <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  waic_iigpp <- WAIC_IIGPP(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res_A, res_B, res_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  # time_6_start <- Sys.time();
  # WAIC_Comp_approx_alt <- WAIC_Competition_Approx2_IS(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75))
  # time_6_end <- Sys.time();
  params <- list("I_A" = I_A, "I_B" = I_B, "sigma_A" = sigma_A, "sigma_B" = sigma_B, "delta" = delta,
                 "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B)
  output <- list("WAIC_Comp" = WAIC_Comp, "WAIC_IIGPP" = waic_iigpp, 
                 "res" = res, "res_A" = res_A, "res_B" = res_B, "res_AB" = res_AB, "data" = dat,
                 "params" = params)
  saveRDS(output, paste0(save_dir, "/Comp_generated/output", iter,".RDS"))
}

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
dir.create(paste0(save_dir, "/Comp_generated")) 
already_ran <- dir(paste0(save_dir, "/Comp_generated"))
to_run <- which(!paste0("output", 1:100, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim2(this_seed))

files <- dir(paste0(save_dir, "/Comp_generated/"))
WAIC <- matrix(0, nrow = length(files), 2)
llpd <- matrix(0, nrow = length(files), 2)
p_theta <- matrix(0, nrow = length(files), 2)
delta <- rep(0, length(files))
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
  output <- readRDS(paste0(save_dir, "/Comp_generated/", files[i]))
  WAIC[i,1] <- output$WAIC_Comp$WAIC
  WAIC[i,2] <- output$WAIC_IIGPP$WAIC
  llpd[i,1] <- output$WAIC_Comp$llpd
  llpd[i,2] <- output$WAIC_IIGPP$llpd
  p_theta[i,1] <- output$WAIC_Comp$Effective_pars
  p_theta[i,2] <- output$WAIC_IIGPP$Effective_pars
  delta[i] <- output$params$delta
  switches[i] <- number_changes(output$data$L_AB, length(output$data$n_AB))
  rel_A_B[i] <- (mean(output$data$n_AB) - mean(output$data$n_A)) / (mean(output$data$n_B) - mean(output$data$n_A))
  AB_dist_scale[i] <- abs(output$params$sigma_A - output$params$sigma_B)
  print(i)
}



time_df <- matrix(0, length(files), 3)
time_df[1:length(files),1] <- WAIC[,1] - WAIC[,2]
time_df[1:length(files),2] <- "WAIC"
time_df[1:length(files),3] <- delta
time_df <- as.data.frame(time_df)
colnames(time_df) <- c("WAIC", "Method", "Delta")
time_df$WAIC <- as.numeric(time_df$WAIC)
time_df$Method <- as.factor(time_df$Method)
time_df$Delta <- as.numeric(time_df$Delta)
p1 <- ggplot(time_df, aes(x=Delta, y=WAIC)) + scale_y_continuous(breaks = c(-2025, -1600, -1225, -900, -625, -400, -225, -100, -25, 0), trans = ssqrt_trans) +  ylab(TeX("WAIC$_{comp}$ - WAIC$_{IIGPP}$"))+
  geom_point() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Delta") + geom_hline(yintercept = 0, col = "red", linetype = "dashed", lwd = 2)

df <- matrix(0, length(files), 3)
df[1:length(files),1] <- WAIC[,1] - WAIC[,2]
df[1:length(files),2] <-  "WAIC"
df[1:length(files),3] <- switches / 25
df <- as.data.frame(df)
colnames(df) <- c("WAIC", "Method", "Switches")
df$WAIC <- as.numeric(df$WAIC)
df$Method <- as.factor(df$Method)
df$Switches <- as.numeric(df$Switches)
p2 <- ggplot(df, aes(x=Switches, y=WAIC)) + scale_y_continuous(breaks = c(-2025, -1600, -1225, -900, -625, -400, -225, -100, -25, 0), trans = ssqrt_trans) + ylab(TeX("WAIC$_{comp}$ - WAIC$_{IIGPP}$"))+
  geom_point() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Number of Switches") + geom_hline(yintercept = 0, col = "red", linetype = "dashed", lwd = 2)


df <- matrix(0, length(files)*2, 3)
df[1:length(files),1] <- p_theta[,1]
df[1:length(files),2] <- "WAIC Comp"
df[1:length(files),3] <- delta
df[(length(files) + 1):(2*length(files)),1] <- p_theta[,2]
df[(length(files) + 1):(2*length(files)),2] <- "WAIC IIGPP"
df[(length(files) + 1):(2*length(files)),3] <- delta
df <- as.data.frame(df)
colnames(df) <- c("Effective_pars", "Method", "Relative_FR")
df$Effective_pars <- as.numeric(df$Effective_pars)
df$Method <- as.factor(df$Method)
df$Relative_FR <- as.numeric(df$Relative_FR)
p5 <- ggplot(df, aes(x=Relative_FR, y=Effective_pars, colour = Method))  + ylab(TeX("Effective Parameters"))+
  geom_point(aes(colour = Method)) +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Delta")

df <- matrix(0, length(files), 3)
df[1:length(files),1] <- llpd[,1] - llpd[,2]
df[1:length(files),2] <- "WAIC Comp Conditional"
df[1:length(files),3] <- delta
df <- as.data.frame(df)
colnames(df) <- c("llpd", "Method", "Relative_FR")
df$llpd <- as.numeric(df$llpd)
df$Method <- as.factor(df$Method)
df$Relative_FR <- as.numeric(df$Relative_FR)
p6 <- ggplot(df, aes(x=Relative_FR, y=llpd)) + ylab(TeX("lppd$_{comp}$ - lppd$_{IIGPP}$"))+
  geom_point() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Delta") 

plts1 <- ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="none")
plts3 <- ggarrange(p5, p6, ncol=2, common.legend = TRUE, legend="bottom")
plot <- ggarrange(plts1,plts3, nrow = 2)
annotate_figure(plot, top = text_grob("WAIC (Competition Generated)", 
                                      color = "Black", size = 20))






################################################################################
################################################################################
############################# Sim from IIGPP Model #############################
################################################################################
################################################################################


### Function to generate data from IIGPP model
generate_data_TI_IIGPP <- function(I_A, I_B, I_AB, basis_coef_A, basis_coef_B, basis_coef_AB,
                                   sigma_A, sigma_B, sigma_AB, N_A, N_B, N_AB, seed, time, basis_degree,
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
    spike <- rinvgauss(1, mean = 1 / (I_AB * exp(spline %*% basis_coef_AB)), shape = (1 / sigma_AB)^2)
    total_time <- spike
    while(total_time < time){
      spline <- GetBSpline(total_time, basis_degree, boundary_knots, internal_knots)
      spike_i <- rinvgauss(1, mean = 1 / (I_AB * exp(spline %*% basis_coef_AB)), shape = (1 / sigma_AB)^2)
      total_time <- total_time + spike_i
      spike <- c(spike, spike_i)
    }
    total_time <- total_time - spike[length(spike)]
    spike <- spike[-length(spike)]
    X_AB[[i]] <- spike
    n_AB[i] <- length(spike)
  }
  
  spline <- GetBSpline(seq(0, 1.0, 0.01), basis_degree, boundary_knots, internal_knots)
  return(list("X_A" = X_A, "X_B" = X_B, "X_AB" = X_AB, "n_A" = n_A,
              "n_B" = n_B, "n_AB" = n_AB, "A_CI" = (I_A * exp(spline %*% basis_coef_A)),
              "B_CI" = (I_B * exp(spline %*% basis_coef_B)),
              "AB_CI" = (I_B * exp(spline %*% basis_coef_AB))))
}


run_sim3 <- function(iter){
  set.seed(iter)
  
  I_A <- rtruncnorm(1, a = 0, mean = 40, sd = 4)
  I_B <- rtruncnorm(1, a = 0, mean = 80, sd = 4)
  I_AB <- rtruncnorm(1, a = 0, mean = 60, sd = 8)
  sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
  sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
  sigma_AB <- rtruncnorm(1, a = 0, mean = sqrt(60), sd = 4)
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  basis_coef_AB <- rnorm(6, 0, 0.3)
  dat <- generate_data_TI_IIGPP(I_A, I_B, I_AB, basis_coef_A, basis_coef_B, basis_coef_AB, sigma_A, sigma_B, sigma_AB, 25, 25, 25, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  
  ## Run Competition Model 
  res <- Sampler_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  ## Run IIGPP Model
  res_A <- Sampler_IIGPP(dat$X_A, dat$n_A, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_B <- Sampler_IIGPP(dat$X_B, dat$n_B, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  res_AB <- Sampler_IIGPP(dat$X_AB, dat$n_AB, 5000, 3, c(0,1), c(0.25, 0.5, 0.75), 1)
  
  WAIC_Comp <- WAIC_Competition(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  waic_iigpp <- WAIC_IIGPP(dat$X_A, dat$X_B, dat$X_AB, dat$n_A, dat$n_B, dat$n_AB, res_A, res_B, res_AB, 3, c(0,1), c(0.25, 0.5, 0.75), 1, burnin_prop = 0.5)
  params <- list("I_A" = I_A, "I_B" = I_B, "I_AB" = I_AB, "sigma_A" = sigma_A, "sigma_B" = sigma_B, "sigma_AB" = sigma_AB,
                 "basis_coef_A" = basis_coef_A, "basis_coef_B" = basis_coef_B, "basis_coef_AB" = basis_coef_AB)
  output <- list("WAIC_Comp" = WAIC_Comp,  "WAIC_IIGPP" = waic_iigpp, 
                 "res" = res, "res_A" = res_A, "res_B" = res_B, "res_AB" = res_AB, "data" = dat,
                 "params" = params)
  saveRDS(output, paste0(save_dir, "/IIGPP_generated/output", iter,".RDS"))
}

library(future.apply)

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
dir.create(paste0(save_dir, "/IIGPP_generated")) 
already_ran <- dir(paste0(save_dir, "/IIGPP_generated"))
to_run <- which(!paste0("output", 1:100, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim3(this_seed))




### Post-Processing
files <- dir(paste0(save_dir, "/IIGPP_generated/"))
WAIC <- matrix(0, nrow = length(files), 2)
llpd <- matrix(0, nrow = length(files), 2)
p_theta <- matrix(0, nrow = length(files), 2)
rel_A_B <- rep(0, length(files))
for(i in 1:length(files)){
  output <- readRDS(paste0(save_dir, "/IIGPP_generated/", files[i]))
  WAIC[i,1] <- output$WAIC_Comp$WAIC
  WAIC[i,2] <- output$WAIC_IIGPP$WAIC
  llpd[i,1] <- output$WAIC_Comp$llpd
  llpd[i,2] <- output$WAIC_IIGPP$llpd
  p_theta[i,1] <- output$WAIC_Comp$Effective_pars
  p_theta[i,2] <- output$WAIC_IIGPP$Effective_pars
  rel_A_B[i] <- (mean(output$data$n_AB) - mean(output$data$n_A)) / (mean(output$data$n_B) - mean(output$data$n_A))
  print(i)
}


df <- matrix(0, length(files), 3)
df[1:length(files),1] <- WAIC[,1] - WAIC[,2]
df[1:length(files),2] <- "WAIC Comp Conditional"
df[1:length(files),3] <- rel_A_B
df <- as.data.frame(df)
colnames(df) <- c("WAIC", "Method", "Relative_FR")
df$WAIC <- as.numeric(df$WAIC)
df$Method <- as.factor(df$Method)
df$Relative_FR <- as.numeric(df$Relative_FR)
p1 <- ggplot(df, aes(x=Relative_FR, y=WAIC)) + scale_y_continuous(breaks = c(-200, 0, 50, 200, 900, 2025, 3600, 5625, 8100), trans = ssqrt_trans) + ylab(TeX("WAIC$_{comp}$ - WAIC$_{IIGPP}$"))+
  geom_point() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Relative AB Firing Rate") + 
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("A", "B", " B + (B-A)")) + geom_hline(yintercept = 0, col = "red", linetype = "dashed", lwd = 2)


df <- matrix(0, length(files)*2, 3)
df[1:length(files),1] <- p_theta[,1]
df[1:length(files),2] <- "WAIC Comp"
df[1:length(files),3] <- rel_A_B
df[(length(files) + 1):(2*length(files)),1] <- p_theta[,2]
df[(length(files) + 1):(2*length(files)),2] <- "WAIC IIGPP"
df[(length(files) + 1):(2*length(files)),3] <- rel_A_B
df <- as.data.frame(df)
colnames(df) <- c("Effective_pars", "Method", "Relative_FR")
df$Effective_pars <- as.numeric(df$Effective_pars)
df$Method <- as.factor(df$Method)
df$Relative_FR <- as.numeric(df$Relative_FR)
p2 <- ggplot(df, aes(x=Relative_FR, y=Effective_pars, colour = Method)) + scale_y_continuous(breaks = c(0, 4, 16, 64, 144, 1024, 3136), trans = ssqrt_trans) + ylab(TeX("Effective Parameters"))+
  geom_point(aes(colour = Method)) +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Relative AB Firing Rate") + 
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("A", "B", " B + (B-A)"))

df <- matrix(0, length(files), 3)
df[1:length(files),1] <- llpd[,1] - llpd[,2]
df[1:length(files),2] <- "WAIC Comp Conditional"
df[1:length(files),3] <- rel_A_B
df <- as.data.frame(df)
colnames(df) <- c("llpd", "Method", "Relative_FR")
df$llpd <- as.numeric(df$llpd)
df$Method <- as.factor(df$Method)
df$Relative_FR <- as.numeric(df$Relative_FR)
p3 <- ggplot(df, aes(x=Relative_FR, y=llpd)) + scale_y_continuous(breaks = c(0, -16, -64, -144, -1024, -3136), trans = ssqrt_trans) + ylab(TeX("lppd$_{comp}$ - lppd$_{IIGPP}$"))+
  geom_point() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Relative AB Firing Rate") + 
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("A", "B", " B + (B-A)")) + geom_hline(yintercept = 0, col = "red", linetype = "dashed", lwd = 2)

ph1 <- ggarrange(p1, common.legend = TRUE, legend="none")
ph2 <- ggarrange(p2, p3, ncol=2, common.legend = TRUE, legend="bottom")
plot <- ggarrange(ph1,ph2, nrow = 2)
annotate_figure(plot, top = text_grob("WAIC (IIGPP Generated)", 
                                      color = "Black", size = 20))






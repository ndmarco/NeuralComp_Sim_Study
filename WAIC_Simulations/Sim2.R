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
  I_AB <- rtruncnorm(1, a = 0, mean = 100, sd = 8)
  sigma_A <- rtruncnorm(1, a = 0, mean = sqrt(40), sd = 2)
  sigma_B <- rtruncnorm(1, a = 0, mean = sqrt(80), sd = 2)
  sigma_AB <- rtruncnorm(1, a = 0, mean = sqrt(100), sd = 4)
  basis_coef_A <- rnorm(6, 0, 0.3)
  basis_coef_B <- rnorm(6, 0, 0.3)
  basis_coef_AB <- rnorm(6, 0, 0.3)
  delta <- rlnorm(1, -2.5, 0.5)
  prop <- runif(1, min = 0.04, max = 0.96)
  n_AB_comp <- round(25 * prop)
  n_AB_IIGPP <- 25 - n_AB_comp
  dat <- generate_data_TI(I_A, I_B, basis_coef_A, basis_coef_B, sigma_A, sigma_B, delta, 25, 25, n_AB_comp, iter, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  dat1 <- generate_data_TI_IIGPP(I_A, I_B, I_AB, basis_coef_A, basis_coef_B, basis_coef_AB, sigma_A, sigma_B, sigma_AB, iter, 1, n_AB_IIGPP, 1, 1, 3, c(0,1), c(0.25, 0.5, 0.75))
  dat$n_AB <- c(dat$n_AB, dat1$n_AB)
  dat$X_AB <- c(dat$X_AB, dat1$X_AB)
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
  saveRDS(output, paste0(save_dir, "/output", iter,".RDS"))
}

library(future.apply)

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)
already_ran <- dir(save_dir)
to_run <- which(!paste0("output", 1:100, ".RDS") %in% already_ran)
future_lapply(to_run, function(this_seed) run_sim3(this_seed))


### Post Processing
files <- dir(save_dir)
WAIC <- matrix(0, nrow = length(files), 2)
prop <- rep(0, length(files))

for(i in 1:length(files)){
  output <- readRDS(paste0(save_dir, "/", files[i]))
  prop[i] <- length(output$data$L_AB) / 25
  WAIC[i,1] <- output$WAIC_Comp$WAIC
  WAIC[i,2] <- output$WAIC_IIGPP$WAIC
  print(i)
}

df <- matrix(0, length(files), 3)
df[1:length(files),1] <- WAIC[,1] - WAIC[,2]
df[1:length(files),2] <-  "WAIC Comp"
df[1:length(files),3] <- prop
df <- as.data.frame(df)
colnames(df) <- c("WAIC", "Method", "Prop_comp")
df$WAIC <- as.numeric(df$WAIC)
df$Method <- as.factor(df$Method)
df$Prop_comp <- as.numeric(df$Prop_comp)
ggplot(df, aes(x=Prop_comp, y=WAIC)) + scale_y_continuous(breaks = c(-900, -625, -400, -225, -100, -25, 0, 25, 100, 225, 400, 625, 900), trans = ssqrt_trans) + ylab(TeX("WAIC$_{comp}$ - WAIC$_{IIGPP}$"))+
  geom_point() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                                         plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + xlab("Proportion Generated from Competition Framework") + 
  geom_hline(yintercept = 0,  color="red", linetype="dashed")


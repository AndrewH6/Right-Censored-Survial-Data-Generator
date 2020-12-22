RC_data_generator_01 <- function(sample_size = 500,
                              outlier_ind = FALSE, c_case = 0){
  
  # 01 Setup Data
  ######
  X <- runif(sample_size, 0, 1)
  err <- rnorm(sample_size, mean = 0, sd = 1)
  logtime <- 6 + 6 * X + err
  time <- round(exp(logtime))
  o_time <- time
  data_g <- data.frame(X,o_time,time)
  
  # 02 Setup Censoring
  ######
  
  if(c_case == 0){
    n_a <- 90000000
    n_b <- 90000001
  } else if (c_case == 1) {
    n_a <- 10
    n_b <- 500000
  } else if (c_case == 2) {
    n_a <- 10
    n_b <- 80000
  } else if (c_case == 3) {
    n_a <- 10
    n_b <- 30000
  }
  
  X <- data_g$X
  time <- data_g$time
  nr <- dim(data_g)[1]
  C_time <- runif(nr,n_a,n_b)
  #cr[i] <- sum(data_g$time > C_time)/nr # Censoring rate
  status <- rep(1,nr)
  status[data_g$time > C_time] <- 0
  data_g$time[data_g$time > C_time] <- C_time[data_g$time > C_time]
  data_i <- data.frame(data_g,C_time,status)
  ######
  cen_list <- which(data_i$status==1)
  
  
  # 03 Setup Outliers
  ######
  
  if (outlier_ind == TRUE){
    otl_1 <- c(0.000001,round(exp(15)),round(exp(15)))
    otl_2 <- c(0.00006,round(exp(12)),round(exp(12)))
    otl_3 <- c(0.5,round(exp(16)),round(exp(16)))
    otl_4 <- c(0.97,20,20)
    otl_5 <- c(0.99,10,10)
    otl_6 <- c(0.98,16,16)
    otl_7 <- c(0.000009,round(exp(12)),round(exp(12)))
    otl_8 <- c(0.92,15,15)
    otl_9 <- c(0.99,26,26)
    data_i[cen_list[1],1:3] <- otl_1
    data_i[cen_list[2],1:3] <- otl_2
    data_i[cen_list[3],1:3] <- otl_3
    data_i[cen_list[4],1:3] <- otl_4
    data_i[cen_list[5],1:3] <- otl_5
    data_i[cen_list[6],1:3] <- otl_6
    data_i[cen_list[7],1:3] <- otl_7
    data_i[cen_list[8],1:3] <- otl_8
    data_i[cen_list[9],1:3] <- otl_9
  }
  
  ######
  fit_bj <- NA
  try(fit_bj <- bj(Surv(log(time), status) ~ X,
                   data = data_i,
                   x = TRUE, y = TRUE, link = "identity"), silent = TRUE)
  b0 <- NA
  b0 <- fit_bj$coefficients[1]
  b1 <- NA
  b1 <- fit_bj$coefficients[2]
  
  pbj_1 <- b0
  pbj_2 <- b1
  
  p_bj <- c(pbj_1,pbj_2)
  names(p_bj) <- NULL
  
  
  bjltime <- NA
  try(bjltime <- as.numeric(fit_bj$y.imputed), silent = TRUE)
  try(data_i$BJTime <- bjltime, silent = TRUE)
  ######
  
  fit_mle <- NA
  try(fit_mle <- survreg(Surv(time, status) ~ X, 
                         data = data_i, dist = "lognormal"), silent = TRUE)
  pmle_1 <- NA
  pmle_2 <- NA
  try(pmle_1 <- fit_mle$coefficients[1], silent = TRUE)
  try(pmle_2 <- fit_mle$coefficients[2], silent = TRUE)
  pmle <- c(pmle_1,pmle_2)
  
  
  Censoring_rate <- sum(data_i$status==0)/dim(data_i)[1]
  N_uc <- round((1-Censoring_rate)*dim(data_i)[1])
  result_list <- list(Cencoring_rate = Censoring_rate,
                      Number_unc = N_uc,
                      B_J = p_bj,
                      MLE = pmle,
                      Data = data_i)
  return(result_list)
}

data_6 <- RC_data_generator_01(c_case = 3,outlier_ind = TRUE)
data_6$Cencoring_rate
names(data_6$Data)

bjerr_try <- data_6$Data$BJTime - 6 - 6 * data_6$Data$X

my_RSE(bjerr_try)

b <- seq(-3,3,by=0.01)
my_epi_c <-my_epi_RSE_1(b,samp=bjerr_try,my_bw=0.9)
plot(x=b,y=my_epi_c,ylim=c(0,0.8))


RC_data_generator_02 <- function(sample_size = 500,
                                 outlier_ind = FALSE, c_case = 0){
  
  # 01 Setup Data
  ######
  X <- runif(sample_size, 0, 1)
  err <- rnorm(sample_size, mean = 0, sd = 1)
  logtime <- 6 + 6 * X + err
  time <- round(exp(logtime))
  o_time <- time
  data_g <- data.frame(X,o_time,time)
  
  # 02 Setup Censoring
  ######
  
  if(c_case == 0){
    n_a <- 90000000
    n_b <- 90000001
  } else if (c_case == 1) {
    n_a <- 10
    n_b <- 500000
  } else if (c_case == 2) {
    n_a <- 10
    n_b <- 80000
  } else if (c_case == 3) {
    n_a <- 10
    n_b <- 30000
  }
  
  X <- data_g$X
  time <- data_g$time
  nr <- dim(data_g)[1]
  C_time <- runif(nr,n_a,n_b)
  #cr[i] <- sum(data_g$time > C_time)/nr # Censoring rate
  status <- rep(1,nr)
  status[data_g$time > C_time] <- 0
  data_g$time[data_g$time > C_time] <- C_time[data_g$time > C_time]
  data_i <- data.frame(data_g,C_time,status)
  ######
  cen_list <- which(data_i$status==1)
  
  
  # 03 Setup Outliers
  ######
  
  if (outlier_ind == TRUE){
    otl_1 <- c(0.000001,round(exp(15)),round(exp(15)))
    otl_2 <- c(0.00006,round(exp(12)),round(exp(12)))
    otl_3 <- c(0.5,round(exp(16)),round(exp(16)))
    otl_4 <- c(0.97,20,20)
    otl_5 <- c(0.000001,round(exp(15)),round(exp(15)))
    otl_6 <- c(0.98,16,16)
    otl_7 <- c(0.000009,round(exp(12)),round(exp(12)))
    otl_8 <- c(0.92,15,15)
    otl_9 <- c(0.99,26,26)
    data_i[cen_list[1],1:3] <- otl_1
    data_i[cen_list[2],1:3] <- otl_2
    data_i[cen_list[3],1:3] <- otl_3
    data_i[cen_list[4],1:3] <- otl_4
    data_i[cen_list[5],1:3] <- otl_5
    data_i[cen_list[6],1:3] <- otl_6
    data_i[cen_list[7],1:3] <- otl_7
    data_i[cen_list[8],1:3] <- otl_8
    data_i[cen_list[9],1:3] <- otl_9
  }
  
  ######
  fit_bj <- NA
  try(fit_bj <- bj(Surv(log(time), status) ~ X,
                   data = data_i,
                   x = TRUE, y = TRUE, link = "identity"), silent = TRUE)
  b0 <- NA
  b0 <- fit_bj$coefficients[1]
  b1 <- NA
  b1 <- fit_bj$coefficients[2]
  
  pbj_1 <- b0
  pbj_2 <- b1
  
  p_bj <- c(pbj_1,pbj_2)
  names(p_bj) <- NULL
  
  
  bjltime <- NA
  try(bjltime <- as.numeric(fit_bj$y.imputed), silent = TRUE)
  try(data_i$BJTime <- bjltime, silent = TRUE)
  ######
  
  fit_mle <- NA
  try(fit_mle <- survreg(Surv(time, status) ~ X, 
                         data = data_i, dist = "lognormal"), silent = TRUE)
  pmle_1 <- NA
  pmle_2 <- NA
  try(pmle_1 <- fit_mle$coefficients[1], silent = TRUE)
  try(pmle_2 <- fit_mle$coefficients[2], silent = TRUE)
  pmle <- c(pmle_1,pmle_2)
  
  
  Censoring_rate <- sum(data_i$status==0)/dim(data_i)[1]
  N_uc <- round((1-Censoring_rate)*dim(data_i)[1])
  result_list <- list(Cencoring_rate = Censoring_rate,
                      Number_unc = N_uc,
                      B_J = p_bj,
                      MLE = pmle,
                      Data = data_i)
  return(result_list)
}
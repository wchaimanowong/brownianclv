################################################################################
# Run demo using CDNOW_sample

run_demo <- function() {
  Data <- read.table(paste(getwd(),'/CDNOW_sample.txt',sep=''), header=FALSE, sep="")
  names(Data) <- c("master_id", "id", "day", "amount", "dollar_amount")
  Data$day <- as.Date(as.character(Data$day), "%Y%m%d")
  Data$day <- as.numeric(Data$day - min(Data$day))/7

  # Get the number of consumers.
  num_ids <- length(unique(Data[,2]))

  # 39 weeks calibration time period
  train_period <- 39
  # length of test period
  test_period <- 39

  # TTT vector stores the amount of time for each consumer from joining the firm to the end of 39 weeks calibration period.
  TTT <- c()

  # In the following we prepare the data by separating them into training (calibration period) and testing set.
  training_Data = list()
  test_period_txns <- c()
  for(i in 1:num_ids){
    id_Data <- Data[Data$id == i,]

    start_date <- id_Data$day[1]
    id_Data$day <- id_Data$day - start_date
    TTT[i] <- time_period - start_date
    training_Data[[i]] <- c(length(id_Data[id_Data$day < TTT[i],]$day), id_Data[id_Data$day < TTT[i],]$day, TTT[i])

    test_period_txns[i] <- length(unique(id_Data[id_Data$day >= TTT[i],]$day))
  }

  # Estimate the model's parameters.
  C <- Estimate_Parameters(training_Data)

  # Predict test period transactions based on the observed transactions in training period.
  res <- Make_Prediction(C, observed_Data=training_Data, testing_period=test_period)
  prediction_B <- res[[1]]
  EXt_B_het_vec <- res[[2]]

  for (i in 1:(length(prediction_B)-1)) {
    print(paste("The expected transactions in test period given", i-1, "observed transactions is", prediction_B[i], collapse=" "))
  }
  print(paste("The expected transactions in test period given", length(prediction_B)-1, "or more observed transactions is",
              prediction_B[length(prediction_B)], collapse=" "))

  # Compute MAE:
  MAE <- sum(abs(test_period_txns - EXt_B_het_vec))/num_ids
}

################################################################################

# Generate NN Brownian motions dv_t = sigma*dW_t starting from zero over the time interval [0,TT] with step size dt.
Gen_Brownian_Motions <- function(NN, sigma, TT, dt) {
  num_steps_ <- as.integer(round(TT/dt))
  v_list <- array(0, dim=c(NN, num_steps_))
  min_v_list <- array(0, dim=c(NN, num_steps_))
  for (j in 1:NN) {
    for (i in 2:num_steps_) {
      v_list[j,i] <- rnorm(1, v_list[j,i-1], sigma*sqrt(dt))
      min_v_list[j,i] <- min(v_list[j,1:i])
    }
  }
  v_df <- as.data.frame(v_list)
  min_v_df <- as.data.frame(min_v_list)
  # We are ignoring the first time frame (t = 0) since all Brownian motions starts at v = 0.
  v_df <- v_df[2:ncol(v_df)]
  min_v_df <- min_v_df[2:ncol(min_v_df)]

  return_df <- list(v_df, min_v_df)
  return(return_df)
}

# # Generate NN Brownian bridges starting from v0 = 0, ending at a fixed vT.
# Gen_Conditional_Brownian_Bridges <- function(vT, v_min, NN, sigma, TT, dt) {
#   v_list <-
# }

Prepare_Simulation <- function(training_Data, NN, sigma, rr, dt, Replicate, test_period) {
  if(missing(test_period)) {
    test_period <- 0
  }
  train_period <- 0

  for (i in 1:length(training_Data)){
    if (tail(training_Data[[i]],1) > train_period) {
      train_period <- tail(training_Data[[i]],1)
    }
  }
  tol_period = train_period + test_period
  # number of steps of Brownian motion on the entire calibration period.
  num_steps <- as.integer(round(tol_period/dt))-1

  # Dataframe to store the simulated Brownian motions for all consumers.
  v_df <- data.frame(matrix(ncol = num_steps, nrow = 0))
  # Dataframe to store the minimum of the simulated Brownian motions for all consumers.
  min_v_df <- data.frame(matrix(ncol = num_steps, nrow = 0))

  # pre-computed vectors to be assembled as lambda_vec later.
  pre_l_vec_1 <- data.frame(matrix(ncol = num_steps, nrow = 0))
  pre_l_vec_2 <- data.frame(matrix(ncol = num_steps, nrow = 0))

  # Data frame and vector which will help determine which simulated path is valid. e.g. Paths that reaches v_t = uv before t = tx is invalid (the consumer couldn't have made a transaction at tx if he is dead).
  select_df <- data.frame(matrix(ncol=num_steps, nrow = 0))
  select_vec <- c()

  # Dataframe to select the testing period.
  select_test <- data.frame(matrix(ncol=num_steps, nrow=0))

  # To construct each dataframe above by concatenating smaller dataframes,is faster to put all the smaller dataframes in a list first before before applying bind_rows.
  v_df_ <- list()
  min_v_df_ <- list()

  pre_l_vec_1_ <- list()
  pre_l_vec_2_ <- list()

  select_df_ <- list()
  select_vec_ <- c()

  select_test_ <- list()

  if (Replicate == TRUE) {
    returned_df <- Gen_Brownian_Motions(NN, sigma, tol_period, dt)
    master_v_df <- returned_df[[1]]
    master_min_v_df <- returned_df[[2]]
  }

  print("Generating random paths...")

  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(training_Data), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "█")   # Character used to create the bar

  for (i in 1:length(training_Data)) {
    TT <- as.integer(round(tail(training_Data[[i]], 1)/dt))-1
    tt_ <- training_Data[[i]][2:(length(training_Data[[i]])-1)]
    tt_ <- tt_[tt_ > 0]
    x <- length(tt_)
    tt <- c(0)
    if (x > 0) {
      # The following for-loop is just to make sure that all the purchases occur before the end of calibration period
      # (ideally, it should do nothing, but just in case...)
      for (j in 1:x) {
        tt[j] <- min(TT, max(1, as.integer(round((1/dt)*tt_[j]))))
      }
    }

    pre_l_vec_1i <- c(rep(1, TT), rep(0, num_steps - TT))
    pre_l_vec_2i <- rep(0, num_steps)
    tx <- tt[length(tt)]
    if (tx > 0) {
      for (ti in tt) {
        pre_l_vec_1i[ti] <- 0
        pre_l_vec_2i[ti] <- 1
      }
    }
    pre_l_vec_1i <- rbind(data.frame(), pre_l_vec_1i)
    pre_l_vec_2i <- rbind(data.frame(), pre_l_vec_2i)
    pre_l_vec_1i <- pre_l_vec_1i[rep(seq_len(nrow(pre_l_vec_1i)), each=NN),]
    pre_l_vec_2i <- pre_l_vec_2i[rep(seq_len(nrow(pre_l_vec_2i)), each=NN),]

    if (Replicate == TRUE) {
      vi_df <- master_v_df
      min_vi_df <- master_min_v_df
    } else {
      returned_df <- Gen_Brownian_Motions(NN, sigma, time_period, dt)
      vi_df <- returned_df[[1]]
      min_vi_df <- returned_df[[2]]
    }

    names(vi_df) <- names(v_df)
    names(min_vi_df) <- names(min_v_df)
    v_df_[[i]] <- vi_df
    min_v_df_[[i]] <- min_vi_df

    selecti_df <- .Machine$double.xmax*pre_l_vec_1i + pre_l_vec_2i*vi_df
    names(selecti_df) <- names(select_df)
    select_df_[[i]] <- selecti_df

    if (tx > 0) {
      select_vec_[[i]] <- min_vi_df[,tx]
    } else {
      select_vec_[[i]] <- rep(.Machine$double.xmax, NN)
    }

    if (test_period > 0) {
      TT_test <- as.integer(round(test_period/dt))
      select_testi <- c(rep(0,TT), rep(1,TT_test), rep(0,num_steps - TT - TT_test))
      select_testi <- rbind(data.frame(), select_testi)
      select_testi <- select_testi[rep(seq_len(nrow(select_testi)), each=NN),]
      names(select_testi) <- names(select_test)
      select_test_[[i]] <- select_testi
    }

    names(pre_l_vec_1i) <- names(pre_l_vec_1)
    names(pre_l_vec_2i) <- names(pre_l_vec_2)
    pre_l_vec_1_[[i]] <- pre_l_vec_1i
    pre_l_vec_2_[[i]] <- pre_l_vec_2i

    setTxtProgressBar(pb, i)
  }
  print("Finish")

  v_df <- dplyr::bind_rows(v_df_)
  rm(v_df_)
  min_v_df <- dplyr::bind_rows(min_v_df_)
  rm(min_v_df_)

  pre_l_vec_1 <- dplyr::bind_rows(pre_l_vec_1_)
  rm(pre_l_vec_1_)
  pre_l_vec_2 <- dplyr::bind_rows(pre_l_vec_2_)
  rm(pre_l_vec_2_)

  select_df <- dplyr::bind_rows(select_df_)
  rm(select_df_)
  select_vec <- Reduce(c, select_vec_)
  rm(select_vec_)

  if(test_period > 0) {
    select_test <- dplyr::bind_rows(select_test_)
    rm(select_test_)

    return(list(v_df, min_v_df, pre_l_vec_1, pre_l_vec_2, select_df, select_vec, select_test))
  } else {
    return(list(v_df, min_v_df, pre_l_vec_1, pre_l_vec_2, select_df, select_vec))
  }
}

# Functions for the homogeneous model.
# Output the big vector of likelihoods of every consumer's purchase history from the given parameters.
L_B <- function(v, uv, lambda, dt, dfs, num_ids, NN) {
  v_df <- dfs[[1]]
  min_v_df <- dfs[[2]]
  pre_l_vec_1 <- dfs[[3]]
  pre_l_vec_2 <- dfs[[4]]
  select_df <- dfs[[5]]
  select_vec <- dfs[[6]]

  if (uv > v) {
    return(rep(0, num_ids))
  }
  lambda_ <- lambda*dt
  lambda_df <- (-lambda_)*pre_l_vec_1 + (lambda_ - 1)*pre_l_vec_2

  # select_rows only select the paths that is positive at the time of purchase and alive at least until the last recorded purchase.
  select_rows <- apply(1*(select_vec > uv - v)*(select_df > -v), MARGIN=1, prod)
  mask <- 1*(min_v_df > uv - v)

  v_vec <- 1*(v_df > -v)
  v_vec <- lambda_df*mask*v_vec
  v_vec <- v_vec + 1
  v_vec <- select_rows*v_vec
  v_vec <- apply(v_vec, MARGIN=1, prod)/NN
  v_vec <- apply(matrix(data=v_vec, nrow=NN), MARGIN=2, sum)

  return(v_vec)
}

# Return probability of 0,1,2,3,4,5,6, and 7+ repeated purchases during the time interval [0,t].
PXt_B <- function(v, uv, lambda, dt, dfs, num_ids, NN) {
  v_df <- dfs[[1]]
  min_v_df <- dfs[[2]]
  pre_l_vec_1 <- dfs[[3]]
  pre_l_vec_2 <- dfs[[4]]
  select_df <- dfs[[5]]
  select_vec <- dfs[[6]]
  select_test <- dfs[[7]]

  lambda_ <- lambda*dt

  # Everything higher than 7 will be grouped as '7+'.
  PXt_B_i <- rep(0,8)

  select_train <- pre_l_vec_1 + pre_l_vec_2
  mask <- 1*(min_v_df > uv - v)
  v_vec <- 1*(v_df > -v)
  v_vec <- select_train*mask*v_vec
  v_vec <- apply(v_vec, MARGIN=1, sum)

  PXt_B_i[8] <- num_ids
  for (i in 0:6) {
    sum_binom <- sum(dbinom(i, v_vec, lambda_))/NN
    PXt_B_i[i+1] <- sum_binom
    PXt_B_i[8] <- PXt_B_i[8] - sum_binom
  }

  return(PXt_B_i)
}

# Expected number of purchases over the period [T, T+t] conditioned on the purchase history over [0,T] of every consumers.
EXt_B <- function(v, uv, lambda, dt, dfs, num_ids, NN) {
  v_df <- dfs[[1]]
  min_v_df <- dfs[[2]]
  pre_l_vec_1 <- dfs[[3]]
  pre_l_vec_2 <- dfs[[4]]
  select_df <- dfs[[5]]
  select_vec <- dfs[[6]]
  select_test <- dfs[[7]]

  if (uv > v) {
    return(rep(0, num_ids))
  }
  lambda_ <- lambda*dt
  lambda_df <- (-lambda_)*pre_l_vec_1 + (lambda_ - 1)*pre_l_vec_2

  select_rows <- apply(1*(select_vec > uv - v)*(select_df > -v), MARGIN=1, prod)
  mask <- 1*(min_v_df > uv - v)

  v_vec <- 1*(v_df > -v)
  v_vec <- mask*v_vec

  prob_vec <- lambda_df*v_vec
  prob_vec <- prob_vec + 1
  prob_vec <- apply(prob_vec, MARGIN=1, prod)
  prob_vec <- select_rows*prob_vec

  # Each path will contribute an expected lambda_ transaction each time it is above zero over the interval [0,T+t].
  v_vec <- select_test*lambda_*v_vec
  v_vec <- apply(v_vec, MARGIN=1, sum) # Summing across the columns

  E_rep_txn <- apply(matrix(data=prob_vec*v_vec, nrow=NN), MARGIN=2, sum)/apply(matrix(data=prob_vec, nrow=NN), MARGIN=2, sum)
  return(E_rep_txn)
}

################################################################################
# Some utilities functions for PDE method.

Np_Mask <- function(uv, dv, num_v_steps) {
  if (uv > 0) {
    return(rep(1,num_v_steps))
  }

  a <- -uv/dv+1
  a_int <- as.integer(floor(a))
  a_frac <- a%%1
  np_mask <- rep(0,max(min(a_int,num_v_steps),0))
  if (length(np_mask) < num_v_steps) {
    np_mask <- c(np_mask, 1-a_frac)
  }
  if (length(np_mask) < num_v_steps) {
    np_mask <- c(np_mask, rep(1,max(min(num_v_steps-length(np_mask),num_v_steps),0)) )
  }
  return(np_mask)
}

################################################################################
# Likelihood function via solving PDE with finite difference method

Initialize_FDM <- function(v, uv, lambda, sigma, dv, num_v_steps) {

  p0 <- rep(0, num_v_steps)
  for (i in 1:num_v_steps) {
    p0[i] <- P_v(uv + dv*(i-1), uv, lambda, v, sigma)
  }

  return(p0)
}

# Evolve the probability distribution of alive and non-purchasing consumers.
Evolve_FDM <- function(p0, num_t_steps, uv, lambda, sigma, dt, dv, num_v_steps) {
  if (num_t_steps <= 0) {
    return(list(p0,0,0))
  }

  coeff <- (1/2)*((sigma^2)*dt/dv^2)

  # collecting dead people on the lower and upper ends.
  dead_L <- 0
  dead_H <- 0

  num_v_steps <- length(p0)
  p1 <- rep(0,num_v_steps)
  for (n in 1:num_t_steps) {
    for (i in 2:(num_v_steps-1)) {
      above <- 0
      if (i*dv > -uv+1) {
        above <- 1
      }
      p1[i] <- coeff*p0[i+1] + (1 - 2*coeff - lambda*dt*above)*p0[i] + coeff*p0[i-1]
    }
    dead_L <- dead_L + coeff*p0[2]
    dead_H <- dead_H + coeff*p0[num_v_steps-1]

    p0 <- p1
  }
  return_list <- list(p1, dead_L, dead_H)
  return(return_list)
}

L_PDE_FDM <- function(X, v, uv, lambda, sigma, dt, dv, num_v_steps) {
  TT <- as.integer(tail(X, 1)/dt)
  tt <- unique(as.integer(X[2:(length(X)-1)]/dt))
  tt[1] <- -1 # For technical convenience (we really want to start from the end of the zeroth day).

  p <- Initialize_FDM(v, uv, lambda, sigma, dv, num_v_steps)

  if (v < uv) {
    return(list(0,0*p))
  }

  np_mask <- Np_Mask(uv, dv, num_v_steps)

  dead_L <- 0
  dead_H <- 0
  for (i in 1:length(tt)) {
    if (i < length(tt)) {
      num_t_steps <- tt[i+1] - tt[i]
      p <- Evolve_FDM(p, num_t_steps, uv, lambda, sigma, dt, dv, num_v_steps)[[1]]
      p <- lambda*dt*np_mask*p
    } else {
      num_t_steps <- TT - tt[i]
      returned_list <- Evolve_FDM(p, num_t_steps, uv, lambda, sigma, dt, dv, num_v_steps)
      p <- returned_list[[1]]
      dead_L <- returned_list[[2]]
      dead_H <- returned_list[[3]]
    }
  }
  val <- (sum(p) + dead_L + dead_H)*dv

  return_list <- list(val, p)
  return(return_list)
}

L_PDE_FDM_vec <- function(training_Data, v, uv, lambda, sigma, dt, dv, num_v_steps) {
  val_vec <- c()
  for (i in 1:length(training_Data)) {
    val_vec[i] <- L_PDE_FDM(training_Data[[i]], v, uv, lambda, sigma, dt, dv, num_v_steps)[[1]]
  }
  return(val_vec)
}

################################################################################
# Likelihood function via solving PDE with finite elements method.

Initialize_FEM <- function(v, uv, lambda, sigma, dv, num_v_steps) {

  U0 <- matrix(rep(0, num_v_steps), nrow=num_v_steps)
  for (i in 1:num_v_steps) {
    U0[i,1] <- P_v(uv + dv*(i-1),uv,lambda,v)
  }

  return(U0)
}

Evolve_FEM <- function(U0, num_t_steps, uv, lambda, sigma, dt, dv, num_v_steps) {
  if (num_t_steps <= 0) {
    return(U0)
  }

  M <- diag(x=(2/3)*dv, nrow=num_v_steps, ncol=num_v_steps)
  M[abs(row(M) - col(M)) == 1] <- (1/6)*dv

  S <- diag(x=2/dv, nrow=num_v_steps, ncol=num_v_steps)
  S[abs(row(S) - col(S)) == 1] <- -1/dv

  a <- -uv/dv+1
  jp <- as.integer(ceiling(a))
  jm <- as.integer(ceiling(a-1))

  mask <- c(rep(0,max(min(num_v_steps,jp),0)), rep(1,min(max(num_v_steps - jp,0), num_v_steps)))

  B <- mask*M
  if (jp > 1) {
    B[jp,jp] <- (2/3)*dv - ((a-jm)^3/3)*dv
    B[jm,jp] <- (((jp - a)^2)*(jp - 3*jm+2*a)/6)*dv
    B[jp,jm] <- B[jm,jp]
    B[jm,jm] <- ((jp - a)^3/3)*dv
  }

  E <- solve(M+(1/2)*dt*(((1/2)*sigma^2)*S+lambda*B)) %*% (M-(1/2)*dt*(((1/2)*sigma^2)*S+lambda*B))

  # collecting dead people on the lower and upper ends.
  dead_L <- 0
  dead_H <- 0

  U <- matrix(data=U0, nrow=num_v_steps)
  for (i in 1:num_t_steps) {
    U <- E %*% U
    dead_L <- dead_L + (1/2)*(sigma^2)*dt*(U[2,1] - U[1,1])/dv
    dead_H <- dead_H + (1/2)*(sigma^2)*dt*(U[nrow(U)-1,1] - U[nrow(U),1])/dv
  }

  return_list <- list(U, dead_L, dead_H)
  return(return_list)
}

L_PDE_FEM <- function(X, v, uv, lambda, sigma, dt, dv, num_v_steps) {
  TT <- as.integer(tail(X, 1)/dt)
  tt <- unique(as.integer(X[2:(length(X)-1)]/dt))
  tt[1] <- -1 # For technical convenience (we really want to start from the end of the zeroth day).

  U <- Initialize_FEM(v, uv, lambda, sigma, dv, num_v_steps)

  if (v < uv) {
    return(list(0, 0*U))
  }

  np_mask <- Np_Mask(uv, dv, num_v_steps)

  dead_H <- 0
  dead_L <- 0
  for (i in 1:length(tt)) {
    if (i < length(tt)) {
      num_t_steps <- tt[i+1] - tt[i]
      U <- Evolve_FEM(U, num_t_steps, uv, lambda, sigma, dt, dv, num_v_steps)[[1]]
      U <- lambda*dt*np_mask*U
    } else {
      num_t_steps <- TT - tt[i]
      returned_list <- Evolve_FEM(U, num_t_steps, uv, lambda, sigma, dt, dv, num_v_steps)
      U <- returned_list[[1]]
      dead_L <- returned_list[[2]]
      dead_H <- returned_list[[3]]
    }
  }

  val <- sum(U)*dv + dead_L + dead_H

  return_list <- list(val, U)

  return(return_list)
}

L_PDE_FEM_vec <- function(training_Data, v, uv, lambda, sigma, dt, dv, num_v_steps) {
  val_vec <- c()
  for (i in 1:length(training_Data)) {
    val_vec[i] <- L_PDE_FEM(training_Data[[i]], v, uv, lambda, sigma, dt, dv, num_v_steps)[[1]]
  }
  return(val_vec)
}

################################################################################

P_v <- function(v, uv, lambda, v0, sigma=1) {
  if ((v < uv) | (v0 < uv)) {
    return(0)
  }
  particular_soln <- (1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v-v0)/sigma)
  a1 <- (exp(2*sqrt(2*lambda)*max(uv,0)/sigma)/(sigma-sqrt(2*lambda)*min(uv,0)))*sqrt(lambda/2)
  a2 <- 2*sqrt(2*lambda)*(v0-uv)*((1-sign(v0))/2)/sigma
  a3 <- (1 + sqrt(2*lambda)*sign(v0)*min(uv,0)/sigma)*exp(-sqrt(2*lambda)*abs(v0)/sigma)
  A <- a1*(a2-a3)
  gen_soln <- A*exp(-sqrt(2*lambda)*v/sigma)
  soln <- particular_soln + gen_soln

  return(soln)
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

UV_fn <- function(uc, lambda, rr, sigma) {
  val_ <- sigma/(sqrt(2*rr))
  if (lambda >= exp(uc/val_)/val_) {
    val <- uc - val_*log(val_*lambda)
  } else {
    val <- exp(uc/val_)/lambda - val_
  }
  return(val)
}

c_From_uc <- function(uc, rr, sigma) {
  val_ <- sigma/(sqrt(2*rr))
  val <- var_*exp(uc/var_)
  return(val)
}

uc_From_c <- function(c, rr, sigma) {
  val_ <- sigma/(sqrt(2*rr))
  val <- var_*log(c/var_)
  return(val)
}

file_sample_v_lambda <- 'sample_v_lambda.csv'

################################################################################
#
# Estimate the model parameters (v0, c, r, alpha) given:
# NN      = number of the simulated random-walks.
# sigma   = diffusion rate (default=1).
# rr      = discounting factor (default=1).
# dt      = time step (default=1)
# Replicate = Option to either simulate new set of Brownian motions for every new id or just copy the same simulation across for all consumers
# training_Data = list of (x,t1,t2,...,tx,T) for each individual customers.
# method  = Finite Element (FEM) or Finite Difference (FDM).
#

Estimate_Parameters_PDE <- function(training_Data, C0=list(v0=0.6973,uc=-2.8896, r=0.3922,alpha=19.8473),
                                    sigma=1, rr=1, dt=1/7, dv=1, max_v=100,
                                    num_sample_points=100, min_sample_lambda=0, max_sample_lambda=1, method="FEM") {

  num_v_steps <- as.integer(max_v/dv)
  sample_lambda <- runif(num_sample_points, min=min_sample_lambda, max=max_sample_lambda)

  ptm <- proc.time()

  epsilon <- 1e-10 # Regularization parameter.
  count <- 0

  nLL_PDE_het <- function(v0, uc, r, alpha) {
    v_het_vec <- rep(0, num_ids)
    P_vec <- rep(0, num_sample_points)
    N1 <- rep(0, num_sample_points)
    for (i in 1:num_sample_points) {
      uv <- UV_fn(uc, sample_lambda[i], rr, sigma)
      P_vec[i] <- P_lambda(sample_lambda[i], r, alpha)
      if (method == "FEM") {
        v_het_vec <- v_het_vec + L_PDE_FEM_vec(training_Data, v0, uv, sample_lambda[i], sigma, dt, dv, num_v_steps)*P_vec[i]
        N1[i] <- sum(Initialize_FEM(v0, uv, sample_lambda[i], sigma, dv, num_v_steps))*dv
      } else if (method == "FDM") {
        v_het_vec <- v_het_vec + L_PDE_FDM_vec(training_Data, v0, uv, sample_lambda[i], sigma, dt, dv, num_v_steps)*P_vec[i]
        N1[i] <- sum(Initialize_FDM(v0, uv, sample_lambda[i], sigma, dv, num_v_steps))*dv
      } else {
        stop("Unknown method specified.")
      }
    }
    P_vec <- P_vec*N1
    if (sum(P_vec) > 0) {
      v_het_vec <- v_het_vec/sum(P_vec)
    }

    LL <- sum(log(v_het_vec + epsilon))
    if (count%%1==0){
      dptm <- proc.time() - ptm
      print(paste("LL = ", as.character(LL), ", theta = ", paste(c(v0, uc, r, alpha), collapse=", "), ', time elapsed = ', as.integer(dptm[[3]]/60)))
    }
    count <<- count + 1
    return(-LL)
  }

  fit <- bbmle::mle2(nLL_PDE_het, start=C0, lower=c(-5, -10, 0.1, 0.1), upper=c(3,0,2,50))

  C[3] <- c_From_uc(C[3], rr, sigma)
  C <- as.numeric(fit@coef)
  return(C)
}

################################################################################
#
# Estimate the model parameters (v0, c, r, alpha) given:
# NN      = number of the simulated random-walks.
# sigma   = diffusion rate (default=1).
# rr      = discounting factor (default=1).
# dt      = time step (default=1)
# Replicate = Option to either simulate new set of Brownian motions for every new id or just copy the same simulation across for all consumers
# training_Data = either list of (x,t1,t2,...,tx,T), or (x,t1,t2,...,tx,T,Tm) for each individual customers.
#                 where Tm is the 'membership' period (observation information tracking period). If Tm is not given, then we assume information tracking is not observable).
#

Estimate_Parameters <- function(training_Data, C0=list(v0=1.5,uc=-3, r=0.1,alpha=20),
                                NN=100, sigma=1, rr=1, dt=1, Replicate=TRUE,
                                num_sample_points=200, min_sample_v=0, max_sample_v=15, min_sample_lambda=0, max_sample_lambda=1) {

  testing_period <- 0
  returned_df <- Prepare_Simulation(training_Data, NN, sigma, rr, dt, Replicate, testing_period)

  sample_v <- runif(num_sample_points, min=min_sample_v, max=max_sample_v)
  sample_lambda <- runif(num_sample_points, min=min_sample_lambda, max=max_sample_lambda)
  write.csv(data.frame(sample_v=sample_v, sample_lambda=sample_lambda), file ='sample_v_lambda.csv', row.names=FALSE)

  ptm <- proc.time()

  epsilon <- 1e-10 # Regularization parameter.
  count <- 0

  nLL_B_het <- function(v0, uc, r, alpha) {
    v_het_vec <- rep(0, length(training_Data))
    P_vec <- rep(0, num_sample_points)
    for (i in 1:num_sample_points) {
      uv <- UV_fn(uc, sample_lambda[i], rr, sigma)
      P_vec[i] <- P_v(sample_v[i], uv, sample_lambda[i], v0, sigma)*P_lambda(sample_lambda[i], r, alpha)*(sample_v[i] > uv)
      v_het_vec <- v_het_vec + L_B(sample_v[i], uv, sample_lambda[i], dt, returned_df, length(training_Data), NN)*P_vec[i]
    }
    tryCatch(
      {
        if (sum(P_vec) > 0) {
          v_het_vec <- v_het_vec/sum(P_vec)
        }
      },
      error = function(e) {
        v_het_vec <- rep(0, length(training_Data))
      }
    )

    LL <- sum(log(v_het_vec + epsilon))
    if (count%%1==0){
      dptm <- proc.time() - ptm
      print(paste("LL = ", as.character(LL), ", theta = ", paste(c(v0, uc, r, alpha), collapse=", "), ', time elapsed = ', as.integer(dptm[[3]]/60)))
    }
    count <<- count + 1
    return(-LL)
  }

  print("Running MLE...")
  fit <- bbmle::mle2(nLL_B_het, start=C0, lower=c(-10, -20, 1e-100, 0.1), upper=c(5,1,10,100))

  C[3] <- c_From_uc(C[3], rr, sigma)
  C <- as.numeric(fit@coef)
  return(C)
}

################################################################################

Make_Prediction <- function(C, observed_Data, testing_period, NN=100, sigma=1, rr=1, dt=1, Replicate=TRUE,
                                num_sample_points=200, min_sample_v=0, max_sample_v=15, min_sample_lambda=0, max_sample_lambda=1, Use_Stored_Sampled_Points) {

  returned_df <- Prepare_Simulation(training_Data, NN, sigma, rr, dt, Replicate, testing_period)

  if(missing(Use_Stored_Sampled_Points)){
    Use_Stored_Sampled_Points <- file.exists(file_sample_v_lambda)
  }
  if (Use_Stored_Sampled_Points) {
    sampled <- read.csv(file_sample_v_lambda)
    sample_v <- sampled$sample_v
    sample_lambda <- sampled$sample_lambda
    num_sample_points <- length(sampled)
  } else {
    sample_v <- runif(num_sample_points, min=min_sample_v, max=max_sample_v)
    sample_lambda <- runif(num_sample_points, min=min_sample_lambda, max=max_sample_lambda)
  }

  sample_uv <- c()
  for (i in 1:length(sample_lambda)) {
    sample_uv[i] <- UV_fn(C[2], sample_lambda[i], rr, sigma)
  }

  print("Pre-compute the homogeneous likelihood function")
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = num_sample_points, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "█")   # Character used to create the bar

  L_B_df_ <- list()
  for (i in 1:num_sample_points) {
    L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i], dt, returned_df, length(training_Data), NN)
    setTxtProgressBar(pb, i)
  }
  print("Finish")
  L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
  rm(L_B_df_)

  print("Pre-compute the probability vector")
  P_vec <- rep(0, num_sample_points)
  for (i in 1:num_sample_points) {
    P_vec[i] <- P_v(sample_v[i], sample_uv[i], sample_lambda[i], C[1], sigma)*P_lambda(sample_lambda[i], C[3], C[4])*(sample_v[i] > sample_uv[i])
  }

  print("Pre-compute the heterogeneous likelihood function")
  L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

  print("Pre-compute the posterior probability vector")
  P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

  ##############################################################################

  prediction_B <- rep(0,8)

  rep_txn_count <- rep(0,8)
  predict_rep_txn_B <- rep(0,8)

  EXt_B_het_vec <- rep(0, num_ids)

  print("Making prediction...")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = num_sample_points, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "█")   # Character used to create the bar

  for (i in 1:num_sample_points) {
    predict_rep_txn_B <- predict_rep_txn_B + PXt_B(sample_v[i], sample_uv[i], sample_lambda[i], dt, returned_df, length(training_Data), NN)*P_vec[i]
    EXt_B_het_vec <- EXt_B_het_vec + EXt_B(sample_v[i], sample_uv[i], sample_lambda[i], dt, returned_df, length(training_Data), NN)*P_pos_df[,i]
    setTxtProgressBar(pb, i)
  }
  print("Finish")
  predict_rep_txn_B <- predict_rep_txn_B/sum(P_vec)
  EXt_B_het_vec <- EXt_B_het_vec/apply(P_pos_df, MARGIN=1, sum)

  for (i in 1:num_ids) {
    rep_txn <- length(unique(id_training_Data[[i]]$day))-1
    if (rep_txn > 7) {
      rep_txn <- 7
    }
    rep_txn_count[rep_txn+1] <- rep_txn_count[rep_txn+1] + 1

    prediction_B[rep_txn+1] <- prediction_B[rep_txn+1] + EXt_B_het_vec[i]
  }

  for (i in 1:length(prediction_B)) {
    prediction_B[i] <- prediction_B[i]/rep_txn_count[i]
  }

  return(list(prediction_B, EXt_B_het_vec, predict_rep_txn_B, rep_txn_count))
}



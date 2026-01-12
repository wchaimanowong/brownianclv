
run_demo <- function(run_estimation=TRUE) {

  Data <- read.table(paste(getwd(),'/CDNOW_membership_aug_smc.txt',sep=''), header=TRUE, sep=",")

  # Get the number of consumers.
  num_ids <- length(unique(Data[,2]))

  # 39 weeks calibration time period
  train_period <- 39
  # length of test period
  test_period <- 39

  # TTT vector stores the amount of time for each consumer from joining the firm to the end of 39 weeks calibration period.
  TTT <- c()

  # In the following we prepare the data by separating them into training (calibration period) and testing set.
  training_Data <- list()
  test_period_txns <- c()

  training_Data_alive <- list() # For evaluation. Since tau is observable, only need to evaluate on the consumers with tau > train_period.
  test_period_txns_alive <- c()
  num_ids_alive <- 1

  for(i in 1:num_ids){
    id_Data <- Data[Data$id == i,]

    tau <- tail(id_Data$day, 1)
    id_Data <- id_Data[1:(nrow(id_Data)-1),]

    start_date <- id_Data$day[1]
    id_Data$day <- id_Data$day - start_date
    TTT[i] <- train_period - start_date
    tau < tau - start_date

    training_Data[[i]] <- c(length(id_Data[id_Data$day < TTT[i],]$day), id_Data[id_Data$day < TTT[i],]$day, tau, TTT[i])

    test_period_txns[i] <- length(unique(id_Data[id_Data$day >= TTT[i],]$day))

    if (tau > TTT[i]) {
      training_Data_alive[[num_ids_alive]] <- training_Data[[i]]
      test_period_txns_alive[num_ids_alive] <- test_period_txns[i]
      num_ids_alive <- num_ids_alive + 1
    }
  }
  num_ids_alive <- length(training_Data_alive)

  ##

  if (run_estimation) {
    print("Running Parameters Estimation...")
    C <- Estimate_Parameters_PDE(training_Data, C0=list(v0=0.1,uc=-3, r=0.8,alpha=12), num_sample_points=100, observe_tau=TRUE)
  } else {
    C <- c(1.0342502390176, -4.04355137720461, 1.04760694706096, 26.0696534069067)
    print("Skipping Estimation. Use the known estimated result for evaluation...")
  }

  ##

  actual_test_txns_alive <- rep(0, 8)
  rep_txns_count_alive <- rep(0, 8)
  for (i in 1:num_ids_alive) {
    x <- training_Data_alive[[i]][1]
    rep_txn <- length(unique(training_Data_alive[[i]][2:(x+1)]))-1
    if (rep_txn > 7) {
      rep_txn <- 7
    }
    rep_txns_count_alive[rep_txn+1] <- rep_txns_count_alive[rep_txn+1] + 1
    actual_test_txns_alive[rep_txn+1] <- actual_test_txns_alive[rep_txn+1] + test_period_txns_alive[i]
  }
  actual_test_txns_alive <- actual_test_txns_alive/rep_txns_count_alive

  ##

  res <- Make_Prediction(C, training_Data_alive, 39, num_sample_points=100, observe_tau=TRUE)
  res <- res$predict_txns

  ##

  prediction <- res[[1]]
  EXt_PDE_het_vec <- res[[2]]

  print(prediction)
  # Compute MAE:
  MAE <- sum(abs(test_period_txns_alive - EXt_PDE_het_vec))/length(test_period_txns_alive)
  print(MAE)

  plot(res[[1]], type='o',ylim=c(0,8))
  points(actual_test_txns_alive, type='o', col='blue')
}

run_demo_noncontractual <- function(run_estimation=TRUE) {

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
    TTT[i] <- train_period - start_date
    training_Data[[i]] <- c(length(id_Data[id_Data$day < TTT[i],]$day), id_Data[id_Data$day < TTT[i],]$day, TTT[i])

    test_period_txns[i] <- length(unique(id_Data[id_Data$day >= TTT[i],]$day))
  }

  ##

  if (run_estimation) {
    print("Running Parameters Estimation...")
    C <- Estimate_Parameters_PDE(training_Data, C0=list(v0=0.1,uc=-3, r=0.8,alpha=8), num_sample_points=100, observe_tau=FALSE)
  } else {
    C <- c(0.174182301150983, -2.61471226189415, 0.889411482884871, 8.07522727474511)
    print("Skipping Estimation. Use the known estimated result for evaluation...")
  }

  ##

  actual_test_txns <- rep(0, 8)
  rep_txns_count <- rep(0, 8)
  for (i in 1:num_ids) {
    x <- training_Data[[i]][1]
    rep_txn <- length(unique(training_Data[[i]][2:(x+1)]))-1
    if (rep_txn > 7) {
      rep_txn <- 7
    }
    rep_txns_count[rep_txn+1] <- rep_txns_count[rep_txn+1] + 1
    actual_test_txns[rep_txn+1] <- actual_test_txns[rep_txn+1] + test_period_txns[i]
  }
  actual_test_txns <- actual_test_txns/rep_txns_count

  ##

  res <- Make_Prediction(C, training_Data, 39, num_sample_points=100, observe_tau=FALSE)
  res <- res$predict_txns

  ##

  prediction <- res[[1]]
  EXt_PDE_het_vec <- res[[2]]

  print(prediction)
  # Compute MAE:
  MAE <- sum(abs(test_period_txns - EXt_PDE_het_vec))/length(test_period_txns)
  print(MAE)

  plot(res[[1]], type='o',ylim=c(0,8))
  points(actual_test_txns, type='o', col='blue')
}

################################################################################

P_v <- function(v, uv, lambda, v0, sigma=1) {
  if ((v < max(uv,0)) | (v0 < uv)) {
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
  if (lambda >= exp(uc/val_)) {
    val <- uc - val_*log(lambda)
  } else {
    val <- val_*exp(uc/val_)/lambda - val_
  }
  return(val)
}

c_From_uc <- function(uc, rr, sigma) {
  val_ <- sigma/(sqrt(2*rr))
  val <- val_*exp(uc/val_)
  return(val)
}

uc_From_c <- function(c, rr, sigma) {
  val_ <- sigma/(sqrt(2*rr))
  val <- val_*log(c/val_)
  return(val)
}

kernel <- function(v1, v0, uv, dt, sigma=1) {
  if ((v1 < uv) | (v0 < uv)) {
    return(0)
  }
  val <- (exp(-(v1-v0)^2/(2*dt*sigma^2))-exp(-(v1-2*uv+v0)^2/(2*dt*sigma^2)))/(sigma*sqrt(2*pi*dt))
  return(val)
}

################################################################################

Generate_Diffusion_Kernel <- function(uv, sigma, dt, dv, num_v_steps) {
  K <- matrix(data=0, nrow=num_v_steps, ncol=num_v_steps)
  for (ii in 1:num_v_steps) {
    for (jj in 1:num_v_steps) {
      v0 <- (ii-1)*dv + uv
      v1 <- (jj-1)*dv + uv
      K[ii, jj] <- kernel(v1, v0, uv, dt, sigma)*dv
    }
  }
  return(K)
}

Generate_Selection_Kernel <- function(uv, lambda, dt, dv, num_v_steps) {
  R <- c()
  for (ii in 1:num_v_steps) {
    R[ii] <- ((ii-1)*dv + uv >= 0)*(1 - exp(-lambda*dt))
    #R[ii] <- ((ii-1)*dv + uv >= 0)*lambda*dt
  }
  return(R)
}

Initialize_Density <- function(v, uv, lambda, sigma, dv, num_v_steps) {

  U0 <- matrix(rep(0, num_v_steps), nrow=num_v_steps)
  for (ii in 1:num_v_steps) {
    U0[ii,1] <- P_v(uv + dv*(ii-1),uv,lambda,v, sigma)
  }

  if (sum(U0) > 0) {
    return(U0/(sum(U0)*dv))
  } else {
    return(0*U0)
  }
}

L_PDE <- function(H, U, K, R, dv) {
  dead_L <- 0
  for (ii in 1:length(H)) {
    if (H[ii]==0) { # No purchase
      U <- (1-R)*(K %*% U)
    }
    if (H[ii]==1) { # A purchase
      U <- R*(K %*% U)
    }
    if (H[ii]==2) { # An observed drop-out
      val <- sum((U - (K %*% U))*dv)
      return(val)
    }
    if (H[ii]==3) { # No purchase, potential non-observable drop-out
      U_ <- K %*% U
      dead_L <- dead_L + sum((U - U_)*dv)
      U <- (1-R)*U_
    }
  }
  val <- sum(U*dv) + dead_L
  return(val)
}

L_PDE_vec <- function(X_Data, U, K, R, dv) {
  val_vec <- c()
  for (i in 1:length(X_Data)) {
    val_vec[i] <- L_PDE(X_Data[[i]], U, K, R, dv)
  }
  return(val_vec)
}

# Expected number of purchases over the period [T, T+t] conditioned on the purchase history over [0,T] of every consumers,
# where t = testing_period.
EXt_PDE <- function(X_Data, U, K, R, dv, testing_period) {

  E_rep_txn <- rep(0, length(X_Data))
  U0 <- U

  for (jj in 1:length(X_Data)) {
    U <- U0
    dead_L <- 0
    H <- X_Data[[jj]]
    if (length(H) > 0) {
      for (ii in 1:length(H)) {
        if (H[ii]==0) {
          U <- (1-R)*(K %*% U)
        }
        if (H[ii]==1) {
          U <- R*(K %*% U)
        }
        if (H[ii]==2) {
          U <- 0*U
        }
        if (H[ii]==3) {
          U_ <- K %*% U
          dead_L <- dead_L + sum((U - U_)*dv)
          U <- (1-R)*U_
        }
      }
    }

    norm <- sum(U*dv) + dead_L
    if (norm > 0) {
      val <- 0
      U <- U/norm
      for (ii in 1:testing_period) {
        U <- K %*% U
        val <- val + sum(R*U*dv)
      }
      E_rep_txn[jj] <- val
    }
  }
  return(E_rep_txn)
}

PSegments <- function(X_Data, U, K, R, dv, testing_period) {

  PBuyers <- matrix(data=0, nrow=length(X_Data), ncol=testing_period)
  PTrackers <- matrix(data=0, nrow=length(X_Data), ncol=testing_period)
  PChurners <- matrix(data=0, nrow=length(X_Data), ncol=testing_period)

  U0 <- U
  buyers_mask <- 1*(R > 0)
  trackers_mask <- 1 - buyers_mask

  for (jj in 1:length(X_Data)) {
    U <- U0
    dead_L <- 0
    H <- X_Data[[jj]]
    if (length(H) > 0) {
      for (ii in 1:length(H)) {
        if (H[ii]==0) {
          U <- (1-R)*(K %*% U)
        }
        if (H[ii]==1) {
          U <- R*(K %*% U)
        }
        if (H[ii]==2) {
          U <- 0*U
        }
        if (H[ii]==3) {
          U_ <- K %*% U
          dead_L <- dead_L + sum((U - U_)*dv)
          U <- (1-R)*U_
        }
      }
    }

    norm <- sum(U*dv) + dead_L
    if (norm > 0) {
      U <- U/norm
      for (ii in 1:testing_period) {
        PBuyers[jj, ii] <- sum(buyers_mask*U*dv)
        PTrackers[jj, ii] <- sum(trackers_mask*U*dv)
        PChurners[jj, ii] <- 1 - sum(U*dv)
        U <- K %*% U
      }
    }
  }

  return(list(PBuyers, PTrackers, PChurners))
}

Estimate_Parameters_PDE <- function(training_Data, C0=list(v0=0.6973,uc=-2.8896, r=0.3922,alpha=19.8473),
                                    Clow=c(-5, -10, 1e-100, 0.1), Chigh=c(5,1,10,100),
                                    sigma=1, rr=1, dt=1, dv=0.3, max_v=30,
                                    num_sample_points=30, min_sample_lambda=0, max_sample_lambda=1.0,
                                    observe_tau=TRUE, fixes_uc=FALSE) {

  X_Data <- list()
  for (ii in 1:length(training_Data)) {
    X <- training_Data[[ii]]
    x <- X[1]
    TT <- as.integer(tail(X, 1)/dt)
    tt <- unique(as.integer(X[2:(x+1)]/dt))

    X_Data[[ii]] <- rep(0, TT)
    for (jj in 2:length(tt)) {
      X_Data[[ii]][tt[jj]] <- 1
    }

    if (observe_tau == TRUE) {
      tau <- as.integer(tail(X,2)[1]/dt)
      if (tau <= TT) {
        X_Data[[ii]][tau] <- 2
      }
    } else {
      if (tt[length(tt)] + 1 <= TT) {
        for (jj in (tt[length(tt)]+1):TT) {
          X_Data[[ii]][jj] <- 3
        }
      }
    }
  }

  num_ids <- length(X_Data)
  num_v_steps <- as.integer(max_v/dv)
  sample_lambda <- seq(from=min_sample_lambda,to=max_sample_lambda,length.out=num_sample_points+1)
  sample_lambda <- sample_lambda[2:(num_sample_points+1)]

  ptm <- proc.time()

  epsilon <- 1e-100 # Regularization parameter.
  count <- 0

  if (fixes_uc) {
    uc <- C0$uc

    sample_uv <- c()
    K_list <- list()
    R_list <- list()
    for (i in 1:num_sample_points) {
      sample_uv[i] <- UV_fn(uc, sample_lambda[i], rr, sigma)
      K_list[[i]] <- Generate_Diffusion_Kernel(sample_uv[i], sigma, dt, dv, num_v_steps)
      R_list[[i]] <- Generate_Selection_Kernel(sample_uv[i], sample_lambda[i], dt, dv, num_v_steps)
    }

    nLL_PDE_het <- function(v0, r, alpha) {
      LL <- num_ids*log(epsilon)

      tryCatch(
        {
          v_het_vec <- rep(0, num_ids)
          P_vec <- rep(0, num_sample_points)
          for (i in 1:num_sample_points) {
            P_vec[i] <- P_lambda(sample_lambda[i], r, alpha)
            U <- Initialize_Density(v0, sample_uv[i], sample_lambda[i], sigma, dv, num_v_steps)

            v_het_vec <- v_het_vec + L_PDE_vec(X_Data, U, K_list[[i]], R_list[[i]], dv)*P_vec[i]
          }

          if (sum(P_vec) > 0) {
            v_het_vec <- v_het_vec/sum(P_vec)
          }
          LL <- sum(log(v_het_vec + epsilon))
          if (is.nan(LL)) {
            LL <- num_ids*log(epsilon)
          }
        },
        error = function(e) {
          LL <- num_ids*log(epsilon)
        }
      )

      if (count%%1==0){
        dptm <- proc.time() - ptm
        print(paste("LL = ", as.character(LL), ", theta = ", paste(c(v0, uc, r, alpha), collapse=", "), ', time elapsed = ', as.integer(dptm[[3]]/60)))
      }
      count <<- count + 1
      return(-LL)
    }

    C0_ <- list(v0=C0$v0, r=C0$r, alpha=C0$alpha)
    Chigh_ <- c(Chigh[1], Chigh[3], Chigh[4])
    Clow_ <- c(Clow[1], Clow[3], Clow[4])
    fit <- bbmle::mle2(nLL_PDE_het, start=C0_, lower=Clow_, upper=Chigh_)

  } else {

    nLL_PDE_het <- function(v0, uc, r, alpha) {
      LL <- num_ids*log(epsilon)
      tryCatch(
        {
          v_het_vec <- rep(0, num_ids)
          P_vec <- rep(0, num_sample_points)
          for (i in 1:num_sample_points) {
            uv <- UV_fn(uc, sample_lambda[i], rr, sigma)
            P_vec[i] <- P_lambda(sample_lambda[i], r, alpha)

            U <- Initialize_Density(v0, uv, sample_lambda[i], sigma, dv, num_v_steps)
            K <- Generate_Diffusion_Kernel(uv, sigma, dt, dv, num_v_steps)
            R <- Generate_Selection_Kernel(uv, sample_lambda[i], dt, dv, num_v_steps)

            v_het_vec <- v_het_vec + L_PDE_vec(X_Data, U, K, R, dv)*P_vec[i]
          }

          if (sum(P_vec) > 0) {
            v_het_vec <- v_het_vec/sum(P_vec)
          }
          LL <- sum(log(v_het_vec + epsilon))
          if (is.nan(LL)) {
            LL <- num_ids*log(epsilon)
          }
        },
        error = function(e) {
          LL <- num_ids*log(epsilon)
        }
      )

      if (count%%1==0){
        dptm <- proc.time() - ptm
        print(paste("LL = ", as.character(LL), ", theta = ", paste(c(v0, uc, r, alpha), collapse=", "), ', time elapsed = ', as.integer(dptm[[3]]/60)))
      }
      count <<- count + 1
      return(-LL)
    }

    fit <- bbmle::mle2(nLL_PDE_het, start=C0, lower=Clow, upper=Chigh)
  }

  C <- as.numeric(fit@coef)
  #C[3] <- c_From_uc(C[3], rr, sigma)
  return(C)
}

Make_Prediction <- function(C, observed_Data, testing_period, sigma=1, rr=1, dt=1, dv=0.3, max_v=30,
                            num_sample_points=30, min_sample_lambda=0, max_sample_lambda=1.0,
                            predict_txns=TRUE, predict_segments=FALSE, observe_tau=TRUE) {

  result_list <- list()

  X_Data <- list()
  for (ii in 1:length(observed_Data)) {
    X <- observed_Data[[ii]]
    x <- X[1]
    TT <- as.integer(tail(X, 1)/dt)
    tt <- unique(as.integer(X[2:(x+1)]/dt))

    X_Data[[ii]] <- rep(0, TT)
    for (jj in 2:length(tt)) {
      X_Data[[ii]][tt[jj]] <- 1
    }

    if (observe_tau == TRUE) {
      tau <- as.integer(tail(X,2)[1]/dt)
      if (tau <= TT) {
        X_Data[[ii]][tau] <- 2
      }
    } else {
      if (tt[length(tt)] + 1 <= TT) {
        for (jj in (tt[length(tt)]+1):TT) {
          X_Data[[ii]][jj] <- 3
        }
      }
    }
  }

  num_ids <- length(X_Data)
  num_v_steps <- as.integer(max_v/dv)
  sample_lambda <- seq(from=min_sample_lambda,to=max_sample_lambda,length.out=num_sample_points+1)
  sample_lambda <- sample_lambda[2:(num_sample_points+1)]

  sample_uv <- c()
  for (i in 1:num_sample_points) {
    sample_uv[i] <- UV_fn(C[2], sample_lambda[i], rr, sigma)
  }

  P_vec <- rep(0, num_sample_points)
  for (i in 1:num_sample_points) {
    P_vec[i] <- P_lambda(sample_lambda[i], C[3], C[4])
  }

  # Initializes the progress bar
  print("Evaluating the model, please wait...")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = num_sample_points, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "█")   # Character used to create the bar

  L_matrix <- matrix(data=0, nrow=num_sample_points, ncol=num_ids)

  if (predict_txns) {
    EXt_matrix <- matrix(data=0, nrow=num_sample_points, ncol=num_ids)
    #PXt_matrix <- matrix(data=0, nrow=num_sample_points, ncol=num_ids)
  }

  if (predict_segments) {
    PBuyers_array <- array(data=0, c(num_sample_points, num_ids, testing_period))
    PTrackers_array <- array(data=0, c(num_sample_points, num_ids, testing_period))
    PChurners_array <- array(data=0, c(num_sample_points, num_ids, testing_period))
  }

  for (i in 1:num_sample_points) {
    U <- Initialize_Density(C[1], sample_uv[i], sample_lambda[i], sigma, dv, num_v_steps)
    K <- Generate_Diffusion_Kernel(sample_uv[i], sigma, dt, dv, num_v_steps)
    R <- Generate_Selection_Kernel(sample_uv[i], sample_lambda[i], dt, dv, num_v_steps)

    L_matrix[i,] <- L_PDE_vec(X_Data, U, K, R, dv)

    if (predict_txns) {
      EXt_matrix[i,] <- EXt_PDE(X_Data, U, K, R, dv, testing_period)
    }

    if (predict_segments) {
      PSegments_list <- PSegments(X_Data, U, K, R, dv, testing_period)
      PBuyers_array[i,,] <- PSegments_list[[1]]
      PTrackers_array[i,,] <- PSegments_list[[2]]
      PChurners_array[i,,] <- PSegments_list[[3]]
    }

    setTxtProgressBar(pb, i)
  }

  L_het_vec <- colSums(L_matrix*P_vec)/sum(P_vec)
  P_pos_matrix <- t(t(L_matrix*P_vec)/L_het_vec)

  if (predict_txns) {
    EXt_het_vec <- colSums(EXt_matrix*P_pos_matrix)/colSums(P_pos_matrix)

    prediction_B <- rep(0,8)
    rep_txn_count <- rep(0,8)
    #predict_rep_txn_B <- rep(0,8)

    #EXt_B_het_vec <- rep(0, num_ids)

    rep_txn_vec <- rep(0, num_ids)
    for (i in 1:num_ids) {
      x <- observed_Data[[i]][1]
      rep_txn_vec[i] <- length(unique(observed_Data[[i]][2:(x+1)]))-1
    }

    for (i in 1:num_ids) {
      rep_txn <- rep_txn_vec[i]
      if (rep_txn > 7) {
        rep_txn <- 7
      }
      rep_txn_count[rep_txn+1] <- rep_txn_count[rep_txn+1] + 1
      prediction_B[rep_txn+1] <- prediction_B[rep_txn+1] + EXt_het_vec[i]
    }

    for (i in 1:length(prediction_B)) {
      prediction_B[i] <- prediction_B[i]/rep_txn_count[i]
    }

    result_list[['predict_txns']] <- list(prediction_B, EXt_het_vec, rep_txn_count)
  }

  if (predict_segments) {

    PBuyers_het_mat <- colSums(sweep(PBuyers_array, c(1,2), P_pos_matrix, `*`))/colSums(P_pos_matrix)
    PTrackers_het_mat <- colSums(sweep(PTrackers_array, c(1,2), P_pos_matrix, `*`))/colSums(P_pos_matrix)
    PChurners_het_mat <- colSums(sweep(PChurners_array, c(1,2), P_pos_matrix, `*`))/colSums(P_pos_matrix)

    result_list[['predict_segments']] <- list(PBuyers_het_mat, PTrackers_het_mat, PChurners_het_mat)
  }

  return(result_list)
}


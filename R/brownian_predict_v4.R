################################################################################
# Importing the data

Data <- read.table(paste(getwd(),'/CDNOW_sample.txt',sep=''), header=FALSE, sep="")
names(Data) <- c("master_id", "id", "day", "amount", "dollar_amount")
Data$day <- as.Date(as.character(Data$day), "%Y%m%d")
Data$day <- as.numeric(Data$day - min(Data$day))/7

# Get the number of consumers.
num_ids <- Data[nrow(Data),2]

################################################################################
# Separating train and test set.

# 39 weeks calibration time period
train_period <- 39
# length of test period
test_period <- 39

# Total recorded period.
tol_period <- train_period + test_period

# TTT vector stores the amount of time for each consumer from joining the firm to the end of 39 weeks calibration period.
TTT <- c()

# In the following we prepare the data by separating them into training (calibration period) and testing set.
id_training_Data = list()
id_testing_Data = list()
for(i in 1:num_ids){
  id_Data <- Data[Data$id == i,]
  start_date <- id_Data$day[1]
  id_Data$day <- id_Data$day - start_date

  TTT[i] <- train_period - start_date
  id_training_Data[[i]] <- id_Data[id_Data$day < TTT[i],]
  id_testing_Data[[i]] <- id_Data[id_Data$day >= TTT[i],]
}

################################################################################
# Generate N Brownian motions dv_t = sigma*dW_t starting from zero over the time interval [0,TT] with step size dt.

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

################################################################################
# We re-adjust data format to be more suitable for the Brownian motion simulation here.

# Parameters for the simulation.
NN = 100
sigma = 1
dt = 1
# the total length of the time period is given by the variable train_period defined above.

# Option to either simulate new set of Brownian motions for every new id or just copy the same simulation across for all consumers
Replicate <- FALSE

# number of steps of Brownian motion on the entire recorded period.
num_steps <- as.integer(round(tol_period/dt))-1

# Dataframe to store the simulated Brownian motions for all consumers.
v_df <- data.frame(matrix(ncol = num_steps, nrow = 0))
# Dataframe to store the minimum of the simulated Brownian motions for all consumers.
min_v_df <- data.frame(matrix(ncol = num_steps, nrow = 0))

# pre-computed vectors to be assembled as lambda_vec later.
pre_l_vec_1 <- data.frame(matrix(ncol = num_steps, nrow = 0))
pre_l_vec_2 <- data.frame(matrix(ncol = num_steps, nrow = 0))

# Dataframe and vector which will help determine which simulated path is valid. e.g. Paths that reaches v_t = uv before t = tx is invalid (the consumer couldn't have made a transaction at tx if he is dead).
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

for (i in 1:num_ids) {
  X <- id_training_Data[[i]]
  TT <- as.integer(round(TTT[i]/dt))-1 # number of steps of Brownian motion from the first purchase until the end of the calibration period
  TT_test <- as.integer(round(test_period/dt))
  X <- X[X$day > 0,]
  x <- nrow(X)

  # Setting up the vector for transaction history during the calibration period
  tt <- c(0)
  if (x > 0) {
    tt <- X$day
    tt <- as.integer(round((1/dt)*tt))
    # The following for-loop is just to make sure that all the purchases occur before the end of calibration period
    # (ideally, it should do nothing, but just in case...)
    for (j in 1:length(X$day)) {
      tt_ <- as.integer(round((1/dt)*X$day[j]))
      tt_ <- min(TT, max(1, tt_))
      tt[j] <- tt_
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
    print(paste("Simulating id = ", as.character(i)))
    returned_df <- Gen_Brownian_Motions(NN, sigma, tol_period, dt)
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

  select_testi <- c(rep(0,TT), rep(1,TT_test), rep(0,num_steps - TT - TT_test))
  select_testi <- rbind(data.frame(), select_testi)
  select_testi <- select_testi[rep(seq_len(nrow(select_testi)), each=NN),]
  names(select_testi) <- names(select_test)
  select_test_[[i]] <- select_testi

  names(pre_l_vec_1i) <- names(pre_l_vec_1)
  names(pre_l_vec_2i) <- names(pre_l_vec_2)
  pre_l_vec_1_[[i]] <- pre_l_vec_1i
  pre_l_vec_2_[[i]] <- pre_l_vec_2i
}

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

select_test <- dplyr::bind_rows(select_test_)
rm(select_test_)

################################################################################

# Output the big vector of likelihoods of every consumer's purchase history from the given parameters.
L_B <- function(v, uv, lambda) {
  if (uv > v) {
    return(rep(0, num_ids))
  }
  lambda_ <- lambda*dt
  lambda_df <- (-lambda_)*pre_l_vec_1 + (lambda_ - 1)*pre_l_vec_2

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
PXt_B <- function(v, uv, lambda) {
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
EXt_B <- function(v, uv, lambda) {
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
# Homogeneous model

num_sample_points <- 1

C <- c(0.50080226613855, -1.00014350672236, 0.0678818130052089) # From Brownian_Homo/slurm-26456038.out

sample_v <- c(C[1])
sample_uv <- c(C[2])
sample_lambda <- c(C[3])

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

print("Pre-compute the probability vector")
P_vec <- rep(1, num_sample_points)

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Full Heterogeneous model (heterogeneous v with stochastic starting point, lambda with gamma distribution, and shiftable lambda sampling grid)

num_sample_points <- 200

#C <- c(0.1,0.508279947943678,0.8,20)
#C <- c(0.01,0.3,0.01,49)
#C <- c(0.922877402462614, 0.530546656576583, 0.567524543464434, 28.6567229840404)
#C <- c(0.0472716773919088, 0.129677888241166, 0.016315498632508, 0.0873215038220269, 52.108539412836) # lucky guess.
#C <- c(0.0472716773919088, 0.209677888241166, 0.026315498632508, 0.0873215038220269, 52.108539412836)

#C <- c(0.0461838490043611, 0.134736818759431, 0.0164024239052434, 0.0893082528283363, 52.0976390694602)
#C <- c(0.0472716773919088, 0.129677888241166, 0.016315498632508, 0.173215038220269, 45.108539412836)
#C <- c(0.0472716773919088, 0.129677888241166, 0.016315498632508, 0.173215038220269, 40.108539412836)
#C <- c(0.202716773919088, 0.129677888241166, 0.016315498632508, 0.173215038220269, 30.108539412836)
#C <- c(0.502716773919088, 0.08177888241166, 0.010315498632508, 0.173215038220269, 20.108539412836)
#C <- c(1.02716773919088, 0.08177888241166, 0.010315498632508, 0.173215038220269, 20.108539412836)
#C <- c(2.02716773919088, 0.08177888241166, 0.010315498632508, 0.173215038220269, 20.108539412836)

#C <- c(2.02716773919088, 0.08177888241166, 0.010315498632508, 0.173215038220269, 15.108539412836)

P_v <- function(v, uv, lambda, v0) {
  if ((v < uv) | (v0 < uv)) {
    return(0)
  }
  particular_soln <- (1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v-v0)/sigma)
  a1 <- (exp(2*sqrt(2*lambda)*max(uv,0)/sigma)/(sigma-sqrt(2*lambda)*min(uv,0)))*sqrt(lambda/2)
  a2 <- 2*sqrt(2*lambda)*(v0-uv)*((1-sign(v0))/2)/sigma
  a3 <- (1 + sqrt(2*lambda)*sign(v0)*min(uv,0)/sigma)*exp(-sqrt(2*lambda)*abs(v0)/sigma)
  A <- a1*(a2-a3)
  gen_soln <- A*exp(-sqrt(2*lambda)*v/sigma)

  #cheap_gen_soln <- -(1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v+v0-2*uv)/sigma)

  soln <- particular_soln + gen_soln

  return(soln)
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

# Sample num_sample_points from the probability distribution of v
#sample_v <- rep(0.01 + (0:19)/(4/3), 10) # 0-15 with 20 steps
sample_v <- rep(0.01 + (0:19)*2, 10)
sample_lambda <- C[3] + (floor((0:199)/20)/10)*0.6 # 0-0.6 with 10 steps

# Numerically solve for uv given c and lambda:
UV_fn <- function(c, lambda) {
  r <- 0.00183
  val_ <- lambda*sigma/((2^(3/2))*sqrt(r))
  if (c <= val_) {
    val <- sigma*log(c/val_)/sqrt(2*r)
  } else {
    val <- optim(fn=function(uv) {
      return(((lambda/r)*max(uv, 0) - c/r + lambda*sigma*exp(-sqrt(2*r)*abs(uv)/sigma)/((2*r)^(3/2)))^2)
    }, par=c(1), method="BFGS")$par
  }
  return(val)
}

sample_uv <- c()
for (i in 1:length(sample_lambda)) {
  sample_uv[i] <- UV_fn(C[2], sample_lambda[i])
}

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_v(sample_v[i], sample_uv[i], sample_lambda[i], C[1])*P_lambda(sample_lambda[i], C[4], C[5])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Full Heterogeneous model (heterogeneous v with stochastic starting point and lambda with gamma distribution)

num_sample_points <- 200

#C <- c(0.1,0.508279947943678,0.8,20)
#C <- c(0.01,0.3,0.01,49)
#C <- c(0.922877402462614, 0.530546656576583, 0.567524543464434, 28.6567229840404)
#C <- c(-0.380493040847759, 0.3903532, 0.1, 50)

#C <- c(0.0465509844483309, 0.13302945707963, 0.174472162057083, 45.102819468649)
#C <- c(2.02716773919088, 0.08177888241166, 0.173215038220269, 20.108539412836)
#C <- c(2.02325124713597, 0.0897965137218653, 0.176004121504952, 20.0894362854947)
#C <- c(2.02325124713597, 0.0853, 0.176004121504952, 20.0894362854947)
#C <- c(1.67816961022208, 0.0823437370977151, 0.384461823250111, 18.7316696673254)
#C <- c(2.02224389034665, 0.0898161334417807, 0.176008457963492, 20.0883971947567)

#C <- c(1.95594666632141, 1e-04, 0.242598179809051, 4.41368831242036) # Fader
#C <- c(1.71745828407736, 1e-04, 0.553, 10.58) # SMC

#C <- c(2.09407770586707, 0.0558690260176269, 0.513830768871799, 13.6147798922145)
C <- c(1.79667508989251, 0.0826114892040001, 0.0521901467652249, 11.0384973513345)

P_v <- function(v, uv, lambda, v0) {
  if ((v < uv) | (v0 < uv)) {
    return(0)
  }
  particular_soln <- (1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v-v0)/sigma)
  a1 <- (exp(2*sqrt(2*lambda)*max(uv,0)/sigma)/(sigma-sqrt(2*lambda)*min(uv,0)))*sqrt(lambda/2)
  a2 <- 2*sqrt(2*lambda)*(v0-uv)*((1-sign(v0))/2)/sigma
  a3 <- (1 + sqrt(2*lambda)*sign(v0)*min(uv,0)/sigma)*exp(-sqrt(2*lambda)*abs(v0)/sigma)
  A <- a1*(a2-a3)
  gen_soln <- A*exp(-sqrt(2*lambda)*v/sigma)

  #cheap_gen_soln <- -(1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v+v0-2*uv)/sigma)

  soln <- particular_soln + gen_soln

  return(soln)
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

#sample_v <- rep(0.01 + (0:19)/2, 5)
#sample_lambda <- 0.0484241946077642 + floor((0:99)/20)/10

# sample_v <- rep(0.01 + (0:19)/(4/3), 10)
# sample_lambda <- 0.016315498632508 + (floor((0:199)/20)/10)*0.6

sample_v <- rep(0.01 + (0:19)*2, 10)
sample_lambda <- 0.010315498632508 + (floor((0:199)/20)/10)*0.6 # 0-0.6 with 10 steps

#sample_v <- rep(0.01 + (0:9)/(2/3), 10)
#sample_lambda <- 0.016315498632508 + 0.6*floor((0:99)/10)/10

# Numerically solve for uv given c and lambda:
UV_fn <- function(c, lambda) {
  r <- 0.00183
  val_ <- lambda*sigma/((2^(3/2))*sqrt(r))
  if (c <= val_) {
    val <- sigma*log(c/val_)/sqrt(2*r)
  } else {
    val <- optim(fn=function(uv) {
      return(((lambda/r)*max(uv, 0) - c/r + lambda*sigma*exp(-sqrt(2*r)*abs(uv)/sigma)/((2*r)^(3/2)))^2)
    }, par=c(1), method="BFGS")$par
  }
  return(val)
}

sample_uv <- c()
for (i in 1:length(sample_lambda)) {
  sample_uv[i] <- UV_fn(C[2], sample_lambda[i])
}

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_v(sample_v[i], sample_uv[i], sample_lambda[i], C[1])*P_lambda(sample_lambda[i], C[3], C[4])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Heterogeneous model (heterogeneous in v only)

num_sample_points <- 20

#C <- c(0.169837695174262, -0.102753503356024, 0.0165482069500899)
#C <- c(-0.380493040847759, -0.412353717894576, 0.0484241946077642)
C <- c(-0.355908529591659, -0.357239390331844, 0.0580481626172486) # Use this result. output file: Brownian_v_het/slurm-27362187.out

#sample_v <- 0.01 + (0:19)/2
sample_v <- 0.01 + (0:19)/(4/3)

sample_lambda <- rep(C[3], num_sample_points)
sample_uv <- rep(C[2], num_sample_points)

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

P_v <- function(v, uv, lambda, v0) {
  if ((v < uv) | (v0 < uv)) {
    return(0)
  }
  particular_soln <- (1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v-v0)/sigma)
  a1 <- (exp(2*sqrt(2*lambda)*max(uv,0)/sigma)/(sigma-sqrt(2*lambda)*uv))*sqrt(lambda/2)
  a2 <- 2*sqrt(2*lambda)*(v0-uv)*((1-sign(v0))/2)/sigma
  a3 <- (1 + sqrt(2*lambda)*sign(v0)*uv/sigma)*exp(-sqrt(2*lambda)*abs(v0)/sigma)
  A <- a1*(a2-a3)
  gen_soln <- A*exp(-sqrt(2*lambda)*v/sigma)

  #cheap_gen_soln <- -(1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v+v0-2*uv)/sigma)

  soln <- particular_soln + gen_soln

  return(soln)
}

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_v(sample_v[i], sample_uv[i], sample_lambda[i], C[1])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Full Heterogeneous model (heterogeneous v with normal distribution and lambda with gamma distribution)

num_sample_points <- 100

#C <- c(0.0467546416113259, 0.266801568308382, 0.510296815611335, 0.797240750736448, 20.0026856708002)
#C <- c(0.144488291265844, 0.474251735014832, 0.509674590661066, 0.773104014358862, 20.0244990854034)
#C <- c(0.139496337251289, 0.462743870231239, 0.512526779937397, 0.771152124031799, 20.0248562350546)
#C <- c(0.264388786689736, 0.733745813071377, 0.532650048035143, 0.727822987238548, 20.056527019651)
#C <- c(0.266363508377535, 0.735789250292922, 0.532681918844418, 0.728651980691742, 20.0567375841424)
#C <- c(0.337510761037386, 0.884858411047251, 0.529037026606364, 0.722222773261275, 20.0709461357287)
#C <- c(0.412055578171408, 1.11394650935233, 0.5340194462943, 0.472289311814273, 20.1319943751643)
#C <- c(0.412055578172972, 1.11394650935658, 0.534019446294555, 0.472289311810907, 20.1319943751653)
#C <- c(0.406904883444302, 1.11394650935243, 0.534019446294306, 0.473289311814193, 20.1319943751644)

#C <- c(0.412055578171445, 1.00255185841719, 0.534019446294306, 0.425960380632774, 20.1319943751644)
#C <- c(0.412055578171445, 1.05824918388481, 0.534019446294306, 0.449624846223484, 20.1319943751644)
#C <- c(0.412055578171445, 1.11394650935243, 0.534012771051228, 0.473289311814193, 20.1319943751644)

#C <- c(-0.691034955831115, 0.204025269675647, 0.500216564018849, 0.358037020038442, 17.1782120798058)

C <- c(0.412055578171445, 1.05824918388481, 1.0659446294306, 0.449624846223484, 20.1319943751644)

P_v <- function(v, mv, sv) {
  return(dnorm(v, mean=mv, sd=sv))
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

sample_v <- rep(0.01 + (0:19)/4, 5)
sample_lambda <- 0.0647020255861506 + floor((0:99)/20)/10

#sample_v <- rep(0.01 + (0:19)/4, 10)
#sample_lambda <- 0.0647020255861506 + floor((0:199)/20)/20

# # Numerically solve for uv given c and lambda:
UV_fn_backup <- function(c, lambda) {
  r <- 0.00183
  val_ <- lambda*sigma/((2^(3/2))*sqrt(r))
  if (c <= val_) {
    val <- sigma*log(c/val_)/sqrt(2*r)
  } else {
    val <- optim(fn=function(uv) {
      return(((lambda/r)*max(uv, 0) - c/r + lambda*sigma*exp(-sqrt(2*r)*abs(uv)/sigma)/((2*r)^(3/2)))^2)
    }, par=c(1), method="BFGS")$par
  }
  return(val)
}

UV_fn <- function(c, lambda) {
  #r <- 0.00183
  r <- 10
  val_ <- sigma/(sqrt(2*r))
  if (c/lambda <= val_) {
    val <- val_*log((c/lambda)/val_)
  } else {
    val <- c/lambda - val_
  }
  return(val)
}

sample_uv <- c()
for (i in 1:length(sample_lambda)) {
  sample_uv[i] <- UV_fn(C[3], sample_lambda[i])
}

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_v(sample_v[i], C[1], C[2])*P_lambda(sample_lambda[i], C[4], C[5])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Heterogeneous model (heterogeneous in lambda only)

num_sample_points <- 10

C <- c(0.219435857342631, 0.120760558607348, 0.637049057932774, 9.9785601099689)

sample_v <- rep(C[1], num_sample_points)
sample_lambda <- 0.016315498632508 + 0.6*(0:9)/10

# Numerically solve for uv given c and lambda:
UV_fn <- function(c, lambda) {
  r <- 0.00183
  val_ <- lambda*sigma/((2^(3/2))*sqrt(r))
  if (c <= val_) {
    val <- sigma*log(c/val_)/sqrt(2*r)
  } else {
    val <- optim(fn=function(uv) {
      return(((lambda/r)*max(uv, 0) - c/r + lambda*sigma*exp(-sqrt(2*r)*abs(uv)/sigma)/((2*r)^(3/2)))^2)
    }, par=c(1), method="BFGS")$par
  }
  return(val)
}

sample_uv <- c()
for (i in 1:length(sample_lambda)) {
  sample_uv[i] <- UV_fn(C[2], sample_lambda[i])
}

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_lambda(sample_v[i], C[3], C[4])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Full Heterogeneous model (heterogeneous v with stochastic starting point, lambda with gamma distribution, and two levels heterogeneous c).

num_sample_points <- 400

C <- c(1.41715459851061, -2.52593950140306, 0.0283697187386233, 0.0135622667517397, 19.1856015715788)

sample_v <- c(1.68175639933906, 14.6588586014695, 4.4386854080949, 0.180352195166051,
              7.30787758016959, 4.83106069033965, 9.04912149882875, 5.91771305655129,
              5.36496828077361, 6.3044267508667, 0.946240414632484, 12.4713125266135,
              14.5305315835867, 0.148589856689796, 6.79520593257621, 2.35032797441818,
              12.4513229134027, 12.4447450961452, 13.9715758874081, 4.09254186437465,
              14.4178418896627, 0.981508480617777, 12.2438996087294, 5.40755015797913,
              2.99821071792394, 6.02869314607233, 10.32597375568, 11.4275448163971,
              4.89102641120553, 8.51200228789821, 12.6759806158952, 1.7274127935525,
              13.4165578836109, 9.36796639929526, 7.67221031128429, 12.0788589853328,
              3.85767260100693, 1.73310055746697, 10.9999717294704, 0.202537203440443,
              14.7608032927383, 11.1253554734867, 7.83911800710484, 3.32871720311232,
              5.84479393204674, 13.5830064676702, 6.74631116562523, 12.5055061234161,
              5.37854732014239, 4.97662026551552, 2.3006677266676, 1.45288348663598,
              8.23106497758999, 1.39285624609329, 9.58308773115277, 6.0480847209692,
              3.67112486623228, 3.40182922314852, 3.72001728625037, 4.13254538783804,
              9.71813543583266, 10.0815505126957, 13.9036998013034, 11.96499493788,
              7.23158593405969, 12.0130914088804, 12.9581530799624, 10.3506857471075,
              6.22453285963275, 0.135304002324119, 2.5669794541318, 14.6217492094729,
              6.04598848731257, 8.58498829300515, 6.71136092045344, 2.95424562762491,
              4.86650783219375, 1.17263372987509, 0.945171463536099, 4.9563392798882,
              13.7436486559454, 5.82052577636205, 2.02084082644433, 1.86940492596477,
              5.39729385171086, 5.84651477518491, 0.629905157256871, 10.0285513582639,
              14.3273848830722, 5.38871180615388, 11.8551964941435, 5.56532601476647,
              4.45374355302192, 13.5709109762684, 10.6145357270725, 9.27137971739285,
              14.4990428374149, 6.795425824821, 5.79916600487195, 1.82444511679932,
              10.3250255528837, 6.51760258362629, 2.94522893498652, 2.03022295958363,
              14.547169636935, 10.4679824283812, 0.141205850522965, 13.2478842814453,
              12.0837348466739, 5.50893392879516, 12.117449051002, 7.65090495580807,
              5.96571061061695, 1.04357680887915, 4.56443091854453, 2.94602528912947,
              5.53992788423784, 10.3177337476518, 5.17763089388609, 14.014276068192,
              9.83365359948948, 0.213837143965065, 6.21003792970441, 4.01309463777579,
              8.71433231164701, 11.3799366471358, 14.6952496189624, 4.59323834977113,
              14.2327320598997, 12.6250618393533, 7.33495730208233, 11.533305037301,
              2.28467576089315, 10.2144341531675, 8.2252348086331, 1.67483045370318,
              0.778937585419044, 1.21399701107293, 11.343251693761, 13.8670885981992,
              10.9229920338839, 8.45869821379893, 8.79408938577399, 8.0550809239503,
              12.9955804080237, 14.8398527945392, 0.270984236849472, 6.61155207781121,
              10.8779728913214, 0.162767440779135, 3.59986743656918, 10.2742229006253,
              10.0499132671393, 7.20462961122394, 5.6100364762824, 9.61975924554281,
              11.1564333876595, 12.7733552106656, 9.78625066461973, 11.2876278231852,
              2.91323354816996, 2.34193512587808, 2.62497638468631, 14.7571719763801,
              4.52229092130437, 8.13451486174017, 0.949066946050152, 11.4620876242407,
              11.9146275299136, 1.26388339558616, 11.9638189394027, 3.11300296336412,
              9.75951808970422, 6.46011732984334, 4.1578357166145, 4.74275907268748,
              4.9968037661165, 6.05411992873996, 11.8639409670141, 12.7385087648872,
              8.98556564818136, 5.68993639200926, 6.9573287141975, 11.2127236684319,
              3.14459359855391, 2.95436383690685, 8.20902776671574, 1.0020282282494,
              14.5150780188851, 12.1932642278261, 11.9664786430076, 5.64989308244549,
              0.276717846281826, 8.55466129258275, 3.65068127168342, 11.1826774408109,
              2.7645058685448, 0.221175451297313, 6.75391786149703, 6.80724637932144)

sample_lambda <- c(0.28877099766396, 0.181045729201287, 0.869207178475335, 0.270777235738933,
                   0.445768818026409, 0.246901531005278, 0.861071688355878, 0.567081271903589,
                   0.0981276091188192, 0.745980692794546, 0.237747552804649, 0.685128286713734,
                   0.283121787477285, 0.326727104606107, 0.816153206164017, 0.899565434781834,
                   0.849676583660766, 0.265099335694686, 0.673374965088442, 0.599258823785931,
                   0.101662459084764, 0.535403144313022, 0.396784826880321, 0.970845782430843,
                   0.726356945699081, 0.220696562202647, 0.694703424116597, 0.789266465930268,
                   0.953594953985885, 0.673619304317981, 0.324313660850748, 0.277806316269562,
                   0.222741171484813, 0.626627802848816, 0.851336323888972, 0.362433645874262,
                   0.0662564372178167, 0.695889139780775, 0.384728528093547, 0.154468357330188,
                   0.46552997129038, 0.750515274005011, 0.916143166134134, 0.222874156665057,
                   0.55096807773225, 0.733305532485247, 0.193973171059042, 0.461730849696323,
                   0.0033594174310565, 0.796322351321578, 0.652022954076529, 0.205127970781177,
                   0.22447640215978, 0.527833426836878, 0.97914626239799, 0.993534402223304,
                   0.491663225926459, 0.514478042488918, 0.535915561253205, 0.361554715316743,
                   0.924915507668629, 0.755945846438408, 0.837975035654381, 0.940069779520854,
                   0.202085100580007, 0.739735312294215, 0.612356108147651, 0.0492368591949344,
                   0.978784467559308, 0.316798771498725, 0.295370659790933, 0.461773524992168,
                   0.0223920028656721, 0.996213541831821, 0.842417206382379, 0.435597244650126,
                   0.482866417383775, 0.800136500038207, 0.0294757904484868, 0.773646077141166,
                   0.959143539657816, 0.836916426662356, 0.825689256424084, 0.70199436834082,
                   0.63952322024852, 0.761368463980034, 0.469117937143892, 0.530816895887256,
                   0.567746959161013, 0.679469482973218, 0.431550114182755, 0.165897716302425,
                   0.204126075608656, 0.879247019765899, 0.222225658129901, 0.583234689896926,
                   0.586412020726129, 0.176084873732179, 0.89218855346553, 0.632343038916588,
                   0.600022586761042, 0.336616350337863, 0.82141621876508, 0.806797077646479,
                   0.068384894169867, 0.773005675757304, 0.0577706703916192, 0.914563352707773,
                   0.0426511731930077, 0.366253996267915, 0.617098974296823, 0.664095965214074,
                   0.853732397779822, 0.425801560282707, 0.0733121901284903, 0.233156747650355,
                   0.950975372456014, 0.00176013936288655, 0.820933631388471, 0.833172924583778,
                   0.379956264747307, 0.563708007568493, 0.186835870612413, 0.395545297302306,
                   0.26978802238591, 0.51101107057184, 0.958305196836591, 0.378640374401584,
                   0.575883569661528, 0.842421066015959, 0.54712823079899, 0.75399328651838,
                   0.130066653946415, 0.987991761183366, 0.910687063122168, 0.0663509571459144,
                   0.0465190082322806, 0.409583583008498, 0.554352622013539, 0.934585153125226,
                   0.330888736061752, 0.904887262964621, 0.434939611935988, 0.330455516465008,
                   0.199151708045974, 0.975528507959098, 0.426458558533341, 0.382877044612542,
                   0.838527224259451, 0.773845428135246, 0.3095755178947, 0.395411380566657,
                   0.799710903083906, 0.122039126697928, 0.226567079080269, 0.127689001848921,
                   0.501281800447032, 0.346394920488819, 0.704341244418174, 0.732859068783,
                   0.0779368800576776, 0.490814453922212, 0.617251590825617, 0.195223369402811,
                   0.65606351941824, 0.724216090748087, 0.931845772778615, 0.60609636688605,
                   0.546191821573302, 0.539795494172722, 0.344922463409603, 0.840645429911092,
                   0.577919841976836, 0.957324761198834, 0.733654704643413, 0.700743275461718,
                   0.530426640063524, 0.892210603458807, 0.527160326717421, 0.961865193909034,
                   0.839904876891524, 0.924858997110277, 0.00760925328359008, 0.381217803107575,
                   0.363098777830601, 0.752417297801003, 0.751745787682012, 0.369844801258296,
                   0.0229710098356009, 0.629751989850774, 0.937836132477969, 0.540559672983363,
                   0.92005866765976, 0.981949874665588, 0.732799916062504, 0.422171018319204,
                   0.393749018898234, 0.0185405674856156, 0.160095678409562, 0.739314859732985)

UV_fn <- function(uc, lambda) {
  #r <- 0.00183
  #r <- 0.6671713
  r <- 1
  #r <- 7
  #r <- 16.01211
  #r <- 960.7266
  val_ <- sigma/(sqrt(2*r))
  if (uc <= val_*log(val_*lambda)) {
    val <- uc - val_*log(val_*lambda)
  } else {
    val <- exp(uc/val_)/lambda - val_
  }
  return(val)
}

sample_uv <- c()
for (i in 1:length(sample_lambda)) {
  sample_uv[i] <- UV_fn(C[2], sample_lambda[i])
}
plot(sample_lambda, sample_uv)

sample_uv <- c(sample_uv, rep(-1/1e-10, length(sample_lambda)))
P_c <- c(rep(C[3], length(sample_lambda)), rep(1-C[3], length(sample_lambda)))

sample_v <- rep(sample_v, 2)
sample_lambda <- rep(sample_lambda, 2)

#sample_uv <- rep(C[2], num_sample_points)

P_v <- function(v, uv, lambda, v0) {
  if ((v < uv) | (v0 < uv)) {
    return(0)
  }
  particular_soln <- (1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v-v0)/sigma)
  a1 <- (exp(2*sqrt(2*lambda)*max(uv,0)/sigma)/(sigma-sqrt(2*lambda)*min(uv,0)))*sqrt(lambda/2)
  a2 <- 2*sqrt(2*lambda)*(v0-uv)*((1-sign(v0))/2)/sigma
  a3 <- (1 + sqrt(2*lambda)*sign(v0)*min(uv,0)/sigma)*exp(-sqrt(2*lambda)*abs(v0)/sigma)
  A <- a1*(a2-a3)
  gen_soln <- A*exp(-sqrt(2*lambda)*v/sigma)

  #cheap_gen_soln <- -(1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v+v0-2*uv)/sigma)

  soln <- particular_soln + gen_soln

  return(soln)
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_c[i]*P_v(sample_v[i], sample_uv[i], sample_lambda[i], C[1])*P_lambda(sample_lambda[i], C[4], C[5])*(sample_v[i] > sample_uv[i])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Full Heterogeneous model (heterogeneous v with stochastic starting point, lambda with gamma distribution)

num_sample_points <- 200

# Format: C <- c(v,uc,gamma,alpha), c = exp(uc/val_)*val_, val := sigma/(sqrt(2*r)).

#C <- c(0.164078997765089, -0.933513496541919, 0.153271457784792, 10.3108946744492) # Best result
#C <- c(-0.155956281948273, -0.925190099847469, 0.276968584905087, 13.6703293554886)

#C <- c(0.909254488191004, -0.931263241769984, 0.0186728434700175, 12.8452899407879)
# C <- c(0.490383634183538, -1.0256549960035, 0.122693802949891, 10.0106216946594) # Best estimation result where c is used instead of uv (with r = 10% per minute)
#C <- c(0.909266922011551, -0.931255556710074, 0.0207296215257242, 12.8457915516191)

#C <- c(0.490383634183538, -2.0256549960035, 0.122693802949891, 10.0106216946594)
# C <- c(0.490383634183538, -4.0256549960035, 0.122693802949891, 20.0106216946594)
#C <- c(0.272054647311783, -1.50627946452175, 0.163270892361543, 18.8129495041847)

C <- c(2.5, -1.5, 0.1, 39.7366500094235)

sample_v <- c(7.13720126193948, 10.2238979400136, 7.49597621499561, 7.50569379539229,
              7.18730436055921, 3.8663922669366, 1.07285407022573, 10.9597547259182,
              12.1122561267111, 2.72072376217693, 6.54472685535438, 3.49708598689176,
              9.19208584818989, 1.49582108482718, 3.76035970752127, 6.68615782284178,
              11.4044801937416, 9.92076971218921, 10.7349454087671, 0.739916724851355,
              5.31245431280695, 9.53822925221175, 10.6766388821416, 0.171052812365815,
              1.52512885630131, 14.9357167014387, 0.351895397761837, 10.2751337690279,
              1.2036546645686, 10.9039172146004, 1.53357302653603, 0.0735185004305094,
              8.90127816237509, 12.0527442754246, 12.990284326952, 3.47540950868279,
              9.99291821848601, 5.09043989237398, 1.03900327929296, 4.49263473507017,
              8.7231546943076, 0.784267511917278, 14.8369204695337, 11.2060338025913,
              1.53958291746676, 11.0273514711298, 11.6367834771518, 13.6870931228623,
              14.7973999043461, 14.528372248169, 13.5573915811256, 3.96596173988655,
              2.9387904109899, 12.6298406906426, 14.9522504454944, 9.3964993476402,
              1.27641429542564, 14.4396325084381, 3.58587417053059, 10.8667156950105,
              5.13850679271854, 13.1689949997235, 1.68923649354838, 0.0433677015826106,
              2.4745819112286, 10.0202490482479, 0.518755627563223, 5.50950618227944,
              10.2066323684994, 14.7106074972544, 2.76616023620591, 7.32604239718057,
              11.023651646683, 4.04849799466319, 0.531205603620037, 10.5719935172237,
              14.6841580769978, 3.3569122350309, 2.90026302449405, 4.32511754217558,
              7.98159257392399, 3.83351788739674, 4.11937577300705, 2.59253523079678,
              6.42042906023562, 3.81501762662083, 14.2969897238072, 12.8209733613767,
              4.34527621488087, 13.7613824056461, 12.1309254330117, 4.26229657139629,
              6.39668103656732, 8.32342522102408, 3.5106775991153, 2.22841741866432,
              6.00406610174105, 9.49397004442289, 2.42312253685668, 13.437351572793)

sample_lambda <- c(0.140555103309453, 0.58494367916137, 0.152061105938628, 0.538590834010392,
                   0.461030945181847, 0.513353197369725, 0.5659702681005, 0.84105627448298,
                   0.710389111191034, 0.390465783886611, 0.465262437006459, 0.39724030578509,
                   0.26102451630868, 0.236271303379908, 0.580386436311528, 0.143074790714309,
                   0.581496251979843, 0.588644138304517, 0.224104268010706, 0.245935294078663,
                   0.490825387416407, 0.821494305506349, 0.0703488651197404, 0.115425863536075,
                   0.925369435222819, 0.739441704703495, 0.600309837376699, 0.53426411212422,
                   0.555669031804428, 0.895066541386768, 0.858686776831746, 0.579684223979712,
                   0.165503344731405, 0.316273862496018, 0.0345178509596735, 0.215697749285027,
                   0.367810460738838, 0.409888820489869, 0.0214630579575896, 0.607413431163877,
                   0.257724018301815, 0.94451712933369, 0.887699117651209, 0.413100291276351,
                   0.831062626326457, 0.12133949669078, 0.398235511267558, 0.693501035217196,
                   0.484576417598873, 0.689768523676321, 0.333759644534439, 0.12977625685744,
                   0.876496590906754, 0.76721339370124, 0.209656664635986, 0.000362426741048694,
                   0.0509776417165995, 0.444497337797657, 0.321338830050081, 0.122982134344056,
                   0.871424932731315, 0.0314744568895549, 0.52174011990428, 0.343335914891213,
                   0.511416089721024, 0.112242203671485, 0.800337044522166, 0.379967038519681,
                   0.935223913285881, 0.745186195941642, 0.310218111379072, 0.288515453692526,
                   0.739404360996559, 0.523198594572023, 0.310040414799005, 0.100403350312263,
                   0.92810521600768, 0.27272280305624, 0.976653260411695, 0.64409016049467,
                   0.972828670637682, 0.5450317079667, 0.966206596232951, 0.0501550044864416,
                   0.447661128127947, 0.694567354628816, 0.228449807036668, 0.333632431225851,
                   0.111259415047243, 0.0615503578446805, 0.451742341043428, 0.0180644183419645,
                   0.982072587590665, 0.811868412885815, 0.195480178110301, 0.246389716863632,
                   0.174116576788947, 0.479228511219844, 0.820034346776083, 0.420130068436265)

#C <- c(0.497983860353072, -1.00275303883075, 0.0792338645003232, 10.0015198044435)
#C <- c(0.416433043941032, -0.999920159859013, 0.0375330172297031, 10.0647683592327)

#C <- c(0.951884145677044, -1.01607479387659, 0.0550803473930457, 13.9581789002256)

# min_sample_v <- 0
# max_sample_v <- 15
# min_sample_lambda <- 0
# max_sample_lambda <- 1
# sample_v <- runif(num_sample_points, min=min_sample_v, max=max_sample_v)
# sample_lambda <- runif(num_sample_points, min=min_sample_lambda, max=max_sample_lambda)

#C <- c(0.71264545466448, -2.89553745054166, 0.0471483133717288, 20.0534078093932) # Ok result with r = 1, output file = Brownian_rc/slurm-28513286.out
#C <- c(0.71264545466448, -2.89553745054166, 0.1, 10.0534078093932)

#C <- c(0.918314699738427, -2.70588185473939, 0.389543831811068, 19.7649369272809) # ok result with r = 1, output file = Browian_lowrc/slurm-28570936.out

# C <- c(0.340723103429562, -5.20055955326476, 0.119592276457191, 16.9716337950085)

#C <- c(0.697301375189045, -2.88964437986363, 0.392193766159119, 19.8472813975261) # good result with r = 1, output file = Brownian_lowrc/slurm-28638327.out
#C <- c(0.743549195653551, -2.87411460914769, 0.292491104410643, 19.8683612837935)
#C <- c(0.716745886359441, -2.8108360466405, 0.193207900422033, 19.8806937644495)
#C <- c(0.709246202300422, -2.81140554895872, 0.21616758019338, 19.8763329419203)
#C <- c(0.706858804576217, -2.81462222884426, 0.224832883293576, 19.8747005420003)

C <- c(0.589697578571773, -2.86689187729988, 0.219100626533936, 19.8596284614822)

sample_v <- c(0.773395114811137, 5.51823433837853, 3.26167005579919, 12.026479630731,
              14.0034884330817, 11.377408298431, 3.53258965420537, 9.7953331656754,
              1.03316146298312, 4.11616688710637, 7.53666533972137, 14.4966303720139,
              14.2144088353962, 4.60875745047815, 5.01729645766318, 6.92830784479156,
              13.679052782245, 10.8299085195176, 9.38310769968666, 14.8868675564881,
              10.8882546715904, 2.11537606897764, 7.12312591960654, 13.9557074778713,
              13.1812609964982, 9.34394070412964, 11.0273661499377, 10.0169758487027,
              2.82914949464612, 10.0968886585906, 0.553244138136506, 3.22195638553239,
              8.24595014797524, 14.9533457960933, 9.64696808019653, 12.9519532585982,
              14.9500529340003, 3.72396443621255, 5.6025135400705, 12.1320319839288,
              6.44604705506936, 1.47788247093558, 10.9675354266074, 3.79665139247663,
              5.01665957272053, 4.15853122947738, 8.09711758629419, 11.9577210315038,
              14.3883298896253, 6.66617801995017, 4.67849762295373, 5.86131123709492,
              0.748267042217776, 11.8694011552725, 13.2860379049089, 3.65411929087713,
              12.4966848653276, 2.17808011570014, 4.13465558085591, 13.0258886783849,
              14.6977156796493, 12.5795014598407, 14.713460089406, 1.01786527549848,
              8.32038767519407, 13.8059474097099, 6.94979692460038, 10.4090430575889,
              8.72176079079509, 2.37931007752195, 1.6717486793641, 7.80689190840349,
              5.07385652046651, 8.82957576191984, 9.50135428458452, 2.93019606499001,
              12.3058738326654, 13.1326375762001, 14.2944467696361, 7.3322074895259,
              8.81834131549112, 5.72378181503154, 5.26949345832691, 9.70092113711871,
              12.0897173357662, 2.54818336805329, 13.7539094255771, 6.58991970238276,
              13.9231240481604, 6.77996269776486, 12.5025409623049, 9.79022440966219,
              12.9535217920784, 12.1869032399263, 4.9643183068838, 14.484330075793,
              8.6389619554393, 14.0252247091848, 13.9669359731488, 10.0636359932832,
              3.96410651039332, 1.11654797801748, 7.1233709785156, 14.8186848289333,
              8.4605671162717, 14.0065615600906, 6.65005139540881, 13.6063331610058,
              1.53921901714057, 9.41459626425058, 9.16335589718074, 11.0721531871241,
              14.2079554894008, 2.47415474848822, 9.19138787314296, 1.75230489461683,
              10.3488673560787, 2.91338818031363, 10.7083081360906, 7.39544197102077,
              13.0655521235894, 1.75683336798102, 2.46553981443867, 4.4240130437538,
              13.3154701441526, 14.5822258398402, 7.15090394136496, 14.3722087412607,
              1.36257692473009, 9.24397574504837, 13.8000195741188, 1.89274551579729,
              14.4797630875837, 1.81933956337161, 0.104722918476909, 0.728089779149741,
              11.3846770499367, 3.81970194750465, 10.0750595214777, 4.51251212507486,
              12.1225641283672, 4.08224310376681, 9.78553559049033, 12.8833274706267,
              0.413870194461197, 9.17689147288911, 0.395853538066149, 4.0140537545085,
              3.12490988988429, 9.44335940759629, 6.12800719449297, 8.70225074118935,
              10.3474508598447, 10.799492856022, 12.6583824085537, 0.70920335361734,
              11.3693516969215, 14.868645707611, 13.9405201084446, 12.3062485666014,
              3.19599001435563, 8.77255637431517, 4.09948291140608, 7.87537869298831,
              7.38911931286566, 2.73658790742047, 5.6559421797283, 11.3627046940383,
              6.93998685688712, 6.28286913852207, 6.96941800648347, 5.80855274456553,
              9.44531097309664, 2.76464223163202, 7.12833961704746, 10.7568159629591,
              11.8344249681104, 1.72759399632923, 12.4578582833055, 0.41388826793991,
              10.5370021169074, 4.48651309241541, 14.3758657691069, 1.11193540273234,
              13.7969714845531, 1.93259864347056, 8.39422021410428, 6.26732889795676,
              9.72123138955794, 14.6175871964078, 6.55017259530723, 12.937250820687,
              8.68708810536191, 9.71256589866243, 8.65891171968542, 8.60729655018076,
              8.27578565455042, 10.4056992265396, 6.80893631419167, 10.977491409285)

sample_lambda <- c(0.685295628150925, 0.524533801712096, 0.513961545191705, 0.205576573032886,
                   0.316045571351424, 0.171596145723015, 0.51223022933118, 0.320223944727331,
                   0.508400801103562, 0.335510130506009, 0.756420037010685, 0.288455079076812,
                   0.328457705676556, 0.110885615693405, 0.514461618382484, 0.883035372709855,
                   0.302886501187459, 0.994387048529461, 0.263317487901077, 0.579933770233765,
                   0.470102443825454, 0.0570952964480966, 0.0139189676847309, 0.0671389130875468,
                   0.84482044656761, 0.131116523873061, 0.295166296185926, 0.419927153736353,
                   0.27563615818508, 0.119739730842412, 0.0300246207043529, 0.682712926529348,
                   0.81089652585797, 0.149418055778369, 0.772070619277656, 0.056529184570536,
                   0.785818117205054, 0.242483457783237, 0.93470539804548, 0.418799663893878,
                   0.889034810243174, 0.329407720826566, 0.911350788781419, 0.692819541553035,
                   0.347276315558702, 0.438890763791278, 0.335350005421788, 0.0296893406193703,
                   0.0465421725530177, 0.280167021788657, 0.581007023341954, 0.989759563235566,
                   0.0327121850568801, 0.462485596071929, 0.834043571958318, 0.194624674506485,
                   0.231662523234263, 0.227399092633277, 0.966397380689159, 0.768063905416057,
                   0.998717583715916, 0.274012678069994, 0.616476257797331, 0.965854706242681,
                   0.226871553109959, 0.629158539697528, 0.970825708471239, 0.491840912727639,
                   0.739244815194979, 0.153735738014802, 0.777293533785269, 0.602523326640949,
                   0.317582216113806, 0.556380774127319, 0.750585563713685, 0.534460620488971,
                   0.223124605370685, 0.561600901884958, 0.3290586466901, 0.132930608000606,
                   0.544134454568848, 0.64592501684092, 0.633772618370131, 0.76788640092127,
                   0.74926318321377, 0.652755832765251, 0.00953536084853113, 0.74790216004476,
                   0.654999587452039, 0.148094350006431, 0.556747068651021, 0.44671824015677,
                   0.374539733165875, 0.0500805000774562, 0.786104690050706, 0.553079142002389,
                   0.381490986561403, 0.481079438468441, 0.835497152525932, 0.626924729673192,
                   0.846243407810107, 0.43796048918739, 0.301781960530207, 0.0856657726690173,
                   0.616910221753642, 0.118005697615445, 0.956903268350288, 0.13703727722168,
                   0.13974245195277, 0.415412486996502, 0.884126174030825, 0.688635450322181,
                   0.818737848429009, 0.390046489657834, 0.885341444518417, 0.171531500993297,
                   0.443103467812762, 0.818357565905899, 0.986961144953966, 0.89328247262165,
                   0.0489113414660096, 0.244614642113447, 0.153146984288469, 0.792325305519626,
                   0.927232511807233, 0.00916835432872176, 0.535857529146597, 0.765383832389489,
                   0.977358356816694, 0.794880585744977, 0.659598926547915, 0.760237926151603,
                   0.477381481090561, 0.621595251373947, 0.929641046095639, 0.165850933641195,
                   0.221124191535637, 0.378087595105171, 0.0613579403143376, 0.909945865860209,
                   0.110222659772262, 0.826997470343485, 0.965124491136521, 0.49126475234516,
                   0.93802160769701, 0.293630329426378, 0.220367837930098, 0.857231108238921,
                   0.0293482125271112, 0.395686159608886, 0.767798001412302, 0.24205255554989,
                   0.518862440250814, 0.0307877347804606, 0.7040888981428, 0.0197018352337182,
                   0.534341987920925, 0.384019812569022, 0.0794651005417109, 0.290586970048025,
                   0.90120127145201, 0.926142692100257, 0.332915576640517, 0.783631427912042,
                   0.742419369053096, 0.45836979104206, 0.640216623432934, 0.868667581118643,
                   0.362292398698628, 0.533720667241141, 0.988822383573279, 0.439149852609262,
                   0.349274213658646, 0.273833075771108, 0.149327145656571, 0.632428423734382,
                   0.894064476713538, 0.418796831276268, 0.0595508343540132, 0.543815910350531,
                   0.901426865486428, 0.211872841697186, 0.557494023581967, 0.341427452163771,
                   0.464103598147631, 0.808015966787934, 0.946670086821541, 0.703273633029312,
                   0.0588582032360137, 0.480189963942394, 0.675346424104646, 0.508347294759005,
                   0.386148748453707, 0.0372101808898151, 0.722868401790038, 0.808446818031371,
                   0.621624490246177, 0.0593952334020287, 0.983567130984738, 0.11030936287716)

UV_fn <- function(uc, lambda) {
  #r <- 0.00183
  #r <- 0.6671713
  r <- 1
  #r <- 7
  #r <- 16.01211
  #r <- 960.7266
  val_ <- sigma/(sqrt(2*r))
  if (uc <= val_*log(val_*lambda)) {
    val <- uc - val_*log(val_*lambda)
  } else {
    val <- exp(uc/val_)/lambda - val_
  }
  return(val)
}

sample_uv <- c()
for (i in 1:length(sample_lambda)) {
  sample_uv[i] <- UV_fn(C[2], sample_lambda[i])
}
plot(sample_lambda, sample_uv)

#sample_uv <- rep(C[2], num_sample_points)

P_v <- function(v, uv, lambda, v0) {
  if ((v < uv) | (v0 < uv)) {
    return(0)
  }
  particular_soln <- (1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v-v0)/sigma)
  a1 <- (exp(2*sqrt(2*lambda)*max(uv,0)/sigma)/(sigma-sqrt(2*lambda)*min(uv,0)))*sqrt(lambda/2)
  a2 <- 2*sqrt(2*lambda)*(v0-uv)*((1-sign(v0))/2)/sigma
  a3 <- (1 + sqrt(2*lambda)*sign(v0)*min(uv,0)/sigma)*exp(-sqrt(2*lambda)*abs(v0)/sigma)
  A <- a1*(a2-a3)
  gen_soln <- A*exp(-sqrt(2*lambda)*v/sigma)

  #cheap_gen_soln <- -(1/sigma)*sqrt(lambda/2)*exp(-sqrt(2*lambda)*abs(v+v0-2*uv)/sigma)

  soln <- particular_soln + gen_soln

  return(soln)
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_v(sample_v[i], sample_uv[i], sample_lambda[i], C[1])*P_lambda(sample_lambda[i], C[3], C[4])*(sample_v[i] > sample_uv[i])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Heterogeneous in v and lambda with uv instead of c, with normal distribution for v.

num_sample_points <- 100

#C <- c(1.75988791967288, -0.369443114630685, 0.646704870719187, 11.746447134459 )
#C <- c(2.60656001205525, -0.338722896520717, 0.930889818416773, 19.0252707054446)
#C <- c(0.916824778544452, 0.177778831793432, -2.57437692865302, 0.618463034441195, 29.5412781880433)
#C <- c(0.0, 2.205692, -0.369443114630685, 0.646704870719187, 11.746447134459)

#C <- c(0.0, 2.205692, 0.0005, 0.646704870719187, 11.746447134459)
#C <- c(4.83817673987683, 4.9939256206626, -0.394202677436942, 1.03685848528259, 15.3086092632073)

#sample_v <- rep(0.01 + ((0:19)/2), 10)
#sample_lambda <- 0.01 + floor((0:99)/10)/10

# min_sample_v <- 0
# max_sample_v <- 5
# min_sample_lambda <- 0
# max_sample_lambda <- 1
# sample_v <- runif(num_sample_points, min=min_sample_v, max=max_sample_v)
# sample_lambda <- runif(num_sample_points, min=min_sample_lambda, max=max_sample_lambda)

#C <- c(0.165811007520969, 2.37521276246873, -0.961571348975179, 0.537168547754793, 11.9336243457009)
C <- c(0.165811007520969, 2.37521276246873, 0.0005, 0.537168547754793, 11.9336243457009)

sample_v <- c(9.76853331550956, 5.06963664898649, 6.28607512684539, 3.3032006630674,
              4.61578730493784, 5.24314971175045, 4.66553797479719, 4.24834446515888,
              3.39309606002644, 8.50095360539854, 0.0757059385068715, 2.32428668998182,
              0.536202404182404, 8.20064031751826, 1.55000230297446, 7.87333103362471,
              9.72146387444809, 7.98630209406838, 0.803718918468803, 5.58928076177835,
              6.05682188645005, 1.01962292334065, 2.44109193095937, 9.63262678356841,
              5.22961625829339, 5.67657563369721, 4.82630877289921, 8.96400656551123,
              6.87015030998737, 3.12800221377984, 0.792783456854522, 3.85982782347128,
              7.47367175063118, 1.48879852611572, 4.68645944492891, 7.56992747541517,
              3.21333840023726, 7.485010556411, 6.59257436171174, 7.48357929522172,
              7.7149643166922, 1.47124944021925, 8.54282982181758, 5.6137865385972,
              2.16383191524073, 3.60345774563029, 0.903663274366409, 8.02030290709808,
              2.43502988480031, 5.80499117262661, 2.14938354212791, 7.38168602110818,
              4.23732332186773, 1.46784445736557, 1.1217476776801, 0.672697997651994,
              1.10718484967947, 1.00436876993626, 3.29712702427059, 9.14727190742269,
              7.59352876571938, 0.0011834385804832, 0.19021765794605, 1.82785857468843,
              6.89051311463118, 4.54415152780712, 4.12735636346042, 6.45460993051529,
              7.28200049605221, 1.82057200232521, 1.33146679727361, 7.4641451286152,
              4.2156188050285, 9.32264673290774, 8.01214536419138, 5.46254583867267,
              8.7639832100831, 0.101939544547349, 9.86967506585643, 8.40062062023208,
              2.49384714057669, 1.4533775462769, 4.07609009882435, 7.6220978749916,
              3.7358191004023, 4.42348951473832, 3.43858133303002, 9.49210463091731,
              9.66960808029398, 2.03581244684756, 1.06307015055791, 1.58924244809896,
              9.31405030656606, 4.67570035718381, 3.73736867913976, 9.36060729902238,
              1.23588900547475, 5.00741207506508, 9.50623262906447, 3.02969945827499)

sample_lambda <- c(0.976162793813273, 0.530264771543443, 0.494442440802231, 0.660209072055295,
                   0.860309258569032, 0.535462349187583, 0.58420978160575, 0.586046514799818,
                   0.397057597292587, 0.782532247947529, 0.749300858238712, 0.446695489110425,
                   0.117832656484097, 0.285324628232047, 0.313062278786674, 0.311176936840639,
                   0.0279256068170071, 0.266172747360542, 0.333466839278117, 0.40180263039656,
                   0.839342754334211, 0.400673164287582, 0.276359601179138, 0.682090745074674,
                   0.808462163666263, 0.937448780518025, 0.404605070594698, 0.063787450781092,
                   0.656649729935452, 0.829899839358404, 0.70440429658629, 0.877657957840711,
                   0.832580404356122, 0.0161103326827288, 0.870701983803883, 0.782410079147667,
                   0.0699446420185268, 0.766662252834067, 0.151806409005076, 0.862086791312322,
                   0.669131957460195, 0.981968907406554, 0.913708637934178, 0.651627411367372,
                   0.187409008154646, 0.467723964015022, 0.901095274602994, 0.0640165661461651,
                   0.271196133922786, 0.870892733568326, 0.0215295713860542, 0.7099663564004,
                   0.499739120714366, 0.97427370143123, 0.864840737776831, 0.0982439909130335,
                   0.127567816525698, 0.483359369914979, 0.484133140416816, 0.688911519944668,
                   0.982656674226746, 0.175129660405219, 0.0900900829583406, 0.195563859771937,
                   0.949094318319112, 0.291959362104535, 0.811371297342703, 0.172585870604962,
                   0.347023356240243, 0.636445988900959, 0.193085166392848, 0.905731934588403,
                   0.563686924986541, 0.546454882482067, 0.61829167464748, 0.123153926339,
                   0.0498951165936887, 0.566465184092522, 0.577949306229129, 0.224360200576484,
                   0.400887280935422, 0.420197926694527, 0.0909187830984592, 0.380277601303533,
                   0.977326963795349, 0.337373004760593, 0.358535267179832, 0.345443163067102,
                   0.18862827308476, 0.0152845662087202, 0.252264607232064, 0.706592309987172,
                   0.675292370375246, 0.596837569726631, 0.213073092745617, 0.907196244457737,
                   0.0941249371971935, 0.375390370609239, 0.927309417398646, 0.155363161116838)

UV_fn <- function(c, lambda) {
  #r <- 0.00183
  r <- 200
  val_ <- sigma/(sqrt(2*r))
  if (c/lambda <= val_) {
    val <- val_*log((c/lambda)/val_)
  } else {
    val <- c/lambda - val_
  }
  return(val)
}

sample_uv <- c()
for (i in 1:length(sample_lambda)) {
  sample_uv[i] <- UV_fn(C[3], sample_lambda[i])
}

#sample_uv <- rep(C[3], num_sample_points)

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

P_v <- function(v, mv, sv) {
  return(dnorm(v, mean=mv, sd=sv))
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_v(sample_v[i], C[1], C[2])*P_lambda(sample_lambda[i], C[4], C[5])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Heterogeneous in v and lambda with uv instead of c

num_sample_points <- 100

# C <- c(1.75988791967288, -0.369443114630685, 0.646704870719187, 11.746447134459 )
# #C <- c(2.60656001205525, -0.338722896520717, 0.930889818416773, 19.0252707054446)
#
# sample_v <- c(4.12999950814992, 3.69633289286867, 3.21354872663505, 0.652489480562508,
#               4.29492570343427, 2.12567484588362, 2.59106015088037, 2.16569548705593,
#               1.53945532278158, 3.59640921698883, 4.37932297121733, 4.87048892537132,
#               1.03409929201007, 4.48572893626988, 4.75324750877917, 0.106454635970294,
#               4.63013596949168, 3.50938543910161, 3.9536379755009, 1.80220969719812,
#               0.0911242549773306, 3.35393951623701, 0.133432896109298, 4.47179717943072,
#               0.1130381366238, 4.34785731835291, 4.19016028987244, 4.33123906841502,
#               3.50365528371185, 4.40776692004874, 3.34024658775888, 4.31588080362417,
#               2.64762821956538, 0.357169786002487, 4.07020782469772, 3.49954359116964,
#               3.90238487627357, 1.88185744336806, 3.18313005845994, 1.66961616487242,
#               3.03655192255974, 0.93425257364288, 1.9997228810098, 0.657894212054089,
#               4.48633232386783, 3.75805961317383, 0.119555647252128, 1.36801708838902,
#               1.14873605896719, 0.387883024523035, 4.26177503191866, 0.157281482825056,
#               1.07788230059668, 3.29599884338677, 3.80763098946773, 1.66512233088724,
#               0.283528257859871, 4.73436967120506, 4.79999975534156, 3.2095903204754,
#               0.524758036481217, 4.05302930041216, 4.05282010440715, 3.28715511714108,
#               1.87387446407229, 4.98448560247198, 2.03728850465268, 3.75573559780605,
#               1.36204213369638, 4.44790036650375, 2.15367441996932, 4.31037819595076,
#               4.70740746241063, 4.03891793452203, 4.87643286935054, 0.552339560817927,
#               3.35052622831427, 3.13704114407301, 3.40779118356295, 1.67753279907629,
#               4.63151719537564, 4.79276527068578, 1.62049875361845, 0.474379421211779,
#               1.37780031072907, 0.476049002027139, 4.63330883765593, 2.554002890829,
#               3.15550035098568, 2.08101065363735, 1.58368552918546, 4.73370558465831,
#               0.267330992501229, 0.967474151402712, 0.795003502862528, 1.34065822814591,
#               3.79445056780241, 0.964817715575919, 3.17544838413596, 0.309615317964926)
#
# sample_lambda <- c(0.850308072054759, 0.556064589880407, 0.85029289056547, 0.997881941264495,
#                    0.285017092712224, 0.390264749526978, 0.423349268967286, 0.674057501601055,
#                    0.66341686132364, 0.700843079481274, 0.862866312032565, 0.21182292024605,
#                    0.670170655474067, 0.565313470549881, 0.806613616878167, 0.184038332430646,
#                    0.0177535119000822, 0.0926143494434655, 0.106683721998706, 0.754633804550394,
#                    0.221237304154783, 0.673755551688373, 0.845521247014403, 0.342713839840144,
#                    0.381399230333045, 0.822355746757239, 0.168950211489573, 0.261941641569138,
#                    0.733049993403256, 0.936028465861455, 0.795731282094494, 0.204213462769985,
#                    0.490877247182652, 0.503583811689168, 0.933895880822092, 0.808955877553672,
#                    0.881635581608862, 0.770337234251201, 0.560144106857479, 0.122357967542484,
#                    0.483258098363876, 0.419826260069385, 0.708167464472353, 0.137961998581886,
#                    0.852174944709986, 0.616065685637295, 0.877141852863133, 0.476675770711154,
#                    0.194821767741814, 0.258901069639251, 0.157527908915654, 0.160536271519959,
#                    0.633640289306641, 0.740542439511046, 0.0726135135628283, 0.990809959359467,
#                    0.541834690840915, 0.483894780278206, 0.829073969041929, 0.98452244582586,
#                    0.686139798955992, 0.417905079433694, 0.587344815023243, 0.530868536792696,
#                    0.0612425173167139, 0.315859517082572, 0.208668164210394, 0.06165567249991,
#                    0.141203071456403, 0.844700614688918, 0.109951876336709, 0.930289818672463,
#                    0.652233351487666, 0.404041068861261, 0.484728938201442, 0.982576236594468,
#                    0.00854553654789925, 0.661161987809464, 0.808924428652972, 0.898620567750186,
#                    0.0568270718213171, 0.542690599570051, 0.312431619735435, 0.697093070950359,
#                    0.348495945334435, 0.724341917783022, 0.922267432091758, 0.273831122554839,
#                    0.166535726049915, 0.510053100064397, 0.0440130960196257, 0.0371885553468019,
#                    0.913370817666873, 0.715710175456479, 0.949244634946808, 0.00754711055196822,
#                    0.515732812229544, 0.582395692588761, 0.820996631169692, 0.544671004405245)

#C <- c(1.19287617394177, -2.22587325395563, 0.89804663567774, 15.9652968367141)
C <- c(1.43032349247673, -1.99211606996115, 1.01882262908409, 19.3346280875697)

sample_v <- c(1.63076030556113, 1.42081345897168, 2.68825458362699, 1.06559204054065,
              3.95796011551283, 1.00841492065229, 3.796194747556, 2.43219825904816,
              2.61813346529379, 3.02485328516923, 3.43551122932695, 1.55090690008365,
              0.838049890007824, 2.56427840562537, 3.49159787292592, 0.608711276436225,
              2.36602609278634, 0.330150100635365, 4.70250926911831, 3.92280973610468,
              1.28950280719437, 1.37559716473334, 3.43849578057416, 3.49009718396701,
              4.9638585082721, 1.31130766123533, 0.0549933465663344, 1.59406274091452,
              3.07059792452492, 0.346484758192673, 1.92892825347371, 3.30411816365086,
              0.961321883369237, 2.93252034112811, 3.36642194888555, 4.37974460772239,
              4.63446655659936, 3.46932601649314, 3.11514775152318, 2.32613084837794,
              3.86331608053297, 4.97527729021385, 0.0898273894563317, 1.45666957017966,
              4.74006540840492, 2.24720311933197, 1.72384787001647, 0.43220670777373,
              2.14526715455577, 0.0480428792070597, 1.78293628618121, 2.67457612324506,
              0.416301870718598, 4.60215834784321, 1.28381164278835, 1.49615809204988,
              2.05636704922654, 2.28984813787974, 3.11342773376964, 3.23056360823102,
              2.30288973194547, 4.56947694416158, 3.87036816799082, 1.40560614992864,
              0.276656143832952, 3.18975281785242, 4.21321097644977, 1.07850864296779,
              0.787914445390925, 3.51926943985745, 4.03823834960349, 2.28094829362817,
              2.83186048036441, 3.18338227691129, 1.22427805210464, 3.84543856955133,
              3.48418980021961, 1.1993874842301, 0.902973493793979, 4.72006158321165,
              1.46367806009948, 4.46121163433418, 0.483547660987824, 3.68800515425391,
              3.53580831433646, 0.716547989286482, 1.55134237371385, 2.44714239030145,
              3.77476818975993, 0.459822666598484, 2.10980731411837, 3.86190376477316,
              0.933758304454386, 2.20210900646634, 1.37067251023836, 3.28628741088323,
              1.10666309134103, 1.7896493722219, 3.82018655305728, 2.05353630939499)

sample_lambda <- c(0.184568778378889, 0.247649414232001, 0.479294002056122, 0.18742071907036,
                   0.626555798109621, 0.965259582269937, 0.125338526908308, 0.284472210565582,
                   0.129373387200758, 0.9712840388529, 0.120655856328085, 0.0504707696381956,
                   0.862425170140341, 0.364076743833721, 0.120689725968987, 0.171556766843423,
                   0.701648371992633, 0.302433780161664, 0.360838222783059, 0.138580456143245,
                   0.34695363859646, 0.316736877895892, 0.0515786726027727, 0.569239754462615,
                   0.0184761493001133, 0.108390801120549, 0.449769010301679, 0.398104301886633,
                   0.106803265400231, 0.488401491660625, 0.518087928183377, 0.00733225722797215,
                   0.295162350637838, 0.302074242150411, 0.688690076582134, 0.930093942210078,
                   0.190792606212199, 0.328962031751871, 0.304114262573421, 0.473321036202833,
                   0.883828061632812, 0.746598243713379, 0.981876422883943, 0.978300262242556,
                   0.0910913571715355, 0.139194839866832, 0.00164437829516828, 0.296292378567159,
                   0.383865190204233, 0.476781149627641, 0.0521963285282254, 0.0636385686229914,
                   0.445114843547344, 0.548564800526947, 0.026791081763804, 0.462110358290374,
                   0.495860942173749, 0.581489593023434, 0.334958199644461, 0.80626315344125,
                   0.769743368960917, 0.581611787201837, 0.973857208155096, 0.614755548071116,
                   0.0202703594695777, 0.508693655719981, 0.0241557527333498, 0.863478661281988,
                   0.237471646629274, 0.747444362379611, 0.375557094812393, 0.46460406598635,
                   0.199026279849932, 0.604309330927208, 0.521087244153023, 0.240768992342055,
                   0.402396495221183, 0.111759891500697, 0.773453035159037, 0.850902524311095,
                   0.229466974502429, 0.464666888816282, 0.121511974139139, 0.219643314601853,
                   0.83743141614832, 0.629812734667212, 0.935145837487653, 0.757250766269863,
                   0.258385443128645, 0.612292520701885, 0.642175579210743, 0.0444105847273022,
                   0.237502862466499, 0.446041788905859, 0.211729276925325, 0.809561612783,
                   0.345832022605464, 0.119239167543128, 0.838425514288247, 0.25203553121537)
# Result (good):
# dput(prediction_B)
# c(0.192560851095246, 0.663376545851364, 1.27671878980907, 1.98652770399018,
#   2.58611739864915, 3.15308032737654, 3.87530822184529, 6.16927385061886
# )
# dput(predict_rep_txn_B)
# c(1432.68509262086, 426.642750587504, 196.649543463404, 108.351146073951,
#   67.9241611405817, 44.5125509612939, 29.0074544837795, 51.2273006686226
# )

sample_uv <- rep(C[2], num_sample_points)

print("Pre-compute the homogeneous likelihood function")
L_B_df_ <- list()
for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))
  L_B_df_[[i]] <- L_B(sample_v[i], sample_uv[i], sample_lambda[i])
}
L_B_df <- as.data.frame(do.call(cbind, L_B_df_))
rm(L_B_df_)

P_v <- function(v, betav) {
  val <- exp(-v/betav)/betav
  return(val)
}

P_lambda <- function(lambda, r, alpha) {
  val <- (alpha^r)*(lambda^(r-1))*exp(-lambda*alpha)/gamma(r)
  return(val)
}

print("Pre-compute the probability vector")
P_vec <- rep(0, num_sample_points)
for (i in 1:num_sample_points) {
  P_vec[i] <- P_v(sample_v[i], C[1])*P_lambda(sample_lambda[i], C[3], C[4])
}

print("Pre-compute the heterogeneous likelihood function")
L_B_het_vec <- apply(sweep(L_B_df, MARGIN=2, P_vec, "*"), MARGIN=1, sum)/sum(P_vec)

print("Pre-compute the posterior probability vector")
P_pos_df <- sweep(L_B_df, MARGIN=2, P_vec, "*")/L_B_het_vec

################################################################################
# Presenting the result

# Make a prediction for the next 39 weeks (78 weeks in total). Compare with actual data.
# For a given number of repeat transactions the prediction will be averaged over all the occurrences of tx.
# Since we don't have a lot of data for 7+ repeat transactions in the first 39 weeks, all the entries with 7+ will be sent to the same bin.

prediction_B <- rep(0,8)

rep_txn_count <- rep(0,8)
predict_rep_txn_B <- rep(0,8)

EXt_B_het_vec <- rep(0, num_ids)

for (i in 1:num_sample_points) {
  print(paste("Computing point i = ", i))

  predict_rep_txn_B <- predict_rep_txn_B + PXt_B(sample_v[i], sample_uv[i], sample_lambda[i])*P_vec[i]
  EXt_B_het_vec <- EXt_B_het_vec + EXt_B(sample_v[i], sample_uv[i], sample_lambda[i])*P_pos_df[,i]
}
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

dput(prediction_B)
dput(predict_rep_txn_B)

################################################################################
# Save logfile:

file_name <- 'prediction_B_record.csv'
#file_name <- 'prediction_B_record_homo.csv'

write.csv(matrix(EXt_B_het_vec, nrow=1), file =file_name, row.names=FALSE)


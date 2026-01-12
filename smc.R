
library(hypergeo)
library(expint)

################################################################################
# Functions for reproducing SMC's model.

# See http://brucehardie.com/notes/008/pareto_nbd_MATLAB.pdf for the likelihood function and conditional expectation.

A0 <- function(r, alpha, s, beta, x, tx, TT) {
  val <- 0
  if (alpha >= beta) {
    F211 <- hypergeo::hypergeo(r+s+x,s+1,r+s+x+1,(alpha-beta)/(alpha+tx))
    F212 <- hypergeo::hypergeo(r+s+x,s+1,r+s+x+1,(alpha-beta)/(alpha+TT))
    val <- F211/((alpha+tx)^(r+s+x)) - F212/((alpha+TT)^(r+s+x))
  } else {
    F211 <- hypergeo::hypergeo(r+s+x,r+x,r+s+x+1,(beta-alpha)/(beta+tx))
    F212 <- hypergeo::hypergeo(r+s+x,r+x,r+s+x+1,(beta-alpha)/(beta+TT))
    val <- F211/((beta+tx)^(r+s+x)) - F212/((beta+TT)^(r+s+x))
  }
  val <- Re(val)
  return(val)
}

L_SMC <- function(r, alpha, s, beta, X, observe_tau=TRUE) {
  x <- X[1]-1
  TT <- X[length(X)]

  if (observe_tau) {
    tau <- X[length(X)-1]
    dead_at_TT <- 1*(tau <= TT)
    val <- (s^dead_at_TT * (alpha^r)*(beta^s)/( (alpha+min(tau,TT))^(r+x) * (beta+min(tau,TT))^(s + dead_at_TT) ))*gamma(r+x)/gamma(r)
  } else {
    tx <- X[x+2]
    val1 <- gamma(r+x)*(alpha^r)*(beta^s)/gamma(r)
    val2 <- 1/(((alpha+TT)^(r+x))*((beta+TT)^s))
    val3 <- s/(r+s+x)
    val <- val1*(val2 + val3*A0(r, alpha, s, beta, x, tx, TT))
  }
  return(val)
}

# In SMC, we need a version of n choose k when n and k are not necessary integers. So we will use gamma function.
choose_SMC <- function(n,k) {
  return(gamma(n+1)/(gamma(k+1)*gamma(n-k+1)))
}

V3 <- function(t, x, r, alpha, s, beta) {
  val <- 0
  if (alpha > beta) {
    for (j in 0:x) {
      F211 <- hypergeo::hypergeo(s+1,r+j+s,r+j+s+1,(alpha-beta)/alpha)
      F212 <- hypergeo::hypergeo(s+1,r+j+s,r+j+s+1,(alpha-beta)/(alpha+t))
      val <- val + choose_SMC(x,j)*((-alpha)^j)*((alpha^(-r-j-s))*F211/(r+j+s) - ((alpha+t)^(-r-j-s))*F212/(r+j+s))
    }
  } else if (alpha < beta) {
    for (j in 0:x) {
      F211 <- hypergeo::hypergeo(r+x,r+j+s,r+j+s+1,(beta-alpha)/beta)
      F212 <- hypergeo::hypergeo(r+x,r+j+s,r+j+s+1,(beta-alpha)/(beta+t))
      val <- val + choose_SMC(x,j)*((-beta)^j)*((beta^(-r-j-s))*F211/(r+j+s) - ((beta+t)^(-r-j-s))*F212/(r+j+s))
    }
  } else {
    for (j in 0:x) {
      val <- val + choose_SMC(x,j)*((-alpha)^j)*((alpha^(-r-j-s))/(r+j+s) - ((alpha+t)^(-r-j-s))/(r+j+s))
    }
  }
  return(val)
}

# (A34) of SMC with V3 given by (A39), (A42), and (A44) of SMC
PXt_SMC <- function(t, x, r, alpha, s, beta) {
  val1 <- choose_SMC(x+r-1,x)*((alpha/(alpha+t))^r)*((t/(alpha+t))^x)*((beta/(beta+t))^s)
  val2 <- choose_SMC(x+r-1,x)*s*(alpha^r)*(beta^s)*V3(t,x,r,alpha,s,beta)
  val <- Re(val1 + val2)
  return(val)
}

# (17) of SMC
EXt_uncond_het_SMC <- function(TT, r,alpha,s,beta) {
  val <- ((r*beta)/(alpha*(s-1)))*(1 - (beta/(beta+TT))^(s-1))
  return(val)
}

# See section 6 of http://brucehardie.com/notes/008/pareto_nbd_MATLAB.pdf.
EXt_SMC <- function(t, X, r, alpha, s, beta, observe_tau=TRUE) {
  x <- X[1]-1
  TT <- X[length(X)]

  EXstar <- EXt_uncond_het_SMC(t, r+x, alpha+TT, s, beta+TT)

  if (observe_tau) {
    tau <- X[length(X)-1]
    P <- 1*(tau > TT)
  } else {
    tx <- X[x+2]
    # See (3) of http://brucehardie.com/notes/008/pareto_nbd_MATLAB.pdf for 'active' probability.
    P <- 1/(1 + (s/(r+x+s))*((alpha+TT)^(r+x))*((beta+TT)^s)*A0(r,alpha,s,beta,x,tx,TT))
  }
  val <- EXstar*P
  return(val)
}

################################################################################

P_mu <- function(mu, X, r, alpha, s, beta) {
  x <- X[1]-1
  TT <- X[length(X)]
  tx <- X[x+2]

  val1 <- ((alpha^r) * (beta^s))*(mu^(s-1))*exp(-mu*(TT+beta)) / (gamma(r)*gamma(s))
  val2 <- (mu^(r+x))*exp(mu*(TT+alpha))*gamma(r+x)*((r+x)*gammainc(-r-x,(TT+alpha)*mu) + gammainc(1-r-x,(tx+alpha)*mu))
  val3 <- L_SMC(r,alpha,s,beta,X,observe_tau=FALSE)

  val <- val1*val2/val3

  return(val)
}

CDNOW_data_augmentation <- function(C_SMC, dt=1/7, tol_period=78) {
  Data <- read.table(paste(getwd(),'/CDNOW_sample.txt',sep=''), header=FALSE, sep="")
  names(Data) <- c("master_id", "id", "day", "amount", "dollar_amount")
  Data$day <- as.Date(as.character(Data$day), "%Y%m%d")
  Data$day <- as.numeric(Data$day - min(Data$day))/7

  # Get the number of consumers.
  num_ids <- length(unique(Data[,2]))

  # In the following we prepare the data by separating them into training (calibration period) and testing set.
  Data_ <- data.frame(matrix(ncol = ncol(Data), nrow = 0))
  names(Data_) <- names(Data)

  print("Synthesizing data from SMC model, please wait...")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = num_ids, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "█")   # Character used to create the bar

  all_ids_Data_ <- list()

  mu_list <- seq(from=0.001,to=2,length.out=1000)
  for(i in 1:num_ids){
    id_Data <- Data[Data$id == i,]

    x <- length(id_Data$day)
    tt <- id_Data$day
    tx <- tail(id_Data$day, 1)
    TT <- tol_period
    X <- c(x, tt, TT)

    r <- runif(1)
    P_mu_vec <- P_mu(mu_list, X, C_SMC[1], C_SMC[2], C_SMC[3], 0.1*C_SMC[4])
    P_mu_vec <- P_mu_vec/sum(P_mu_vec)
    j <- 1
    accum <- 0
    while((accum < r) && (j < length(mu_list))) {
      accum <- accum + P_mu_vec[j]
      j <- j + 1
    }
    mu <- mu_list[j]

    t <- tx
    alive <- TRUE
    while((t < TT) && (alive)) {
      t <- t + dt
      r <- runif(1)
      alive <- (r < exp(-mu*dt))
    }

    id_Data[nrow(id_Data)+1, ] <- c(id_Data$master_id[1], i, t, 0, 0)
    all_ids_Data_[[i]] <- id_Data

    setTxtProgressBar(pb, i)
  }
  Data_ <- dplyr::bind_rows(all_ids_Data_)
  row.names(Data_) <- NULL

  write.table(Data_, file="CDNOW_membership_aug.txt", row.names=FALSE, sep=",")
}

C_SMC_ <- c(0.553, 10.580, 0.606, 11.656)
#CDNOW_data_augmentation(C_SMC=C_SMC_)

################################################################################

Data <- read.table(paste(getwd(),'/CDNOW_membership_aug_unif.txt',sep=''), header=TRUE, sep=",")

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


training_Data_alive <- list()
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

################################################################################
# SMC's MLE

observe_tau <- TRUE

count <- 0
nLL_SMC <- function(r,alpha,s,beta) {
  LL <- 0
  for(i in 1:length(training_Data)) {
    # Ignore people who make too many (> 8) repeated transactions, they causes overflow in calculation. There are a few of them anyway.
    #if (nrow(id_training_Data[[i]]) <= 8) {
    LL <- LL + log(L_SMC(r,alpha,s,beta, training_Data[[i]], observe_tau))
    #print(paste("LL = ", as.character(LL), ", i = ", as.character(i)))
    #}
  }
  if (count%%10==0){
    print(paste("LL = ", as.character(LL), ", theta = ", paste(c(r,alpha,s,beta), collapse=", ")))
  }
  count <- count + 1
  return(-LL)
}

fit_SMC <- stats4::mle(nLL_SMC, start=list(r=0.553,alpha=10.580,s=0.606,beta=11.656), lower=c(0.01,0.01,0.01,0.01), upper=c(10,100,50,1000))
C_SMC <- as.numeric(fit_SMC@coef)

################################################################################

observe_tau <- FALSE
Input_Data <- training_Data
Output_GT <- test_period_txns
actual <- actual_test_txns

num_ids <- length(Output_GT)

prediction_SMC <- rep(0,8)
rep_txn_count <- rep(0,8)

rep_txn_vec <- rep(0, num_ids)
EXt_SMC_vec <- rep(0, num_ids)
for (i in 1:num_ids) {
  rep_txn_vec[i] <- Input_Data[[i]][1]-1
  EXt_SMC_vec[i] <- EXt_SMC(test_period, Input_Data[[i]], C_SMC[1], C_SMC[2], C_SMC[3], C_SMC[4], observe_tau)
}

for (i in 1:num_ids) {
  rep_txn <- rep_txn_vec[i]
  if (rep_txn > 7) {
    rep_txn <- 7
  }
  rep_txn_count[rep_txn+1] <- rep_txn_count[rep_txn+1] + 1
  prediction_SMC[rep_txn+1] <- prediction_SMC[rep_txn+1] + EXt_SMC_vec[i]
}

for (i in 1:length(prediction_SMC)) {
  prediction_SMC[i] <- prediction_SMC[i]/rep_txn_count[i]
}

MAE_SMC <- sum(abs(Output_GT - EXt_SMC_vec))/num_ids

print(MAE_SMC)

plot(prediction_SMC, type='o',ylim=c(0,8))
points(actual, type='o', col='blue')

################################################################################

C_SMC_ <- c(0.553, 10.580, 0.606, 11.656)

label_vec = 0:7

actual <- c(0.2367116, 0.6970387, 1.392523, 1.560000, 2.532258, 2.947368, 3.862069, 6.359375)
prediction_SMC <- c(0.138, 0.600, 1.196, 1.714, 2.399, 2.907, 3.819, 6.403)
prediction_PDE <- c(0.1818177, 0.7657967, 1.3354794, 1.7826441, 2.4946101, 2.9251404, 3.8692302, 5.8741765)

plot(label_vec, type='o', prediction_PDE, col='brown', ylim=c(0,8), xlab="# Repeat Transactions in Weeks 1-39", ylab="Expected # Transactions in Weeks 40-78", xaxt='n', lwd=2, lty=1, pch=2,cex.lab=1.3, cex.axis=1.3)
points(label_vec, type='o', prediction_SMC, col='blue', ylim=c(0,8), lwd=2, lty=5,pch=0)
points(label_vec, type='o', actual, col='black', ylim=c(0,8), lwd=2, lty=2,pch=1)
axis(1, at=label_vec, labels=c('0', '1', '2','3','4','5','6','7+'), cex.axis=1)
op <- par(cex=1)
legend(x="topleft", legend=c('Actual', 'Pareto/NBD', 'Our Model'), col=c('black', 'blue', 'brown'), lwd=c(2,2,2), lty=c(2,6,1), pch=c(1, 0, 2))

# AUG_UNIF
actual_contractual <- c(1.730570, 1.366071, 1.922581, 2.136986, 2.854545, 3.393939, 4.480000, 6.672131)
prediction_SMC_contractual <- c(0.5451542, 1.1328202, 1.7237065, 2.3562227, 2.9050291, 3.4379793, 4.0547169, 6.7643780)
prediction_PDE_contractual <- c(0.7158068, 1.1811701, 1.7909049, 2.3289669, 3.0875035, 3.6764880, 4.4919265, 6.5889622)

plot(label_vec, type='o', prediction_PDE_contractual, col='brown', ylim=c(0,8), xlab="# Repeat Transactions in Weeks 1-39", ylab="Expected # Transactions in Weeks 40-78", xaxt='n', lwd=2, lty=1, pch=2,cex.lab=1.3, cex.axis=1.3)
points(label_vec, type='o', prediction_SMC_contractual, col='blue', ylim=c(0,8), lwd=2, lty=6,pch=0)
points(label_vec, type='o', actual_contractual, col='black', ylim=c(0,8), lwd=2, lty=2,pch=1)
axis(1, at=label_vec, labels=c('0', '1', '2','3','4','5','6','7+'), cex.axis=1)
op <- par(cex=1)
legend(x="topleft", legend=c('Actual', 'Pareto/NBD', 'Our Model'), col=c('black', 'blue', 'brown'), lwd=c(2,2,2), lty=c(2,6,1), pch=c(1, 0, 2))

# AUG_SMC
actual_contractual <- c(0.6162362, 0.9683544, 1.6108108, 1.9746835, 2.7543860, 3.1111111, 4.1481481, 6.5645161)
prediction_SMC_contractual <- c(0.3641299, 1.0143559, 1.6661777, 2.3270000, 2.9772387, 3.5935483, 4.2299678, 7.2251521)
prediction_PDE_contractual <- c(0.7417928, 1.1149617, 1.5782882, 1.9792096, 2.4988476, 2.9108850, 3.4889297, 5.1015773)

plot(label_vec, type='o', prediction_PDE_contractual, col='brown', ylim=c(0,8), xlab="# Repeat Transactions in Weeks 1-39", ylab="Expected # Transactions in Weeks 40-78", xaxt='n', lwd=2, lty=1, pch=2,cex.lab=1.3, cex.axis=1.3)
points(label_vec, type='o', prediction_SMC_contractual, col='blue', ylim=c(0,8), lwd=2, lty=6,pch=0)
points(label_vec, type='o', actual_contractual, col='black', ylim=c(0,8), lwd=2, lty=2,pch=1)
axis(1, at=label_vec, labels=c('0', '1', '2','3','4','5','6','7+'), cex.axis=1)
op <- par(cex=1)
legend(x="topleft", legend=c('Actual', 'Pareto/NBD', 'Our Model'), col=c('black', 'blue', 'brown'), lwd=c(2,2,2), lty=c(2,6,1), pch=c(1, 0, 2))

################################################################################

# dt <- 1/7
# time_period <- train_period + test_period
# num_steps_d <- as.integer(time_period/dt)
# num_customers_incr_d <- rep(0, num_steps_d)
#
# num_ids <- length(training_Data)
# daily_txns_d <- rep(0, num_steps_d)
# for (i in 1:num_ids) {
#   id_txns <- as.integer(unique(Data[Data$id == i,]$day)/dt)
#   start_day <- id_txns[1]
#   num_customers_incr_d[start_day+1] <- num_customers_incr_d[start_day+1] + 1
#
#   for (j in 1:(length(id_txns)-1)) {
#     t <- id_txns[j]
#     if (t != start_day) {
#       daily_txns_d[t+1] <- daily_txns_d[t+1] + 1
#     }
#   }
# }
#
# predict_txns_SMC_d <- c()
# for (i in 1:num_steps_d) {
#   predict_txns_SMC_d[i] <- EXt_uncond_het_SMC(i*dt, C_SMC[1], C_SMC[2], C_SMC[3], C_SMC[4]) - EXt_uncond_het_SMC((i-1)*dt, C_SMC[1], C_SMC[2], C_SMC[3], C_SMC[4])
# }
#
# tol_predict_txns_SMC_d <- rep(0, num_steps_d)
# for (i in 1:num_steps_d) {
#   tol_predict_txns_SMC_d <- tol_predict_txns_SMC_d + num_customers_incr_d[i]*c(rep(0,i-1), predict_txns_SMC_d[i:num_steps_d - i +1])
# }
#
# num_steps <- time_period
#
# daily_txns <- c()
#
# tol_predict_txns_SMC <- c()
#
# for (i in 1:num_steps) {
#   daily_txns[i] <- sum(daily_txns_d[as.integer((i-1)/dt + 1):as.integer(i/dt)])
#   tol_predict_txns_SMC[i] <- sum(tol_predict_txns_SMC_d[as.integer((i-1)/dt + 1):as.integer(i/dt)])
# }
#
# label_vec <- 1:78
#
# plot(label_vec, type='l', daily_txns, col='black', ylim=c(0,130), xlab="Week", ylab="Transactions", xaxt='n', cex.lab=1.3, cex.axis=1.3, lwd=2, lty=2)
# points(label_vec, type='l', tol_predict_txns_SMC, col='red', lwd=2, lty=5)
# #points(label_vec, type='l', tol_predict_txns_B, col='brown', lwd=2)
# axis(1, at=1+7*(0:11), labels=as.character(1+7*(0:11)), cex.axis=1.3)
# op <- par(cex=1.3)
# legend(x="topright", legend=c('Actual', 'Pareto/NBD'), col=c('black', 'red'), lwd=c(2,2), lty=c(2,5))

# dmu <- 1e-13
# val <- 0
# mu <- dmu
# X <- training_Data[[1]]
# for (ii in 1:30000) {
#   val <- val + P_mu(mu, X,C_SMC_[1],C_SMC_[2], C_SMC_[3], C_SMC_[4])*dmu
#   mu <- mu + dmu
# }
# val <- val/L_SMC(C_SMC_[1],C_SMC_[2], C_SMC_[3], C_SMC_[4], X, TRUE)

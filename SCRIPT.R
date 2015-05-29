# crr.vs.crprep.1: compare time needed to estimate a Fine & Gray model using "crr" function vs the time needed using "crprep" + "coxph"
# each procedure runs M times, and results are averaged

crr_vs_crprep_1 <- function(n, M = 1000, seed = 1234){
   require(mstate)
   require(cmprsk)

   #seed for reproducibility
   set.seed <- seed

   results_matrix <- matrix(0, nrow = length(n), ncol = 4)
   colnames(results_matrix) = c("n", "CRR.ELAP", "CRPREP.ELAP", "ratio")

   for(j in 1:length(n)){
      results_matrix[j, 1] <- n[j]
      # generate data
      nn <- n[j]
      ev1times <- rexp(nn, rate = 0.1)
      ev2times <- rexp(nn, rate = 0.2)
      censtimes <- runif(n = nn, min = 0, max = .5 * max(ev1times, ev2times))
      time <- pmin(ev1times, ev2times, censtimes)
      status <- rep(0, nn)
      status[(ev1times < censtimes) & (ev1times < ev2times)] <- 1
      status[(ev2times < censtimes) & (ev2times < ev1times)] <- 2
      trt <- sample(c("Y", "N"), size = nn, replace = TRUE)
      dataset <- data.frame(time, status, trt)

      #crr
      crr_elap_tmp = rep(0, M)
      for(i in 1:M){
         tstart <- proc.time()
         fit_tmp <- with(dataset, crr(ftime = time, fstatus = status, cov1 = c(trt), failcode = 2, cencode = 0))
         tstop <- proc.time()
         telapsed <- tstop - tstart
         crr_elap_tmp[i] <- telapsed[3]
      }
      results_matrix[j, 2] <- mean(crr_elap_tmp)

      #crprep
      crprep_elap_tmp <- rep(0,M)
      for(i in 1:M){
         tstart <- proc.time()
         dataset_wt <- with(dataset, crprep(Tstop = time, status = status, data = dataset, trans = c(1, 2), cens = 0, keep = dataset))
         fit_tmp <- with(dataset_wt, coxph(Surv(Tstart, Tstop, status == 1) ~ trt, weigh = weight.cens, subset = (failcode == 1)))
         tstop <- proc.time()
         telapsed <- tstop - tstart
         crprep_elap_tmp[i] <- telapsed[3]
      }
      results_matrix[j, 3] <- mean(crprep_elap_tmp)

   }
   results_matrix[, 4] <- results_matrix[, 3] / results_matrix[, 2]
   return(results_matrix)
}

# crr.vs.crprep.2: compare time needed to estimate M Fine & Gray models using "crr" function vs the time needed using "crprep" once + "coxph" M times
# total elapsed time is returned

crr_vs_crprep_2 <- function(n, M, seed = 1234){
   require(mstate)
   require(cmprsk)

   #seed for reproducibility
   set.seed <- seed

   results_matrix <- matrix(0, nrow = length(n), ncol = 4)
   colnames(results_matrix) = c("n", "CRR.ELAP", "CRPREP.ELAP", "ratio")

   for(j in 1:length(n)){
      results_matrix[j, 1] <- n[j]
      # generate data
      nn <- n[j]
      ev1times <- rexp(nn, rate = 0.1)
      ev2times <- rexp(nn, rate = 0.2)
      censtimes <- runif(n = nn, min = 0, max = .5 * max(ev1times, ev2times))
      time <- pmin(ev1times, ev2times, censtimes)
      status <- rep(0, nn)
      status[(ev1times < censtimes) & (ev1times < ev2times)] <- 1
      status[(ev2times < censtimes) & (ev2times < ev1times)] <- 2
      trt <- sample(c("Y", "N"), size = nn, replace = TRUE)
      dataset <- data_frame(time, status, trt)

      #crr
      crr_elap_tmp <- rep(0, M)
      for(i in 1:M){
         tstart <- proc.time()
         fit_tmp <- with(dataset, crr(ftime = time, fstatus = status, cov1 = c(trt), failcode = 2, cencode = 0))
         tstop <- proc.time()
         telapsed <- tstop - tstart
         crr_elap_tmp[i] <- telapsed[3]
      }
      results_matrix[j, 2] <- sum(crr_elap_tmp)

      #crprep
      crprep_elap_tmp <- rep(0, M)
      tstart <- proc.time()
      dataset_wt <- with(dataset, crprep(Tstop = time, status = status, data = dataset, trans = c(1, 2), cens = 0, keep = dataset))
      tstop <- proc.time()
      trw <- tstop - tstart
      for(i in 1:M){
         tstart <- proc.time()
         fit_tmp <- with(dataset_wt, coxph(Surv(Tstart, Tstop, status == 1) ~ trt, weigh = weight.cens, subset = (failcode == 1)))
         tstop <- proc.time()
         telapsed <- tstop - tstart
         crprep_elap_tmp[i] <- telapsed[3]
      }
      results_matrix[j, 3] <- sum(crprep_elap_tmp) + trw[3]

   }
   results_matrix[, 4] <- results_matrix[, 3] / results_matrix[, 2]
   return(results_matrix)
}
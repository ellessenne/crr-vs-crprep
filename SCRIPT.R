# gen_cr_data: function to generate a competing risks dataset, with customizable sample size and just one competing event for simplicity
gen_cr_data <- function(n = 1000, seed = 1234){
  set.seed(seed)
  require(dplyr)
  data <- data_frame(time0 = rexp(n, 0.1),
                     time1 = rexp(n, 0.5),
                     time2 = rexp(n, 1)) %>%
    mutate(time = round(1 + 100 * pmin(time0, time1, time2)),
           status = ifelse(time1 == pmin(time0, time1, time2),
                           1,
                           ifelse(time2 == pmin(time0, time1, time2),
                                  2,
                                  0)),
           trt = sample(0:1, n, replace = TRUE),
           age = abs(rnorm(n, 50, 10))) %>%
    select(-time0, -time1, -time2)
  return(tbl_df(data))
}

gen_cr_data()

# crr_vs_crprep_1: compare time needed to estimate a Fine & Gray model using "crr" function vs the time needed using "crprep" + "coxph"

crr_vs_crprep1 <- function(ssize, B = 1000){
  if (length(ssize) != 1) stop("ssize cannot be a vector.")
  require(pacman)
  p_load("cmprsk", "mstate", "microbenchmark")
  d <- gen_cr_data(round(ssize))
  mb <- microbenchmark(
    "CRR" = with(d, crr(ftime = time, fstatus = status, cov1 = cbind(age, trt), tf = function(t) t, failcode = 1, cencode = 0)),
    "CRPREP + COXPH" = with(with(d, crprep(Tstop = time, status = status, data = d, trans = 1:2, cens = 0, keep = cbind(age, trt))),
                            coxph(Surv(Tstart, Tstop, status == 1) ~ age + trt, weights = weight.cens, subset = (failcode == 1))),
    times = B
  )
  return(mb)
}

# crr_vs_crprep2: compare time needed to estimate M Fine & Gray models using "crr" function vs the time needed using "crprep" once + "coxph" M times

# regular
f1 <- function(d, M){
  with(d, lapply(1:M, function(x) crr(ftime = time, fstatus = status, cov1 = cbind(age, trt), tf = function(t) t, failcode = 1, cencode = 0)))
}
# weighted
f2 <- function(d, M){
  d_wt <- with(d, crprep(Tstop = time, status = status, data = d, trans = 1:2, cens = 0, keep = cbind(age, trt)))
  with(d_wt, lapply(1:M, function(x) coxph(Surv(Tstart, Tstop, status == 1) ~ age + trt, weights = weight.cens, subset = (failcode == 1))))
}

crr_vs_crprep2 <- function(M, ssize = 100, B = 1000){
  require(pacman)
  p_load("cmprsk", "mstate", "microbenchmark", "parallel")
  d <- gen_cr_data(round(ssize))
  mb <- microbenchmark(
    "CRR" = f1(d, M),
    "CRPREP + COXPH" = f2(d, M),
    times = B
  )
  return(mb)
}
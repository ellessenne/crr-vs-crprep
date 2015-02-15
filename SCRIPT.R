# crr.vs.crprep.1: compare time needed to estimate a Fine & Gray model using "crr" function vs the time needed using "crprep" + "coxph"
# each procedure runs M times, and results are averaged

crr.vs.crprep.1 = function(n, M=1000, seed=1234){
   require(mstate)
   require(cmprsk)
   
   #seed for reproducibility
   set.seed=seed
   
   results.matrix = matrix(0, nrow=length(n), ncol=4)
   colnames(results.matrix)=c("n", "CRR.ELAP", "CRPREP.ELAP", "ratio")
   
   for(j in 1:length(n)){
      results.matrix[j,1] = n[j]
      # generate data
      nn=n[j]
      ev1times = rexp(nn, rate = 0.1)
      ev2times = rexp(nn, rate = 0.2)
      censtimes = runif(n=nn, min=0, max=.5*max(ev1times, ev2times))
      time = pmin(ev1times, ev2times, censtimes)
      status = rep(0, nn)
      status[(ev1times<censtimes)&(ev1times<ev2times)] = 1
      status[(ev2times<censtimes)&(ev2times<ev1times)] = 2
      trt = sample(c("Y","N"),size=nn, replace=TRUE)
      dataset = data.frame(time, status, trt)
      
      #crr
      crr.elap.tmp = rep(0,M)
      for(i in 1:M){
         tstart = proc.time()
         fit.tmp = with(dataset, crr(ftime=time, fstatus=status, cov1=c(trt), failcode=2, cencode=0))
         tstop = proc.time()
         telapsed = tstop - tstart
         crr.elap.tmp[i] = telapsed[3]
      }
      results.matrix[j,2] = mean(crr.elap.tmp)
      
      #crprep
      crprep.elap.tmp = rep(0,M)
      for(i in 1:M){
         tstart = proc.time()
         dataset.wt = with(dataset, crprep(Tstop=time, status=status, data=dataset, trans = c(1,2), cens = 0, keep=dataset))
         fit.tmp = with(dataset.wt, coxph(Surv(Tstart, Tstop, status==1)~trt, weigh=weight.cens, subset=failcode==1))
         tstop = proc.time()
         telapsed = tstop - tstart
         crprep.elap.tmp[i] = telapsed[3]
      }
      results.matrix[j,3] = mean(crprep.elap.tmp)
      
   }
   results.matrix[,4] = results.matrix[,3]/results.matrix[,2]
   return(results.matrix)
}

# crr.vs.crprep.2: compare time needed to estimate M Fine & Gray models using "crr" function vs the time needed using "crprep" once + "coxph" M times
# total elapsed time is returned

crr.vs.crprep.2 = function(n, M, seed=1234){
   require(mstate)
   require(cmprsk)
   
   #seed for reproducibility
   set.seed=seed
   
   results.matrix = matrix(0, nrow=length(n), ncol=4)
   colnames(results.matrix)=c("n","CRR.ELAP", "CRPREP.ELAP", "ratio")
   
   for(j in 1:length(n)){
      results.matrix[j,1] = n[j]
      # generate data
      nn=n[j]
      ev1times = rexp(nn, rate = 0.1)
      ev2times = rexp(nn, rate = 0.2)
      censtimes = runif(n=nn, min=0, max=.5*max(ev1times, ev2times))
      time = pmin(ev1times, ev2times, censtimes)
      status = rep(0, nn)
      status[(ev1times<censtimes)&(ev1times<ev2times)] = 1
      status[(ev2times<censtimes)&(ev2times<ev1times)] = 2
      trt = sample(c("Y","N"),size=nn, replace=TRUE)
      dataset = data.frame(time, status, trt)
      
      #crr
      crr.elap.tmp = rep(0,M)
      for(i in 1:M){
         tstart = proc.time()
         fit.tmp = with(dataset, crr(ftime=time, fstatus=status, cov1=c(trt), failcode=2, cencode=0))
         tstop = proc.time()
         telapsed = tstop - tstart
         crr.elap.tmp[i] = telapsed[3]
      }
      results.matrix[j,2] = sum(crr.elap.tmp)
      
      #crprep
      crprep.elap.tmp = rep(0,M)
      tstart=proc.time()
      dataset.wt = with(dataset, crprep(Tstop=time, status=status, data=dataset, trans = c(1,2), cens = 0, keep=dataset))
      tstop=proc.time()
      trw = tstop - tstart
      for(i in 1:M){
         tstart = proc.time()
         fit.tmp = with(dataset.wt, coxph(Surv(Tstart, Tstop, status==1)~trt, weigh=weight.cens, subset=failcode==1))
         tstop = proc.time()
         telapsed = tstop - tstart
         crprep.elap.tmp[i] = telapsed[3]
      }
      results.matrix[j,3] = sum(crprep.elap.tmp) + trw[3]
      
   }
   results.matrix[,4] = results.matrix[,3]/results.matrix[,2]
   return(results.matrix)
}
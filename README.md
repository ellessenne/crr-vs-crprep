# Comparison of crr and crprep + coxph performances
This code has been used to compare the perfomance of `R` code for comparing the performance of `crr` function (`cmprsk` package) and `crprep` + `coxph` (`mstate` and `survival` packages respectively).
Both procedures can be used to estimate a Fine & Gray model for the subdistribution hazard in a competing risks framework.

Recall that the time required to run these procedures is machine-dependent.

***

## Basic documentation

### crr.vs.crprep.1()

#### Usage
`crr.vs.crprep.1(n, M=1000, seed=1234)`

#### Arguments
* `n` vector with sample sizes
* `M` number of times to repeat both procedures. Default value is 1000
* `seed` for reproducibility. Default value is 1234

#### Examples
```
sample.sizes = seq(100, 1000, by=100)
crr.vs.crprep.1(n = sample.sizes)
```

***

### crr.vs.crprep.2()

#### Usage
`crr.vs.crprep.2(n, M, seed=1234)`

#### Arguments 
* `n` vector with sample sizes
* `M` number of models to estimate
* `seed` for reproducibility. Default value is 1234

#### Examples
```
# assuming the need to estimate 10 models:
sample.sizes = seq(100, 1000, by=100)
crr.vs.crprep.2(n = sample.sizes, M=10)
```

```
# create a list with different number of models required:
results.list = list()
sample.sizes = seq(100, 1000, by=100)
required.models = c(10, 50, 100)
for(i in 1:length(required.models)){
  results.list[[i]] = crr.vs.crprep.2(n = sample.sizes, M=required.models[i])
  names(results.list)[i] = paste(i, "models")
}
results.list
```

***
# FloodFreqPlot
## A package to plot flood quantiles and their probabilities
FloodFreqPlot is a package including some functions to plot flood quantiles and their probabilities on
    the probability papers. The package is useful for the flood frequency analysis studies.

Here is a simple example:

```
# Probability Plotting for the floods obseved at Harricana River at Amos (Quebec, Canada)
## Loading Harricana dataset
data(Harricana)
## Flood Probability Plotting
ProbPlot(data_obs = Harricana, PP = 'Cunnane', dist = 'LPea3', T_rp = c(100, 1000))
```
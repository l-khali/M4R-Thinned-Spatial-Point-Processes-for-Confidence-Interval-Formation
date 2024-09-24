# M4R: Thinned Spatial Point Processes for Confidence Interval Formation

## Report can be found [here](https://github.com/l-khali/M4R-Thinned-Spatial-Point-Processes-for-Confidence-Interval-Formation/blob/main/01717289_LK.pdf)

This page contains the code used to carry out experiments and generate plots for my MSci Mathematics research project. Code is organised into folders corresponding to chapters of the report.

Throughout this project, we present a novel resampling method, which we call thinning. We use this method to form confidence intervals for statistics estimated from a set of observations, and focus on spatial point processes. A thinned sample is a sample where each point in the observed sample is retained with a predetermined probability, $p$, and otherwise discarded. We prove that for some desired statistic, $\theta$, with estimate $\hat{\theta}$ and thinned replication $\hat{\theta}_s$, that the normalised distribution of $\hat{\theta}_s$ approaches that of $\hat{\theta}$ assuming the underlying data is i.i.d and $\theta$ is a linear statistical functional. This means that with a significance level of $\alpha$, the studentized confidence interval contains the true value about $(1-\alpha)$ percent of the time. We demonstrate this result numerically in the setting of probability distributions. We see that the theory extends to the spatial setting, in the case that the statistic is a homogeneous or inhomogeneous intensity function of a Poisson point process. However, the desired result does not hold when estimating Ripley's $K$, which is a measure of the amount of clustering in a process. We conclude that similarly to the jackknife and bootstrap, thinning does not adequately preserve the interactions between dependent points, and therefore is not useful for all spatial statistics, although it does provide accurate results for some.

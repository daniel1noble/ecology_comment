---
  title: "Tutorial on how to implement simple corrections to deal with non-independent effect sizes in multi-level meta-analysis"
  author: Shinichi Nakagawa, Alistair M. Senior, Wolfgang Viechtbauer and Daniel W. A. Noble
  date: "2021-03-01"
  bibliography: refs.bib
  csl: ecology.csl
  output: 
    bookdown::html_document2:
      code_folding: hide
      number_sections: no
      toc: yes
      toc_depth: 6
      toc_float: yes
---



<!--html_preserve--><script>
  addClassKlippyTo("pre.r, pre.markdown");
  addKlippy('right', 'top', '#FF0000', '1', 'Click to Copy Code', 'Done');
</script><!--/html_preserve-->

# Introduction
In this tutorial, we demonstrate how meta-analyst's can implement approaches for correcting inflated Type I error rates when using multilevel meta-analytic (MLMA) models with non-independent effect size data. For ease of presentation, we use simulated data to demonstrate the implementation of a few easy solutions. 

# R Packages Required

First, we'll load some of the packages that we'll need.

```{.r .klippy}
# Clean workspace
rm(list = ls())

# Loading packages & Functions

pacman::p_load(tidyverse, MASS, kableExtra, gridExtra, MCMCglmm, brms, metafor, robumeta, clubSandwich, 
    pander, tidyverse)
```

# Simulating Non-independent Effect Size Data

Here, we will simulate some meta-analytic data. We will keep this very simple, just for demonstration purposes. Hence, we will assume that we have collected data from a total of 20 studies, and we'll assume that we were able to extract n = 3 effect sizes from each of these 20 studies. In total, we have a data set that contains n = 60 effect sizes.


```{.r .klippy}
# Simulate a dataset composed of 20 papers, each having 3 effects from the same study.

set.seed(87)

# Parameters
no.paper = 20  # Numbers of unique papers
n.effects = 3  # Number of effects per paper. We will keep simple and balanced
rho.e = 0.8  # Correlation among sampling variances
study_id = rep(1:no.paper, each = n.effects)  # Study ID's
var_paper = 1  # Between-study variance (i.e., tau2)
var_effect = 0.8  # Effect size (or within study) variance
mu = 0.4  # Average, or overall, effect size
rho = 0.1  # Correlation among effect sizes within study; could vary

# Add sampling variance First, sample
n.sample <- rnbinom(no.paper, mu = 40, size = 0.5) + 4

# Assume logRR. So, sampling variance approximated by the CV.
cv <- sample(seq(0.3, 0.35, by = 0.001), no.paper, replace = TRUE)

# Assuming sample size is the same for all effects within study (pretty sensible)
rep_sample <- rep(n.sample, each = n.effects)

# Sampling variances
vi <- (2 * cv^2)/rep_sample

# Create the sampling (co)variance matrix
cov_vi <- as.matrix(Matrix::bdiag(impute_covariance_matrix(vi = vi, cluster = study_id, r = rho.e)))

# Now create the effect size level covariance matrix. This would simulate a situation where, for
# example, you have effect sizes on reading and maths that are correlated, but to varying extents
# across studies. Of course, this could vary within study too.

# Create list of matrices, 1 for each study
matrices <- list()
for (i in 1:no.paper) {
    r <- sample(rho, size = 1, replace = TRUE)
    tmp <- matrix(r, nrow = n.effects, ncol = n.effects)
    diag(tmp) <- 1
    matrices[[i]] <- tmp
}

# Build the full correlation matrix
cor_yi <- as.matrix(Matrix::bdiag(matrices))

# Calculate the full covariance matrix
cov_yi <- cor_yi * sqrt(var_effect) * sqrt(var_effect)

# Now simulate effect sizes, assuming that the average effect size and all the relevant within study
# correlations.

yi <- mu + rep(rnorm(no.paper, 0, sqrt(var_paper)), each = n.effects) + mvrnorm(n = 1, mu = rep(0, n.effects * 
    no.paper), cov_yi) + mvrnorm(n = 1, mu = rep(0, n.effects * no.paper), cov_vi)

# Create the data. We'll just assume that meta-analysts have already derived their effect size and
# sampling variance
data <- data.frame(study_id = study_id, yi = yi, vi = vi, obs = c(1:(n.effects * no.paper)))
```

Now that we have our simulated data, we can demonstrate a few corrections that can be applied to MLMA models that will offset any possible inflated Type I error rates.

# Step 1: Fit the Multi-level Meta-analytic (MLMA) Model

First, lets just fit our multilevel meta-analytic (MLMA) model. We can do that using our simulated data as follows:


```{.r .klippy}
mod_multilevel = metafor::rma.mv(yi = yi, V = vi, mods = ~1, random = list(~1 | study_id, ~1 | obs), 
    data = data, test = "t")
summary(mod_multilevel)
```

```
## 
## Multivariate Meta-Analysis Model (k = 60; method: REML)
## 
##   logLik  Deviance       AIC       BIC      AICc 
## -95.5412  191.0824  197.0824  203.3150  197.5188   
## 
## Variance Components:
## 
##             estim    sqrt  nlvls  fixed    factor 
## sigma^2.1  0.9791  0.9895     20     no  study_id 
## sigma^2.2  0.9304  0.9646     60     no       obs 
## 
## Test for Heterogeneity:
## Q(df = 59) = 21024.3536, p-val < .0001
## 
## Model Results:
## 
## estimate      se    tval    pval    ci.lb   ci.ub 
##   0.4262  0.2545  1.6747  0.0993  -0.0830  0.9354  . 
## 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
We have fit a simple MLMA model that estimates the overall meta-analytic mean. We can see that our model estimates this correctly. Remember that the true mean is 0.4, and we are pretty close to this value (i.e., 0.43). This is also true of our random effect variance estimates. In this case, we know that the MLMA model is not completely dealing with the dependence among effect sizes. We expect that this 'should' (at least on average) inflate Type I error rates. Assuming we did not know any better we would want to account for this dependence. Below, we describe a few corrections that meta-analysts can apply that should overcome the problems associated with not accounting for effect size dependence. 

# Correction 1: Using a Robust Variance Estimator (RVE) with Saitterwaite Degrees of Freedom Correction

A very simple solution is to make use of robust variance estimators (RVE). This can be done in a few packages, but very easily using the `clubSandwich` package [@Pustejovsky2020; @Hedges2010; @Tipon2015], which also works well with `metafor` [@Wolfgang2010] models. This approach also makes use of a Saitterwaite degrees of freedom correction [@SW; @Tipon2015]. This works with `metafor` objects quite elegantly. The benefit of such an approach is simply that we need not make any assumptions about what the correlation between effect sizes actually is (assuming we didn't know the true correlation) [@Hedges2010; @Tipon2015]. In addition, it also will account for possible heteroscedascity. This solution can be implemented as follows using our MLMA model we fit in the above section.


```{.r .klippy}
mod_RVE <- coef_test(mod_multilevel, vcov = "CR2", cluster = data$study_id)
print(mod_RVE)
```

```
##     Coef. Estimate    SE t-stat d.f. p-val (Satt) Sig.
## 1 intrcpt    0.426 0.254   1.67   19         0.11
```

A better, but slightly more restricted RVE can be implemented in the `robumeta` package in R. It is better at dealing with non-independence, but is currently limited to a single random effect level. Nonetheless, with our simple model we can fit a RVE model that completely deals with non-independence as follows:


```{.r .klippy}
mod <- robumeta::robu(formula = yi ~ 1, data = data, studynum = study_id, var.eff.size = vi, method = "HIER", 
    small = FALSE)
print(mod)
```

```
## RVE: Correlated Effects Model  
## 
## Model: yi ~ 1 
## 
## Number of studies = 20 
## Number of outcomes = 60 (min = 3 , mean = 3 , median = 3 , max = 3 )
## Rho = 0.8 
## I.sq = 100 
## Tau.sq = 1.8 
## 
##                Estimate StdErr t-value dfs P(|t|>) 95% CI.L 95% CI.U Sig
## 1 X.Intercept.    0.425  0.255    1.67  19   0.112   -0.108    0.958    
## ---
## Signif. codes: < .01 *** < .05 ** < .10 *
## ---
```

With this simple (and well balanced) data, our RVE approaches don't change the results much, but this won't always be the case. 

# Correction 2:  Modeling the Entire Sampling Covariance Matrix

Of course, we can also take an approach proposed by @Noble2017, where we fit the covariance matrix directly by simply assuming that effects that come from the same study are correlated by r = 0.5. Ultimately, one could change this correlation, depending on the situation and context, but r = 0.5 will probably suffice in many situations. This assumes, however, that the degree of correlation among effect sizes within a study is the same across studies. This assumption is relaxed in the RVE approaches described above. We can also test whether this is a safe assumption by combining it with a ClubSandwich estimator. We can build the matrix an implement this approach as follows:


```{.r .klippy}
vcv <- impute_covariance_matrix(vi = data$vi, cluster = data$study_id, r = 0.5)

mod_multilevel_vcv <- metafor::rma.mv(yi = yi, V = vcv, mods = ~1, random = list(~1 | study_id, ~1 | 
    obs), data = data, test = "t")

mod_multilevel_vcv <- coef_test(mod_multilevel_vcv, vcov = "CR2", cluster = data$study_id)
print(mod_multilevel_vcv)
```

```
##     Coef. Estimate    SE t-stat d.f. p-val (Satt) Sig.
## 1 intrcpt    0.425 0.255   1.67   19        0.112
```

# Correction 3: Applying Bayesian Multi-level Meta-analytic Model

As we describe in our comment, Bayesian approaches, assuming one has a good sample size, do a very good job correcting for inflated Type I errors across a variety of situations. Bayesian MLMA models can be fit in various packages. Probably the most flexible for meta-analyst's are `MCMCglmm` [@Hadfield2010] and `brms` [@Brkner2017; @Brkner2018]. Here, we demonstrate how to fit the same MLMA model using `MCMCglmm` which has a syntax that is different from the typical one used in packages such as `metafor` and `lme4` [@Bates2015], which meta-analysts might be more accustomed too. 


```{.r .klippy}
prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)))

bayes_multilevel <- MCMCglmm(yi ~ 1, mev = data$vi, random = ~study_id, data = data, prior = prior, verbose = FALSE)
summary(bayes_multilevel)
```

```
## 
##  Iterations = 3001:12991
##  Thinning interval  = 10
##  Sample size  = 1000 
## 
##  DIC: 185 
## 
##  G-structure:  ~study_id
## 
##          post.mean l-95% CI u-95% CI eff.samp
## study_id      1.17    0.353     2.41      578
## 
##  R-structure:  ~units
## 
##       post.mean l-95% CI u-95% CI eff.samp
## units     0.984    0.605     1.49     1000
## 
##  Location effects: yi ~ 1 
## 
##             post.mean l-95% CI u-95% CI eff.samp pMCMC
## (Intercept)     0.423   -0.104    0.953     1000  0.12
```

As expected, our credible intervals get a little bit wider. Bayesian models are the most conservative here given this is not a large data set.

# Conclusions
The goal of our short tutorial was to dispel the idea that overcoming, and implementing, solutions to deal with non-independent effect sizes when working with multi-level meta-analytic models is challenging. Our simulations show [@Nakagawa2021] that there are a number of very easily implemented solutions. As such, meta-analyst's can harness the power of MLMA models without the need to average effect sizes, as suggested by @Song2020. 

# Session Info
**R version 3.6.0 (2019-04-26)**

**Platform:** x86_64-apple-darwin15.6.0 (64-bit) 

**locale:**
en_AU.UTF-8||en_AU.UTF-8||en_AU.UTF-8||C||en_AU.UTF-8||en_AU.UTF-8

**attached base packages:** 
_stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_pander(v.0.6.3)_, _clubSandwich(v.0.5.1)_, _robumeta(v.2.0)_, _metafor(v.2.4-0)_, _brms(v.2.12.0)_, _Rcpp(v.1.0.5)_, _MCMCglmm(v.2.29)_, _ape(v.5.4-1)_, _coda(v.0.19-3)_, _Matrix(v.1.2-18)_, _gridExtra(v.2.3)_, _kableExtra(v.1.1.0)_, _MASS(v.7.3-51.6)_, _forcats(v.0.5.0)_, _stringr(v.1.4.0)_, _dplyr(v.1.0.2)_, _purrr(v.0.3.4)_, _readr(v.1.3.1)_, _tidyr(v.1.1.1)_, _tibble(v.3.0.3)_, _ggplot2(v.3.3.2)_, _tidyverse(v.1.3.0)_, _klippy(v.0.0.0.9500)_ and _knitr(v.1.29)_

**loaded via a namespace (and not attached):** 
_cubature(v.2.0.4.1)_, _colorspace(v.1.4-1)_, _ellipsis(v.0.3.1)_, _ggridges(v.0.5.2)_, _rsconnect(v.0.8.16)_, _corpcor(v.1.6.9)_, _markdown(v.1.1)_, _base64enc(v.0.1-3)_, _fs(v.1.5.0)_, _rstudioapi(v.0.11)_, _rstan(v.2.19.3)_, _DT(v.0.17)_, _fansi(v.0.4.1)_, _mvtnorm(v.1.1-0)_, _lubridate(v.1.7.9)_, _xml2(v.1.3.2)_, _bridgesampling(v.1.0-0)_, _shinythemes(v.1.1.2)_, _bayesplot(v.1.7.1)_, _jsonlite(v.1.7.2)_, _broom(v.0.7.0)_, _dbplyr(v.1.4.4)_, _shiny(v.1.4.0.2)_, _compiler(v.3.6.0)_, _httr(v.1.4.2)_, _backports(v.1.1.9)_, _assertthat(v.0.2.1)_, _fastmap(v.1.0.1)_, _cli(v.2.0.2)_, _later(v.1.1.0.1)_, _formatR(v.1.7)_, _prettyunits(v.1.1.1)_, _htmltools(v.0.5.1.1)_, _tools(v.3.6.0)_, _igraph(v.1.2.5)_, _gtable(v.0.3.0)_, _glue(v.1.4.1)_, _reshape2(v.1.4.4)_, _cellranger(v.1.1.0)_, _vctrs(v.0.3.2)_, _nlme(v.3.1-147)_, _crosstalk(v.1.1.0.1)_, _tensorA(v.0.36.1)_, _xfun(v.0.16)_, _ps(v.1.3.4)_, _rvest(v.0.3.6)_, _mime(v.0.9)_, _miniUI(v.0.1.1.1)_, _lifecycle(v.0.2.0)_, _pacman(v.0.5.1)_, _gtools(v.3.8.2)_, _zoo(v.1.8-8)_, _scales(v.1.1.1)_, _colourpicker(v.1.0)_, _hms(v.0.5.3)_, _promises(v.1.1.1)_, _Brobdingnag(v.1.2-6)_, _sandwich(v.2.5-1)_, _parallel(v.3.6.0)_, _inline(v.0.3.15)_, _shinystan(v.2.5.0)_, _StanHeaders(v.2.21.0-1)_, _loo(v.2.2.0)_, _stringi(v.1.5.3)_, _dygraphs(v.1.1.1.6)_, _pkgbuild(v.1.1.0)_, _rlang(v.0.4.10)_, _pkgconfig(v.2.0.3)_, _matrixStats(v.0.56.0)_, _evaluate(v.0.14)_, _lattice(v.0.20-41)_, _rstantools(v.2.0.0)_, _htmlwidgets(v.1.5.1)_, _processx(v.3.4.5)_, _tidyselect(v.1.1.0)_, _plyr(v.1.8.6)_, _magrittr(v.1.5)_, _R6(v.2.4.1)_, _generics(v.0.0.2)_, _DBI(v.1.1.0)_, _pillar(v.1.4.6)_, _haven(v.2.3.1)_, _withr(v.2.4.1)_, _xts(v.0.12-0)_, _abind(v.1.4-5)_, _modelr(v.0.1.8)_, _crayon(v.1.3.4)_, _rmarkdown(v.2.3)_, _grid(v.3.6.0)_, _readxl(v.1.3.1)_, _callr(v.3.5.1)_, _blob(v.1.2.1)_, _threejs(v.0.3.3)_, _reprex(v.0.3.0)_, _digest(v.0.6.27)_, _webshot(v.0.5.2)_, _xtable(v.1.8-4)_, _httpuv(v.1.5.2)_, _stats4(v.3.6.0)_, _munsell(v.0.5.0)_, _viridisLite(v.0.3.0)_ and _shinyjs(v.1.1)_

# References

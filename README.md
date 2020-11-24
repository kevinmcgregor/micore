[![Travis-CI Build Status](https://travis-ci.org/kevinmcgregor/micore.svg?branch=master)](https://travis-ci.org/kevinmcgregor/micore)

# micore
**M**icrobiome **Co**variance **Re**gression (**micore**) performs covariance regression in a multinomial logistic-normal model to estimate how microbiome co-occurrence networks vary with respect to covariates.  This work was developed in the [Greenwood Lab](https://www.mcgill.ca/statisticalgenetics/) at McGill University.

## Installation
You can install micore directly from Github:
```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("kevinmcgregor/micore", dependencies=TRUE, build_vignettes = TRUE)
```
## Vignette
Once you've successfully installed **micore**, you can access the vignette by running:
```r
vignette("micore-vignette")
```

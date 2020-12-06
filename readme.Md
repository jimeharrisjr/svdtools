svdtools
=====
### SVDTOOLS - tools for reducing and filtering matrices
<!-- badges: start -->
[![Travis build status](https://travis-ci.com/jimeharrisjr/svdtools.svg?branch=master)](https://travis-ci.com/jimeharrisjr/svdtools)
<!-- badges: end -->

Single Value Decomposition is a useful tool for the analysis and filtering of matrices. I regularly use these SVD functions for things like:

- Reducing noise in images (with the imager package, for example)
- Filtering collections of signal trace spectra to remove noise
- Preventing overfitting in tensorflow models with limited training data

### Usage Examples

Single Value Decomposition can remove noise from matrices in the process of reducing them. One example is image filtering. The included data object _noisymatrix_ is one such example. This is a grayscale image of Linux Penguin on Wikimedia Commons from Larry Ewing <lewing@isc.tamu.edu>. At random, some pixels have been deleted, creating noise.

``` r
library(svdtools)
# show the noisy image
image(noisymatrix, col=gray.colors(65536))
# Plot the explanation
plot_explanation(noisymatrix, limit=.99)
```

There are a lot of components to play with. One brute force way of doing this would be to try specifying some percentage of components to keep. We could try 98%:

``` r
m<-reduce_percentage(noisymatrix, .98)
image(m, col=gray.colors(65536))
```
This image is obviously overfiltered. Instead, we could simply eliminate all the components above some threshold - for example, discard 51-411:
```{r}
# Keep the first 50 components
m<-reduce_components(noisymatrix, 50)
image(m, col=gray.colors(65536))
```

See the vignettes for more detail, and for signal filter examples.

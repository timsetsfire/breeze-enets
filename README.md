# en4s (Elastic Net 4 Scala)

Simple framework for running generalized linear models or robust linear models with elastic net regularization in scala.  Optimization is completed via coordinate descent as laid out in [Regularization Paths for Generalized Linear Models via Coordinate Descent](http://web.stanford.edu/~hastie/Papers/glmnet.pdf)

## Families

This currently supports 
* Gaussian with identity link
* Poisson with log link
* Binomial with log link 
* Negative Binomial with logit link

## Robust Norms

When specifying the Gaussian family with identity link, you can also specify a robust norm for a regularized Robust Linear Model (RLM).  These classes are based on the implementation in statsmodels for Python.  

Availabe robust norms include
* Least Squares - default norm
* Huber T
* Ramsey E
* Tukey Biweight
* Cauchy
* Trimmed Mean
* Approximate Huber (smooth)
* L1

### Deprecation warnings

Element-wise multipilication in breeze via `:*` has been deprecated in favor of `*:*`.  Please be aware you will receive many 
deprecation warnings.  

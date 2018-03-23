package com.github.timsetsfire.enets

import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import com.github.timsetsfire.enets.utils._
import com.github.timsetsfire.enets.robust.norms._
import com.github.timsetsfire.enets.robust.scale.Mad._
import com.github.timsetsfire.enets.families._
import com.github.timsetsfire.enets.optimize.CoordinateDescent
import breeze.plot._


/** Estimate Robust Linear Model with Elastic Net Regularization 
  * using methods decribed in Regularization Paths for Generalized 
  * linear model via Coordinate Descent 
  * and Elements of Statistical Learning 
  * @constructor create a RlmNet model 
  * @param feature DenseMatrix of doubles containing features 
  * @param target DenseVector of doubles containing the target 
  * @param rnorm robust norm, default is HuberT()
  * @param alpha mixing parameter for L1 and L2 Regularization
  * @param tolerance for coordinate descent 
  * @param standardizeFeatures boolean value, whether to standardize
  * the feature matrix to mean 0 and unit variance, default is true 
  * @param standardizeTarget boolean value, whether to standardize 
  * the target vector to mean 0 and unit variance, default is false 
  * @param intercept boolean value, whether or not to include a constant, 
  * default is true 
  */ 

class RlmNet (
               features: DenseMatrix[Double],
               target: DenseVector[Double],
               rnorm: RobustNorm = LeastSquares,
               lambdaSeq: DenseVector[Double] = DenseVector.zeros[Double](100),
               alpha: Double = 1d,
               tolerance: Double = 1e-3,
               standardizeFeatures: Boolean = true,
               standardizeTarget: Boolean = false,
               intercept: Boolean = true
              ) extends ElasticNet(features, target, DenseVector(1d), "gaussian", "identity", rnorm, lambdaSeq, alpha, tolerance, standardizeFeatures, standardizeTarget, intercept)
              
/** Estimate Generalized Linear Model with Elastic Net Regularization 
  * using methods decribed in Regularization Paths for Generalized 
  * linear model via Coordinate Descent 
  * and Elements of Statistical Learning 
  * @constructor create a RlmNet model 
  * @param feature DenseMatrix of doubles containing features 
  * @param target DenseVector of doubles containing the target
  * @param offset DenseVector of doubles containing the offset
  * @param family for the model, values include gaussian, poisson,
  * negbin, binomial.  Default is Gaussian
  * @param link for the expected value of the target.  Doesn't really 
  * need to be included because it is set depending on the family choosen
  * @param alpha mixing parameter for L1 and L2 Regularization
  * @param tolerance for coordinate descent 
  * @param standardizeFeatures boolean value, whether to standardize
  * the feature matrix to mean 0 and unit variance, default is true 
  * @param standardizeTarget boolean value, whether to standardize 
  * the target vector to mean 0 and unit variance, default is false 
  * @param intercept boolean value, whether or not to include a constant, 
  * default is true 
  */ 
class GlmNet(
               features: DenseMatrix[Double],
               target: DenseVector[Double],
               offset: DenseVector[Double] = DenseVector(1d),
               family: String = "gaussian",
               link: String = "identity",
               lambdaSeq: DenseVector[Double] = DenseVector.zeros[Double](100),
               alpha: Double = 1d,
               tolerance: Double = 1e-3,
               standardizeFeatures: Boolean = true,
               standardizeTarget: Boolean = false,
               intercept: Boolean = true 
             ) extends ElasticNet(features, target, offset, family, link, LeastSquares, lambdaSeq, alpha, tolerance, standardizeFeatures, standardizeTarget, intercept)
             
/** create an ElasticNet instance.  Shouldn't be used directly
  */ 
class ElasticNet (
                   val features: DenseMatrix[Double],
                   val target: DenseVector[Double],
                   val offset: DenseVector[Double] = DenseVector(1d),
                   val family: String = "gaussian",
                   val link: String = "identity",
                   val rnorm: RobustNorm = LeastSquares,
                   val lambdaSeq: DenseVector[Double] = DenseVector.zeros[Double](100),
                   val alpha: Double = 1d,
                   val tolerance: Double = 1e-3,
                   val standardizeFeatures: Boolean = true,
                   val standardizeTarget: Boolean = false,
                   val intercept: Boolean = false
                 ) {

  // standardize features
  val (x, xm, xsig) = {
    val z = stdizeMatrix(features)
    if(standardizeFeatures) z else (features.copy, z._2, z._3)
  }
  // standardize target
  val (y, ym, ysig) = {
    val z = stdizeVector(target)
    if(standardizeTarget) z else (target.copy, z._2, z._3)
  }
  // nobs and nfeatures
  val (m, n) = (features.rows, features.cols)

  // initialize weight 
  val weight = y - mean(y)

  // apply robust norm
  weight := rnorm.w( weight / mad(weight) )

  // initialize parameters
  val parms = Params( DenseVector.zeros[Double](n), DenseVector(0d), intercept)
  parms.scale(0) = 1e-4
  parms.b0(0) = if(link=="log" & (offset.length == y.length)) {
    log(mean(y /:/ offset))
  } else if(link=="log") {
    log(mean(y))
  } else if(link=="logit") {
    log( mean(y) / (1d - mean(y)))
  } else {
    mean(y)
  }

  val exposure = if(offset.length == 1) DenseVector.ones[Double](y.length) else offset.copy

  // set distribution
  val dist: Family = {
    if (family == "negbin") new NegativeBinomial(x,y,weight,exposure)
    else if(family == "poisson") new Poisson(x,y,weight,exposure)
    else if(family == "binomial") new Binomial(x,y,weight, exposure)
    else Gaussian(x,y,weight, exposure)
  }

  /** weight function that will be passed to coordinate descent 
    * @param r DenseVector of doubles.  Will typically be the residual
    */ 
  def weightFunction(r: DenseVector[Double]) = {
    rnorm.w( r / mad(r) )
  }
  /** regularized cost function.  curried
    * @param dist family for the cost function - will exposue the deviance for 
    * the calculation 
    * @param parms Params of the linear Model
    * @param lambda regularization term 
    * @param alpha mixing of the L1 and L2 regularization 
    */ 
  def costFunc(dist: Family)(parms: Params, lambda: Double, alpha: Double) = {
    //deviance based cost function
    val d = -1d*dist.ll(parms)
    val l1 = norm(parms.b, 1)
    val l2 = norm(parms.b, 2)
    d + lambda * ( alpha * 1d + (1d - alpha)* 1d)
  }

  val b = DenseMatrix.zeros[Double](n, lambdaSeq.length)
  val b0 = DenseVector.zeros[Double](lambdaSeq.length)

  val nz = (0 until n).toArray
  
  // create the lambda sequence 
  if ( sum(abs(lambdaSeq)) == 0d ) {
    //  val yt = if(family == "negbin") y / ( mean(y) + mean(y) * mean(y) ) else y.copy
    val (xt, xs, xm) = stdizeMatrix(x)
    val wx = xt(::, *).map{ _ *:* weight}
    val xTy = 1 / (m).toDouble * wx.t*(y)
    val lambdaMax = max( abs( xTy ) ) / (if(alpha==0) 0.001 else alpha)
    val lambdaMin = if(x.rows < x.cols) lambdaMax * 1e-2 else lambdaMax * 1e-4
    lambdaSeq := linspace(log(lambdaMax), log(lambdaMin), 100)
  }

  // initialize vector to store the deviance
  val deviance = DenseVector.zeros[Double](lambdaSeq.length)

  /** fit the elastic net 
    */ 
  def fit: Unit = {
    for(iter <- 0 until lambdaSeq.length) {
      val lambda = 2*exp(lambdaSeq(iter)) - exp(lambdaSeq(iter - 1))
      val z = dist.z(x, y, exposure, parms)
      val s = dist.w(x, y, exposure, parms)
      val r = (y - dist.yhat(x, exposure, parms) )
      weight := rnorm.w(r / mad(r))
      parms.nzbv := (- abs( x.t * r )/x.rows.toDouble <:< - lambda * alpha) |:| parms.nzbv
      val nz = parms.nzbv.toArray.zipWithIndex.filter{ _._1 == true}.map{ _._2}
      val cd = new CoordinateDescent(costFunc(dist), x, z, weight, s, weightFunction)
      cd.optimize( alpha, lambda, parms, tolerance, nz = nz)
      parms.nzbv := ( parms.b :!= 0d )
      if(family=="negbin") {
        val df = (x.rows - parms.nzbv.activeSize).toDouble
        val deviance = dist.dev(parms)
        parms.scale(0) = if (df > 0) deviance / df else 0.0001
      }
      b(::, iter) := parms.b
      b0(iter) = parms.b0(0)
      deviance(iter) = dist.dev(parms)
    }
  }
  /** plot the coordinate path 
    */ 
  def plotCoordinatePath = {
    val nonZero = (b(::, -1) :!= 0d).activeKeysIterator.toArray
    val f = Figure()
    val p = f.subplot(0)
    for(i <- nonZero) p += plot(lambdaSeq, b.t(::, i))
    f.refresh
  }
}

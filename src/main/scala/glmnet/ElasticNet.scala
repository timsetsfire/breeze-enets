package com.github.timsetsfire.en4s

import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import com.github.timsetsfire.en4s.utils._
import com.github.timsetsfire.en4s.robust.norms._
import com.github.timsetsfire.en4s.robust.scale.Mad._
import com.github.timsetsfire.en4s.families._
import com.github.timsetsfire.en4s.optimize.CoordinateDescent
import breeze.plot._

/**
  * Created by WhittakerT on 10/26/2017.
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
               intercept: Boolean = false
              ) extends ElasticNet(features, target, DenseVector(1d), "gaussian", "identity", rnorm, lambdaSeq, alpha, tolerance, standardizeFeatures, standardizeTarget, intercept)

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

  val (x, xm, xsig) = {
    val z = stdizeMatrix(features)
    if(standardizeFeatures) z else (features.copy, z._2, z._3)
  }
  val (y, ym, ysig) = {
    val z = stdizeVector(target)
    if(standardizeTarget) z else (target.copy, z._2, z._3)
  }
  val (m, n) = (features.rows, features.cols)

  val weight = y - mean(y)

  weight := rnorm.w( weight / mad(weight) )

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

  val dist: Family = {
    if (family == "negbin") new NegativeBinomial(x,y,weight,exposure)
    else if(family == "poisson") new Poisson(x,y,weight,exposure)
    else if(family == "binomial") new Binomial(x,y,weight, exposure)
    else Gaussian(x,y,weight, exposure)
  }

  //if(family=="negbin") parms.scale(0) = dist.dev(x,y,w)(parms) / (m).toDouble
  def weightFunction(r: DenseVector[Double]) = {
    rnorm.w( r / mad(r) )
  }
  def costFunc(dist: Family)(parms: Params, lambda: Double, alpha: Double) = {
    //deviance based cost function
    val d = -1d*dist.ll(parms)
    val l1 = norm(parms.b, 1)
    val l2 = norm(parms.b, 2)
    d + lambda * ( alpha * 1d + (1d - alpha)* 1d)
  }

  val b = DenseMatrix.zeros[Double](n, lambdaSeq.length)
  val b0 = DenseVector.zeros[Double](lambdaSeq.length)

  //val nullCost = costFunc(x, y)(parms)

  val nz = (0 until n).toArray

  // def c = costFunc(x,y)_


  if ( sum(abs(lambdaSeq)) == 0d ) {
    //  val yt = if(family == "negbin") y / ( mean(y) + mean(y) * mean(y) ) else y.copy
    val (xt, xs, xm) = stdizeMatrix(x)
    val wx = xt(::, *).map{ _ *:* weight}
    val xTy = 1 / (m).toDouble * wx.t*(y)
    val lambdaMax = max( abs( xTy ) ) / (if(alpha==0) 0.001 else alpha)
    val lambdaMin = if(x.rows < x.cols) lambdaMax * 1e-2 else lambdaMax * 1e-4
    lambdaSeq := linspace(log(lambdaMax), log(lambdaMin), 100)
  }

  val deviance = DenseVector.zeros[Double](lambdaSeq.length)

  def fit: Unit = {
    // don't touch fast3 this is it!!!!!
    for(iter <- 0 until lambdaSeq.length) {
      val lambda = 2*exp(lambdaSeq(iter)) - exp(lambdaSeq(iter - 1))
      val z = dist.z(x, y, exposure, parms)
      val s = dist.w(x, y, exposure, parms)
      val r = (y - dist.yhat(x, exposure, parms) )
      weight := rnorm.w(r / mad(r))
      parms.nzbv := (- abs( x.t * r )/x.rows.toDouble <:< - lambda * alpha) |:| parms.nzbv
      val nz = parms.nzbv.toArray.zipWithIndex.filter{ _._1 == true}.map{ _._2}
      val cd = new CoordinateDescent(costFunc(dist), x, z, weight, s, weightFunction)
	  //cd.optimize( alpha, lambdaSeq(iter), parms, tolerance, nz = nz)
      cd.optimize( alpha, lambda, parms, tolerance, nz = nz)
      parms.nzbv := ( parms.b :!= 0d )
      if(family=="negbin") {
        //val info = nbest.fitAndReturn
        val df = (x.rows - parms.nzbv.activeSize).toDouble
        //parms.scale(0) = if( abs(info.grad(0)) < 1e-4 & df > 0) info.x(0) else 1e-4
        val deviance = dist.dev(parms)
        parms.scale(0) = if (df > 0) deviance / df else 0.0001
      }
      b(::, iter) := parms.b
      b0(iter) = parms.b0(0)
      deviance(iter) = dist.dev(parms)
    }
  }

  def plotCoordinatePath = {
    val nonZero = (b(::, -1) :!= 0d).activeKeysIterator.toArray
    val f = Figure()
    val p = f.subplot(0)
    for(i <- nonZero) p += plot(lambdaSeq, b.t(::, i))
    f.refresh
  }
}

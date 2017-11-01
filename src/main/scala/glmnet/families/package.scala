package com.github.timsetsfire.en4s

/**
  * Created by WhittakerT on 10/26/2017.
  */
package object families {

  import breeze.numerics._
  import breeze.linalg._
  import com.github.timsetsfire.en4s.utils._

  trait Family {
    def ll(b: Params): Double
    def dev(b: Params): Double
    def yhat(x: DenseMatrix[Double], offset: DenseVector[Double], b: Params): DenseVector[Double]
    def w(x: DenseMatrix[Double], y: DenseVector[Double],offset: DenseVector[Double], b: Params): DenseVector[Double]
    def z(x: DenseMatrix[Double], y: DenseVector[Double],offset: DenseVector[Double], b: Params): DenseVector[Double]
  }

  object Poisson  {
    def deviance(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val mu = if(offset.length == y.length) exp( g(x,b) + log(offset)) else exp( g(x,b))
      val t1 = y.map{ v => if(v==0) 0 else v*log(v)} - (y*:*log(mu))
      2d * sum( w*:* ( t1 - (y-mu)))
    }
    def logLikelihood(x: DenseMatrix[Double], y: DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val mu = if(offset.length == y.length) exp( g(x,b) + log(offset)) else exp( g(x,b))
      val ll = w*:*( (y*:*log(mu)) - mu - lgamma(y + 1d))
      sum(ll)
    }
  }
  case class Poisson(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double] ) extends Family {
    import Poisson._
    def ll(b: Params) = logLikelihood(x,y,w,offset)(b)/ y.length.toDouble
    def dev(b: Params) = deviance(x,y,w,offset)(b)
    def yhat(x: DenseMatrix[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = exp( g(x,b) + log(offset))
    def w(x: DenseMatrix[Double], y:DenseVector[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = exp(g(x,b) + log(offset))
    def z(x: DenseMatrix[Double], y: DenseVector[Double],offset: DenseVector[Double], b: Params): DenseVector[Double] = g(x,b) + ((y - yhat(x, offset, b)) *:* (1d / w(x,y,offset, b)))
  }

  object NegativeBinomial {
    def deviance(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val k = b.scale(0)
      val mu = if(offset.length == y.length) exp( g(x,b) + log(offset)) else exp( g(x,b))
      val t1 = y.map{ v => if(v==0) 0 else v*log(v)} - (y*:*log(mu))
      2d*sum(  t1 - ((y + w / k)*:*(log( y + w/k)-log(mu + w /k))))
    }
    def logLikelihood(x: DenseMatrix[Double], y: DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val k = b.scale(0)
      val mu = if(offset.length == y.length) exp( g(x,b) + log(offset)) else exp( g(x,b))
      val ll = (y*:*log( k*mu /:/ w)) - ((y + w * 1d/k)*:*log( (k*mu /:/ w) + 1d)) + lgamma(y + w*1d/k) - lgamma(y + 1d) - lgamma(w*1d/k)
      sum(ll)
    }
  }
  case class NegativeBinomial(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double] )  extends Family {
    import NegativeBinomial._
    def ll(b: Params) = logLikelihood(x,y,w,offset)(b)/ y.length.toDouble
    def dev(b: Params) = deviance(x,y,w,offset)(b)
    def yhat(x: DenseMatrix[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = exp( g(x,b) + log(offset))
    def w(x: DenseMatrix[Double], y:DenseVector[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = {
      val p = yhat(x, offset, b)
      val den = ( b.scale(0)*p + 1d).map{ r => math.pow(r, 2)}
      val num = (b.scale(0)*y + 1d) *:* p
      num /:/ den
    }
    def z(x: DenseMatrix[Double], y: DenseVector[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = {
      val p = yhat(x, offset, b)
      val num = ( b.scale(0)*p + 1d)
      val den = (b.scale(0)*y + 1d) *:* p
      val q = num /:/ den
      g(x,b) + ((y - p) *:* q)
    }
  }
  object Binomial {
    def deviance(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val mu = sigmoid( g(x,b))
      2d *sum( w *:* ( y*:*log(y/:/mu) + (-y + 1d)*:*log( (-y + 1d)/:/(-mu + 1d))))
    }
    def logLikelihood(x: DenseMatrix[Double], y: DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val mu = sigmoid( g(x,b))
      //val ll = (w*:*( (y *:* log(mu)) + ((-y + 1d)*:*log(-mu + 1d))))
      val ll = (y *:* ( g(x,b))) - log( exp( g(x, b)) + 1d)
      sum(ll)
    }

  }

  case class Binomial(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double] )  extends Family {
    import Binomial._
    def ll(b: Params) = logLikelihood(x, y, w, offset)(b)/ y.length.toDouble
    def dev(b: Params) = -2*ll(b) // deviance(x,y,w,offset)(b)
    def yhat(x: DenseMatrix[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = sigmoid( g(x,b))
    def w(x: DenseMatrix[Double], y:DenseVector[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = {
      val p = yhat(x, offset, b)
      p *:* (1d - p)
    }
    def z(x: DenseMatrix[Double], y: DenseVector[Double],offset: DenseVector[Double], b: Params): DenseVector[Double] = {
      g(x,b) + ((y - yhat(x, offset, b)) *:* (1d /:/ w(x, y,offset, b)))
    }
  }
  object Gaussian {
    def logLikelihood(x: DenseMatrix[Double], y: DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val e = y - g(x, b)
      -(w*:*e).dot(e) / 2d
    }
    def deviance(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double])(b: Params) = {
      val e = y - g(x, b)
      (w*:*e).dot(e)
    }
  }
  case class Gaussian(x: DenseMatrix[Double], y:DenseVector[Double], w: DenseVector[Double], offset: DenseVector[Double] )  extends Family {
    import Gaussian._
    def ll(b: Params) = logLikelihood(x, y, w, offset)(b)/ y.length.toDouble
    def dev(b: Params) = deviance(x,y,w, offset)(b)
    def yhat(x: DenseMatrix[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = g(x,b)
    def w(x: DenseMatrix[Double], y: DenseVector[Double], offset: DenseVector[Double], b: Params): DenseVector[Double] = DenseVector.ones[Double](x.rows)
    def z(x: DenseMatrix[Double], y: DenseVector[Double],offset: DenseVector[Double], b: Params): DenseVector[Double] = y
  }
}

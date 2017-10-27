package com.github.timsetsfire.glmnet.robust

/**
  * Created by WhittakerT on 10/26/2017.
  */
package object norms {

  import breeze.linalg._
  import breeze.numerics._

  trait RobustNorm {
    def rho(z: DenseVector[Double]): DenseVector[Double]
    // the robust ciriterion estiator function
    def psi(z: DenseVector[Double]): DenseVector[Double]
    // derivative of rho.  sometimes referred to as the influence function
    def psiDeriv(z: DenseVector[Double]): DenseVector[Double]
    // derivative of psi.
    // this can be used to obtain robust covariance matrix
    // for the purpose of enet, will not implement
    def w(z: DenseVector[Double]): DenseVector[Double]
    // returns the value of psi(z) / z
  }

  case object LeastSquares extends RobustNorm{
    def rho(z: DenseVector[Double]) = pow(z, 2) /2d
    def psi(z: DenseVector[Double]) = z
    def psiDeriv(z: DenseVector[Double]) = DenseVector.ones[Double](z.length)
    def w(z: DenseVector[Double]) = DenseVector.ones[Double](z.length)
  }

  case class HuberT(k: Double = 1.345) extends RobustNorm {
    def rho(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) < k) pow(elem,2) / 2d else k*(abs(elem) - k/2d) }
    def psi(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) < k) elem else k*signum(elem)}
    def psiDeriv(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) < k) 1d else 0d}
    def w(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) < k) 1 else k / abs (elem) }
  }

  case class RamsayE(a: Double = 0.3) extends RobustNorm {
    def rho(z: DenseVector[Double]): DenseVector[Double] = {
      val a2inv = pow(a, -2)
      z.map{ elem =>  a2inv * (1 - exp( -a * abs(elem)))*(1d + a * abs(elem)) }
    }
    def psi(z: DenseVector[Double]): DenseVector[Double] = z.map{ elem => elem * exp( - a * abs( elem) )}
    def psiDeriv(z: DenseVector[Double]): DenseVector[Double] = ???
    def w(z: DenseVector[Double]): DenseVector[Double] = exp(abs(z) * (-a))

  }

  case class TukeyBiweight(c: Double  = 4.685) extends RobustNorm {
    def rho(z: DenseVector[Double]) = z.map{ elem => if(abs(elem) <= c) pow(c,2)/6d * (1d - pow(1d - pow(elem / c,2), 3)) else pow(c,2) / 6d }
    def psi(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) <= c) elem * pow( 1d - pow(elem/c, 2), 2) else 0d}
    def psiDeriv(z: DenseVector[Double]) = ???
    def w(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) <= c) pow( 1d - pow(elem/c, 2), 2) else 0d}
  }

  case class Cauchy(c: Double = 2.3849) extends RobustNorm {
    val c2 = c*c
    def rho(z: DenseVector[Double]) = z.map{ elem => c2 / 2d * log( 1d + pow(elem/c,2))}
    def psi(z: DenseVector[Double]) = z.map{ elem => elem / (1d + pow(elem/c,2))}
    def psiDeriv(z: DenseVector[Double]) = ???
    def w(z: DenseVector[Double]) = z.map{ elem => 1d / (1d + pow(elem/c,2))}
  }

  case class TrimmedMean(c: Double = 2d) extends RobustNorm {
    def rho(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) <= c) 1/2d * pow(elem,2) else 0d}
    def psi(z: DenseVector[Double]) = z.map{ elem => if( abs(elem) <= c) elem else 0d}
    def psiDeriv(z: DenseVector[Double]) = ???
    def w(z: DenseVector[Double]) = z.map{ elem => if(abs(elem) <= c) 1 else 0d}
  }

  case object ApproxHuber extends RobustNorm {
    def rho(z: DenseVector[Double]) = z.map{ elem => 2d* ( sqrt( 1d + pow(elem,2)/2d) - 1d)}
    def psi(z: DenseVector[Double]) = z.map{ elem => elem / sqrt( 1d + pow(elem,2)/2d) }
    def psiDeriv(z: DenseVector[Double]) = ???
    def w(z: DenseVector[Double]) = z.map{ elem => 1d / sqrt( 1d + pow(elem,2)/2d) }
  }

  case object L1 extends RobustNorm {
    def rho(z: DenseVector[Double]) = abs(z)
    def psi(z: DenseVector[Double]) = signum(z)
    def psiDeriv(z: DenseVector[Double]) = ???
    def w(z: DenseVector[Double]) = z.map{ elem => if(elem == 0d) 1e6 else 1 / abs(elem) }
  }


}

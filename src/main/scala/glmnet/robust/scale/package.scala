package com.github.timsetsfire.en4s.robust


/**
  * Created by WhittakerT on 10/26/2017.
  */
package object scale {
  import breeze.linalg.DenseVector
  import breeze.numerics._
  import breeze.stats._
  import breeze.stats.distributions.Gaussian

  object Mad  {
    def mad(
             x: DenseVector[Double],
             c: Double = Gaussian(0,1).inverseCdf(3/4d),
             center: (DenseVector[Double] => Double) = (x: DenseVector[Double]) => 0d
    ) = {
      val m = center(x)
      median( abs( x - m) / c)
    }
  }

  class HuberScale(d: Double = 2.5, tol: Double = 1e-8, maxiter:Double = 30)
}

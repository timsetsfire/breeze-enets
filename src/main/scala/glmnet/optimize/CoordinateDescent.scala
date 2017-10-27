package com.github.timsetsfire.glmnet.optimize

import breeze.linalg._
import breeze.numerics._
import com.github.timsetsfire.glmnet.utils._

/**
  * Created by WhittakerT on 10/26/2017.
  */
class CoordinateDescent( cost: (Params, Double, Double) => Double,
                         x: DenseMatrix[Double],
                         y: DenseVector[Double],
                         w: DenseVector[Double],
                         s: DenseVector[Double],
                         weightFunction: (DenseVector[Double]) => DenseVector[Double]) {

  def optimize(alpha: Double, lambda: Double, params: Params, tolerance: Double = 1e-6, nz: Array[Int] = Array.empty) = {
    //println("optimizing")
    def descend(params: Params, ind: Int): Unit = {
      val b = params.b(ind)
      val wx = (x(::, ind) :* w) :* s
      //val wx = (x(::, ind) :* w)
      val t = 1/ (x.rows.toDouble) *( wx.t * ( (y - g(x,params) + (params.b(ind)*x(::, ind)) )))
      params.b(ind) = softThresh(t, alpha * lambda) / ( 1/x.rows.toDouble * (wx.t*x(::, ind)) + (1 - alpha)*lambda )
      if( abs( b - params.b(ind)) < tolerance) Unit
      else descend(params, ind)
    }
    def interceptDescend(params: Params): Unit = {
      val b0 = params.b0(0)
      params.b0(0) = 1 / sum(w:* s) * ((w:* s).t *( ( y - x*params.b)))
      //params.b0(0) = 1 / sum(w) * ((w).t *( ( y - x*params.b)))
      if( abs(b0 - params.b0(0)) < tolerance) Unit
      else interceptDescend(params)
    }

    def update(params: Params, tolerance: Double = 1e-6): Unit = {
      val j = cost(params, lambda, alpha)
      //val randind = DenseVector.rand(params.nparams, Binomial(1, scd)).toArray.zipWithIndex.filter{ tup => tup._1 == 1 & tup._2 > 0}.map{ _._2}
      if(params.intercept) interceptDescend(params)
      for( ind <- nz) {
        descend(params, ind)
      }

      if( abs(j - cost(params, lambda, alpha)) < tolerance) {
        //  println(s"lambda => $lambda\nalpha => $alpha\nk => ${params.scale(0)}")
        //  println( cost(params, lambda, alpha, params.scale(0)) )
        Unit
      }
      else {
        // println(s"lambda => $lambda\nalpha => $alpha\nk => ${params.scale(0)}")
        // println( cost(params, lambda, alpha, params.scale(0)) )
        val r = (y - g(x, params))
        w := weightFunction(r)
        update(params, tolerance)
      }
    }
    update(params, tolerance)
  }
}

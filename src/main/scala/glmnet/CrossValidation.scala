//package glmnet
//
//import breeze.linalg._
//import breeze.numerics._
//import glmnet.utils._
//
//
///**
//  * Created by WhittakerT on 10/26/2017.
//  */
//class CrossValidationBlocked(
//                              //costFunc: (DenseMatrix[Double], DenseVector[Double]) => Params => Double,
//                              features: DenseMatrix[Double],
//                              target: DenseVector[Double],
//                              lossType: String = "square-loss",
//                              link: String = "identity",
//                              alpha: Double = 1d,
//                              lambdaDepth: Int = 5,
//                              tolerance: Double = 1e-3,
//                              stepTolerance: Double = 1e-4,
//                              nahead: Int = 4,
//                              crossValidations: Int = 16,
//                              standardizeFeatures: Boolean = true,
//                              standardizeTarget: Boolean = true,
//                              intercept: Boolean = false,
//                              fastSolve: Boolean = true,
//                              arOrder: Int = 0,
//                              names: Array[String] = Array()
//                            )  {
//  val (m, n) = (features.rows, features.cols)
//  //val p = (2 * math.floor(percent * m)).toInt
//  val (x, xm, xsig) = {
//    val z = stdizeMatrix(features)
//    if(standardizeFeatures) z else (features.copy, z._2, z._3)
//  }
//  val (y, ym, ysig) = {
//    val z = stdizeVector(target)
//    if(standardizeTarget) z else (target.copy, z._2, z._3)
//  }
//
//  val xTy = 1 / m.toDouble * x.t*y
//  val lambdaMax = max( abs( xTy ) ) / (if(alpha==0) 0.001 else alpha)
//  val lambdaMin = if(x.rows < x.cols) lambdaMax * 1e-2 else lambdaMax * 1e-4
//  val lambdaSeq = linspace(log(lambdaMax), log(lambdaMin), 100).apply(0 to lambdaDepth)
//
//
//  // val ind = new DenseVector((0 until x.rows).toArray)
//  // val indices = (0 to p - 1).map{ i =>
//  //   val training = ((ind :> i) :& (ind :< (m - p + i + 1))).toArray.zipWithIndex.filter{ _._1 == true}.map{ _._2}.toIndexedSeq
//  //   val test = ((ind :> (m - p + i)) & (ind :<= (m - p + i + nahead) )).toArray.zipWithIndex.filter{_._1 == true}.map{_._2}.toIndexedSeq
//  //   (training, test)
//  // }.filter{ _._2.length == nahead}
//  // val crossValdiations = 20
//  // val nahead = 10
//
//  val indices = for{i <- 0 until crossValidations; mpi = m - crossValidations + i;
//                    trainStart = i;
//                    trainStop = mpi - nahead ;
//                    cvStart = mpi - nahead +1;
//                    cvStop = mpi
//  } yield (trainStart to trainStop,cvStart to cvStop)
//
//  val enets: IndexedSeq[(ElasticNet, DenseMatrix[Double], DenseVector[Double])] = indices.map{
//    case(training, test) => {
//      val (x1, y1) = ( x(training, ::).toDenseMatrix, y(training).toDenseVector )
//      val (x2, y2) = ( x(test, ::).toDenseMatrix, y(test).toDenseVector)
//      val en = new ElasticNet(x1, y1, lossType, link, lambdaSeq, alpha, tolerance, stepTolerance, false, false, true )
//      (en, x2, y2)
//      }
//    }
//
//
//  def g(x: DenseMatrix[Double], parms: Params): DenseVector[Double] = {
//    val z = if(parms.intercept) (x * parms.b).map{ _ + parms.b0(0)} else x * parms.b
//    if(link == "logit") sigmoid(z)
//    else if (link == "log") exp(z)
//    else z
//  }
//
//  def costFunc(x: DenseMatrix[Double], y:DenseVector[Double])(p: Params) = {
//    val yhat = g(x, p)
//    val e = y - yhat;
//    if (lossType == "absolute-loss") e.dot(e) / e.length.toDouble
//    else if (lossType == "huber-loss") CoordinateDescentH.huber(x,y)(p)
//    else sum( abs( e )) / e.length.toDouble
//  }
//
//  lazy val arrayCvMatrix = enets.map{ case(enet, x, y) =>
//    val s = enet.b.cols
//    val errors = DenseVector.range(0,s).map{ i =>
//      //val xb = x * enet.b
//      //val yhat = xb(*, ::).map{ _ + enet.b0}
//      //val err = pow( yhat(::, *).map{ _ - y}, 2);
//      //mean(err(::, *)).t
//      val b = Params( enet.b(::, i), DenseVector(enet.b0(i)), true)
//      val e2 = if(lossType == "huber-loss") huber(x,y)(b) else costFunc(x,y)(b)
//      //val yhat = g(x , b)
//      //val e = yhat - y;
//      //val e2 = e.t * e / (2 * e.length)
//      DenseVector(e2)
//    }
//    errors.toArray
//  }.map{ set => DenseMatrix( set:_*)}
//
//  lazy val cvErrorMatrix = DenseMatrix.horzcat(arrayCvMatrix:_*)
//
//  lazy val (cvMean, cvStddev) = ( mean( cvErrorMatrix(*, ::)), stddev(cvErrorMatrix(*, ::)))
//
//  // def plotCvError = {
//  //   val f = Figure()
//  //   val p = f.subplot(0)
//  //   val n = cvMean.t.length
//  //   p += plot( linspace(0,n,n), cvMean.t)
//  //   p += plot( linspace(0,n,n), cvMean.t + 2d*cvStddev.t)
//  //   p += plot( linspace(0,n,n), cvMean.t - 2d*cvStddev.t)
//  // }
//  // def plotCvErrorJupyter(plotOut: String) = {
//  //   val f = Figure()
//  //   f.visible = false
//  //   val p = f.subplot(0)
//  //   val n = cvMean.t.length
//  //   p += plot( linspace(0,n,n), cvMean.t)
//  //   p += plot( linspace(0,n,n), cvMean.t + 2d*cvStddev.t)
//  //   p += plot( linspace(0,n,n), cvMean.t - 2d*cvStddev.t)
//  //   f.saveas(plotOut)
//  // }
//
//  lazy val optLambdaIndex = cvMean.toArray.zipWithIndex.sortWith( _._1 < _._1).head._2
//
//  def fieldVote(lambdaIndex: Int, cutOff: Int): Array[(Int, String)] = {
//    val votes = DenseMatrix.zeros[Int](n, lambdaSeq.length)
//    enets.foreach{ case(enet,x,y) =>
//      val nzb = (enet.b :!= 0d)
//      votes += nzb(::, *).map{ _.map{ i => if(i==true) 1 else 0 }}
//    }
//    val sumVotes = votes(::, lambdaIndex).toArray.zip(names)
//    sumVotes.filter{ _._1 >= cutOff}.foreach{ case(value, label) => println( s"$label => $value")}
//    sumVotes.filter{ _._1 >= cutOff}
//    if( sumVotes.filter{ _._1 >= cutOff}.length == 0) fieldVote(lambdaIndex + 1, cutOff)
//    else sumVotes.filter{ _._1 >= cutOff}
//  }
//
//}

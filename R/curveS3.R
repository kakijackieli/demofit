#' Common S3 methods for mortality curves
#'
#' @description
#' Common S3 methods for objects returned by \code{MC()}:
#'
#' coef(object) gives estimated parameter values.
#'
#' fitted(object) gives fitted values.
#'
#' predict(object,x) gives estimated mortality rate(s) for given age(s) x.
#'
#' plot(object) plots data and fitted mortality curve.
#'
#' deviance(object) gives weighted sum of squared residuals on \eqn{ln(m_x)}.
#'
#' residuals(object) gives residuals on \eqn{ln(m_x)}.
#'
#' @name curveS3
#'
#' @method coef fit1
#' @method coef fit2
#' @method coef fit3
#' @method coef fit4
#' @method coef fit5
#' @method coef fit6
#' @method coef fit7
#' @method coef fit8
#' @method coef fit9
#' @method coef fit10
#' @method coef fit11
#' @method coef fit12
#' @method coef fit13
#' @method coef fit14
#' @method coef cfit1
#' @method coef cfit2
#' @method coef cfit3
#' @method coef cfit4
#' @method coef cfit5
#' @method coef cfit6
#' @method coef cfit7
#' @method coef cfit8
#' @method coef cfit9
#' @method coef cfit10
#' @method coef cfit11
#' @method coef cfit12
#' @method coef cfit13
#' @method coef cfit14
#'
#' @method fitted fit1
#' @method fitted fit2
#' @method fitted fit3
#' @method fitted fit4
#' @method fitted fit5
#' @method fitted fit6
#' @method fitted fit7
#' @method fitted fit8
#' @method fitted fit9
#' @method fitted fit10
#' @method fitted fit11
#' @method fitted fit12
#' @method fitted fit13
#' @method fitted fit14
#' @method fitted cfit1
#' @method fitted cfit2
#' @method fitted cfit3
#' @method fitted cfit4
#' @method fitted cfit5
#' @method fitted cfit6
#' @method fitted cfit7
#' @method fitted cfit8
#' @method fitted cfit9
#' @method fitted cfit10
#' @method fitted cfit11
#' @method fitted cfit12
#' @method fitted cfit13
#' @method fitted cfit14
#'
#' @method predict fit1
#' @method predict fit2
#' @method predict fit3
#' @method predict fit4
#' @method predict fit5
#' @method predict fit6
#' @method predict fit7
#' @method predict fit8
#' @method predict fit9
#' @method predict fit10
#' @method predict fit11
#' @method predict fit12
#' @method predict fit13
#' @method predict fit14
#' @method predict cfit1
#' @method predict cfit2
#' @method predict cfit3
#' @method predict cfit4
#' @method predict cfit5
#' @method predict cfit6
#' @method predict cfit7
#' @method predict cfit8
#' @method predict cfit9
#' @method predict cfit10
#' @method predict cfit11
#' @method predict cfit12
#' @method predict cfit13
#' @method predict cfit14
#'
#' @method plot fit1
#' @method plot fit2
#' @method plot fit3
#' @method plot fit4
#' @method plot fit5
#' @method plot fit6
#' @method plot fit7
#' @method plot fit8
#' @method plot fit9
#' @method plot fit10
#' @method plot fit11
#' @method plot fit12
#' @method plot fit13
#' @method plot fit14
#' @method plot cfit1
#' @method plot cfit2
#' @method plot cfit3
#' @method plot cfit4
#' @method plot cfit5
#' @method plot cfit6
#' @method plot cfit7
#' @method plot cfit8
#' @method plot cfit9
#' @method plot cfit10
#' @method plot cfit11
#' @method plot cfit12
#' @method plot cfit13
#' @method plot cfit14
#'
#' @method deviance fit1
#' @method deviance fit2
#' @method deviance fit3
#' @method deviance fit4
#' @method deviance fit5
#' @method deviance fit6
#' @method deviance fit7
#' @method deviance fit8
#' @method deviance fit9
#' @method deviance fit10
#' @method deviance fit11
#' @method deviance fit12
#' @method deviance fit13
#' @method deviance fit14
#' @method deviance cfit1
#' @method deviance cfit2
#' @method deviance cfit3
#' @method deviance cfit4
#' @method deviance cfit5
#' @method deviance cfit6
#' @method deviance cfit7
#' @method deviance cfit8
#' @method deviance cfit9
#' @method deviance cfit10
#' @method deviance cfit11
#' @method deviance cfit12
#' @method deviance cfit13
#' @method deviance cfit14
#'
#' @method residuals fit1
#' @method residuals fit2
#' @method residuals fit3
#' @method residuals fit4
#' @method residuals fit5
#' @method residuals fit6
#' @method residuals fit7
#' @method residuals fit8
#' @method residuals fit9
#' @method residuals fit10
#' @method residuals fit11
#' @method residuals fit12
#' @method residuals fit13
#' @method residuals fit14
#' @method residuals cfit1
#' @method residuals cfit2
#' @method residuals cfit3
#' @method residuals cfit4
#' @method residuals cfit5
#' @method residuals cfit6
#' @method residuals cfit7
#' @method residuals cfit8
#' @method residuals cfit9
#' @method residuals cfit10
#' @method residuals cfit11
#' @method residuals cfit12
#' @method residuals cfit13
#' @method residuals cfit14
#'
NULL

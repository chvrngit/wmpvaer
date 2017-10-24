#' additive_eigmECM matrix dataset for wmpvaer
#'
#' The Plastic dataset is from an industrial experiment found in a textbook by Johnson and Wichern,
#' 1992. Data is also found in R documentation as an example of the function 
#' manova.The dependent variables are taken to be "tear","gloss" and "opacity". The
#' independent variables are taken to be "rate" (actually "rate of extrusian" orginally)
#' and "additive". This dataset was segregated by the variable additive into low_additive and 
#' high_additive. Then covariance matricies were computed on the dependent variables. The covariance
#' matrix from low_additive was then mutiplied by the inverse of the covariance matrix of high_additive,
#' resulting in this matrix dataset.
#'
#' @docType data
#' @usage data(additive_eigmECM)
#' @format typeof(additive_eigmECM)=="matrix"
#' @keywords datasets
#' @references Johnson,R.A. and Wichern, D.W. 1992, Applied Multivariate Statistical Analysis  
 "additive_eigmECM"
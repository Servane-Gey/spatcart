#' Load Paracou data set
#'
#' @description Load the Paracou data set as a bivaraite spatial point process.
#'
#' @return the bivaraite spatial point process of class 'ppp'.
#' @export
#'
#' @examples
#' Paracoudata()
Paracoudata = function(){
  ypp = Paracou.pp
  return(ypp)
}

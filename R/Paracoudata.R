#' Load Paracou data set
#'
#' @description Load the Paracou data set as a spatial bivaraite point process.
#'
#' @return the spatial bivaraite point process of class 'ppp'.
#' @export
#'
#' @examples
#' Paracoudata()
Paracoudata = function(){
  ypp = Paracou.pp
  return(ypp)
}


#' Use Single Value Decomposition to reduce information in a matrix by components
#'
#' Take a numerical matrix and return a matrix produced only from the number of components specified
#'
#' @param mat A numeric matrix or a data frame coercible into a numeric matrix
#' @param n An integer number of components to keep. Default is 1
#' @return A matrix or an error if the inputs are illogical
#' @examples
#'
#'# Uses matrix made from the Linux Penguin: Attribution: Larry Ewing <lewing@isc.tamu.edu>
#'\dontrun{
#' image(noisymatrix, col=gray.colors(65536))
#' cleanmatrix<-reduce_components(noisymatrix,50)
#' image(cleanmatrix, col=gray.colors(65536))
#'}
#' @export
reduce_components<-function(mat, n=1){
  # First check the class of the mat object to see if it is a data frame, attempt coercion if so
  if (grepl('data.frame', class(mat))) {
    message('Object is a data frame. Attempting coersion')
    mat<-as.matrix(mat)
  }
  # if the resulting object is not a matrix, stop
  if (class(mat)!='matrix') stop('Input must be a numerical matrix, or coercible to a numerical matrix')
  # Check to see if it is numeric, or if coercion to numeric results in NAs
  isnum<-mean(apply(mat,c(1,2),is.numeric))
  if (isnum<1){
    message('Matrix is not numeric. Attempting to coerce.')
    mat<-suppressWarnings(apply(mat,c(1,2), as.numeric))
    isna<-mean(apply(mat,c(1,2), is.na))
    if(isna>0) stop('Object not coercible to a numeric matrix.')
  }
  matsvd<-svd(mat)
  # get vector of percentage of variance explained by each component
  sv<-(matsvd$d^2)/sum(matsvd$d^2)
  if (length(sv)<n) stop('More components requested than are possible with this matrix size')
  if (n==1){
    m<-matsvd$u[,1]%*%t(matsvd$v[,1])*matsvd$d[1]
  } else {
    m<-matsvd$u[,1:n]%*%diag(matsvd$d[1:n])%*%t(matsvd$v[,1:n])
  }
  return(m)

}

#' Use Single Value Decomposition to reduce information in a matrix by percentage
#'
#' Take a numerical matrix and return a matrix produced only from the percentage of variance explained
#'
#' @param mat A numeric matrix or a data frame coercible into a numeric matrix
#' @param p Numeric value such that 0 < p <= 1. Defaults to .9
#' @return A matrix or an error if the inputs are illogical
#' @examples
#' # Uses matrix made from the Linux Penguin: Attribution: Larry Ewing <lewing@isc.tamu.edu>
#'\dontrun{
#'
#' image(noisymatrix, col=gray.colors(65536))
#' cleanmatrix<-reduce_percentage(noisymatrix,.988)
#' image(cleanmatrix, col=gray.colors(65536))
#'}
#' @export
reduce_percentage<-function(mat, p=.9){
  # is 0 < p <= 1?
  if (!p >0) stop("Invalid value for p. Use a value 0 < p <= 1")
  if (p>1) {
    warning('p > 1. Resetting to 1')
    p<-1
  }
  # First check the class of the mat object to see if it is a data frame, attempt coercion if so
  if (grepl('data.frame', class(mat))) {
    message('Object is a data frame. Attempting coersion')
    mat<-as.matrix(mat)
  }
  # if the resulting object is not a matrix, stop
  if (class(mat)!='matrix') stop('Input must be a numerical matrix, or coercible to a numerical matrix')
  # Check to see if it is numeric, or if coercion to numeric results in NAs
  isnum<-mean(apply(mat,c(1,2),is.numeric))
  if (isnum<1){
    warning('Matrix is not numeric. Attempting to coerce.')
    mat<-suppressWarnings(apply(mat,c(1,2), as.numeric))
    isna<-mean(apply(mat,c(1,2), is.na))
    if(isna>0) stop('Object not coercible to a numeric matrix.')
  }
  matsvd<-svd(mat)
  # get vector of percentage of variance explained by each component
  sv<-(matsvd$d^2)/sum(matsvd$d^2)
  # find the first component which adds up to, or exceeds the percentage specified
  n<-which(cumsum(sv)>=p)[1]
  if(is.na(n)) n<-length(sv)
  if (n==1){
    m<-matsvd$u[,1]%*%t(matsvd$v[,1])*matsvd$d[1]
  } else {
    m<-matsvd$u[,1:n]%*%diag(matsvd$d[1:n])%*%t(matsvd$v[,1:n])
  }
  return(m)
}


#' Use return the number of SVD components available in a matrix
#'
#' Take a numerical matrix and return the number of SVD components available to reduce
#'
#' @param mat A numeric matrix or a data frame coercible into a numeric matrix
#' @return An integer number of components contained in the matrix
#' @examples
#'
#'# Uses matrix made from the Linux Penguin: Attribution: Larry Ewing <lewing@isc.tamu.edu>
#'\dontrun{
#' print(num_components(noisymatrix))
#'}
#'
#' @export
num_components<-function(mat){
  if (grepl('data.frame', class(mat))) {
    message('Object is a data frame. Attempting coersion')
    mat<-as.matrix(mat)
  }
  # if the resulting object is not a matrix, stop
  if (class(mat)!='matrix') stop('Input must be a numerical matrix, or coercible to a numerical matrix')
  # Check to see if it is numeric, or if coercion to numeric results in NAs
  isnum<-mean(apply(mat,c(1,2),is.numeric))
  if (isnum<1){
    warning('Matrix is not numeric. Attempting to coerce.')
    mat<-apply(mat,c(1,2), as.numeric)
    isna<-mean(apply(mat,c(1,2), is.na))
    if(isna>0) stop('Object not coercible to a numeric matrix.')
  }
  matsvd<-svd(mat)
  # get vector of percentage of variance explained by each component
  sv<-(matsvd$d^2)/sum(matsvd$d^2)
  return(length(sv))
}

#' Use Single Value Decomposition to reduce information in a matrix by specific components
#'
#' Take a numerical matrix and return a matrix produced by leaving out specified components
#'
#' @param mat A numeric matrix or a data frame coercible into a numeric matrix
#' @param exclude A vector of integers specifying the components to exclude, defaults to all but the first
#' @return A matrix or an error if the inputs are illogical
#' @examples
#'
#'# Uses matrix made from the Linux Penguin: Attribution: Larry Ewing <lewing@isc.tamu.edu>
#'\dontrun{
#' image(noisymatrix, col=gray.colors(65536))
#' cleanmatrix<-exclude_components(noisymatrix,50:250)
#' image(cleanmatrix, col=gray.colors(65536))
#'}
#' @export
exclude_components<-function(mat, exclude=2:num_components(mat)){
  numComp<-num_components(mat)
  excl<-exclude[exclude %in% 1:numComp]
  excllen<-length(excl)
  if (excllen<1) {
    stop('No valid exclude components selected')
  }
  if (excllen<length(exclude)){
    warning('Some excluded components specified are invalid. Ignoring')
  }

  # First check the class of the mat object to see if it is a data frame, attempt coercion if so
  if (grepl('data.frame', class(mat))) {
    message('Object is a data frame. Attempting coersion')
    mat<-as.matrix(mat)
  }
  # if the resulting object is not a matrix, stop
  if (class(mat)!='matrix') stop('Input must be a numerical matrix, or coercible to a numerical matrix')
  # Check to see if it is numeric, or if coercion to numeric results in NAs
  isnum<-mean(apply(mat,c(1,2),is.numeric))
  if (isnum<1){
    warning('Matrix is not numeric. Attempting to coerce.')
    mat<-apply(mat,c(1,2), as.numeric)
    isna<-mean(apply(mat,c(1,2), is.na))
    if(isna>0) stop('Object not coercible to a numeric matrix.')
  }
  matsvd<-svd(mat)
  sv<-(matsvd$d^2)/sum(matsvd$d^2)

  # get components remaining
  remain<-1:numComp
  remain<-remain[!remain %in% excl]
  n<-length(remain)
  if (n==1){
    m<-matsvd$u[,remain]%*%t(matsvd$v[,remain])*matsvd$d[remain]
  } else {
    m<-matsvd$u[,remain]%*%diag(matsvd$d[remain])%*%t(matsvd$v[,remain])
  }
  return(m)

}


#' Plot the variance explained by each component of the SVD
#'
#' Take a numerical matrix and plot the variance explained by the components
#'
#' @param mat A numeric matrix or a data frame coercible into a numeric matrix
#' @param limit (Optional) A numeric value specifying the maximum percentage to display.
#' @return A double-plot of the
#' @examples
#'
#'# Uses matrix made from the Linux Penguin: Attribution: Larry Ewing <lewing@isc.tamu.edu>
#'\dontrun{
#' plot_explanation(noisymatrix, limit=.99)
#'}
#'
#' @importFrom graphics par plot
#' @export
plot_explanation<-function(mat, limit=NULL){
  #Check to see if limit is there & makes sense
  if (!is.null(limit)) {
    if (limit >1 | limit <=0) {
      warning('Limit > 1 (100%) or limit <= 0, Setting to 1')
      limit<-1
    }
  } else {
    limit<-1
  }
  # First check the class of the mat object to see if it is a data frame, attempt coercion if so
  if (grepl('data.frame', class(mat))) {
    message('Object is a data frame. Attempting coersion')
    mat<-as.matrix(mat)
  }
  # if the resulting object is not a matrix, stop
  if (class(mat)!='matrix') stop('Input must be a numerical matrix, or coercible to a numerical matrix')
  # Check to see if it is numeric, or if coercion to numeric results in NAs
  isnum<-mean(apply(mat,c(1,2),is.numeric))
  if (isnum<1){
    warning('Matrix is not numeric. Attempting to coerce.')
    mat<-apply(mat,c(1,2), as.numeric)
    isna<-mean(apply(mat,c(1,2), is.na))
    if(isna>0) stop('Object not coercible to a numeric matrix.')
  }
  matsvd<-svd(mat)
  # get vector of percentage of variance explained by each component
  percent_explained<-(matsvd$d^2)/sum(matsvd$d^2)
  total=length(percent_explained)
  cumulative<-cumsum(percent_explained)
  n<-which(cumulative>=limit)[1]
  if (is.na(n)) n<-total
  num_components<-1:n
  cumulative<-cumulative[1:n]
  percent_explained<-percent_explained[1:n]
  op<-graphics::par(mfrow=c(2,1), mar=c(1,1,1,1))
  graphics::plot(x=num_components,y=percent_explained, type='b', pch=16, main = sprintf('Percent by component. Total components=%d',total))
  graphics::plot(x=num_components, y=cumulative, type='b',pch=16, main = sprintf('Cumulative Sum Percent by component. Total components=%d',total))
  graphics::par(op)
}

#' A numerical matrix of the Linux logo with noise added
#'
#' A dataset containing the grayscale image of the linux logo with noise added from Wikimedia Commons Attribution: Larry Ewing <lewing@isc.tamu.edu>
#'
#' @format A Matrix with 421 rows and 500 columns:
#
#' @source \url{https://isc.tamu.edu/~lewing/linux/}
"noisymatrix"

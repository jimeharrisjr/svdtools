# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

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
    warning('Matrix is not numeric. Attempting to coerce.')
    mat<-apply(mat,c(1,2), as.numeric)
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

reduce_percentage<-function(mat, p=.9){
  # is 0 < p <= 1?
  if (!p >0 | p>1) stop("Invalid value for p. Use a value 0 < p <= 1")
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
  sv<-(matsvd$d^2)/sum(matsvd$d^2)
  # find the first component which adds up to, or exceeds the percentage specified
  n<-which(cumsum(sv)>=p)[1]
  if (n==1){
    m<-matsvd$u[,1]%*%t(matsvd$v[,1])*matsvd$d[1]
  } else {
    m<-matsvd$u[,1:n]%*%diag(matsvd$d[1:n])%*%t(matsvd$v[,1:n])
  }
  return(m)
}

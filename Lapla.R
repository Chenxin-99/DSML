Lapla_nor <- function(S) {
  rL <- rowSums(S)
  for (i in 1:length(rL)) {
    rL[i] <- rL[i]^(-1/2)
  }
  DS <- diag(rL)
  L <- DS%*%S%*%DS
  return(L)
}


Lapla_matrix_I <- function(S) {
  I <- diag(dim(S)[1])
  L <- I-S
  return(L)
}
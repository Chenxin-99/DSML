PathSim <- function(RGG, len) {
  gene_num <- dim(RGG)[1]
  sim_mat <- diag(gene_num)
  Path_count <- diag(gene_num)
  for (i in 1:len) {
    Path_count <- Path_count %*% RGG
  }
  
  for (i in 1:gene_num) {
    for (j in 1:gene_num) {
      if (Path_count[i,i]+Path_count[j,j]==0) {
        sim_mat[i,j] <- 0
      }
      else {
        sim_mat[i,j] <- 2*Path_count[i,j]/(Path_count[i,i]+Path_count[j,j])
      }
      
    }
  }
  return(sim_mat)
}




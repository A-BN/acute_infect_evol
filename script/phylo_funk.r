# Phylo-funk


make_mut_matrix <- function(mut_table, subclon_names, letter = FALSE){
  # Transforms a mutation table in a presence(T/1) absence(A/0) matrix.
  mut_absence <- ifelse(test = letter, yes = "A", no = 0)
  mut_presence <- ifelse(test = letter, yes = "T", no = 1)
  out_mat <- 
    matrix(nrow = nrow(mut_table), ncol = length(subclon_names), data = mut_absence)
  colnames(out_mat) <- subclon_names
  for (i in 1:nrow(out_mat)) {
      curr_mut <-
      mut_table[i,]
    out_mat[i, str_detect(string = curr_mut$subclon,
                          pattern = subclon_names) ] <- mut_presence
  }
  return(out_mat)
}



my_clearwater <- function(mut_matrix){
  # Watterson Theta's calculation
  # Theta W ~ Number of segregating sites
  k <- nrow(mut_matrix)
  n_indiv <- ncol(mut_matrix)
  an <- sum(1/1:(n_indiv - 1)) # harmonic number!
  # print(paste( "k =" ,k))
  # print(paste( "an =" ,an))
  return(k/an)
}

my_apple_pie <- function(mut_matrix){
  # Theta Pi calculation
  # Theta Pi ~ Number of pairwise differences
  n <- ncol(mut_matrix)
  n_mut <- nrow(mut_matrix)
  my_pair <- 2
  my_mat <- mut_matrix
  my_mat[my_mat == "A"] <- 0
  my_mat[my_mat == "T"] <- 1
  my_mat <-
    matrix(data = sapply(my_mat, FUN = as.numeric),
           nrow = n_mut)
  pairwise_comb <- 
    combn(x = 1:n, 
          m = my_pair, 
          simplify = FALSE)
  pair_dist <-
    sum(sapply(X = pairwise_comb, 
               FUN = function(x) sum(my_mat[, x[[1]] ] != my_mat[,x[[2]] ], simplify = FALSE) ))
  
  denominator <- n*(n - 1) / 2 # equivalent to length(pairwise_com)
  
  theta_pi <- pair_dist / denominator
  return(theta_pi)
}

codon_compare <- function(cod_ref, cod_mut){
  # compare two codons, find the positions that differs 
  # returns a list containing the divergent positions from ref and mut.
  cod_ref_ <- str_split(string = cod_ref, pattern = "")[[1]]
  cod_mut_ <- str_split(string = cod_mut, pattern = "")[[1]]
  diff_base <- !(cod_ref_ == cod_mut_)
  nuc_ref <- cod_ref_[diff_base]
  nuc_mut <- cod_mut_[diff_base]
  return(list(nuc_ref, nuc_mut))
}

# TREE SHAPE CALC R redundancy
tree_shape_R <- function(phylo){
  TL <- sum(phylo$edge.length) #total length of tree
  d <- max(ape::node.depth(phylo)) #Depth of tree 
  n <- phylo$Nnode #+ length(phylo$tip.label)
  R <- 1 - ((TL-d) / (d*n-d))
  print(TL)
  print(d)
  print(n)
  return(R)
}



# codon_comparison <- function(cod1 = "ATG", cod2 = "TTT"){
#   library(stringr)
#   cod1 <- str_split(string = cod1, pattern = "")[[1]]
#   cod2 <- str_split(string = cod2, pattern = "")[[1]]
#   diff_1 <- cod1[!(cod1 == cod2)]
#   diff_2 <- cod2[!(cod1 == cod2)]
#   return(list(diff_1, diff_2))
# }

# ## theta pi test value should be 1.8
# test_mat <-
#   matrix(data = c(0,1,1,1,0,1,0,0,1,1,0,0,0,1,0),
#          nrow = 3)
# my_apple_pie(test_mat)

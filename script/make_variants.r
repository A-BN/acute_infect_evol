# Make variants
#

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(DT)
library(ape)

# setwd("~/Documents/2018_three_patients/script/")
# Some unassigned mutations have been manually added to 2018_all_mutations.tsv
# all_mut <-
  # read_tsv(file = "../data/variants/tables/all_mutations_final.tsv")

source(file = "phylo_funk.r")

###### Patient 11
p_11_mutator <- "4254"
p_11_mut <-
  all_mut %>% 
  filter(patient == "patient_11") %>%
  filter(type == "SNP") %>% 
  filter(subclon != UQ(p_11_mutator)) %>% 
  identity()

p_11_subclon <-
  str_split(string = p_11_mut$subclon, pattern = "\\|") %>% 
  unlist() %>% 
  unique() %>% 
  c("4250") # 4250 is the ref, hence has to be manually added
p_11_subclon <-
  p_11_subclon[p_11_subclon != "4254"]

p_11_matrix <- 
  make_mut_matrix(mut_table = p_11_mut, subclon_names = p_11_subclon, letter = TRUE)

p_11_seq <-
  as.list(apply(X = p_11_matrix, MARGIN = 2, FUN = paste0, collapse = ""))

p_11_align <- 
  ape::as.DNAbin(x = t(p_11_matrix))

p_11_watt <-
  my_clearwater(mut_matrix = p_11_matrix)
  
p_11_pi <-
  my_apple_pie(mut_matrix = p_11_matrix)

p_11_watt <-
  pegas::theta.s(x = nrow(p_11_matrix), n = ncol(p_11_matrix))

p_11_D <-
  pegas::tajima.test(x = p_11_align)


####### Patient 13
p_13_mutator <- "4433"
p_13_mut <-
  all_mut %>% 
  filter(patient == "patient_13") %>%
  filter(type == "SNP") %>% 
  filter(subclon != UQ(p_13_mutator)) %>% 
  identity()
  
p_13_n <- nrow(p_13_mut)

p_13_subclon <-
  str_split(string = p_13_mut$subclon, pattern = "\\|") %>% 
  unlist() %>% 
  unique() %>%
  c("4429") # 4429 is the ref, hence has to be manually added
p_13_subclon <-
  p_13_subclon[p_13_subclon != "4433"]

p_13_matrix <- 
  make_mut_matrix(mut_table = p_13_mut, subclon_names = p_13_subclon, letter = TRUE)

p_13_seq <-
  as.list(apply(X = p_13_matrix, MARGIN = 2, FUN = paste0, collapse = ""))

p_13_align <- 
  ape::as.DNAbin(x = t(p_13_matrix))

p_13_watt <-
my_clearwater(mut_matrix = p_13_matrix)

p_13_pi <-
my_apple_pie(mut_matrix = p_13_matrix)

p_13_watt <-
pegas::theta.s(x = nrow(p_13_matrix), n = ncol(p_13_matrix))

p_13_D <-
pegas::tajima.test(x = p_13_align)


###### Patient 17

p_17_mutator <- "5-36"

p_17_mut <-
  all_mut %>% 
  filter(patient == "patient_17") %>%
  filter(type == "SNP") %>% 
  filter(subclon != UQ(p_17_mutator)) %>% 
  identity()

p_17_n <- nrow(p_17_mut)

p_17_subclon <-
  str_split(string = p_17_mut$subclon, pattern = "\\|") %>% 
  unlist() %>% 
  unique() %>% 
  c("5-35", "5-32") # p5-35 is the ref, hence has to be manually added

p_17_subclon <-
  p_17_subclon[p_17_subclon != "5-36"]

p_17_matrix <- 
  make_mut_matrix(mut_table = p_17_mut, subclon_names = p_17_subclon, letter = TRUE)

p_17_seq <-
  as.list(apply(X = p_17_matrix, MARGIN = 2, FUN = paste0, collapse = ""))

p_17_align <- 
  ape::as.DNAbin(x = t(p_17_matrix))

p_17_watt <-
my_clearwater(mut_matrix = p_17_matrix)

p_17_pi <-
my_apple_pie(mut_matrix = p_17_matrix)

p_17_watt <-
pegas::theta.s(x = nrow(p_17_matrix), n = ncol(p_17_matrix))

p_17_D <-
pegas::tajima.test(x = p_17_align)




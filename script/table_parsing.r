#
#

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(DT)
# library(ape)
setwd("~/longterm_abn/2019_three_patients/script/")


# PATIENT 13
patient_13 <-
  read_tsv(file = "../data/variants/tables/patient_13/01_COMPARE_13.tsv", col_names = TRUE)
# datatable(patient_13)
# How many shared mutations (by at least 2 quasi-clone)
nrow(patient_13)
patient_13 %>%
  select(-title) %>%
  distinct() %>%
  nrow()
table(patient_13$title)

# Create uniq_ID:
patient_13_compact <-
patient_13 %>%
  unite(col = uniq_ID, aa_new_seq:time, remove = FALSE) %>%
  group_by(uniq_ID) %>%
  mutate(subclon = paste0(title, collapse = "|")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-uniq_ID, everything())

# Saving this compact table for manual curation.
write_tsv(x = patient_13_compact, path = "../data/variants/tables/patient_13/02_COMPARE_13_compact.tsv")
##### 
# 
# 
# # The mapping of the reference on itself generated no false positive.
# # We'll have to manually add it later
# patient_13_names <- c(unique(patient_13$title), "patient_13_4429")
# 
# for(namus in patient_13_names){
#   patient_13 <-
#   patient_13 %>%
#     # filter(title == namus) %>%
#     mutate(!!namus := ifelse(test = str_detect(string = namus, subclon), yes = 1, no = 0))
# }
# 
# patient_13_mat <- as.matrix(patient_13[ ,42:48])
# patient_13_mat <- as.matrix(dist(t(patient_13_mat), diag = TRUE, upper = TRUE))
# plot(nj(patient_13_mat), type = "rad", edge.width = 5, cex = 1.2)
# plot(nj(patient_13_mat), type = "fan", edge.width = 5, cex = 1.2)
#-------------------------------------------------------------------------------

# Patient 11
patient_11 <-
  read_tsv(file = "../data/variants/tables/patient_11/01_COMPARE_11.tsv", col_names = TRUE)
# datatable(patient_11)


# How many shared mutations (by at least 2 quasi-clone)
nrow(patient_11)
patient_11 %>%
  select(-title) %>%
  distinct() %>%
  nrow()
table(patient_11$title)

# Create uniq_ID:
patient_11_compact <-
  patient_11 %>%
  unite(col = uniq_ID, aa_new_seq:time, remove = FALSE) %>%
  group_by(uniq_ID) %>%
  mutate(subclon = paste0(title, collapse = "|")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-uniq_ID, everything(), -title)

# Saving this compact table for manual curation.
write_tsv(x = patient_11_compact, path = "../data/variants/tables/patient_11/02_COMPARE_11_compact.tsv")
#-------------------------------------------------------------------------------

# Patient 17
patient_17 <-
  read_tsv(file = "../data/variants/tables/patient_17/01_COMPARE_17.tsv", col_names = TRUE)
# datatable(patient_17)

# How many shared mutations (by at least 2 quasi-clone)
nrow(patient_17)
patient_17 %>%
  select(-title) %>%
  distinct() %>%
  nrow()
table(patient_17$title)

# Create uniq_ID:
patient_17_compact_pre <-
  patient_17 %>%
  unite(col = uniq_ID, aa_new_seq:time, remove = FALSE) %>%
  group_by(uniq_ID) %>%
  mutate(subclon = paste0(title, collapse = "|")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-uniq_ID, everything(), -title)

# Remove bad genes (coverage problem, or mutations called for all subclones)
bad_genes <- # Pay attention to ysgA_2: one bad mut, one good mut.
  c("cpsB", "insEF‑1_1", "FOBKHAMO_00153", "rsxC", "insC‑1_2","gadB_1", "insB‑1_1",
    "FOBKHAMO_00675", "stfR_1", "FOBKHAMO_00708", "FOBKHAMO_00710", "FOBKHAMO_00745",
    "lomR_2", "FOBKHAMO_00861", "ydfR_1", "ydaT_2", "FOBKHAMO_00917", "FOBKHAMO_00949", 
    "FOBKHAMO_01411", "ybfD", "ybfQ", "rhsC", "yhhI_2", "rhsD_3", "copA", "FOBKHAMO_01798",
    "FOBKHAMO_01838", "FOBKHAMO_01839", "yjjW", "mdtM", "FOBKHAMO_02354", "intB_3", 
    "FOBKHAMO_02470", "argE", "rhaA", "FOBKHAMO_02736", "corA_1", "rhsA_2", "bcsZ", 
    "bcsC_1", "yibA_2", "FOBKHAMO_03342", "fadH", "FOBKHAMO_03617", "yeeT_2", "yfjX_2",
    "insC‑1_6", "insCD‑1_6", "FOBKHAMO_04951", "fbaA", "insO‑2_4", "ypjA_1", "tfaE_2",
    "yfeA_1", "doc_1", "FOBKHAMO_04822")

patient_17_compact <-
  patient_17_compact_pre %>%
  filter(!(gene_name %in% bad_genes))
  

# Saving this compact table for manual curation.
write_tsv(x = patient_17_compact, path = "../data/variants/tables/patient_17/02_COMPARE_17_compact.tsv")


# Loading Manually curated tables
curr_11 <- 
  read_tsv(file = "../data/variants/tables/patient_11/04_COMPARE_11_compact_annot_curated.csv", col_names = TRUE)
curr_13 <- 
  read_tsv(file = "../data/variants/tables/patient_13/04_COMPARE_13_compact_annot_curated.csv", col_names = TRUE)
curr_17 <-
  read_tsv(file = "../data/variants/tables/patient_17/04_COMPARE_17_compact_annot_curated.csv", col_names = TRUE)
# the column title is still present in curr_13 and should be removed
curr_13$title <- NULL

curr_all <-
  curr_11 %>%
  rbind(curr_13) %>%
  rbind(curr_17)



curr_all <-
  curr_all %>%
  mutate(patient = str_extract(string = subclon, pattern = "patient_[0-9]{2}")) %>%
  mutate(subclon = str_remove_all(string = subclon, pattern = "patient_[0-9]{2}_")) %>%
  select(patient, subclon, aa_new_seq, aa_position, aa_ref_seq, codon_new_seq, codon_ref_seq, nuc_new_seq, nuc_ref_seq,
         seq_id, position, gene_name, gene_position, gene_product, snp_type, type) %>% 
  mutate(synonymous = ifelse(test = aa_new_seq == aa_ref_seq, yes = TRUE, no = FALSE)) %>%
  identity()
  # write_tsv(x = curr_all, path = "../data/variants/tables/2018_all_mutations.tsv")
  # DT::datatable(data = filter = "top", options = list("pageLength = 100")) %>%
  # DT::saveWidget(file = "~/Desktop/20180908_mutations.html") %>% 


curated_all <- 
  read_delim(file = "../data/variants/tables/2019_all_mutations_CURATED.csv", delim = "\t")

# Loading reference genomes and extract reference nucleotide for each mutation
ref_11 <- Biostrings::readDNAStringSet(filepath = "../data/reference/patient_11/4250_filtered.fasta")
ref_13 <- Biostrings::readDNAStringSet(filepath = "../data/reference/patient_13/4429_filtered.fasta")
ref_17 <- Biostrings::readDNAStringSet(filepath = "../data/reference/patient_17/p5_35_filtered.fasta")

curated_all <-
  curated_all %>% 
  # SNP size set to 1
  mutate(mut_size = if_else(condition = type == "SNP", true = 1, false = mut_size)) %>%
  # Synonymous
  mutate(synonymous = case_when(
    snp_type == "synonymous" ~ TRUE,
    snp_type == "nonsynonymous" ~ FALSE,
    snp_type == "nonsynonymous;nonsynonymous" ~ FALSE,
    aa_new_seq == "*" ~ FALSE,
    TRUE ~ NA)) %>%
  # Phase disruption
  mutate(phase_disrupt = case_when(
    aa_new_seq == "*" ~ TRUE,
    type %in% c("DEL", "INS", "SUB") & mut_size %% 3 != 0 & ! str_detect(string = str_replace_na(string = gene_position, replacement = ""), pattern = "intergenic") ~ TRUE,
    type %in% c("DEL", "INS", "SUB") & mut_size %% 3 == 0 & ! str_detect(string = str_replace_na(string = gene_position, replacement = ""), pattern = "intergenic") ~ FALSE,
    TRUE ~ NA)) %>%
  # Inversion by parcimony
  # group_by(patient) %>%
  # mutate(names_subclon = paste(unique(unlist(str_split(string = subclon, pattern = "\\|"))), collapse = "/")) %>% 
  # mutate(names_subclon = case_when(
  #   patient == "patient_11" ~ str_split(paste(unique(names_subclon), "4250", collapse = "/"), pattern = "/"),
  #   patient == "patient_13" ~ str_split(paste(unique(names_subclon), "4429", collapse = "/"), pattern = "/"),
  #   patient == "patient_17" ~ str_split(paste(unique(names_subclon), "5-35", collapse = "/"), pattern = "/")
  #   )) %>% 
  # mutate(total_subclon = str_count(string = names_subclon, pattern = " ")+ 1 + 1) %>% # plus 1 for the absence of | when unique subclone , +1 for the reference.
  # 
  # mutate(shared_by = 1 + str_count(string = subclon, pattern = "\\|")) %>% # separator counting, not the classiest move
  # mutate(to_invert = ifelse(test = total_subclon %/% shared_by < 2, yes = TRUE, no = FALSE)) %>%
  # ungroup() %>% 
  # mutate(shared_by = 1 + str_count(string = subclon, pattern = "\\|")) %>% # separator counting, not the classiest move
  arrange(patient) %>% 
  identity()

# Add the plasmid vs chromosome information
p_11_plasmid <- c(4,5)
p_13_plasmid <- c(2)
p_17_plasmid <- c(3,4,5,7,8)

curated_all %>% 
  mutate(plasmidic = 
           case_when(
             patient == "patient_11" & seq_id %in% p_11_plasmid ~ TRUE,
             patient == "patient_13" & seq_id %in% p_13_plasmid ~ TRUE,
             patient == "patient_17" & seq_id %in% p_17_plasmid ~ TRUE
           )) %>% 
  identity() -> curated_all


write_tsv(x = curated_all, path = "../data/variants/tables/all_mutations_final.tsv")
  

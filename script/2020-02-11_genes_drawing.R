library(tidyverse)
library(gggenes)
library(ggnewscale)

parse_gff <- function(gff, n_records, contig, from, to, isolate_ref = ""){
  gff_df <- read_tsv(file = gff, 
                     comment = "#", 
                     n_max = n_records, 
                     col_names = FALSE)
  
  gff_df %>% 
    mutate(chr = X1, start = X4, end = X5, strand = X7, name = X9) %>% 
    select(-starts_with("X")) %>% 
    mutate(strand = if_else(condition = strand == "+", 
                            true = 1, 
                            false = -1)) %>% 
    mutate(name = if_else(condition = str_detect(string = name, pattern = "Name="),
                          true = str_replace(string = name, 
                                             pattern = ".*Name=(.*?);.*", 
                                             replacement = "\\1"), 
                          false = "")) %>% 
    filter(chr == contig & start > from & end < to) %>% 
    mutate(isolate = isolate_ref) -> gff_df
  return(gff_df)
}




plot_genes <- function(genes_df, mutation_df, show_legend = F){

  fill_asso <- c("#619CFF", "#00BA38","#F3383C")
  # fill_asso <- c("blue", "green","red")
  names(fill_asso) <- c("SNP", "INS", "DEL")
  genes_plot <-
  ggplot(data = genes_df, mapping = aes(xmin = end,
                                          xmax = start, 
                                          y = isolate, 
                                          fill = name, 
                                          label = name,
                                          forward = -strand)) +
    geom_gene_arrow(arrow_body_height = grid::unit(5, "mm"),
                    arrowhead_width = grid::unit(3, "mm"),
                    arrowhead_height = grid::unit(5, "mm")) +
    # geom_gene_label(grow = F, 
    #                 reflow = F,
    #                 align = "center", height = grid::unit(1000, "mm"),
    #                 ) +
    geom_text(aes(x=(start+end) / 2, y = isolate, label = name), 
              nudge_y = -1.5, angle = 65, fontface = "italic", size = 6) +
    ylab("") +
    xlab("") +
    scale_y_discrete(expand = c(0,3)) +
    scale_x_continuous(breaks = NULL) +
    theme_genes() +
    scale_fill_discrete(guide = F) 
  
  genes_plot <-
  genes_plot +
    ggnewscale::new_scale_fill() +
    geom_gene_arrow(data = mutation_df[abs(mutation_df$start - mutation_df$end) >= 30, ], 
                    aes(xmin = end, 
                        xmax = start, 
                        fill = type),
                        arrowhead_height = grid::unit(0, "mm"),
                        arrowhead_width = grid::unit(0, "mm")) +
    geom_point(data = mutation_df[abs(mutation_df$start - mutation_df$end) <= 30, ],
               aes(x = start, y = isolate, fill = type), 
               colour = rgb(0,0,0,0), 
               size = 5, 
               shape = 23) +
    scale_fill_manual(values = fill_asso) +
    theme(legend.title = element_blank(), axis.text.y = element_text(face = "bold", size = 16, hjust = 0))
  if (! show_legend){
    genes_plot <-
      genes_plot +
      theme(legend.position = "none")
  } else {
    genes_plot <-
      genes_plot +
      theme(legend.position = "bottom")
  }
  return(genes_plot)
} 
plot_genes(genes_df = p_13_genes, mutation_df = p_13_mut, show_legend = F)


p_11_genes <- 
  parse_gff(gff = "data/reference/patient_11/annotation/patient_11.gff", 
            n_records = 4742-14, 
            contig = 6, 
            from = 15e3,
            to =   29e3, 
            isolate_ref = "")

p_11_mut <- 
  read_tsv(file = "data/variants/tables/genes_plot/patient_11_mut_envir.tsv", 
           col_types = cols(isolate = "c")) %>% 
  mutate(type = factor(type, levels = c("SNP", "DEL")))

p_13_genes <- parse_gff(gff = "data/reference/patient_13/annotation/patient_13.gff", 
                        n_records = 4810-3, 
                        contig = 1, 
                        from = 1973800,
                        to =   1980000, 
                        isolate_ref = "")
p_13_mut <- 
  read_tsv(file = "data/variants/tables/genes_plot/patient_13_mut_envir.tsv", 
           col_types = cols(isolate = "c")) %>% 
  mutate(type = factor(type, levels = c("SNP", "DEL", "INS")))

p_17_genes <- parse_gff(gff = "data/reference/patient_17/annotation/patient_17.gff", 
                        n_records = 4968-15,
                        contig = 2, 
                        from = 244000,
                        to =   260000, 
                        isolate_ref = "")
p_17_mut <- 
  read_tsv(file = "data/variants/tables/genes_plot/patient_17_mut_envir.tsv", 
           col_types = cols(isolate = "c")) %>% 
  mutate(type = factor(type, levels = c("SNP", "DEL")))

p_11_plot <- plot_genes(genes_df = p_11_genes, mutation_df = p_11_mut)
p_13_leg_plot <- plot_genes(genes_df = p_13_genes, mutation_df = p_13_mut, show_legend = T)
p_13_plot <- plot_genes(genes_df = p_13_genes, mutation_df = p_13_mut, show_legend = F)
p_17_plot <- plot_genes(genes_df = p_17_genes, mutation_df = p_17_mut)

# plot_grid(
#   plot_genes(genes_df = p_11_genes, mutation_df = p_11_mut),
#   plot_genes(genes_df = p_13_genes, mutation_df = p_13_mut),
#   plot_genes(genes_df = p_17_genes, mutation_df = p_17_mut), 
#   ncol = 1)

cowplot::save_plot(plot = p_11_plot, filename = "docs/Figs/trees/p_11_genes.png")
cowplot::save_plot(plot = p_13_leg_plot, filename = "docs/Figs/trees/p_13_leg_genes.png")
cowplot::save_plot(plot = p_13_plot, filename = "docs/Figs/trees/p_13_genes.png")
cowplot::save_plot(plot = p_17_plot, filename = "docs/Figs/trees/p_17_genes.png")


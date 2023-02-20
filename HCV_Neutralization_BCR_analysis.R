##### Title: Convergent antibody responses are associated with broad neutralization of hepatitis C virus
##### Methods: Performed Bulk BCR-seq on E2-reactive and non-reactive B cells from 5 subjects with highly neutralizing plasma and 5 subjects with poorly neutralizing plasma
##### Software versions: R version 4.1.2 (2021-11-01); R Studio version 2022.12.0+353 (2022.12.0+353)

##### BEFORE ANALYSIS #####
### Load packages
library(immunarch)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(pheatmap)
library(forcats)
library(Hmisc)
library(ggrepel)
library(ff)
library(RecordLinkage)
library(readxl)
library(Biostrings)
library(rstatix)
library(vegan)
library(gridExtra)
library(ggplotify)
library(grid)
library(eulerr)
library(ggpattern)
library(ggtext)

### Functions
## function for italicizing
make_italics <- function(x) {
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}

### Load data 
## IGH
immdataG3 = repLoad("/Volumes/nes002/Projects_Ongoing/HCV_clear_pers/Bulk_BCR/Raw_data/mixcr_HCV/full_deID/IgG")
## IGK
immdataK3 = repLoad("/Volumes/nes002/Projects_Ongoing/HCV_clear_pers/Bulk_BCR/Raw_data/mixcr_HCV/full_deID/IgK")
## IGL
immdataL3 = repLoad("/Volumes/nes002/Projects_Ongoing/HCV_clear_pers/Bulk_BCR/Raw_data/mixcr_HCV/full_deID/IgL")

## Sample names
samples <- c("C117n", "C117p", "C172n", "C172p", "C176n", "C176p", "C429n", "C429p", "C48n", "C48p", "P111n", "P111p", "P157n", "P157p", "P49n", "P49p", "P53n", "P53p", "P54n", "P54p" )


##### ANALYSIS #####
##### GENERATE CLONOTYPE DATAFRAMES
#### IgG
### clonotype data
## all data
all_dfs_VJ <- list()
for(i in 1:20){
  df_sample <- immdataG3$data[[i]]
  df_sample$V.name <- gsub("\\*.*", "", df_sample$V.name)
  df_sample$J.name <- gsub("\\*.*", "", df_sample$J.name)
  df_agg <- aggregate(Clones ~ CDR3.aa + V.name + J.name, df_sample, sum)
  df_named <- cbind.data.frame(Sample=rep(samples[i], nrow(df_agg)), df_agg, stringsAsFactors=FALSE)
  df_named$V.name[df_named$V.name == "IGHV1-69" | df_named$V.name == "IGHV1-69D" | 
                    df_named$V.name == "IGHV1-69, IGHV1-69D" |
                    df_named$V.name == "IGHV1-69-2"] <- "IGHV1-69"
  all_dfs_VJ[[i]] <- df_named
}

cdr3G_seq_VJ_raw <- all_dfs_VJ %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

## downsampled to get equal number of clonotypes per subject
all_dfs_VJds <- list()
for(i in 1:20){
  full_df <- all_dfs_VJ[[i]]
  set.seed(20)
  ifelse(grepl("p", full_df$Sample[1]), df_ds <- sample_n(full_df, 537, replace=FALSE), df_ds <- sample_n(full_df, 220, replace=FALSE))
  all_dfs_VJds[[i]] <- df_ds
}

cdr3G_seq_VJ <- all_dfs_VJds %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

### Full sequence data
## all data
all_dfs_VJseq <- list()
for(i in 1:20){
  df_sample <- immdataG3$data[[i]]
  df_sample$V.name <- gsub("\\*.*", "", df_sample$V.name)
  df_sample$J.name <- gsub("\\*.*", "", df_sample$J.name)
  df_agg <- aggregate(Clones ~ CDR3.aa + V.name + J.name + Sequence, df_sample, sum)
  df_named <- cbind.data.frame(Sample=rep(samples[i], nrow(df_agg)), df_agg, stringsAsFactors=FALSE) %>%
    ungroup() %>%
    mutate(clon_add = ifelse(i==1, 0, max(all_dfs_VJseq[[i-1]]$clon_num))) %>%
    mutate(clon_num = 1:n()) %>%
    mutate(clon_num = clon_num + clon_add) %>%
    select(-clon_add)
  df_named$V.name[df_named$V.name == "IGHV1-69" | df_named$V.name == "IGHV1-69D" | 
                    df_named$V.name == "IGHV1-69, IGHV1-69D" |
                    df_named$V.name == "IGHV1-69-2"] <- "IGHV1-69"
  all_dfs_VJseq[[i]] <- df_named
}

cdr3G_seq_VJseq_raw <- all_dfs_VJseq %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name", "Sequence", "clon_num")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

## downsampled to get equal number of clonotypes per subject
all_dfs_VJseq_ds <- list()
for(i in 1:20){
  full_df <-all_dfs_VJseq[[i]]
  set.seed(25)
  ifelse(grepl("p", full_df$Sample[1]), df_ds <- sample_n(full_df, 579, replace=FALSE), df_ds <- sample_n(full_df, 227, replace=FALSE))
  all_dfs_VJseq_ds[[i]] <- df_ds
}

cdr3G_seq_VJseq <- all_dfs_VJseq_ds %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name", "Sequence", "clon_num")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

#### IGK
### clonotype data
## all data
all_dfsK_VJ <- list()
for(i in 1:20){
  df_sample <- immdataK3$data[[i]]
  df_sample$V.name <- gsub("\\*.*", "", df_sample$V.name)
  df_sample$J.name <- gsub("\\*.*", "", df_sample$J.name)
  df_agg <- aggregate(Clones ~ CDR3.aa + V.name + J.name, df_sample, sum)
  df_named <- cbind.data.frame(Sample=rep(samples[i], nrow(df_agg)), df_agg, stringsAsFactors=FALSE)
  all_dfsK_VJ[[i]] <- df_named
}

cdr3K_seq_VJ_raw <- all_dfsK_VJ %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

## downsampled to get equal number of clonotypes per subject
all_dfsK_VJds <- list()
for(i in 1:20){
  full_df <- all_dfsK_VJ[[i]]
  set.seed(20)
  ifelse(grepl("p", full_df$Sample[1]), df_ds <- sample_n(full_df, 601, replace=FALSE), df_ds <- sample_n(full_df, 250, replace=FALSE))
  all_dfsK_VJds[[i]] <- df_ds
}

cdr3K_seq_VJ <- all_dfsK_VJds %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

### Full sequence data
## all data
all_dfsK_VJseq <- list()
for(i in 1:20){
  df_sample <- immdataK3$data[[i]]
  df_sample$V.name <- gsub("\\*.*", "", df_sample$V.name)
  df_sample$J.name <- gsub("\\*.*", "", df_sample$J.name)
  df_agg <- aggregate(Clones ~ CDR3.aa + V.name + J.name + Sequence, df_sample, sum)
  df_named <- cbind.data.frame(Sample=rep(samples[i], nrow(df_agg)), df_agg, stringsAsFactors=FALSE) %>%
    ungroup() %>%
    mutate(clon_add = ifelse(i==1, 0, max(all_dfsK_VJseq[[i-1]]$clon_num))) %>%
    mutate(clon_num = 1:n()) %>%
    mutate(clon_num = clon_num + clon_add) %>%
    select(-clon_add)
  all_dfsK_VJseq[[i]] <- df_named
}

cdr3K_seq_VJseq_raw <- all_dfsK_VJseq %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name", "Sequence", "clon_num")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

## downsampled to get equal number of clonotypes per subject
all_dfsK_VJseq_ds <- list()
for(i in 1:20){
  full_df <-all_dfsK_VJseq[[i]]
  set.seed(25)
  ifelse(grepl("p", full_df$Sample[1]), df_ds <- sample_n(full_df, 579, replace=FALSE), df_ds <- sample_n(full_df, 227, replace=FALSE))
  all_dfsK_VJseq_ds[[i]] <- df_ds
}

cdr3K_seq_VJseq <- all_dfsK_VJseq_ds %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name", "Sequence", "clon_num")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

#### IGL
### clonotype data
## all data
all_dfsL_VJ <- list()
for(i in 1:20){
  df_sample <- immdataL3$data[[i]]
  df_sample$V.name <- gsub("\\*.*", "", df_sample$V.name)
  df_sample$J.name <- gsub("\\*.*", "", df_sample$J.name)
  df_agg <- aggregate(Clones ~ CDR3.aa + V.name + J.name, df_sample, sum)
  df_named <- cbind.data.frame(Sample=rep(samples[i], nrow(df_agg)), df_agg, stringsAsFactors=FALSE)
  all_dfsL_VJ[[i]] <- df_named
}

cdr3L_seq_VJ_raw <- all_dfsL_VJ %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

## downsampled to get equal number of clonotypes per subject
all_dfsL_VJds <- list()
for(i in 1:20){
  full_df <- all_dfsL_VJ[[i]]
  set.seed(20)
  ifelse(grepl("p", full_df$Sample[1]), df_ds <- sample_n(full_df, 334, replace=FALSE), df_ds <- sample_n(full_df, 197, replace=FALSE))
  all_dfsL_VJds[[i]] <- df_ds
}

cdr3L_seq_VJ <- all_dfsL_VJds %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

### Full sequence data
## all data
all_dfsL_VJseq <- list()
for(i in 1:20){
  df_sample <- immdataL3$data[[i]]
  df_sample$V.name <- gsub("\\*.*", "", df_sample$V.name)
  df_sample$J.name <- gsub("\\*.*", "", df_sample$J.name)
  df_agg <- aggregate(Clones ~ CDR3.aa + V.name + J.name + Sequence, df_sample, sum)
  df_named <- cbind.data.frame(Sample=rep(samples[i], nrow(df_agg)), df_agg, stringsAsFactors=FALSE) %>%
    ungroup() %>%
    mutate(clon_add = ifelse(i==1, 0, max(all_dfsL_VJseq[[i-1]]$clon_num))) %>%
    mutate(clon_num = 1:n()) %>%
    mutate(clon_num = clon_num + clon_add) %>%
    select(-clon_add)
  all_dfsL_VJseq[[i]] <- df_named
}

cdr3L_seq_VJseq_raw <- all_dfsL_VJseq %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name", "Sequence", "clon_num")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

## downsampled to get equal number of clonotypes per subject
all_dfsL_VJseq_ds <- list()
for(i in 1:20){
  full_df <-all_dfsL_VJseq[[i]]
  set.seed(25)
  ifelse(grepl("p", full_df$Sample[1]), df_ds <- sample_n(full_df, 346, replace=FALSE), df_ds <- sample_n(full_df, 206, replace=FALSE))
  all_dfsL_VJseq_ds[[i]] <- df_ds
}

cdr3L_seq_VJseq <- all_dfsL_VJseq_ds %>%
  purrr::reduce(full_join, by = c("CDR3.aa", "V.name", "J.name", "Sequence", "clon_num")) %>%
  rename_at(vars(starts_with("Clones.")), ~ samples) %>%
  select(-contains("Sample")) %>%
  replace(is.na(.), 0)

#### IGKL
cdr3L_seq_VJseq2 <- cdr3L_seq_VJseq %>%
  mutate(clon_num = 31427+clon_num)

cdr3L_seq_VJseq_raw2 <- cdr3L_seq_VJseq_raw %>%
  mutate(clon_num = 31427+clon_num)

cdr3KL_seq_VJseq <- rbind.data.frame(cdr3K_seq_VJseq, cdr3L_seq_VJseq2, stringsAsFactors = FALSE)
cdr3KL_seq_VJseq_raw <- rbind.data.frame(cdr3K_seq_VJseq_raw, cdr3L_seq_VJseq_raw2, stringsAsFactors = FALSE)


##### V-GENE USAGE
#### IGH
### Quantify V-gene usage for each group
cdr3G_VJ_grp <- cdr3G_seq_VJ %>%
  rowwise() %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  mutate(E2neg = sum(C48n, C117n, C172n, C176n, C429n, P49n, P157n, P53n, P54n, P111n, na.rm = T))

hneut_pos_VgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, hneut_pos) %>%
  filter(hneut_pos > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "hneut_pos") %>%
  mutate(V.name=str_remove(V.name, "IGHV")) 

hneut_neg_VgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, hneut_neg) %>%
  filter(hneut_neg > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "hneut_neg") %>%
  mutate(V.name=str_remove(V.name, "IGHV"))

lneut_pos_VgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, lneut_pos) %>%
  filter(lneut_pos > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "lneut_pos") %>%
  mutate(V.name=str_remove(V.name, "IGHV"))

lneut_neg_VgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, lneut_neg) %>%
  filter(lneut_neg > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "lneut_neg") %>%
  mutate(V.name=str_remove(V.name, "IGHV"))

### Statistical analysis of differences between groups
VstatG <- rbind.data.frame(hneut_pos_VgenesG, hneut_neg_VgenesG, lneut_pos_VgenesG, lneut_neg_VgenesG) %>%
  complete(V.name, Group, fill = list(V_count = 0)) %>%
  dplyr::select(V.name, Group, V_count) %>%
  group_by(Group) %>%
  mutate(Total = sum(V_count)) 

listVG <- unique(VstatG$V.name)

Vgenes_statsG <- data.frame()
for(i in 1:length(listVG)){
  VG <- filter(VstatG, V.name == listVG[i])
  hnp_hnn <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[1], VG$Total[2]-VG$V_count[2], VG$Total[1]-VG$V_count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[4], VG$Total[2]-VG$V_count[2], VG$Total[4]-VG$V_count[4]), ncol=2))$p.value
  lnp_lnn <- fisher.test(matrix(c(VG$V_count[4], VG$V_count[3], VG$Total[4]-VG$V_count[4], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  hnp_lnn <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[3], VG$Total[2]-VG$V_count[2], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  lnp_hnn <- fisher.test(matrix(c(VG$V_count[4], VG$V_count[1], VG$Total[4]-VG$V_count[4], VG$Total[1]-VG$V_count[1]), ncol=2))$p.value
  hnn_lnn <- fisher.test(matrix(c(VG$V_count[1], VG$V_count[3], VG$Total[1]-VG$V_count[1], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  VG_df <- cbind.data.frame(V.name = listVG[i], hnp_hnn, hnp_lnp, lnp_lnn, hnp_lnn, lnp_hnn, hnn_lnn)
  Vgenes_statsG <- rbind.data.frame(Vgenes_statsG, VG_df, stringsAsFactors = FALSE)
}

Vgenes_sigG <-data.frame(Vgenes_statsG$V.name, matrix(p.adjust(as.vector(as.matrix(Vgenes_statsG[,-1])), method="BH"),ncol=6))
colnames(Vgenes_sigG) <- colnames(Vgenes_statsG)

Vgenes_sig_filtG <- Vgenes_sigG %>%
  filter_all(any_vars(. <0.05)) %>%
  pivot_longer(., c("hnp_hnn", "hnp_lnp", "lnp_lnn", "hnp_lnn", "lnp_hnn", "hnn_lnn"),
               names_to = "Comp",
               values_to = "p.adj") %>%
  mutate(.y = "Prop") %>%
  mutate(group1 = gsub("_.*","", Comp)) %>%
  mutate(group2 = gsub("..._","", Comp)) %>%
  mutate(group1 = ifelse(group1=="hnp", "hneut_pos", 
                         ifelse(group1=="hnn", "hneut_neg",
                                ifelse(group1=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(group2 = ifelse(group2=="hnp", "hneut_pos", 
                         ifelse(group2=="hnn", "hneut_neg",
                                ifelse(group2=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(p.adj.signif = ifelse(p.adj < 0.0001, "****", 
                               ifelse(p.adj < 0.001, "***",
                                      ifelse(p.adj < 0.01, "**",
                                             ifelse(p.adj <0.05, "*", "ns"))))) %>%
  select(V.name, .y, group1, group2, p.adj, p.adj.signif) 

### Plots
## Heatmap of V-gene expression for all 4 groups
VgenesG <- rbind.data.frame(hneut_pos_VgenesG, hneut_neg_VgenesG, lneut_pos_VgenesG, lneut_neg_VgenesG) %>%
  complete(V.name, Group, fill = list(V_count = 0)) %>%
  dplyr::select(V.name, Group, V_count, Prop) %>%
  replace(is.na(.), 0) %>%
  reshape2::dcast(., Group ~ V.name, value.var = "Prop")
row.names(VgenesG) <- c("Hneut E2-", "Hneut E2+",  "Lneut E2-", "Lneut E2+" )

VgenesG_hm <- VgenesG[c("Hneut E2+", "Hneut E2-", "Lneut E2+", "Lneut E2-"),-1]

Vgenes_G_annot <- data.frame(row.names = rownames(VgenesG_hm), Neut = c("High", "High", "Low", "Low"))
Vgenes_G_annot_cols <- list(Neut = c(High = "#009999", Low= "#b66dff"))

reorder_VG <- data.frame(vname = colnames(VgenesG_hm)) %>%
  mutate(numbering = sub("[0-9]-", "", vname)) %>%
  mutate(numbering = sub("-[0-9]", "", numbering)) %>%
  mutate(fam = sub("-.*", "", vname)) %>%
  group_by(fam) %>%
  arrange(as.numeric(numbering), .by_group = TRUE)

VgenesG_hm_arr <- t(VgenesG_hm)[reorder_VG$vname,]

pheatmap(VgenesG_hm_arr, scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 13, labels_row = make_italics(rownames(VgenesG_hm_arr)), labels_col = c("E2+", "E2-", "E2+", "E2-"), annotation_col = Vgenes_G_annot, annotation_names_col = TRUE, annotation_colors = Vgenes_G_annot_cols, annotation_legend = FALSE, main = "")
grid.text("High", x=0.402, y= 0.976, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("Low", x=0.458, y=0.976, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("A", x=0.35, y=0.98, gp=gpar(fontsize=30))

## Volcano plot comparing V-gene expression between high and low neut E2-reactive B cells
FC_preG <- cbind.data.frame(Vgenes = Vgenes_sigG$V.name, HNP =t(VgenesG_hm)[,1], LNP = t(VgenesG_hm)[,3]) 
FC_preG[FC_preG == 0] <- 0.0002

FC_G <- FC_preG %>%
  mutate(FC = log2(HNP/LNP))

volcG <- data.frame(Vgenes = Vgenes_sigG$V.name,
                    Log2foldchange = FC_G$FC,
                    Pvalues = Vgenes_sigG$hnp_lnp) %>%
  mutate(Change = ifelse(Log2foldchange >= 0.6 & Pvalues < 0.05, "inc", 
                         ifelse(Log2foldchange <= -0.6 & Pvalues < 0.05, "dec", "unchanged"))) %>%
  mutate(Label = ifelse(Change == "inc" | Change == "dec", Vgenes, NA)) 

volcG$Change = factor(volcG$Change, levels=c("inc", "dec", "unchanged"))

volcG_p <- ggplot(data=volcG , aes(x=Log2foldchange, y=-log10(Pvalues), color=Change, label=Label))+
  geom_point(size=4) +
  scale_color_manual(labels = c("High Neut", "Low Neut", "No Change"), values=c("#009999", "#b66dff", "dark gray")) +
  geom_text_repel(size=7, show.legend = FALSE, fontface="bold.italic", color="black", max.overlaps=50) +
  geom_vline(xintercept=c(-0.6, 0.6), color="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), color="black", linetype="dashed") +
  scale_x_continuous(breaks = seq(from = -5, to = 5, by = 1), limits = c(-5.25,5.25))+
  scale_y_continuous(breaks = seq(from = -5, to = 30, by = 5), limits = c(-6,30))+
  ggtitle("A")+
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size=22),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, color = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        legend.key = element_rect(fill = "white", color = NA),
        legend.position = "right",
        plot.title = element_text(size = 45, hjust = -0.135, vjust = -1),
        aspect.ratio=0.9)+
  annotate("rect", xmin = -0.6, xmax = 0.6, ymin = -Inf, ymax = Inf,alpha = 0.1)+
  annotate("rect", xmin = -Inf, xmax = -0.6, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  annotate("rect", xmin = 0.6, xmax = Inf, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  labs(x= expression("Log"[2]~"Fold Change"),
       y= expression("-Log"[10]~"P-value"))
volcG_p

#### IGKL
### Quantify V-gene usage for each group
## IGK
cdr3K_VJ_grp <- cdr3K_seq_VJ %>%
  rowwise() %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  mutate(E2neg = sum(C48n, C117n, C172n, C176n, C429n, P49n, P157n, P53n, P54n, P111n, na.rm = T))

# IGL
cdr3L_VJ_grp <- cdr3L_seq_VJ %>%
  rowwise() %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  mutate(E2neg = sum(C48n, C117n, C172n, C176n, C429n, P49n, P157n, P53n, P54n, P111n, na.rm = T))

# Combine IGK and IGL
cdr3KL_VJ_grp <- rbind.data.frame(cdr3K_VJ_grp, cdr3L_VJ_grp, stringsAsFactors = FALSE)

hneut_pos_VgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, hneut_pos) %>%
  filter(hneut_pos > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "hneut_pos") 

hneut_neg_VgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, hneut_neg) %>%
  filter(hneut_neg > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "hneut_neg") 

lneut_pos_VgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, lneut_pos) %>%
  filter(lneut_pos > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "lneut_pos") 

lneut_neg_VgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, lneut_neg) %>%
  filter(lneut_neg > 0) %>%
  filter(., !grepl(",", V.name)) %>%
  group_by(V.name) %>%
  dplyr::summarise(V_count = n()) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total) %>%
  mutate(Group = "lneut_neg") 


### Statistical analysis of differences between groups
VstatKL <- rbind.data.frame(hneut_pos_VgenesKL, hneut_neg_VgenesKL, lneut_pos_VgenesKL, lneut_neg_VgenesKL) %>%
  complete(V.name, Group, fill = list(V_count = 0)) %>%
  dplyr::select(V.name, Group, V_count) %>%
  group_by(Group) %>%
  mutate(Total = sum(V_count)) 

listVKL <- unique(VstatKL$V.name)

Vgenes_statsKL <- data.frame()
for(i in 1:length(listVKL)){
  VG <- filter(VstatKL, V.name == listVKL[i])
  hnp_hnn <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[1], VG$Total[2]-VG$V_count[2], VG$Total[1]-VG$V_count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[4], VG$Total[2]-VG$V_count[2], VG$Total[4]-VG$V_count[4]), ncol=2))$p.value
  lnp_lnn <- fisher.test(matrix(c(VG$V_count[4], VG$V_count[3], VG$Total[4]-VG$V_count[4], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  hnp_lnn <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[3], VG$Total[2]-VG$V_count[2], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  lnp_hnn <- fisher.test(matrix(c(VG$V_count[4], VG$V_count[1], VG$Total[4]-VG$V_count[4], VG$Total[1]-VG$V_count[1]), ncol=2))$p.value
  hnn_lnn <- fisher.test(matrix(c(VG$V_count[1], VG$V_count[3], VG$Total[1]-VG$V_count[1], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  VG_df <- cbind.data.frame(V.name = listVKL[i], hnp_hnn, hnp_lnp, lnp_lnn, hnp_lnn, lnp_hnn, hnn_lnn)
  Vgenes_statsKL <- rbind.data.frame(Vgenes_statsKL, VG_df, stringsAsFactors = FALSE)
}

Vgenes_sigKL <-data.frame(Vgenes_statsKL$V.name, matrix(p.adjust(as.vector(as.matrix(Vgenes_statsKL[,-1])), method='BH'),ncol=6))
colnames(Vgenes_sigKL) <- colnames(Vgenes_statsKL)

Vgenes_sig_filtKL <- Vgenes_sigKL %>%
  filter_all(any_vars(. <0.05)) %>%
  pivot_longer(., c("hnp_hnn", "hnp_lnp", "lnp_lnn", "hnp_lnn", "lnp_hnn", "hnn_lnn"),
               names_to = "Comp",
               values_to = "p.adj") %>%
  mutate(.y = "Prop") %>%
  mutate(group1 = gsub("_.*","", Comp)) %>%
  mutate(group2 = gsub("..._","", Comp)) %>%
  mutate(group1 = ifelse(group1=="hnp", "hneut_pos", 
                         ifelse(group1=="hnn", "hneut_neg",
                                ifelse(group1=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(group2 = ifelse(group2=="hnp", "hneut_pos", 
                         ifelse(group2=="hnn", "hneut_neg",
                                ifelse(group2=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(p.adj.signif = ifelse(p.adj < 0.0001, "****", 
                               ifelse(p.adj < 0.001, "***",
                                      ifelse(p.adj < 0.01, "**",
                                             ifelse(p.adj <0.05, "*", "ns"))))) %>%
  select(V.name, .y, group1, group2, p.adj, p.adj.signif) 

### Plots
## Heatmap of V-gene expression for all 4 groups
VgenesKL <- rbind.data.frame(hneut_pos_VgenesKL, hneut_neg_VgenesKL, lneut_pos_VgenesKL, lneut_neg_VgenesKL) %>%
  complete(V.name, Group, fill = list(V_count = 0)) %>%
  dplyr::select(V.name, Group, V_count, Prop) %>%
  replace(is.na(.), 0) %>%
  reshape2::dcast(., Group ~ V.name, value.var = "Prop")
row.names(VgenesKL) <- c("Hneut E2-", "Hneut E2+",  "Lneut E2-", "Lneut E2+" )

VgenesKL_hm <- VgenesKL[c("Hneut E2+", "Hneut E2-", "Lneut E2+", "Lneut E2-"),-1]

Vgenes_KL_annot <- data.frame(row.names = rownames(VgenesKL_hm), Neut = c("High Neut.", "High Neut.", "Low Neut.", "Low Neut."))
Vgenes_KL_annot_cols <- list(Neut = c("High Neut."= "#009999", "Low Neut."= "#b66dff"))

Vgenes_KL_colnames <- gsub("IGKV", "", colnames(VgenesKL_hm))
Vgenes_KL_colnames <- gsub("IGLV", "", Vgenes_KL_colnames)

reorder_VKL <- data.frame(vname_full = colnames(VgenesKL_hm), vname = sub("IG[A-Z]V", "", colnames(VgenesKL_hm))) %>%
  mutate(numbering = sub(".{1,2}-", "", vname)) %>%
  mutate(fam = sub("-.*", "", vname)) %>%
  mutate(chain = substr(vname_full, 1,4)) %>%
  group_by(chain, fam) %>%
  arrange(as.numeric(numbering), .by_group = TRUE)

VgenesKL_hm_arr <- t(VgenesKL_hm)[reorder_VKL$vname_full,]

#K only
pheatmap(VgenesKL_hm_arr[1:47,], scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 13, labels_row = make_italics(reorder_VKL$vname[1:47]), labels_col = c("E2+", "E2-", "E2+", "E2-"), annotation_col = Vgenes_KL_annot, annotation_names_col = TRUE, annotation_colors = Vgenes_KL_annot_cols, annotation_legend = FALSE)
grid.text("High", x=0.407, y= 0.961, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("Low", x=0.463, y=0.961, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("B", x=0.36, y=0.965, gp=gpar(fontsize=30))

#L only
pheatmap(VgenesKL_hm_arr[48:86,], scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 13, labels_row = make_italics(reorder_VKL$vname[48:86]), labels_col = c("E2+", "E2-", "E2+", "E2-"), annotation_col = Vgenes_KL_annot, annotation_names_col = TRUE, annotation_colors = Vgenes_KL_annot_cols, annotation_legend = TRUE)
grid.text("High", x=0.407, y= 0.925, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("Low", x=0.463, y=0.925, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("C", x=0.36, y=0.93, gp=gpar(fontsize=30))


## volcano
FC_preKL <- cbind.data.frame(Vgenes = Vgenes_statsKL$V.name, HNP =t(VgenesKL_hm)[,1], LNP = t(VgenesKL_hm)[,3]) 
FC_preKL[FC_preKL == 0] <- 0.0001

FC_KL <- FC_preKL %>%
  mutate(FC = log2(HNP/LNP))

volcKL <- data.frame(Vgenes = Vgenes_statsKL$V.name,
                     Log2foldchange = FC_KL$FC,
                     Pvalues = Vgenes_sigKL$hnp_lnp) %>%
  mutate(Change = ifelse(Log2foldchange >= 0.6 & Pvalues < 0.05, "inc", 
                         ifelse(Log2foldchange <= -0.6 & Pvalues < 0.05, "dec", "unchanged"))) %>%
  mutate(Label = ifelse(Change == "inc" | Change == "dec", Vgenes, NA)) 


volcKL$Change = factor(volcKL$Change, levels=c("inc", "dec", "unchanged"))

volcKL_p <- ggplot(data=volcKL , aes(x=Log2foldchange, y=-log10(Pvalues), color=Change, label=Label))+
  geom_point(size=4) +
  scale_color_manual(labels = c("High Neut", "Low Neut", "No Change"), values=c( "#009999", "#b66dff", "dark gray")) +
  geom_text_repel(size=7, show.legend = FALSE, fontface="bold.italic", color="black", max.overlaps=50) +
  geom_vline(xintercept=c(-0.6, 0.6), color="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), color="black", linetype="dashed") +
  scale_x_continuous(breaks = seq(from = -5, to = 5, by = 1), limits = c(-5,5))+
  scale_y_continuous(breaks = seq(from = -2, to = 6, by = 2), limits = c(-2,6))+
  ggtitle(" C")+
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, color = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        legend.key = element_rect(fill = "white", color = NA),
        legend.position = "right",
        plot.title = element_text(size = 45, hjust = -0.15, vjust = -1),
        aspect.ratio=0.9)+
  annotate("rect", xmin = -0.6, xmax = 0.6, ymin = -Inf, ymax = Inf,alpha = 0.1)+
  annotate("rect", xmin = -Inf, xmax = -0.6, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  annotate("rect", xmin = 0.6, xmax = Inf, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  labs(x= expression("Log"[2]~"Fold Change"),
       y= expression("-Log"[10]~"P-value"))
volcKL_p


##### VJ-GENE ANALYSIS
#### IGH
### Quantify VJ-gene usage for each group
hneut_pos_VJgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, hneut_pos) %>%
  filter(hneut_pos > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "hneut_pos") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IGH"))

hneut_neg_VJgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, hneut_neg) %>%
  filter(hneut_neg > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "hneut_neg") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IGH"))

lneut_pos_VJgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, lneut_pos) %>%
  filter(lneut_pos > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "lneut_pos") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IGH"))

lneut_neg_VJgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, lneut_neg) %>%
  filter(lneut_neg > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "lneut_neg") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IGH"))

E2n_VJgenesG <- cdr3G_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, E2neg) %>%
  filter(E2neg > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "E2neg") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IGH"))

### Statistical analysis of differences between groups
VJstatG <- rbind.data.frame(hneut_pos_VJgenesG, hneut_neg_VJgenesG, lneut_pos_VJgenesG, lneut_neg_VJgenesG) %>%
  complete(VJ.name, Group, fill = list(VJ_count = 0)) %>%
  dplyr::select(VJ.name, Group, VJ_count) %>%
  group_by(Group) %>%
  mutate(Total = sum(VJ_count)) 

listVJG <- unique(VJstatG$VJ.name)

VJgenes_statsG <- data.frame()
for(i in 1:length(listVJG)){
  VJG <- filter(VJstatG, VJ.name == listVJG[i])
  hnp_hnn <- fisher.test(matrix(c(VJG$VJ_count[2], VJG$VJ_count[1], VJG$Total[2]-VJG$VJ_count[2], VJG$Total[1]-VJG$VJ_count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(VJG$VJ_count[2], VJG$VJ_count[4], VJG$Total[2]-VJG$VJ_count[2], VJG$Total[4]-VJG$VJ_count[4]), ncol=2))$p.value
  lnp_lnn <- fisher.test(matrix(c(VJG$VJ_count[4], VJG$VJ_count[3], VJG$Total[4]-VJG$VJ_count[4], VJG$Total[3]-VJG$VJ_count[3]), ncol=2))$p.value
  hnp_lnn <- fisher.test(matrix(c(VJG$VJ_count[2], VJG$VJ_count[3], VJG$Total[2]-VJG$VJ_count[2], VJG$Total[3]-VJG$VJ_count[3]), ncol=2))$p.value
  lnp_hnn <- fisher.test(matrix(c(VJG$VJ_count[4], VJG$VJ_count[1], VJG$Total[4]-VJG$VJ_count[4], VJG$Total[1]-VJG$VJ_count[1]), ncol=2))$p.value
  hnn_lnn <- fisher.test(matrix(c(VJG$VJ_count[1], VJG$VJ_count[3], VJG$Total[1]-VJG$VJ_count[1], VJG$Total[3]-VJG$VJ_count[3]), ncol=2))$p.value
  VJG_df <- cbind.data.frame(VJ.name = listVJG[i], hnp_hnn, hnp_lnp, lnp_lnn, hnp_lnn, lnp_hnn, hnn_lnn)
  VJgenes_statsG <- rbind.data.frame(VJgenes_statsG, VJG_df, stringsAsFactors = FALSE)
}

VJgenes_sigG <-data.frame(VJgenes_statsG$VJ.name, matrix(p.adjust(as.vector(as.matrix(VJgenes_statsG[,-1])), method='BH'),ncol=6))
colnames(VJgenes_sigG) <- colnames(VJgenes_statsG)

VJgenes_sig_filtG <- VJgenes_sigG %>%
  filter_all(any_vars(. <0.05)) %>%
  pivot_longer(., c("hnp_hnn", "hnp_lnp", "lnp_lnn", "hnp_lnn", "lnp_hnn", "hnn_lnn"),
               names_to = "Comp",
               values_to = "p.adj") %>%
  mutate(.y = "Prop") %>%
  mutate(group1 = gsub("_.*","", Comp)) %>%
  mutate(group2 = gsub("..._","", Comp)) %>%
  mutate(group1 = ifelse(group1=="hnp", "hneut_pos", 
                         ifelse(group1=="hnn", "hneut_neg",
                                ifelse(group1=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(group2 = ifelse(group2=="hnp", "hneut_pos", 
                         ifelse(group2=="hnn", "hneut_neg",
                                ifelse(group2=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(p.adj.signif = ifelse(p.adj < 0.0001, "****", 
                               ifelse(p.adj < 0.001, "***",
                                      ifelse(p.adj < 0.01, "**",
                                             ifelse(p.adj <0.05, "*", "ns"))))) %>%
  select(VJ.name, .y, group1, group2, p.adj, p.adj.signif) 

### Plots
## Heatmap of V-gene expression for all 4 groups
VJgenesG <- rbind.data.frame(hneut_pos_VJgenesG, hneut_neg_VJgenesG, lneut_pos_VJgenesG, lneut_neg_VJgenesG) %>%
  complete(VJ.name, Group, fill = list(VJ_count = 0)) %>%
  dplyr::select(VJ.name, Group, VJ_count, Prop) %>%
  replace(is.na(.), 0) %>%
  reshape2::dcast(., Group ~ VJ.name, value.var = "Prop")
row.names(VJgenesG) <- c("Hneut E2-", "Hneut E2+",  "Lneut E2-", "Lneut E2+" )

VJgenesG_hm <- VJgenesG[c("Hneut E2+", "Hneut E2-", "Lneut E2+", "Lneut E2-"),-1]

VJgenes_G_annot <- data.frame(row.names = rownames(VJgenesG_hm), Neut = c("High", "High", "Low", "Low"))
VJgenes_G_annot_cols <- list(Neut = c(High = "#009999", Low= "#b66dff"))

pheatmap(t(VJgenesG_hm), scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 14, labels_row = make_italics(colnames(VJgenesG_hm)), labels_col = c("E2+", "E2-", "E2+", "E2-"), annotation_col = VJgenes_G_annot, annotation_names_col = FALSE, annotation_colors = VJgenes_G_annot_cols, annotation_legend = FALSE)
grid.text("High\nNeut", x=0.367, y=.983, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("Low\nNeut", x=0.443, y=.983, gp=gpar(fontsize=12, fontface="bold.italic"))

## Volcano plot comparing V-gene expression between high and low neut E2-reactive B cells
FC_preVJG <- cbind.data.frame(VJgenes = VJgenes_statsG$VJ.name, HNP =t(VJgenesG_hm)[,1], LNP = t(VJgenesG_hm)[,3]) 
FC_preVJG[FC_preVJG == 0] <- 0.0002

FC_VJG <- FC_preVJG %>%
  mutate(FC = log2(HNP/LNP))

volcVJG <- data.frame(VJgenes = VJgenes_statsG$VJ.name,
                      Log2foldchange = FC_VJG$FC,
                      Pvalues = VJgenes_sigG$hnp_lnp) %>%
  mutate(Change = ifelse(Log2foldchange >= 0.6 & Pvalues < 0.05, "inc", 
                         ifelse(Log2foldchange <= -0.6 & Pvalues < 0.05, "dec", "unchanged"))) %>%
  mutate(Label = ifelse(Change == "inc" | Change == "dec", VJgenes, NA)) 


volcVJG$Change = factor(volcVJG$Change, levels = c("inc", "dec", "unchanged"))

volcVJG_p <- ggplot(data=volcVJG , aes(x=Log2foldchange, y=-log10(Pvalues), color=Change, label=Label))+
  geom_point(size=4) +
  scale_color_manual(labels = c("High Neut", "Low Neut", "No Change"), values=c("#009999", "#b66dff", "dark gray")) +
  geom_text_repel(size=7, show.legend = FALSE, fontface="bold.italic", color="black", max.overlaps=50) +
  geom_vline(xintercept=c(-0.6, 0.6), color="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), color="black", linetype="dashed") +
  scale_x_continuous(breaks = seq(from = -5, to = 5, by = 1), limits = c(-5,5))+
  scale_y_continuous(breaks = seq(from = -2, to = 10, by = 2), limits = c(-3,10))+
  ggtitle("B")+
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title = element_text(size=22),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, color = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        legend.key = element_rect(fill = "white", color = NA),
        legend.position = "right",
        plot.title = element_text(size = 45, hjust = -0.135, vjust = -1),
        aspect.ratio=0.9)+
  annotate("rect", xmin = -0.6, xmax = 0.6, ymin = -Inf, ymax = Inf,alpha = 0.1)+
  annotate("rect", xmin = -Inf, xmax = -0.6, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  annotate("rect", xmin = 0.6, xmax = Inf, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  labs(x= expression("Log"[2]~"Fold Change"),
       y= expression("-Log"[10]~"P-value"))
volcVJG_p

#### IGKL
### Quantify V-gene usage for each group
hneut_pos_VJgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, hneut_pos) %>%
  filter(hneut_pos > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "hneut_pos") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IG"))

hneut_neg_VJgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, hneut_neg) %>%
  filter(hneut_neg > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "hneut_neg") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IG"))

lneut_pos_VJgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, lneut_pos) %>%
  filter(lneut_pos > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "lneut_pos") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IG"))

lneut_neg_VJgenesKL <- cdr3KL_VJ_grp %>%
  dplyr::select(CDR3.aa, V.name, J.name, lneut_neg) %>%
  filter(lneut_neg > 0) %>%
  mutate(VJ.name = paste(V.name, J.name, sep = "/")) %>%
  group_by(VJ.name) %>%
  dplyr::summarise(VJ_count = n()) %>%
  mutate(Total = sum(VJ_count)) %>%
  mutate(Prop = VJ_count/Total) %>%
  mutate(Group = "lneut_neg") %>%
  mutate(VJ.name=str_remove_all(VJ.name, "IG"))

### Statistical analysis of differences between groups
VJstatKL <- rbind.data.frame(hneut_pos_VJgenesKL, hneut_neg_VJgenesKL, lneut_pos_VJgenesKL, lneut_neg_VJgenesKL) %>%
  complete(VJ.name, Group, fill = list(VJ_count = 0)) %>%
  dplyr::select(VJ.name, Group, VJ_count) %>%
  group_by(Group) %>%
  mutate(Total = sum(VJ_count)) 

listVJKL <- unique(VJstatKL$VJ.name)

VJgenes_statsKL <- data.frame()
for(i in 1:length(listVJKL)){
  VJG <- filter(VJstatKL, VJ.name == listVJKL[i])
  hnp_hnn <- fisher.test(matrix(c(VJG$VJ_count[2], VJG$VJ_count[1], VJG$Total[2]-VJG$VJ_count[2], VJG$Total[1]-VJG$VJ_count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(VJG$VJ_count[2], VJG$VJ_count[4], VJG$Total[2]-VJG$VJ_count[2], VJG$Total[4]-VJG$VJ_count[4]), ncol=2))$p.value
  lnp_lnn <- fisher.test(matrix(c(VJG$VJ_count[4], VJG$VJ_count[3], VJG$Total[4]-VJG$VJ_count[4], VJG$Total[3]-VJG$VJ_count[3]), ncol=2))$p.value
  hnp_lnn <- fisher.test(matrix(c(VJG$VJ_count[2], VJG$VJ_count[3], VJG$Total[2]-VJG$VJ_count[2], VJG$Total[3]-VJG$VJ_count[3]), ncol=2))$p.value
  lnp_hnn <- fisher.test(matrix(c(VJG$VJ_count[4], VJG$VJ_count[1], VJG$Total[4]-VJG$VJ_count[4], VJG$Total[1]-VJG$VJ_count[1]), ncol=2))$p.value
  hnn_lnn <- fisher.test(matrix(c(VJG$VJ_count[1], VJG$VJ_count[3], VJG$Total[1]-VJG$VJ_count[1], VJG$Total[3]-VJG$VJ_count[3]), ncol=2))$p.value
  VJG_df <- cbind.data.frame(VJ.name = listVJKL[i], hnp_hnn, hnp_lnp, lnp_lnn, hnp_lnn, lnp_hnn, hnn_lnn)
  VJgenes_statsKL <- rbind.data.frame(VJgenes_statsKL, VJG_df, stringsAsFactors = FALSE)
}

VJgenes_sigKL <-data.frame(VJgenes_statsKL$VJ.name, matrix(p.adjust(as.vector(as.matrix(VJgenes_statsKL[,-1])), method='BH'),ncol=6))
colnames(VJgenes_sigKL) <- colnames(VJgenes_statsKL)

VJgenes_sig_filtKL <- VJgenes_sigKL %>%
  filter_all(any_vars(. <0.05)) %>%
  pivot_longer(., c("hnp_hnn", "hnp_lnp", "lnp_lnn", "hnp_lnn", "lnp_hnn", "hnn_lnn"),
               names_to = "Comp",
               values_to = "p.adj") %>%
  mutate(.y = "Prop") %>%
  mutate(group1 = gsub("_.*","", Comp)) %>%
  mutate(group2 = gsub("..._","", Comp)) %>%
  mutate(group1 = ifelse(group1=="hnp", "hneut_pos", 
                         ifelse(group1=="hnn", "hneut_neg",
                                ifelse(group1=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(group2 = ifelse(group2=="hnp", "hneut_pos", 
                         ifelse(group2=="hnn", "hneut_neg",
                                ifelse(group2=="lnp", "lneut_pos", "lneut_neg")))) %>%
  mutate(p.adj.signif = ifelse(p.adj < 0.0001, "****", 
                               ifelse(p.adj < 0.001, "***",
                                      ifelse(p.adj < 0.01, "**",
                                             ifelse(p.adj <0.05, "*", "ns"))))) %>%
  select(VJ.name, .y, group1, group2, p.adj, p.adj.signif) 

### Plots
## Heatmap of V-gene expression for all 4 groups
VJgenesKL <- rbind.data.frame(hneut_pos_VJgenesKL, hneut_neg_VJgenesKL, lneut_pos_VJgenesKL, lneut_neg_VJgenesKL) %>%
  complete(VJ.name, Group, fill = list(VJ_count = 0)) %>%
  dplyr::select(VJ.name, Group, VJ_count, Prop) %>%
  replace(is.na(.), 0) %>%
  reshape2::dcast(., Group ~ VJ.name, value.var = "Prop")
row.names(VJgenesKL) <- c("Hneut E2-", "Hneut E2+",  "Lneut E2-", "Lneut E2+" )

VJgenesKL_hm <- VJgenesKL[c("Hneut E2+", "Hneut E2-", "Lneut E2+", "Lneut E2-"),-1]

VJgenes_KL_annot <- data.frame(row.names = rownames(VJgenesKL_hm), Neut = c("High", "High", "Low", "Low"))
VJgenes_KL_annot_cols <- list(Neut = c(High = "#009999", Low= "#b66dff"))

pheatmap(t(VJgenesKL_hm), scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 14, labels_row = make_italics(colnames(VJgenesKL_hm)), labels_col = c("E2+", "E2-", "E2+", "E2-"), annotation_col = VJgenes_KL_annot, annotation_names_col = FALSE, annotation_colors = VJgenes_KL_annot_cols, annotation_legend = FALSE)
grid.text("High\nNeut", x=0.367, y=.983, gp=gpar(fontsize=12, fontface="bold.italic"))
grid.text("Low\nNeut", x=0.443, y=.983, gp=gpar(fontsize=12, fontface="bold.italic"))

## Volcano plot comparing V-gene expression between high and low neut E2-reactive B cells
FC_preVJKL <- cbind.data.frame(VJgenes = VJgenes_statsKL$VJ.name, HNP =t(VJgenesKL_hm)[,1], LNP = t(VJgenesKL_hm)[,3]) 
FC_preVJKL[FC_preVJKL == 0] <- 0.0001

FC_VJKL <- FC_preVJKL %>%
  mutate(FC = log2(HNP/LNP))

volcVJKL <- data.frame(VJgenes = VJgenes_statsKL$VJ.name,
                       Log2foldchange = FC_VJKL$FC,
                       Pvalues = VJgenes_sigKL$hnp_lnp) %>%
  mutate(Change = ifelse(Log2foldchange >= 0.6 & Pvalues < 0.05, "inc", 
                         ifelse(Log2foldchange <= -0.6 & Pvalues < 0.05, "dec", "unchanged"))) %>%
  mutate(Label = ifelse(Change == "inc" | Change == "dec", VJgenes, NA)) 

volcVJKL_p <- ggplot(data=volcVJKL , aes(x=Log2foldchange, y=-log10(Pvalues), color=Change, label=Label))+
  geom_point(size=4) +
  scale_color_manual(labels = c("Low Neut", "High Neut", "No change"), values=c("#b66dff", "#009999", "dark gray")) +
  geom_text_repel(size=7, show.legend = FALSE, fontface="bold.italic", color="black", max.overlaps=50) +
  geom_vline(xintercept=c(-0.6, 0.6), color="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), color="black", linetype="dashed") +
  scale_x_continuous(breaks = seq(from = -5, to = 5, by = 1), limits = c(-5,5))+
  scale_y_continuous(breaks = seq(from = -2, to = 10, by = 2), limits = c(-3,10))+
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(size=22),
        axis.text.y=element_text(size=22),
        axis.title = element_text(size=24),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, color = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        legend.key = element_rect(fill = "white", color = NA),
        legend.position = "bottom")+
  annotate("rect", xmin = -0.6, xmax = 0.6, ymin = -Inf, ymax = Inf,alpha = 0.1)+
  annotate("rect", xmin = -Inf, xmax = -0.6, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  annotate("rect", xmin = 0.6, xmax = Inf, ymin = -Inf, ymax = -log10(0.05),alpha = 0.1)+
  labs(x= expression("Log"[2]~"Fold Change"),
       y= expression("-Log"[10]~"P-value"))
volcVJKL_p


##### CDR3 LENGTH
#### IGH
### Quantify CDRH3 length for each group
HCDR3_len <- cdr3G_seq_VJ %>%
  mutate_if(is.numeric, function(x)ifelse(x>0, 1, 0)) %>%
  rowwise() %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, C176p, C429p, P111p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, C176n, C429n, P111n, na.rm = T))

HCDR3_len_grps <- rbind.data.frame(
  mutate(filter(HCDR3_len, hneut_pos>0), E2 = "pos", Neut = "hneut", Group = "hneut_pos"),
  mutate(filter(HCDR3_len, lneut_pos>0), E2 = "pos", Neut = "lneut", Group = "lneut_pos"),
  mutate(filter(HCDR3_len, hneut_neg>0), E2 = "neg", Neut = "hneut", Group = "hneut_neg"),
  mutate(filter(HCDR3_len, lneut_neg>0), E2 = "neg", Neut = "lneut", Group = "lneut_neg"),
  stringsAsFactors = FALSE) %>%
  rowwise() %>%
  dplyr::select(CDR3.aa, E2, Neut, Group) %>%
  mutate(len = nchar(CDR3.aa))

### Statistical analysis of differences between groups
kruskal_test(HCDR3_len_grps, len ~ Group)

HCDR3_len_stat_neut <- as.data.frame(HCDR3_len_grps) %>%
  dunn_test(len ~ Group, 
            p.adjust.method = "BH") %>%
  add_significance() 

### Plot
## Specify order of groups
HCDR3_len_grps$Group <- factor(HCDR3_len_grps$Group, levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

# Make plot
HCDR3_len_p <- ggplot(HCDR3_len_grps, aes(x= Group, y = len)) +
  geom_point(aes(color=Group), shape=1, stroke = 0.75, position = position_jitter(width=0.35), size=1.5)+
  geom_violin(color= "#676767", show.legend = FALSE, lwd=1.5, alpha=0.4) +
  geom_boxplot(color= "#676767", outlier.color = NA, show.legend = FALSE, lwd=1, width = 0.2, alpha=0.4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray")) +
  scale_y_continuous(breaks = seq(from = 5, to = 38.95, by = 5), limits = c(5,38.95)) +
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  stat_pvalue_manual(HCDR3_len_stat_neut, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(36.5, 35, 38.9 ,37.4), tip.length = 0.005,vjust = 0.7, size=10, bracket.size = 0.75) +
  ggtitle("A")+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=24),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.165, vjust = -1),
        aspect.ratio=1.25)+
  labs(x= "",
       y= "CDRH3 AA Length")
HCDR3_len_p

#### IGKL
### Quantify CDRL3 length for each group
## IGK
LCDR3_lenK <- cdr3K_seq_VJ %>%
  mutate_if(is.numeric, function(x)ifelse(x>0, 1, 0)) %>%
  rowwise() %>%
  mutate(E2neg = sum(C48n, C117n, C172n, C176n, C429n, P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, C176p, C429p, P111p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, C176n, C429n, P111n, na.rm = T))

LCDR3_len_grpsK <- rbind.data.frame(
  mutate(filter(LCDR3_lenK, hneut_pos>0), E2 = "pos", Neut = "hneut", Group = "hneut_pos"),
  mutate(filter(LCDR3_lenK, lneut_pos>0), E2 = "pos", Neut = "lneut", Group = "lneut_pos"),
  mutate(filter(LCDR3_lenK, hneut_neg>0), E2 = "neg", Neut = "hneut", Group = "hneut_neg"),
  mutate(filter(LCDR3_lenK, lneut_neg>0), E2 = "neg", Neut = "lneut", Group = "lneut_neg"),
  stringsAsFactors = FALSE) %>%
  dplyr::select(CDR3.aa, E2, Neut, Group) %>%
  mutate(len = nchar(CDR3.aa))

## IGL
LCDR3_lenL <- cdr3L_seq_VJ %>%
  mutate_if(is.numeric, function(x)ifelse(x>0, 1, 0)) %>%
  rowwise() %>%
  mutate(E2neg = sum(C48n, C117n, C172n, C176n, C429n, P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, C176p, C429p, P111p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, C176n, C429n, P111n, na.rm = T))

LCDR3_len_grpsL <- rbind.data.frame(
  mutate(filter(LCDR3_lenL, hneut_pos>0), E2 = "pos", Neut = "hneut", Group = "hneut_pos"),
  mutate(filter(LCDR3_lenL, lneut_pos>0), E2 = "pos", Neut = "lneut", Group = "lneut_pos"),
  mutate(filter(LCDR3_lenL, hneut_neg>0), E2 = "neg", Neut = "hneut", Group = "hneut_neg"),
  mutate(filter(LCDR3_lenL, lneut_neg>0), E2 = "neg", Neut = "lneut", Group = "lneut_neg"),
  stringsAsFactors = FALSE) %>%
  dplyr::select(CDR3.aa, E2, Neut, Group) %>%
  mutate(len = nchar(CDR3.aa))

## Combine IGK and IGL
LCDR3_len_neut <-rbind.data.frame(LCDR3_len_grpsK, LCDR3_len_grpsL, stringsAsFactors = FALSE)

### Statistical analysis of differences between groups
kruskal_test(LCDR3_len_neut, len ~ Group)

LCDR3_len_stat_neut <- as.data.frame(LCDR3_len_neut) %>%
  dunn_test(len ~ Group,
            p.adjust.method = "BH") %>%
  add_significance() 

### Plot
## Specify order of groups
LCDR3_len_neut$Group <- factor(LCDR3_len_neut$Group, levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

# Make plot
LCDR3_len_p <- ggplot(LCDR3_len_neut, aes(x= Group, y = len)) +
  geom_point(aes(color= Group), shape=1, stroke=0.75, position = position_jitter(width = 0.35), size=1.5)+
  geom_violin(color= "#676767", show.legend = FALSE, lwd=1.5, alpha=0.4) +
  geom_boxplot(color= "#676767", outlier.color = NA, show.legend = FALSE, lwd=1, width = 0.2, alpha=0.4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray")) +
  scale_y_continuous(breaks = seq(from = 5, to = 21, by = 5), limits = c(5,21))+
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  stat_pvalue_manual(LCDR3_len_stat_neut, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(19.5, 20.25, 21, 19.5), tip.length = 0.005,vjust = 0.7, size=10, bracket.size = 0.75) +
  ggtitle("B")+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=24),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.165, vjust = -1),
        aspect.ratio=1.25)+
  labs(x= "",
       y= "CDRL3 AA Length")
LCDR3_len_p


##### SOMATIC HYPERMUTATION OF THE V-GENE
#### IGH
### Read in excel spreadsheet with IGHV mutation rates (obtained from IMGT), remove nonproductive or out-of-frame sequences, re-format V/J/D genes to remove allele information, and filter out overly short CDR1/CDR2 sequences
IMGT_IGH_all <- read_excel("/Volumes/nes002/Projects_Ongoing/HCV_clear_pers/Bulk_BCR/Sequences/For_revision/Deidentified/IGHV_SHM.xlsx") %>%
  select(clon_num = "Sequence number...1", V_function = "V-DOMAIN Functionality...3", V_gene_IMGT = "V-GENE and allele...4", V_identity= "V-REGION identity %", J_gene_IMGT = "J-GENE and allele", CDR3aa_IMGT = "AA JUNCTION", CDR3_frame = "JUNCTION frame", V_mutation = "V-REGION", CDR1_mut = "CDR1-IMGT", CDR2_mut = "CDR2-IMGT", CDR1_len = "CDR1-IMGT length", CDR2_len = "CDR2-IMGT length", Sequence) %>%
  mutate(Sequence = toupper(Sequence)) %>%
  filter(V_function == "productive" | V_function ==	"productive (see comment)") %>%
  filter(CDR3_frame == "in-frame") %>%
  filter(V_identity > 85) %>%
  replace(is.na(.), "none") %>%
  select(-V_function, -CDR3_frame) %>%
  mutate(V_gene_IMGT = gsub("(\\*+[0-9]+\\s+[A-Z]*)|(\\*+[0-9]+\\s\\()","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("[A-Z]+\\)","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("\\(see comment\\)","",V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("\\,.*","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("D","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("Homsap ","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub(" ","", V_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("(\\*+[0-9]+\\s+[A-Z]*)|(\\*+[0-9]+\\s\\()","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("[A-Z]+\\)","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("\\(see comment\\)","",J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("\\,.*","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("D","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("Homsap ","", J_gene_IMGT)) %>%
  mutate(CDR1_len = ifelse(CDR1_len == "X", 0, CDR1_len)) %>%
  mutate(CDR2_len = ifelse(CDR2_len == "X", 0, CDR2_len)) %>%
  rowwise() %>%
  mutate(CDR12_len = sum(as.numeric(CDR1_len), as.numeric(CDR2_len))) %>%
  filter(CDR12_len >=14) %>%
  select(-CDR12_len)

## Change format of mutations
IMGT_IGH_all$V_mutation <- sapply(IMGT_IGH_all$V_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))

## Add subject and group information
IMGT_IGH_all_grp <- IMGT_IGH_all %>%
  inner_join(cdr3G_seq_VJseq, by="clon_num") %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, C176p, C429p, P111p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, C176n, C429n, P111n, na.rm = T)) %>%
  mutate(E2neg = sum(C48n, C117n, C172n, C176n, C429n, P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
  ungroup() 

### Quantify SHM frequencies for each group
SHM_IGH_neut <- rbind.data.frame(
  mutate(filter(IMGT_IGH_all_grp, hneut_pos>0), E2 = "pos", Neut = "hneut", Group = "hneut_pos"),
  mutate(filter(IMGT_IGH_all_grp, lneut_pos>0), E2 = "pos", Neut = "lneut", Group = "lneut_pos"),
  mutate(filter(IMGT_IGH_all_grp, hneut_neg>0), E2 = "neg", Neut = "hneut", Group = "hneut_neg"),
  mutate(filter(IMGT_IGH_all_grp, lneut_neg>0), E2 = "neg", Neut = "lneut", Group = "lneut_neg"), 
  stringsAsFactors = FALSE) %>%
  dplyr::select(CDR3_aa = CDR3.aa, V_gene = V.name, J_gene = J.name, V_gene_IMGT, V_identity, V_mutation, E2, Neut, Group) 

### Statistical analysis of differences between groups
kruskal_test(SHM_IGH_neut, V_identity ~ Group)

SHM_IGH_stat_neut <- as.data.frame(SHM_IGH_neut) %>%
  dunn_test(V_identity ~ Group, 
            p.adjust.method = "BH") %>%
  add_significance() 



### Plot
## Specify order of groups
SHM_IGH_neut$Group <- factor(SHM_IGH_neut$Group,levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

## Make Plot
SHM_IGH_p <- ggplot(SHM_IGH_neut, aes(x= Group, y = V_identity)) +
  geom_point(aes(color = Group), shape =1, stroke = 0.75, position = position_jitter(width=0.35), size=1.5)+
  geom_violin(color="#676767", show.legend = FALSE, lwd=1.2, alpha = 0.4) +
  geom_boxplot(outlier.color = NA, color="#676767", show.legend = FALSE, lwd=1, alpha = 0.5, width = 0.4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray"))+
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  ggtitle("A")+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=24),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.18, vjust = -1),
        aspect.ratio=1.25)+
  stat_pvalue_manual(SHM_IGH_stat_neut, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(100.5, 102.5, 103.25, 101, 101.75), tip.length = 0.005, vjust = 0.7, size=10, bracket.size = 0.75) +
  labs(x= "",
       y= expression("%"~italic(V[H])~ " Identity"))
SHM_IGH_p

#### IGKLV
### Read in excel spreadsheet with IGHV mutation rates (obtained from IMGT), remove nonproductive or out-of-frame sequences, re-format V/J/D genes to remove allele information, and filter out overly short CDR1/CDR2 sequences
IMGT_IGKL_all <- read_excel("/Volumes/nes002/Projects_Ongoing/HCV_clear_pers/Bulk_BCR/Sequences/For_revision/Deidentified/IGKLV_SHM.xlsx") %>%
  select(clon_num = "Sequence number...1", V_function = "V-DOMAIN Functionality...3", V_gene_IMGT = "V-GENE and allele...4", V_identity= "V-REGION identity %", J_gene_IMGT = "J-GENE and allele", CDR3aa_IMGT = "AA JUNCTION", CDR3_frame = "JUNCTION frame", V_mutation = "V-REGION", CDR1_mut = "CDR1-IMGT", CDR2_mut = "CDR2-IMGT", CDR1_len = "CDR1-IMGT length", CDR2_len = "CDR2-IMGT length", Sequence) %>%
  mutate(Sequence = toupper(Sequence)) %>%
  filter(V_function == "productive" | V_function ==	"productive (see comment)") %>%
  filter(CDR3_frame == "in-frame") %>%
  filter(V_identity > 85) %>%
  replace(is.na(.), "none") %>%
  select(-V_function, -CDR3_frame) %>%
  mutate(V_gene_IMGT = gsub("(\\*+[0-9]+\\s+[A-Z]*)|(\\*+[0-9]+\\s\\()","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("[A-Z]+\\)","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("\\(see comment\\)","",V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("\\,.*","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("D","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub("Homsap ","", V_gene_IMGT)) %>%
  mutate(V_gene_IMGT = gsub(" ","", V_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("(\\*+[0-9]+\\s+[A-Z]*)|(\\*+[0-9]+\\s\\()","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("[A-Z]+\\)","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("\\(see comment\\)","",J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("\\,.*","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("D","", J_gene_IMGT)) %>%
  mutate(J_gene_IMGT = gsub("Homsap ","", J_gene_IMGT)) %>%
  mutate(CDR1_len = ifelse(CDR1_len == "X", 0, CDR1_len)) %>%
  mutate(CDR2_len = ifelse(CDR2_len == "X", 0, CDR2_len)) %>%
  rowwise() %>%
  mutate(CDR12_len = sum(as.numeric(CDR1_len), as.numeric(CDR2_len))) %>%
  filter(CDR12_len >=8) %>%
  select(-CDR12_len)

## Change format of mutations
IMGT_IGKL_all$V_mutation <- sapply(IMGT_IGKL_all$V_mutation, function(x) paste(unique(str_extract_all(x,"[A-Z]+[0-9]+>+[A-Z]")[[1]]), collapse = " "))

## Add subject and group information
IMGT_IGKL_all_grp <- IMGT_IGKL_all %>%
  inner_join(cdr3KL_seq_VJseq, by="clon_num") %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, C176p, C429p, P111p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, C176n, C429n, P111n, na.rm = T)) %>%
  mutate(E2neg = sum(C48n, C117n, C172n, C176n, C429n, P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
  ungroup() 

### Quantify SHM frequencies for each group
SHM_IGKL_neut <- rbind.data.frame(
  mutate(filter(IMGT_IGKL_all_grp, hneut_pos>0), E2 = "pos", Neut = "hneut", Group = "hneut_pos"),
  mutate(filter(IMGT_IGKL_all_grp, lneut_pos>0), E2 = "pos", Neut = "lneut", Group = "lneut_pos"),
  mutate(filter(IMGT_IGKL_all_grp, hneut_neg>0), E2 = "neg", Neut = "hneut", Group = "hneut_neg"),
  mutate(filter(IMGT_IGKL_all_grp, lneut_neg>0), E2 = "neg", Neut = "lneut", Group = "lneut_neg"), 
  stringsAsFactors = FALSE) %>%
  dplyr::select(CDR3_aa = CDR3.aa, V_gene = V.name, J_gene = J.name, V_gene_IMGT, V_identity, V_mutation, E2, Neut, Group) 

### Statistical analysis of differences between groups
kruskal_test(SHM_IGKL_neut, V_identity ~ Group)

SHM_IGKL_stat_neut <- as.data.frame(SHM_IGKL_neut) %>%
  dunn_test(V_identity ~ Group, 
            p.adjust.method = "BH") %>%
  add_significance() 

### Plot
## Specify order of groups
SHM_IGKL_neut$Group <- factor(SHM_IGKL_neut$Group,levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

## Make plot
SHM_IGKL_p <- ggplot(SHM_IGKL_neut, aes(x= Group, y = V_identity)) +
  geom_point(aes(color = Group), shape = 1, stroke=0.75, position = position_jitter(width=0.35), size=1.5)+
  geom_violin(color="#676767", show.legend = FALSE, lwd=1.2, alpha = 0.4) +
  geom_boxplot(outlier.color = NA, color="#676767", show.legend = FALSE, lwd=1, alpha = 0.5, width = 0.4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray"))+
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  ggtitle("B")+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=24),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.185, vjust = -1),
        aspect.ratio=1.25)+
  stat_pvalue_manual(SHM_IGKL_stat_neut, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(100.5, 102.5, 101, 101.75, 100.5), tip.length = 0.005, vjust = 0.7, size=10, bracket.size = 0.75) +
  labs(x= "",
       y= expression("%"~italic(V[K])~"/"~italic(V[L])~" Identity"))
SHM_IGKL_p


##### ASSOCIATION OF V-GENE SHM AND DURATION OF INFECTION
#### IGH (did not perform for IGK/IGL)
### Input duration of infection for each subject
dur_inf <- data.frame(samples = samples, DOI = c(379,379,156,156,314,314,318,318,343.5,343.5,285,285, 376,376,91,91,87,87,81.5,81.5))

### Quantify median SHM frequencies for each subject
SHM_subjG <- c(median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$C48p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$P49p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$P157p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$P53p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$P54p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$C117p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$C172p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$C176p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$P111p > 0,]$V_identity),
               median(IMGT_IGH_all_grp[IMGT_IGH_all_grp$C429p > 0,]$V_identity))

### Combine median SHM and duration of infection data
SHM_durG <- cbind.data.frame(dur_inf[grepl("p", dur_inf$samples),], V_ident_med = SHM_subjG, Group = c("hneut", "lneut","lneut", "hneut", "hneut", "hneut", "hneut", "lneut", "lneut", "lneut"))

### Plot
SHM_durG_p <- ggplot(SHM_durG, aes(x=DOI, y=V_ident_med)) +
  stat_cor(method = "kendall", size=7, label.x.npc="middle", label.y=99, digits=1)+
  geom_smooth(method = "lm", se = TRUE, color = "black", size=1, fill="light gray")+
  geom_point(size=6, aes(color = Group)) +
  scale_color_manual(values = c("#009999", "#b66dff"), labels = c("High Neut", "Low Neut"))+
  scale_y_continuous(limits=c(93.5,100))+
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text=element_text(size=18),
        legend.position="right",
        axis.line = element_line(size = 0.75, color = "black"),
        plot.title = element_text(size = 45, hjust = -0.1, vjust = -1),
        aspect.ratio = 0.85)+
  labs(x= "Duration of Infection (Days)",
       y= expression("Median "~italic(V[H])~ " Identity (%)"))
SHM_durG_p


##### V MUTATION QUANTITATION
#### base function quantitating mutation rates for plotting or stats:
# Vgene can be the name of a Vgene (e.g., "IGHV1-69") or "all" if all V genes are being analyzed
# mutation_type is "position" (e.g., 36>), "subsitution" (e.g., N36>S), or "delta" (e.g., 36>S)
# mutation_subset is "all" or a vector containing positions to analyze (e.g c(27:38) for CDR1)

cdr1cdr2 = c(27:38, 56:65) 

mutation_quant <- function(df1, df2, Vgene, mutation_type, mutation_subset) {
  ##extract mutations relevant for specific V gene of interest ("Vgene") or for all V genes
  if(Vgene != "all"){                     ##filtering for a specific V gene
    df1_V <- filter(df1, V_gene == Vgene) %>%
      filter(V_gene_IMGT == Vgene) %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%   #format mutations to reveal IMGT position only
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))  #format mutations to reveal IMGT position + mutated amino acid (but not germline amino acid)
    df2_V <- filter(df2, V_gene == Vgene) %>%
      filter(V_gene_IMGT == Vgene) %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))
  } else {                                ##if analyzing all V genes
    df1_V <- df1 %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))
    df2_V <- df2 %>%
      mutate(V_mutation_pos = gsub("[A-Z]","", V_mutation)) %>%
      mutate(V_mutation_delta = sapply(V_mutation, function(x) paste(sub(".","", strsplit(x," ")[[1]]), collapse = " ")))
  }
  
  ###extract positions intended for analysis
  if(mutation_subset[1] != "all"){      ##if analyzing specified IMGT positions
    filt_Vmut1 <- c()
    filt_Vmut_pos1 <- c()
    filt_Vmut_delta1 <- c()
    for(i in 1:nrow(df1_V)){      # filtering IMGT positions for each df of V mutations
      positions = paste(mutation_subset, collapse="|")
      
      V_muts <- unlist(strsplit(df1_V$V_mutation[i]," "))
      filt <- paste(V_muts[grepl(positions, V_muts)], collapse = " ")
      filt_Vmut1 <- c(filt_Vmut1, filt)
      
      V_muts_pos <- unlist(strsplit(df1_V$V_mutation_pos[i]," "))
      filt_pos <- paste(V_muts_pos[grepl(positions, V_muts_pos)], collapse = " ")
      filt_Vmut_pos1 <- c(filt_Vmut_pos1, filt_pos)
      
      V_muts_delta <- unlist(strsplit(df1_V$V_mutation_delta[i]," "))
      filt_delta <- paste(V_muts_delta[grepl(positions, V_muts_delta)], collapse = " ")
      filt_Vmut_delta1 <- c(filt_Vmut_delta1, filt_delta)
    }
    
    df1_V <- df1_V %>%          #inputing filtered positions into dfs
      mutate(V_mutation = filt_Vmut1) %>%
      mutate(V_mutation_pos = filt_Vmut_pos1) %>%
      mutate(V_mutation_delta = filt_Vmut_delta1)
    df1_V[df1_V==""] <- "none"
    
    
    # repeating the process for df2
    filt_Vmut2 <- c()
    filt_Vmut_pos2 <- c()
    filt_Vmut_delta2 <- c()
    for(i in 1:nrow(df2_V)){
      positions = paste(mutation_subset, collapse="|")
      V_muts <- unlist(strsplit(df2_V$V_mutation[i]," "))
      filt <- paste(V_muts[grepl(positions, V_muts)], collapse = " ")
      filt_Vmut2 <- c(filt_Vmut2, filt)
      
      V_muts_pos <- unlist(strsplit(df2_V$V_mutation_pos[i]," "))
      filt_pos <- paste(V_muts_pos[grepl(positions, V_muts_pos)], collapse = " ")
      filt_Vmut_pos2 <- c(filt_Vmut_pos2, filt_pos)
      
      V_muts_delta <- unlist(strsplit(df2_V$V_mutation_delta[i]," "))
      filt_delta <- paste(V_muts_delta[grepl(positions, V_muts_delta)], collapse = " ")
      filt_Vmut_delta2 <- c(filt_Vmut_delta2, filt_delta)
    }
    df2_V <- df2_V %>%
      mutate(V_mutation = filt_Vmut2) %>%
      mutate(V_mutation_pos = filt_Vmut_pos2) %>%
      mutate(V_mutation_delta = filt_Vmut_delta2)
    df2_V[df2_V==""] <- "none"
  } else{}              ## no filtering if analyzing all positions
  
  ## generate list of mutations (modified by whether analyzing substitution, position, or just the delta substitution)
  df1_subs <- unique(unlist(str_split(df1_V$V_mutation, " "))) 
  df1_pos <- unique(gsub("[A-Z]","", df1_subs))
  df1_delta <- unique(sub(".","", df1_subs))
  
  df2_subs <- unique(unlist(str_split(df2_V$V_mutation, " "))) 
  df2_pos <- unique(gsub("[A-Z]","", df2_subs))
  df2_delta <- unique(sub(".","", df2_subs))
  
  ## quantitate df1 and df2 mutations
  if(mutation_type == "substitution"){      # if interested in substitutions
    mut_list1 = df1_subs
    df_mut_list1 = df1_V$V_mutation
    mut_list2 = df2_subs
    df_mut_list2 = df2_V$V_mutation
  } else{if(mutation_type == "position"){   # if interested in position
    mut_list1 = df1_pos
    df_mut_list1 = df1_V$V_mutation_pos
    mut_list2 = df2_pos
    df_mut_list2 = df2_V$V_mutation_pos
  } else{if(mutation_type == "delta"){     # if interested in delta substitution
    mut_list1 = df1_delta
    df_mut_list1 = df1_V$V_mutation_delta
    mut_list2 = df2_delta
    df_mut_list2 = df2_V$V_mutation_delta
  } else{print("ERROR: You must enter an accepted mutation type: substitution (e.g., A36>G), position (e.g., 36>), or delta (e.g., 36>G)")}
  }}
  
  # generate df with mutation counts and proportions for df1 and df2
  df1_pos_counts <- data.frame()
  for(i in 1:length(mut_list1)){
    subs_count <- sum(str_count(df_mut_list1, mut_list1[i]))
    subs_prop <- subs_count/nrow(df1_V)
    df_count <- cbind.data.frame(Substitution = mut_list1[i], Count = subs_count, Total = nrow(df1_V), Prop = subs_prop, Group = "df1")
    df1_pos_counts <- rbind.data.frame(df1_pos_counts, df_count)
  }
  
  df2_pos_counts <- data.frame()
  for(i in 1:length(mut_list2)){
    subs_count <- sum(str_count(df_mut_list2, mut_list2[i]))
    subs_prop <- subs_count/nrow(df2_V)
    df_count <- cbind.data.frame(Substitution = mut_list2[i], Count = subs_count, Total = nrow(df2_V), Prop = subs_prop, Group = "df2")
    df2_pos_counts <- rbind.data.frame(df2_pos_counts, df_count)
  }
  
  # combine df1 and df2 into one df for plotting
  df_pos_counts <- rbind.data.frame(df1_pos_counts, df2_pos_counts) %>%
    complete(Substitution, Group, fill = list(Count = 0)) %>%
    dplyr::select(Substitution, Group, Count, Prop) %>%
    replace(is.na(.), 0) %>%
    dcast(., Group ~ Substitution, value.var = "Prop")
  row.names(df_pos_counts) <- c("df1", "df2")
  
  # format df for pheatmap
  df_pos_hm <- t(df_pos_counts[c("df1", "df2"),-1])
  
  return(list(df_pos_hm, df1_pos_counts, df2_pos_counts))   ## generate a list with the df structured for pheatmap [[1]], df1 mutation counts/props [[2]], df2 mutation counts/props [[3]] 
}

#### function for generating heatmaps using output from the function above; it therefore requires the same inputs as the base function (mutation_quant)
pheatmap_mut <- function(df1, df2, Vgene, mutation_type, mutation_subset) {
  pheatmap(mutation_quant(df1, df2, Vgene, mutation_type, mutation_subset)[[1]], , scale="none", cluster_rows = FALSE, cluster_cols=FALSE, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 11, cellwidth = 11, fontsize_row = 10)
}

#### function for generating stats using output from the mutation_quant function
stat_mut <- function(df1, df2, Vgene, mutation_type, mutation_subset){
  
  ## generate one dataframe from df1 df2
  df_stat <- rbind.data.frame(mutation_quant(df1, df2, Vgene, mutation_type, mutation_subset)[[2]], 
                              mutation_quant(df1, df2, Vgene, mutation_type, mutation_subset)[[3]]) %>%
    complete(Substitution, Group, fill = list(Count = 0)) %>%   # for each mutation that is only present in one group (df1 or df2), add a row for the other group with count = 0
    mutate(Total = rep(c(unique(filter(., Group == "df1")$Total[1]), 
                         unique(filter(., Group == "df2")$Total[1])), 
                       nrow(.)/2)) %>%    #fill in totals for added rows
    replace(is.na(.), 0)
  
  ## obtain uncorrected P values using Fisher test to compare proportion of mutation usage (# times used/total # cloonotypes) between the 2 groups
  df_pvalue <- data.frame()
  for(i in 1:length(unique(df_stat$Substitution))){
    subst <- filter(df_stat, Substitution == unique(df_stat$Substitution)[i])
    fisherP <- fisher.test(matrix(c(subst$Count[1], subst$Total[1], subst$Count[2], subst$Total[2]), ncol=2))$p.value
    subs_df <- cbind.data.frame(Substitution = unique(df_stat$Substitution)[i], 
                                df1_df2 = fisherP)
    df_pvalue  <- rbind.data.frame(df_pvalue, subs_df, stringsAsFactors = FALSE)
  }
  
  ## correct P values using Benjamini-Hochberg correction (aka FDR)
  df_pvalue_corr <- cbind.data.frame(df_pvalue, df_pvalue_corr = p.adjust(df_pvalue$df1_df2, method = "fdr")) 
  
  ## extract mutations with significant P values (p < 0.01 is used; for ~ 20 comparisons, this would mean 0.2 false positives; for 100 comparisons, it would mean 1 false positive) 
  df_sig_pvalue <- df_pvalue_corr %>%
    filter(df_pvalue_corr < 0.01) %>%
    dplyr::select(Substitution, df_pvalue_corr)
  
  return(list(df_stat, df_pvalue_corr, df_sig_pvalue))    ## generate a list with the df containing mutation counts/proportions [[1]], P values and corrected P values for each mutation [[2]], and just the mutations with significant P values [[3]] 
}


#### IGH
### Quantify V mutations for each group
Vmut_hneut_pos <- filter(IMGT_IGH_all_grp, hneut_pos > 0)
Vmut_hneut_neg <- filter(IMGT_IGH_all_grp, hneut_neg > 0)
Vmut_lneut_pos <- filter(IMGT_IGH_all_grp, lneut_pos > 0)
Vmut_lneut_neg <- filter(IMGT_IGH_all_grp, lneut_neg > 0)

Vmut_hnp_hnn <- mutation_quant(Vmut_hneut_pos, Vmut_hneut_neg, "all", "delta", cdr1cdr2)
Vmut_lnp_lnn <- mutation_quant(Vmut_lneut_pos, Vmut_lneut_neg, "all", "delta", cdr1cdr2)


Vmut_neut_dfG <- rbind.data.frame(mutate(Vmut_hnp_hnn[[2]], Group = "hneut_pos"),
                                  mutate(Vmut_hnp_hnn[[3]], Group = "hneut_neg"),
                                  mutate(Vmut_lnp_lnn[[2]], Group = "lneut_pos"),
                                  mutate(Vmut_lnp_lnn[[3]], Group = "lneut_neg"), stringsAsFactors = FALSE) %>%
  complete(Substitution, Group, fill = list(Count = 0)) %>%
  dplyr::select(-Prop) %>%
  group_by(Group) %>%
  fill(Total, .direction = "down") %>%
  mutate(Substitution = gsub("one", "None", Substitution))


### Statistical analysis of differences between groups
subs_neutG <- unique(Vmut_neut_dfG$Substitution)

Vmut_neut_statG <- data.frame()
for(i in 1:length(subs_neutG)){ 
  sub_df <- filter(Vmut_neut_dfG, Substitution == subs_neutG[i])
  hnp_hnn <- fisher.test(matrix(c(sub_df$Count[2], sub_df$Count[1], sub_df$Total[2]-sub_df$Count[2], sub_df$Total[1]-sub_df$Count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(sub_df$Count[2], sub_df$Count[4], sub_df$Total[2]-sub_df$Count[2], sub_df$Total[4]-sub_df$Count[4]), ncol=2))$p.value
  lnp_lnn <- fisher.test(matrix(c(sub_df$Count[4], sub_df$Count[3], sub_df$Total[4]-sub_df$Count[4], sub_df$Total[3]-sub_df$Count[3]), ncol=2))$p.value
  hnp_lnn <- fisher.test(matrix(c(sub_df$Count[2], sub_df$Count[3], sub_df$Total[2]-sub_df$Count[2], sub_df$Total[3]-sub_df$Count[3]), ncol=2))$p.value
  lnp_hnn <- fisher.test(matrix(c(sub_df$Count[4], sub_df$Count[1], sub_df$Total[4]-sub_df$Count[4], sub_df$Total[1]-sub_df$Count[1]), ncol=2))$p.value
  hnn_lnn <- fisher.test(matrix(c(sub_df$Count[1], sub_df$Count[3], sub_df$Total[1]-sub_df$Count[1], sub_df$Total[3]-sub_df$Count[3]), ncol=2))$p.value
  fish_df <- cbind.data.frame(Substitution = subs_neutG[i], hnp_hnn, hnp_lnp, lnp_lnn, hnp_lnn, lnp_hnn, hnn_lnn)
  Vmut_neut_statG <- rbind.data.frame(Vmut_neut_statG, fish_df, stringsAsFactors = FALSE)
}

Vmut_neut_statG <- Vmut_neut_statG[Vmut_neut_statG$Substitution != "None",]

Vmut_neut_sigG <- Vmut_neut_statG[,c(1,3)] %>%
  mutate(p_adj = p.adjust(hnp_lnp, "BH"))

## Convert IMGT numbering to Kabat numbering
mutation_convG <- data.frame(IMGT = c(27,28,29,30,31,32,33,34,35,36,37,38, 56,57,58,59,60,61,62,63,64,65), Kabat = c("26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "35A ", "35B ", "51", "52", "52A ", "52B ", "52C ", "53", "54", "55", "56", "57"))

Vmut_neut_sig_filtG <- Vmut_neut_sigG %>%
  filter(p_adj < 0.01) %>%
  mutate(IMGT = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(AA_change = gsub("[0-9][0-9]>", "", Substitution)) %>%
  left_join(mutation_convG, by="IMGT") %>%
  mutate(Kabat_subs = paste0(Kabat,AA_change)) %>%
  dplyr::select(-IMGT, -AA_change, -Kabat)

Vmut_prop_sig <- Vmut_neut_dfG %>%
  filter(Substitution %in% Vmut_neut_sig_filtG$Substitution) %>%
  filter(!grepl("neg", Group)) %>%
  group_by(Substitution) %>%
  mutate(Prop = Count/Total) %>%
  group_by(Group) %>%
  group_split()

Vmut_sig_filt_FCG <-cbind.data.frame(Vmut_neut_sig_filtG, FC = Vmut_prop_sig[[1]]$Prop/Vmut_prop_sig[[2]]$Prop) %>%
  filter(FC >= 1.5 | FC <= -1.5)


### Identify which V-genes contained each mutation
VgermG <- data.frame()
for(i in 1:nrow(Vmut_neut_sig_filtG)){
  vmut = Vmut_neut_sig_filtG$Substitution[i]
  
  Vgerm <- IMGT_IGH_all_grp %>%
    filter(grepl(vmut, V_mutation)) %>%
    filter(hneut_pos > 0 | lneut_pos > 0) %>%
    mutate(group = ifelse(hneut_pos > 0 & lneut_pos > 0, "both", ifelse(hneut_pos > 0, "hneut", "lneut"))) %>%
    group_by(group, V.name) %>%
    mutate(numV = n()) %>%
    group_by(group) %>%
    mutate(total = n()) %>%
    distinct(group, V.name, numV, total) %>%
    mutate(prop = numV/total) %>%
    mutate(Substitution = vmut)
  
  VgermG <- rbind.data.frame(VgermG, Vgerm, stringsAsFactors = FALSE)
}


#### IGKL
### Quantify V mutations for each group
Vmut_hneut_posKL <- filter(IMGT_IGKL_all_grp, hneut_pos > 0)
Vmut_hneut_negKL <- filter(IMGT_IGKL_all_grp, hneut_neg > 0)
Vmut_lneut_posKL <- filter(IMGT_IGKL_all_grp, lneut_pos > 0)
Vmut_lneut_negKL <- filter(IMGT_IGKL_all_grp, lneut_neg > 0)

Vmut_hnp_hnnKL <- mutation_quant(Vmut_hneut_posKL, Vmut_hneut_negKL, "all", "delta", cdr1cdr2)
Vmut_lnp_lnnKL <- mutation_quant(Vmut_lneut_posKL, Vmut_lneut_negKL, "all", "delta", cdr1cdr2)

Vmut_neut_dfKL <- rbind.data.frame(mutate(Vmut_hnp_hnnKL[[2]], Group = "hneut_pos"),
                                   mutate(Vmut_hnp_hnnKL[[3]], Group = "hneut_neg"),
                                   mutate(Vmut_lnp_lnnKL[[2]], Group = "lneut_pos"),
                                   mutate(Vmut_lnp_lnnKL[[3]], Group = "lneut_neg"), stringsAsFactors = FALSE) %>%
  complete(Substitution, Group, fill = list(Count = 0)) %>%
  dplyr::select(-Prop) %>%
  group_by(Group) %>%
  fill(Total, .direction = "down") %>%
  replace(is.na(.), 1132) %>%
  mutate(Substitution = gsub("one", "None", Substitution))

### Statistical analysis of differences between groups
subs_neutKL <- unique(Vmut_neut_dfKL$Substitution)

Vmut_neut_statKL <- data.frame()
for(i in 1:length(subs_neutKL)){ 
  sub_df <- filter(Vmut_neut_dfKL, Substitution == subs_neutKL[i])
  hnp_hnn <- fisher.test(matrix(c(sub_df$Count[2], sub_df$Count[1], sub_df$Total[2]-sub_df$Count[2], sub_df$Total[1]-sub_df$Count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(sub_df$Count[2], sub_df$Count[4], sub_df$Total[2]-sub_df$Count[2], sub_df$Total[4]-sub_df$Count[4]), ncol=2))$p.value
  lnp_lnn <- fisher.test(matrix(c(sub_df$Count[4], sub_df$Count[3], sub_df$Total[4]-sub_df$Count[4], sub_df$Total[3]-sub_df$Count[3]), ncol=2))$p.value
  hnp_lnn <- fisher.test(matrix(c(sub_df$Count[2], sub_df$Count[3], sub_df$Total[2]-sub_df$Count[2], sub_df$Total[3]-sub_df$Count[3]), ncol=2))$p.value
  lnp_hnn <- fisher.test(matrix(c(sub_df$Count[4], sub_df$Count[1], sub_df$Total[4]-sub_df$Count[4], sub_df$Total[1]-sub_df$Count[1]), ncol=2))$p.value
  hnn_lnn <- fisher.test(matrix(c(sub_df$Count[1], sub_df$Count[3], sub_df$Total[1]-sub_df$Count[1], sub_df$Total[3]-sub_df$Count[3]), ncol=2))$p.value
  fish_df <- cbind.data.frame(Substitution = subs_neutKL[i], hnp_hnn, hnp_lnp, lnp_lnn, hnp_lnn, lnp_hnn, hnn_lnn)
  Vmut_neut_statKL <- rbind.data.frame(Vmut_neut_statKL, fish_df, stringsAsFactors = FALSE)
}

Vmut_neut_statKL <- Vmut_neut_statKL[Vmut_neut_statKL$Substitution != "None",]

Vmut_neut_sigKL <- Vmut_neut_statKL[,c(1,3)] %>%
  mutate(p_adj = p.adjust(hnp_lnp, "BH"))

## Convert IMGT numbering to Kabat numbering
mutation_convKL <- data.frame(IMGT = c(27,28,29,30,31,32,33,34,35,36,37,38, 56,57,58,59,62,63,64,65), Kabat = c("27", "27A ", "27B ", "27C ", "27D ", "27E ", "27F ", "28", "29", "30", "31", "32", "50", "50A ", "50B ", "50C ", "50D ", "50E ", "51", "52"))

Vmut_neut_sig_filtKL <- Vmut_neut_sigKL %>%
  filter(p_adj < 0.05) %>%
  mutate(IMGT = as.numeric(gsub(">[A-Z]", "", Substitution))) %>%
  mutate(AA_change = gsub("[0-9][0-9]>", "", Substitution)) %>%
  left_join(mutation_convKL, by="IMGT") %>%
  mutate(Kabat_subs = paste0(Kabat,AA_change)) %>%
  dplyr::select(-IMGT, -AA_change, -Kabat)

Vmut_prop_sigKL <- Vmut_neut_dfKL %>%
  filter(Substitution %in% Vmut_neut_sig_filtKL$Substitution) %>%
  filter(!grepl("neg", Group)) %>%
  group_by(Substitution) %>%
  mutate(Prop = Count/Total) %>%
  group_by(Group) %>%
  group_split()

Vmut_sig_filt_FCKL <-cbind.data.frame(Vmut_neut_sig_filtKL, FC = Vmut_prop_sigKL[[1]]$Prop/Vmut_prop_sigKL[[2]]$Prop) %>%
  filter(FC>=1.5 | FC<=-1.5)

### Identify which V-genes contained each mutation
VgermKL <- data.frame()
for(i in 1:nrow(Vmut_sig_filt_FCKL)){
  vmut = Vmut_sig_filt_FCKL$Substitution[i]
  
  Vgerm <- IMGT_IGKL_all_grp %>%
    filter(grepl(vmut, V_mutation)) %>%
    filter(hneut_pos > 0 | lneut_pos > 0) %>%
    mutate(group = ifelse(hneut_pos > 0 & lneut_pos > 0, "both", ifelse(hneut_pos > 0, "hneut", "lneut"))) %>%
    group_by(group, V.name) %>%
    mutate(numV = n()) %>%
    group_by(group) %>%
    mutate(total = n()) %>%
    distinct(group, V.name, numV, total) %>%
    mutate(prop = numV/total) %>%
    mutate(Substitution = vmut)
  
  VgermKL <- rbind.data.frame(VgermKL, Vgerm, stringsAsFactors = FALSE)
}


##### PUBLIC CLONOTYPES
#### IGH
### Quantify # shared clonotypes within each group (and not present in the other groups)
hneut_pos100VG <- cdr3G_seq_VJ_raw[rowSums(cdr3G_seq_VJ_raw[, c("C48p", "C117p", "C172p", "P53p", "P54p")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, lneut_pos==0, lneut_neg ==0)

hneut_neg100VG <- cdr3G_seq_VJ_raw[rowSums(cdr3G_seq_VJ_raw[, c("C48n", "C117n", "C172n", "P53n", "P54n")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48n, C117n, C172n, P53n, P54n, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_pos == 0, lneut_pos==0, lneut_neg ==0)

lneut_pos100VG <- cdr3G_seq_VJ_raw[rowSums(cdr3G_seq_VJ_raw[, c("P49p", "P157p", "P111p", "C176p", "C429p")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, P49p, P157p, P111p, C176p, C429p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, hneut_pos==0, lneut_neg ==0)

lneut_neg100VG <- cdr3G_seq_VJ_raw[rowSums(cdr3G_seq_VJ_raw[, c("P49n", "P157n", "P111n", "C176n", "C429n")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, P49p, P157n, P111n, C176n, C429n, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, hneut_pos==0, lneut_pos ==0)

cdr3G_seqVJ_grp <- cdr3G_seq_VJ_raw %>%
  rowwise() %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) 

## Quantify total number of clonotypes per group
hnp_total_seqVJG <- filter(cdr3G_seqVJ_grp, hneut_pos > 0)
hnn_total_seqVJG <- filter(cdr3G_seqVJ_grp, hneut_neg > 0)
lnp_total_seqVJG <- filter(cdr3G_seqVJ_grp, lneut_pos > 0)
lnn_total_seqVJG <- filter(cdr3G_seqVJ_grp, lneut_neg > 0)

### Statistical analysis of differences between groups
neut100G_comps <-cbind.data.frame(group= c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                                  shared = c(nrow(hneut_pos100VG),
                                             nrow(hneut_neg100VG),
                                             nrow(lneut_pos100VG),
                                             nrow(lneut_neg100VG)),
                                  total = c(nrow(hnp_total_seqVJG),
                                            nrow(hnn_total_seqVJG),
                                            nrow(lnp_total_seqVJG),
                                            nrow(lnn_total_seqVJG)))

neut100G_comps <- data.frame(Comp = c("hnp_hnn", "hnp_lnp", "lnp_lnn", "hnp_lnn", "lnp_hnn", "hnn_lnn"),
                             P_values = c(fisher.test(matrix(c(neut100G_comps$shared[1], 
                                                               neut100G_comps$shared[2], 
                                                               neut100G_comps$total[1]-neut100G_comps$shared[1],
                                                               neut100G_comps$total[2]-neut100G_comps$shared[2]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100G_comps$shared[1], 
                                                               neut100G_comps$shared[3], 
                                                               neut100G_comps$total[1]-neut100G_comps$shared[1],
                                                               neut100G_comps$total[3]-neut100G_comps$shared[3]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100G_comps$shared[3], 
                                                               neut100G_comps$shared[4], 
                                                               neut100G_comps$total[3]-neut100G_comps$shared[3],
                                                               neut100G_comps$total[4]-neut100G_comps$shared[4]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100G_comps$shared[1], 
                                                               neut100G_comps$shared[4], 
                                                               neut100G_comps$total[1]-neut100G_comps$shared[1],
                                                               neut100G_comps$total[4]-neut100G_comps$shared[4]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100G_comps$shared[3], 
                                                               neut100G_comps$shared[2], 
                                                               neut100G_comps$total[3]-neut100G_comps$shared[3],
                                                               neut100G_comps$total[2]-neut100G_comps$shared[2]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100G_comps$shared[2], 
                                                               neut100G_comps$shared[4], 
                                                               neut100G_comps$total[2]-neut100G_comps$shared[2],
                                                               neut100G_comps$total[4]-neut100G_comps$shared[4]), ncol=2))$p.value)) %>%
  mutate(P_adj = p.adjust(P_values, "BH"))


seqVJgrpsG_stat <- data.frame(group1 = c("hneut_pos", "hneut_pos", "lneut_pos", "hneut_pos", "lneut_pos", "hneut_neg"),
                              group2 = c("hneut_neg", "lneut_pos", "lneut_neg", "lneut_neg", "hneut_neg", "lneut_neg"),
                              p = neut100G_comps$P_values,
                              p.adj = neut100G_comps$P_adj,
                              p.adj.signif = c("****", "****", "****", "****", "****", "ns"))



### Plot
neut100G_seqVJgrps <- cbind.data.frame(group= c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                                       prop = c(nrow(hneut_pos100VG)/nrow(hnp_total_seqVJG),
                                                nrow(hneut_neg100VG)/nrow(hnn_total_seqVJG),
                                                nrow(lneut_pos100VG)/nrow(lnp_total_seqVJG),
                                                nrow(lneut_neg100VG)/nrow(lnn_total_seqVJG)),
                                       stringsAsFactors=FALSE)

neut100G_seqVJgrps$group <- factor(neut100G_seqVJgrps$group,levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

seqVJgrpsG_p <- ggplot(neut100G_seqVJgrps) +
  geom_bar_pattern(aes(x=group, y=prop, fill = group, color = group, pattern = group, pattern_color = group), size = 3, position="dodge", stat="identity",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.035,
                   pattern_size = 1,
                   pattern_key_scale_factor = 0.6, inherit.aes=FALSE)+
  scale_fill_manual(values = c("#009999", "white", "#b66dff", "white"))+
  scale_pattern_manual(values = c("none", "stripe", "none", "stripe"))+
  scale_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_pattern_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_x_discrete(breaks = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"), labels = c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  ggtitle("A")+
  stat_pvalue_manual(seqVJgrpsG_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.0405, 0.0425, 0.0325, 0.0445, 0.034), tip.length = 0.005, vjust= 0.7, bracket.size=0.75, size=11.5) +
  theme(axis.title.x = element_blank(), 
        axis.text=element_text(color="black", size = 20),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.17, vjust = -1.3),
        aspect.ratio=1.1)+
  labs(x= "",
       y= "Prop. Public Clonotypes")
seqVJgrpsG_p

#### IGK
### Quantify # shared clonotypes within each group (and not present in the other groups)
hneut_pos100VK <- cdr3K_seq_VJ_raw[rowSums(cdr3K_seq_VJ_raw[, c("C48p", "C117p", "C172p", "P53p", "P54p")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, lneut_pos==0, lneut_neg ==0)

hneut_neg100VK <- cdr3K_seq_VJ_raw[rowSums(cdr3K_seq_VJ_raw[, c("C48n", "C117n", "C172n", "P53n", "P54n")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_pos == 0, lneut_pos==0, lneut_neg ==0)

lneut_pos100VK <- cdr3K_seq_VJ_raw[rowSums(cdr3K_seq_VJ_raw[, c("P49p", "P157p", "P111p", "C176p", "C429p")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, hneut_pos==0, lneut_neg ==0)

lneut_neg100VK <- cdr3K_seq_VJ_raw[rowSums(cdr3K_seq_VJ_raw[, c("P49n", "P157n", "P111n", "C176n", "C429n")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, hneut_pos==0, lneut_pos ==0)

cdr3K_seqVJ_grp <- cdr3K_seq_VJ_raw %>%
  rowwise() %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) 

## Quantify total number of clonotypes per group
hnp_total_seqVJK <- filter(cdr3K_seqVJ_grp, hneut_pos > 0)
hnn_total_seqVJK <- filter(cdr3K_seqVJ_grp, hneut_neg > 0)
lnp_total_seqVJK <- filter(cdr3K_seqVJ_grp, lneut_pos > 0)
lnn_total_seqVJK <- filter(cdr3K_seqVJ_grp, lneut_neg > 0)

### Statistical analysis of differences between groups
neut100K_comps <-cbind.data.frame(group= c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                                  shared = c(nrow(hneut_pos100VK),
                                             nrow(hneut_neg100VK),
                                             nrow(lneut_pos100VK),
                                             nrow(lneut_neg100VK)),
                                  total = c(nrow(hnp_total_seqVJK),
                                            nrow(hnn_total_seqVJK),
                                            nrow(lnp_total_seqVJK),
                                            nrow(lnn_total_seqVJK)))

neut100K_comps <- data.frame(Comp = c("hnp_hnn", "hnp_lnp", "lnp_lnn", "hnp_lnn", "lnp_hnn", "hnn_lnn"),
                             P_values = c(fisher.test(matrix(c(neut100K_comps$shared[1], 
                                                               neut100K_comps$shared[2], 
                                                               neut100K_comps$total[1]-neut100K_comps$shared[1],
                                                               neut100K_comps$total[2]-neut100K_comps$shared[2]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100K_comps$shared[1], 
                                                               neut100K_comps$shared[3], 
                                                               neut100K_comps$total[1]-neut100K_comps$shared[1],
                                                               neut100K_comps$total[3]-neut100K_comps$shared[3]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100K_comps$shared[3], 
                                                               neut100K_comps$shared[4], 
                                                               neut100K_comps$total[3]-neut100K_comps$shared[3],
                                                               neut100K_comps$total[4]-neut100K_comps$shared[4]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100K_comps$shared[1], 
                                                               neut100K_comps$shared[4], 
                                                               neut100K_comps$total[1]-neut100K_comps$shared[1],
                                                               neut100K_comps$total[4]-neut100K_comps$shared[4]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100K_comps$shared[3], 
                                                               neut100K_comps$shared[2], 
                                                               neut100K_comps$total[3]-neut100K_comps$shared[3],
                                                               neut100K_comps$total[2]-neut100K_comps$shared[2]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100K_comps$shared[2], 
                                                               neut100K_comps$shared[4], 
                                                               neut100K_comps$total[2]-neut100K_comps$shared[2],
                                                               neut100K_comps$total[4]-neut100K_comps$shared[4]), ncol=2))$p.value)) %>%
  mutate(P_adj = p.adjust(P_values, "BH"))


seqVJgrpsK_stat <- data.frame(group1 = c("hneut_pos", "hneut_pos", "lneut_pos", "hneut_pos", "lneut_pos", "hneut_neg"),
                              group2 = c("hneut_neg", "lneut_pos", "lneut_neg", "lneut_neg", "hneut_neg", "lneut_neg"),
                              p = neut100K_comps$P_values,
                              p.adj = neut100K_comps$P_adj,
                              p.adj.signif = c("****", "****", "**", "****", "****", "****"))



### Plot
neut100K_seqVJgrps <- cbind.data.frame(group= c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                                       prop = c(nrow(hneut_pos100VK)/nrow(hnp_total_seqVJK),
                                                nrow(hneut_neg100VK)/nrow(hnn_total_seqVJK),
                                                nrow(lneut_pos100VK)/nrow(lnp_total_seqVJK),
                                                nrow(lneut_neg100VK)/nrow(lnn_total_seqVJK)),
                                       stringsAsFactors=FALSE)

neut100K_seqVJgrps$group <- factor(neut100K_seqVJgrps$group,levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

seqVJgrpsK_p <- ggplot(neut100K_seqVJgrps) +
  geom_bar_pattern(aes(x=group, y=prop, fill = group, color = group, pattern = group, pattern_color = group), size = 3, position="dodge", stat="identity",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.035,
                   pattern_size = 1,
                   pattern_key_scale_factor = 0.6, inherit.aes=FALSE)+
  scale_fill_manual(values = c("#009999", "white", "#b66dff", "white"))+
  scale_pattern_manual(values = c("none", "stripe", "none", "stripe"))+
  scale_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_pattern_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_x_discrete(breaks = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"), labels = c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  stat_pvalue_manual(seqVJgrpsK_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.0405, 0.0425, 0.0325, 0.0445, 0.034, 0.036), tip.length = 0.005, vjust= 0.7, bracket.size=0.75, size=11.5) +
  ggtitle("B")+
  theme(axis.title.x = element_blank(), 
        axis.text=element_text(color="black", size = 20),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.18, vjust = -1.3),
        aspect.ratio=1.1)+
  labs(x= "",
       y= "Prop. Public Clonotypes")
seqVJgrpsK_p

#### IGL
### Quantify # shared clonotypes within each group (and not present in the other groups)
hneut_pos100VL <- cdr3L_seq_VJ_raw[rowSums(cdr3L_seq_VJ_raw[, c("C48p", "C117p", "C172p", "P53p", "P54p")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, lneut_pos==0, lneut_neg ==0)

hneut_neg100VL <- cdr3L_seq_VJ_raw[rowSums(cdr3L_seq_VJ_raw[, c("C48n", "C117n", "C172n", "P53n", "P54n")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_pos == 0, lneut_pos==0, lneut_neg ==0)

lneut_pos100VL <- cdr3L_seq_VJ_raw[rowSums(cdr3L_seq_VJ_raw[, c("P49p", "P157p", "P111p", "C176p", "C429p")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, hneut_pos==0, lneut_neg ==0)

lneut_neg100VL <- cdr3L_seq_VJ_raw[rowSums(cdr3L_seq_VJ_raw[, c("P49n", "P157n", "P111n", "C176n", "C429n")] >= 1 ) >=2,] %>%
  rowwise() %>% 
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) %>%
  dplyr::select(CDR3.aa, V.name, J.name, C48p, C117p, C172p, P53p, P54p, hneut_pos, hneut_neg, lneut_pos, lneut_neg) %>%
  filter(hneut_neg == 0, hneut_pos==0, lneut_pos ==0)

cdr3L_seqVJ_grp <- cdr3L_seq_VJ_raw %>%
  rowwise() %>%
  mutate(hneut_pos=sum(C48p, C117p, C172p, P53p, P54p, na.rm = T)) %>%
  mutate(hneut_neg=sum(C48n, C117n, C172n, P53n, P54n, na.rm = T)) %>%
  mutate(lneut_pos=sum(P49p, P157p, P111p, C176p, C429p, na.rm = T)) %>%
  mutate(lneut_neg=sum(P49n, P157n, P111n, C176n, C429n, na.rm = T)) 

## Quantify total number of clonotypes per group
hnp_total_seqVJL <- filter(cdr3K_seqVJ_grp, hneut_pos > 0)
hnn_total_seqVJL <- filter(cdr3K_seqVJ_grp, hneut_neg > 0)
lnp_total_seqVJL <- filter(cdr3K_seqVJ_grp, lneut_pos > 0)
lnn_total_seqVJL <- filter(cdr3K_seqVJ_grp, lneut_neg > 0)

### Statistical analysis of differences between groups
neut100L_comps <-cbind.data.frame(group= c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                                  shared = c(nrow(hneut_pos100VL),
                                             nrow(hneut_neg100VL),
                                             nrow(lneut_pos100VL),
                                             nrow(lneut_neg100VL)),
                                  total = c(nrow(hnp_total_seqVJL),
                                            nrow(hnn_total_seqVJL),
                                            nrow(lnp_total_seqVJL),
                                            nrow(lnn_total_seqVJL)))

neut100L_comps <- data.frame(Comp = c("hnp_hnn", "hnp_lnp", "lnp_lnn", "hnp_lnn", "lnp_hnn", "hnn_lnn"),
                             P_values = c(fisher.test(matrix(c(neut100L_comps$shared[1], 
                                                               neut100L_comps$shared[2], 
                                                               neut100L_comps$total[1]-neut100L_comps$shared[1],
                                                               neut100L_comps$total[2]-neut100L_comps$shared[2]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100L_comps$shared[1], 
                                                               neut100L_comps$shared[3], 
                                                               neut100L_comps$total[1]-neut100L_comps$shared[1],
                                                               neut100L_comps$total[3]-neut100L_comps$shared[3]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100L_comps$shared[3], 
                                                               neut100L_comps$shared[4], 
                                                               neut100L_comps$total[3]-neut100L_comps$shared[3],
                                                               neut100L_comps$total[4]-neut100L_comps$shared[4]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100L_comps$shared[1], 
                                                               neut100L_comps$shared[4], 
                                                               neut100L_comps$total[1]-neut100L_comps$shared[1],
                                                               neut100L_comps$total[4]-neut100L_comps$shared[4]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100L_comps$shared[3], 
                                                               neut100L_comps$shared[2], 
                                                               neut100L_comps$total[3]-neut100L_comps$shared[3],
                                                               neut100L_comps$total[2]-neut100L_comps$shared[2]), ncol=2))$p.value,
                                          fisher.test(matrix(c(neut100L_comps$shared[2], 
                                                               neut100L_comps$shared[4], 
                                                               neut100L_comps$total[2]-neut100L_comps$shared[2],
                                                               neut100L_comps$total[4]-neut100L_comps$shared[4]), ncol=2))$p.value)) %>%
  mutate(P_adj = p.adjust(P_values, "BH"))


seqVJgrpsL_stat <- data.frame(group1 = c("hneut_pos", "hneut_pos", "lneut_pos", "hneut_pos", "lneut_pos", "hneut_neg"),
                              group2 = c("hneut_neg", "lneut_pos", "lneut_neg", "lneut_neg", "hneut_neg", "lneut_neg"),
                              p = neut100L_comps$P_values,
                              p.adj = neut100L_comps$P_adj,
                              p.adj.signif = c("****", "****", "**", "****", "****", "***"))



### Plot
neut100L_seqVJgrps <- cbind.data.frame(group= c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                                       prop = c(nrow(hneut_pos100VL)/nrow(hnp_total_seqVJL),
                                                nrow(hneut_neg100VL)/nrow(hnn_total_seqVJL),
                                                nrow(lneut_pos100VL)/nrow(lnp_total_seqVJL),
                                                nrow(lneut_neg100VL)/nrow(lnn_total_seqVJL)),
                                       stringsAsFactors=FALSE)

neut100L_seqVJgrps$group <- factor(neut100L_seqVJgrps$group,levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

seqVJgrpsL_p <- ggplot(neut100L_seqVJgrps) +
  geom_bar_pattern(aes(x=group, y=prop, fill = group, color = group, pattern = group, pattern_color = group), size = 3, position="dodge", stat="identity",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.035,
                   pattern_size = 1,
                   pattern_key_scale_factor = 0.6, inherit.aes=FALSE)+
  scale_fill_manual(values = c("#009999", "white", "#b66dff", "white"))+
  scale_pattern_manual(values = c("none", "stripe", "none", "stripe"))+
  scale_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_pattern_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_x_discrete(breaks = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"), labels = c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  stat_pvalue_manual(seqVJgrpsL_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.0405, 0.0425, 0.0325, 0.0445, 0.034, 0.036), tip.length = 0.005, vjust= 0.7, bracket.size=0.75, size=11.5) +
  ggtitle("C")+
  theme(axis.title.x = element_blank(), 
        axis.text=element_text(color="black", size = 20),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.175, vjust = -1.3),
        aspect.ratio = 1.1)+
  labs(x= "",
       y= "Prop. Public Clonotypes")
seqVJgrpsL_p


##### RANDOM PERMUTATION TO DETERMINE EXPECTED SHARING
#### IGH
### Randomly permute clonotype/subject data 1000x
sharedVJ100G_perm_rep <- data.frame()
for(i in 1:1000){
  set.seed(i)
  print(i)
  perm_repVJ <- t(apply(cdr3G_seq_VJ_raw[,4:23], 1, sample))
  perm_cdr3G_seq_repVJ <- cbind.data.frame(cdr3G_seq_VJ_raw[,1:3], perm_repVJ)
  colnames(perm_cdr3G_seq_repVJ) <- colnames(cdr3G_seq_VJ_raw)
  
  # clear/pos
  cpos100G_perm_repVJ <- perm_cdr3G_seq_repVJ[rowSums(perm_cdr3G_seq_repVJ[, c("C48p", "C117p", "C172p", "C176p", "C429p")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, C48p, C117p, C172p, C176p, C429p, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_neg == 0, pers_pos==0, pers_neg ==0)
  
  # clear/neg
  cneg100G_perm_repVJ <- perm_cdr3G_seq_repVJ[rowSums(perm_cdr3G_seq_repVJ[, c("C48n", "C117n", "C172n", "C176n", "C429n")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, C48n, C117n, C172n, C176n, C429n, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, pers_pos==0, pers_neg ==0)
  
  # pers/pos
  ppos100G_perm_repVJ <- perm_cdr3G_seq_repVJ[rowSums(perm_cdr3G_seq_repVJ[, c("P49p", "P157p", "P53p", "P54p", "P111p")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, P49p, P157p, P53p, P54p, P111p, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, clear_neg ==0, pers_neg ==0)
  
  # pers/neg
  pneg100G_perm_repVJ <- perm_cdr3G_seq_repVJ[rowSums(perm_cdr3G_seq_repVJ[, c("P49n", "P157n", "P53n", "P54n", "P111n")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, P49n, P157n, P53n, P54n, P111n, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, clear_neg ==0, pers_pos ==0)
  
  #totals
  cdr3G_VJ_grp_perm_rep <- perm_cdr3G_seq_repVJ %>%
    rowwise() %>%
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T))
  
  cp_total_VJG_perm_rep <- filter(cdr3G_VJ_grp_perm_rep, clear_pos > 0)
  cn_total_VJG_perm_rep <- filter(cdr3G_VJ_grp_perm_rep, clear_neg > 0)
  pp_total_VJG_perm_rep <- filter(cdr3G_VJ_grp_perm_rep, pers_pos > 0)
  pn_total_VJG_perm_rep <- filter(cdr3G_VJ_grp_perm_rep, pers_neg > 0)
  
  #df
  shared100G_perm_repVJ <- cbind.data.frame(run = rep(i, 4),
                                            group= c("cl_pos", "cl_neg", "pers_pos", "pers_neg"),
                                            shared = c(nrow(cpos100G_perm_repVJ),
                                                       nrow(cneg100G_perm_repVJ),
                                                       nrow(ppos100G_perm_repVJ),
                                                       nrow(pneg100G_perm_repVJ)),
                                            total = c(nrow(cp_total_VJG_perm_rep),
                                                      nrow(cn_total_VJG_perm_rep),
                                                      nrow(pp_total_VJG_perm_rep),
                                                      nrow(pn_total_VJG_perm_rep)),
                                            prop = c(nrow(cpos100G_perm_repVJ)/nrow(cp_total_VJG_perm_rep),
                                                     nrow(cneg100G_perm_repVJ)/nrow(cn_total_VJG_perm_rep),
                                                     nrow(ppos100G_perm_repVJ)/nrow(pp_total_VJG_perm_rep),
                                                     nrow(pneg100G_perm_repVJ)/nrow(pn_total_VJG_perm_rep)),
                                            stringsAsFactors=FALSE)
  
  # all dfs
  sharedVJ100G_perm_rep <- rbind.data.frame(sharedVJ100G_perm_rep, shared100G_perm_repVJ, stringsAsFactors = FALSE)
}

### Quantify proportion of shared clonotypes in permuted data
histo_perm100VJ_cl <- filter(sharedVJ100G_perm_rep, group=="cl_pos")

## Median
median(histo_perm100VJ_cl$prop) #0.015

## 99.9th percentile
perc_prop_cl_permVJ <- histo_perm100VJ_cl[order(histo_perm100VJ_cl$prop),]
perc_prop_cl_permVJ[999,] #0.02

### Plot
## Show # shared clonotypes for each group (not permuted)
lines_data <- data.frame(Group = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"), Xpos = c(0.028, 0.0046, 0.018, 0.011)) %>%
  mutate(Xmin = Xpos - 0.0004) %>%
  mutate(Xmax = Xpos + 0.0004)
lines_data$Group <- factor(lines_data$Group, levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg")) 

## Make plot
perm_cl_histoVJ_prop_p <- ggplot(histo_perm100VJ_cl, aes(x=prop)) + 
  geom_histogram(binwidth = 0.001, color="#c1c1c1", fill="light gray", size=0.35) +
  geom_density(size =0.75, color="#c1c1c1", fill="light gray", alpha = 0.25)+
  geom_rect_pattern(data=lines_data,aes(xmin=Xmin,ymin=0,xmax=Xmax,ymax=255,fill=Group, pattern = Group, pattern_color = Group, color = Group),
                    pattern_angle = 45,
                    pattern_density = 0.1,
                    pattern_spacing = 0.015,
                    pattern_size = 0.7,
                    pattern_key_scale_factor = 0.6, size=1, inherit.aes=FALSE)+
  scale_fill_manual(values = c("#009999", "white", "#b66dff", "white"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_pattern_manual(values = c("hneut_pos" = "none", "hneut_neg" = "stripe", "lneut_pos" = "none", "lneut_neg" = "stripe"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_pattern_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_x_continuous(breaks = seq(from = 0, to = 0.034, by = .01), limits = c(0,0.034))+
  scale_y_continuous(breaks = seq(from = 0, to = 262, by = 50), limits = c(0, 262))+
  geom_segment(x=0.015, y=0, xend=0.015, yend=255,
               color="black", linetype="dashed", size=1) +
  geom_segment(x=0.02, y=0, xend=0.02, yend=255,
               color="black", linetype="dashed", size=1) +
  annotate(geom = "text",
           label = c("M", "99.9%"),
           x = c(0.015, 0.02),
           y = c(262, 262),
           size = 8, fontface="bold.italic")+
  ggtitle("D") +
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.text= element_text(size=20),
        legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 45, hjust = -0.13, vjust = -1.3),
        aspect.ratio = 1.1)+
  labs(x= "Prop. Public Clonotypes",
       y= "Count")
perm_cl_histoVJ_prop_p 

#### IGK
### Randomly permute clonotype/subject data 1000x
sharedVJ100K_perm_rep <- data.frame()
for(i in 1:1000){
  set.seed(i)
  print(i)
  perm_repVJK <- t(apply(cdr3K_seq_VJ_raw[,4:23], 1, sample))
  perm_cdr3K_seq_repVJ <- cbind.data.frame(cdr3K_seq_VJ_raw[,1:3], perm_repVJK)
  colnames(perm_cdr3K_seq_repVJ) <- colnames(cdr3K_seq_VJ_raw)
  
  # clear/pos
  cpos100K_perm_repVJ <- perm_cdr3K_seq_repVJ[rowSums(perm_cdr3K_seq_repVJ[, c("C48p", "C117p", "C172p", "C176p", "C429p")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, C48p, C117p, C172p, C176p, C429p, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_neg == 0, pers_pos==0, pers_neg ==0)
  
  # clear/neg
  cneg100K_perm_repVJ <- perm_cdr3K_seq_repVJ[rowSums(perm_cdr3K_seq_repVJ[, c("C48n", "C117n", "C172n", "C176n", "C429n")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, C48n, C117n, C172n, C176n, C429n, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, pers_pos==0, pers_neg ==0)
  
  # pers/pos
  ppos100K_perm_repVJ <- perm_cdr3K_seq_repVJ[rowSums(perm_cdr3K_seq_repVJ[, c("P49p", "P157p", "P53p", "P54p", "P111p")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, P49p, P157p, P53p, P54p, P111p, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, clear_neg ==0, pers_neg ==0)
  
  # pers/neg
  pneg100K_perm_repVJ <- perm_cdr3K_seq_repVJ[rowSums(perm_cdr3K_seq_repVJ[, c("P49n", "P157n", "P53n", "P54n", "P111n")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, P49n, P157n, P53n, P54n, P111n, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, clear_neg ==0, pers_pos ==0)
  
  #totals
  cdr3K_VJ_grp_perm_rep <- perm_cdr3K_seq_repVJ %>%
    rowwise() %>%
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T))
  
  cp_total_VJK_perm_rep <- filter(cdr3K_VJ_grp_perm_rep, clear_pos > 0)
  cn_total_VJK_perm_rep <- filter(cdr3K_VJ_grp_perm_rep, clear_neg > 0)
  pp_total_VJK_perm_rep <- filter(cdr3K_VJ_grp_perm_rep, pers_pos > 0)
  pn_total_VJK_perm_rep <- filter(cdr3K_VJ_grp_perm_rep, pers_neg > 0)
  
  #df
  shared100K_perm_repVJ <- cbind.data.frame(run = rep(i, 4),
                                            group= c("cl_pos", "cl_neg", "pers_pos", "pers_neg"),
                                            shared = c(nrow(cpos100K_perm_repVJ),
                                                       nrow(cneg100K_perm_repVJ),
                                                       nrow(ppos100K_perm_repVJ),
                                                       nrow(pneg100K_perm_repVJ)),
                                            total = c(nrow(cp_total_VJK_perm_rep),
                                                      nrow(cn_total_VJK_perm_rep),
                                                      nrow(pp_total_VJK_perm_rep),
                                                      nrow(pn_total_VJK_perm_rep)),
                                            prop = c(nrow(cpos100K_perm_repVJ)/nrow(cp_total_VJK_perm_rep),
                                                     nrow(cneg100K_perm_repVJ)/nrow(cn_total_VJK_perm_rep),
                                                     nrow(ppos100K_perm_repVJ)/nrow(pp_total_VJK_perm_rep),
                                                     nrow(pneg100K_perm_repVJ)/nrow(pn_total_VJK_perm_rep)),
                                            stringsAsFactors=FALSE)
  
  # all dfs
  sharedVJ100K_perm_rep <- rbind.data.frame(sharedVJ100K_perm_rep, shared100K_perm_repVJ, stringsAsFactors = FALSE)
}

### Quantify proportion of shared clonotypes in permuted data
histo_perm100VJ_clK <- filter(sharedVJ100K_perm_rep, group=="cl_pos")

# Median
median(histo_perm100VJ_clK$prop) #0.019

# 99.9th percentile
perc_prop_clK_permVJ <- histo_perm100VJ_clK[order(histo_perm100VJ_clK$prop),]
perc_prop_clK_permVJ[999,] #0.024

### Plot
## Show # shared clonotypes for each group (not permuted)
lines_dataK <- data.frame(Group = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"), Xpos = c(0.04, 0.0036, 0.023, 0.016)) %>%
  mutate(Xmin = Xpos - 0.0004) %>%
  mutate(Xmax = Xpos + 0.0004)
lines_dataK$Group <- factor(lines_dataK$Group, levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg")) 

## Make plot
perm_cl_histoVJ_propK_p <- ggplot(histo_perm100VJ_clK, aes(x=prop)) + 
  geom_histogram(binwidth = 0.001, color="#c1c1c1", fill="light gray", size=0.35) +
  geom_density(size =0.75, color="#c1c1c1", fill="light gray", alpha = 0.25)+
  geom_rect_pattern(data=lines_dataK,aes(xmin=Xmin,ymin=0,xmax=Xmax,ymax=255,fill=Group, pattern = Group, pattern_color = Group, color = Group),
                    pattern_angle = 45,
                    pattern_density = 0.1,
                    pattern_spacing = 0.015,
                    pattern_size = 0.7,
                    pattern_key_scale_factor = 0.6, size=1, inherit.aes=FALSE)+
  scale_fill_manual(values = c("#009999", "white", "#b66dff", "white"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_pattern_manual(values = c("hneut_pos" = "none", "hneut_neg" = "stripe", "lneut_pos" = "none", "lneut_neg" = "stripe"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_pattern_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_x_continuous(breaks = seq(from = 0, to = 0.04, by = .01), limits = c(0,0.0425))+
  scale_y_continuous(breaks = seq(from = 0, to = 262, by = 50), limits = c(0, 262))+
  geom_segment(x=0.018, y=0, xend=0.018, yend=255,
               color="black", linetype="dashed", size=1) +
  geom_segment(x=0.024, y=0, xend=0.024, yend=255,
               color="black", linetype="dashed", size=1) +
  annotate(geom = "text",
           label = c("M", "99.9%"),
           x = c(0.018, 0.024),
           y = c(262, 262),
           size = 8, fontface="bold.italic")+
  ggtitle("E")+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.text= element_text(size=20),
        legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 45, hjust = -0.13, vjust = -1.3),
        aspect.ratio = 1.1)+
  labs(x= "Prop. Public Clonotypes",
       y= "Count")
perm_cl_histoVJ_propK_p 


#### IGL
### Randomly permute clonotype/subject data
sharedVJ100L_perm_rep <- data.frame()
for(i in 1:1000){
  set.seed(i)
  print(i)
  perm_repVJL <- t(apply(cdr3L_seq_VJ_raw[,4:23], 1, sample))
  perm_cdr3L_seq_repVJ <- cbind.data.frame(cdr3L_seq_VJ_raw[,1:3], perm_repVJL)
  colnames(perm_cdr3L_seq_repVJ) <- colnames(cdr3L_seq_VJ_raw)
  
  # clear/pos
  cpos100L_perm_repVJ <- perm_cdr3L_seq_repVJ[rowSums(perm_cdr3L_seq_repVJ[, c("C48p", "C117p", "C172p", "C176p", "C429p")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, C48p, C117p, C172p, C176p, C429p, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_neg == 0, pers_pos==0, pers_neg ==0)
  
  # clear/neg
  cneg100L_perm_repVJ <- perm_cdr3L_seq_repVJ[rowSums(perm_cdr3L_seq_repVJ[, c("C48n", "C117n", "C172n", "C176n", "C429n")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, C48n, C117n, C172n, C176n, C429n, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, pers_pos==0, pers_neg ==0)
  
  # pers/pos
  ppos100L_perm_repVJ <- perm_cdr3L_seq_repVJ[rowSums(perm_cdr3L_seq_repVJ[, c("P49p", "P157p", "P53p", "P54p", "P111p")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, P49p, P157p, P53p, P54p, P111p, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, clear_neg ==0, pers_neg ==0)
  
  # pers/neg
  pneg100L_perm_repVJ <- perm_cdr3L_seq_repVJ[rowSums(perm_cdr3L_seq_repVJ[, c("P49n", "P157n", "P53n", "P54n", "P111n")] >= 1 ) >=2,] %>%
    rowwise() %>% 
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T)) %>%
    dplyr::select(CDR3.aa, P49n, P157n, P53n, P54n, P111n, clear_pos, clear_neg, pers_pos, pers_neg) %>%
    filter(clear_pos == 0, clear_neg ==0, pers_pos ==0)
  
  #totals
  cdr3L_VJ_grp_perm_rep <- perm_cdr3L_seq_repVJ %>%
    rowwise() %>%
    mutate(clear_pos=sum(C48p, C117p, C172p, C176p, C429p, na.rm = T)) %>%
    mutate(clear_neg=sum(C48n, C117n, C172n, C176n, C429n, na.rm = T)) %>%
    mutate(pers_pos=sum(P49p, P157p, P53p, P54p, P111p, na.rm = T)) %>%
    mutate(pers_neg=sum(P49n, P157n, P53n, P54n, P111n, na.rm = T))
  
  cp_total_VJL_perm_rep <- filter(cdr3L_VJ_grp_perm_rep, clear_pos > 0)
  cn_total_VJL_perm_rep <- filter(cdr3L_VJ_grp_perm_rep, clear_neg > 0)
  pp_total_VJL_perm_rep <- filter(cdr3L_VJ_grp_perm_rep, pers_pos > 0)
  pn_total_VJL_perm_rep <- filter(cdr3L_VJ_grp_perm_rep, pers_neg > 0)
  
  #df
  shared100L_perm_repVJ <- cbind.data.frame(run = rep(i, 4),
                                            group= c("cl_pos", "cl_neg", "pers_pos", "pers_neg"),
                                            shared = c(nrow(cpos100L_perm_repVJ),
                                                       nrow(cneg100L_perm_repVJ),
                                                       nrow(ppos100L_perm_repVJ),
                                                       nrow(pneg100L_perm_repVJ)),
                                            total = c(nrow(cp_total_VJL_perm_rep),
                                                      nrow(cn_total_VJL_perm_rep),
                                                      nrow(pp_total_VJL_perm_rep),
                                                      nrow(pn_total_VJL_perm_rep)),
                                            prop = c(nrow(cpos100L_perm_repVJ)/nrow(cp_total_VJL_perm_rep),
                                                     nrow(cneg100L_perm_repVJ)/nrow(cn_total_VJL_perm_rep),
                                                     nrow(ppos100L_perm_repVJ)/nrow(pp_total_VJL_perm_rep),
                                                     nrow(pneg100L_perm_repVJ)/nrow(pn_total_VJL_perm_rep)),
                                            stringsAsFactors=FALSE)
  
  # all dfs
  sharedVJ100L_perm_rep <- rbind.data.frame(sharedVJ100L_perm_rep, shared100L_perm_repVJ, stringsAsFactors = FALSE)
}


### Quantify proportion of shared clonotypes in permuted data
histo_perm100VJ_clL <- filter(sharedVJ100L_perm_rep, group=="cl_pos")

## Median
median(histo_perm100VJ_clL$prop) #0.017

## 99.9th percentile
perc_prop_clL_permVJ <- histo_perm100VJ_clL[order(histo_perm100VJ_clL$prop),]
perc_prop_clL_permVJ[999,] #0.022

### Plot
## Show # shared clonotypes for each group (not permuted)
lines_dataL <- data.frame(Group = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"), Xpos = c(0.028, 0.01046, 0.018, 0.011)) %>%
  mutate(Xmin = Xpos - 0.0005) %>%
  mutate(Xmax = Xpos + 0.0005)
lines_dataL$Group <- factor(lines_dataL$Group, levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg")) 

## Make plot
perm_cl_histoVJ_propL_p <- ggplot(histo_perm100VJ_clL, aes(x=prop)) + 
  geom_histogram(binwidth = 0.001, color="#c1c1c1", fill="light gray", size=0.35) +
  geom_density(size =0.75, color="#c1c1c1", fill="light gray", alpha = 0.25)+
  geom_rect_pattern(data=lines_dataK,aes(xmin=Xmin,ymin=0,xmax=Xmax,ymax=255,fill=Group, pattern = Group, pattern_color = Group, color = Group),
                    pattern_angle = 45,
                    pattern_density = 0.1,
                    pattern_spacing = 0.015,
                    pattern_size = 0.7,
                    pattern_key_scale_factor = 0.6, size=1, inherit.aes=FALSE)+
  scale_fill_manual(values = c("#009999", "white", "#b66dff", "white"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_pattern_manual(values = c("hneut_pos" = "none", "hneut_neg" = "stripe", "lneut_pos" = "none", "lneut_neg" = "stripe"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_pattern_color_manual(values = c("#009999", "#009999", "#b66dff", "#b66dff"), labels=c("High Neut E2+", "High Neut E2-", "Low Neut E2+", "Low Neut E2-"))+
  scale_x_continuous(breaks = seq(from = 0, to = 0.042, by = .01), limits = c(0,0.042))+
  scale_y_continuous(breaks = seq(from = 0, to = 260, by = 50), limits = c(0, 260))+
  geom_segment(x=0.017, y=0, xend=0.017, yend=255,
               color="black", linetype="dashed", size=0.7) +
  geom_segment(x=0.022, y=0, xend=0.022, yend=255,
               color="black", linetype="dashed", size=0.7) +
  annotate(geom = "text",
           label = c("M", "99.9%"),
           x = c(0.017, 0.022),
           y = c(260, 260),
           size = 8, fontface="bold.italic")+
  ggtitle("F")+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.text= element_text(size=20),
        legend.title = element_blank(),
        legend.position="right",
        plot.title = element_text(size = 45, hjust = -0.13, vjust = -1.3),
        aspect.ratio = 1.1)+
  labs(x= "Prop. Public Clonotypes",
       y= "Count")
perm_cl_histoVJ_propL_p 


##### CLONOTYPE SHARING BETWEEN EACH PAIR OF SUBJECTS 
#### IGH
### Quantify total number of clonotypes for each sample
seq_totalsG <- data.frame()
for(i in 4:23){
  ct_clonoG <- filter(cdr3G_seq_VJ_raw, cdr3G_seq_VJ_raw[,i] >0)
  add_seqG <- cbind.data.frame(subject=samples[i-3], num_clonos = nrow(ct_clonoG))
  seq_totalsG <- rbind.data.frame(seq_totalsG, add_seqG)
}

seq_totalsG_pos <- seq_totalsG %>%
  filter(grepl("p",subject))

seq_totalsG_neg <- seq_totalsG %>%
  filter(grepl("n",subject))

### Quantify sharing between each subject (normalize by dividing shared by nonshared)
shareG <- data.frame()
for(i in 4:23){
  filt_subjG <- filter(cdr3G_seq_VJ_raw, cdr3G_seq_VJ_raw[,i] > 0)
  share_quantG <- colSums(filt_subjG[,4:23] != 0)
  denomG <- seq_totalsG[,2] + seq_totalsG[i-3,2] - share_quantG
  norm_shareG <- share_quantG/denomG %>%
    replace(., i-3, NA)
  shareG <- rbind.data.frame(shareG, norm_shareG)
}

colnames(shareG) <- colnames(cdr3G_seq_VJ_raw[,4:23])
rownames(shareG) <- colnames(cdr3G_seq_VJ_raw[,4:23])

### Plot
## E2-reactive B cells
share_posG <- shareG[c(2,4,6,8,10,12,14,16,18,20),c(2,4,6,8,10,12,14,16,18,20)]

colnames(share_posG) <- gsub("p", "", colnames(share_posG))
rownames(share_posG) <- gsub("p", "", rownames(share_posG))

# Specify genotype and group
hm_groups <- data.frame(Genotype = c("Gen1a", "Gen1b", "Gen1a", "Gen3a", "Gen3a", "Gen1a", "Gen2b", "Gen1a", "Gen1a", "Gen1b"), Group = c("High Neut", "High Neut", "Low Neut", "Low Neut", "High Neut", "Low Neut", "Low Neut", "Low Neut", "High Neut", "High Neut"))
rownames(hm_groups) <- rownames(share_posG)

# Color code groups
hm_cols <- list(Genotype = c(Gen1a="black", Gen1b="#484848", Gen2b="#727272", Gen3a="#d4d4d4"), Group = c("High Neut"="#009999", "Low Neut"="#b66dff"))

# Make heatmap
share_posG_hm <- as.ggplot(pheatmap(share_posG, scale="none", cluster_rows = TRUE, cluster_cols=TRUE, show_rownames = TRUE, show_colnames = TRUE, border_color = NA, breaks = seq(0, 0.14, length.out = 100), cellheight = 14, cellwidth = 14, fontsize = 14, fontsize_row = 14, fontsize_col = 14, treeheight_row = 0, clustering_method = "average", annotation_col = hm_groups, annotation_names_col = FALSE, annotation_colors = hm_cols, annotation_legend = FALSE, legend=TRUE, main = "IGH"))
grid.text("A", x=0.24, y=0.73, gp=gpar(fontsize=30))

## E2-nonreactive
share_negG <- shareG[c(1,3,5,7,9,11,13,15,17,19),c(1,3,5,7,9,11,13,15,17,19)]

colnames(share_negG) <- gsub("n", "", colnames(share_negG))
rownames(share_negG) <- gsub("n", "", rownames(share_negG))

# Make heatmap
share_negG_hm <- as.ggplot(pheatmap(share_negG, scale="none", cluster_rows = TRUE, cluster_cols=TRUE, show_rownames = TRUE, show_colnames = TRUE, breaks = seq(0, 0.14, length.out = 100), border_color = NA, cellheight = 14, cellwidth = 14, fontsize = 14, fontsize_row = 14, fontsize_col = 14, treeheight_row = 0, clustering_method = "average", annotation_col = hm_groups, annotation_names_col = FALSE, annotation_colors = hm_cols, annotation_legend = FALSE, legend=FALSE))

#### IGK
### Quantify total number of clonotypes for each sample
seq_totalsK <- data.frame()
for(i in 4:23){
  ct_clonoK <- filter(cdr3K_seq_VJ, cdr3K_seq_VJ[,i] >0)
  add_seqK <- cbind.data.frame(subject=samples[i-3], num_clonos = nrow(ct_clonoK))
  seq_totalsK <- rbind.data.frame(seq_totalsK, add_seqK)
}

### Quantify sharing between each subject (normalize by dividing shared by nonshared)
shareK <- data.frame()
for(i in 4:23){
  filt_subjK <- filter(cdr3K_seq_VJ, cdr3K_seq_VJ[,i] > 0)
  share_quantK <- colSums(filt_subjK[,4:23] != 0)
  denomK <- seq_totalsK[,2] + seq_totalsK[i-3,2] - share_quantK
  norm_shareK <- share_quantK/denomK %>%
    replace(., i-3, NA)
  shareK <- rbind.data.frame(shareK, norm_shareK)
}

colnames(shareK) <- colnames(cdr3K_seq_VJ_raw[,4:23])
rownames(shareK) <- colnames(cdr3K_seq_VJ_raw[,4:23])

### Plot
## E2-reactive B cells
share_posK <- shareK[c(2,4,6,8,10,12,14,16,18,20),c(2,4,6,8,10,12,14,16,18,20)]

colnames(share_posK) <- gsub("p", "", colnames(share_posK))
rownames(share_posK) <- gsub("p", "", rownames(share_posK))

# Make heatmap
share_posK_hm <- as.ggplot(pheatmap(share_posK, scale="none", cluster_rows = TRUE, cluster_cols=TRUE, show_rownames = TRUE, show_colnames = TRUE, breaks = seq(0, 0.08, length.out = 100), border_color = NA, cellheight = 14, cellwidth = 14, fontsize = 14, fontsize_row = 14, fontsize_col = 14, treeheight_row = 0, clustering_method = "average", annotation_col = hm_groups, annotation_names_col = FALSE, annotation_colors = hm_cols, annotation_legend = FALSE, legend=TRUE, main = "IGK"))
grid.text("B", x=0.24, y=0.73, gp=gpar(fontsize=30))

## E2-nonreactive B cells
share_negK <- shareK[c(1,3,5,7,9,11,13,15,17,19),c(1,3,5,7,9,11,13,15,17,19)]

colnames(share_negK) <- gsub("n", "", colnames(share_negK))
rownames(share_negK) <- gsub("n", "", rownames(share_negK))

# Make heatmap
share_negK_hm <- as.ggplot(pheatmap(share_negK, scale="none", cluster_rows = TRUE, cluster_cols=TRUE, show_rownames = TRUE, show_colnames = TRUE,  breaks = seq(0, 0.08, length.out = 100), border_color = NA, cellheight = 14, cellwidth = 14, fontsize = 14, fontsize_row = 14, fontsize_col = 14, treeheight_row = 0, clustering_method = "average", annotation_col = hm_groups, annotation_names_col = FALSE, annotation_colors = hm_cols, annotation_legend = FALSE, legend=FALSE))

#### IGL
### Quantify total number of clonotypes for each sample
seq_totalsL <- data.frame()
for(i in 4:23){
  ct_clonoL <- filter(cdr3L_seq_VJ, cdr3L_seq_VJ[,i] >0)
  add_seqL <- cbind.data.frame(subject=samples[i-3], num_clonos = nrow(ct_clonoL))
  seq_totalsL <- rbind.data.frame(seq_totalsL, add_seqL)
}

### Quantify sharing between each subject (normalize by dividing shared by nonshared)
shareL <- data.frame()
for(i in 4:23){  
  filt_subjL <- filter(cdr3L_seq_VJ, cdr3L_seq_VJ[,i] > 0)
  share_quantL <- colSums(filt_subjL[,4:23] != 0)
  denomL <- seq_totalsL[,2] + seq_totalsL[i-3,2] - share_quantL
  norm_shareL <- share_quantL/denomL %>%
    replace(., i-3, NA)
  shareL <- rbind.data.frame(shareL, norm_shareL)
}

colnames(shareL) <- colnames(cdr3L_seq_VJ_raw[,4:23])
rownames(shareL) <- colnames(cdr3L_seq_VJ_raw[,4:23])

### Plot
## E2-reactive B cells
share_posL <- shareL[c(2,4,6,8,10,12,14,16,18,20),c(2,4,6,8,10,12,14,16,18,20)]

colnames(share_posL) <- gsub("p", "", colnames(share_posL))
rownames(share_posL) <- gsub("p", "", rownames(share_posL))

 # Make heatmap
share_posL_hm <- as.ggplot(pheatmap(share_posL, scale="none", cluster_rows = TRUE, cluster_cols=TRUE, show_rownames = TRUE, show_colnames = TRUE, breaks = seq(0, 0.09, length.out = 100), border_color = NA, cellheight = 14, cellwidth = 14, fontsize = 14, fontsize_row = 14, fontsize_col = 14, treeheight_row = 0, clustering_method = "average", annotation_col = hm_groups, annotation_names_col = FALSE, annotation_colors = hm_cols, annotation_legend = FALSE, legend=TRUE, main = "IGL"))
grid.text("C", x=0.24, y=0.73, gp=gpar(fontsize=30))

## E2-nonreactive B cells
share_negL <- shareL[c(1,3,5,7,9,11,13,15,17,19),c(1,3,5,7,9,11,13,15,17,19)]

colnames(share_negL) <- gsub("n", "", colnames(share_negL))
rownames(share_negL) <- gsub("n", "", rownames(share_negL))

# Make heatmap
share_negL_hm <- as.ggplot(pheatmap(share_negL, scale="none", cluster_rows = TRUE, cluster_cols=TRUE, show_rownames = TRUE, show_colnames = TRUE, breaks = seq(0, 0.09, length.out = 100), border_color = NA, cellheight = 14, cellwidth = 14, fontsize = 14, fontsize_row = 14, fontsize_col = 14, treeheight_row = 0, clustering_method = "average", annotation_col = hm_groups, annotation_names_col = FALSE, annotation_colors = hm_cols, annotation_legend = FALSE, legend=FALSE))

  # with legend
share_negL_hm2 <- as.ggplot(pheatmap(share_negL, scale="none", cluster_rows = TRUE, cluster_cols=TRUE, show_rownames = TRUE, show_colnames = TRUE, breaks = seq(0, 0.09, length.out = 100), border_color = NA, cellheight = 14, cellwidth = 14, fontsize = 14, fontsize_row = 14, fontsize_col = 14, treeheight_row = 0, clustering_method = "average", annotation_col = hm_groups, annotation_names_col = FALSE, annotation_colors = hm_cols, annotation_legend = TRUE, legend=TRUE, main = "E2+  E2-"))

#### Arrange all heatmaps (IGH, IGK, and IGL, both E2-reactive and nonreactive) in one plot
grid.arrange(share_posG_hm, share_posK_hm, share_posL_hm, share_negG_hm, share_negK_hm, share_negL_hm, ncol=3)


##### DIVERSITY
#### IGH
### Quantify Shannon entropy and species richness for each subject (normalized by max possible diversity for each subject)
immdataG3_df <- data.frame()
for(i in 1:20){
  subj_data <- immdataG3$data[[i]] %>%
    mutate(V.name = gsub("\\*.*", "", V.name)) %>%
    mutate(J.name = gsub("\\*.*", "", J.name)) %>%
    mutate(clonotype = paste(CDR3.aa, V.name, J.name)) %>%
    group_by(clonotype) %>%
    mutate(clono_sum = n()) 
  subj_name <- sub("\\..*", "", names(immdataG3$data[i]))
  immdataG3_df <- rbind.data.frame(immdataG3_df, mutate(subj_data, subj = subj_name), stringsAsFactors = FALSE)
}

HCDR3_div <- immdataG3_df %>%
  group_by(subj) %>%
  mutate(clono_max =1) %>%
  mutate(max_div = diversity(clono_max)) %>%
  distinct(subj, clonotype, clono_sum, .keep_all = TRUE) %>%
  group_by(subj) %>%
  mutate(shan_div = diversity(clono_sum)) %>%
  mutate(norm_div = shan_div/max_div) %>%
  mutate(num_clon = n()) %>%
  mutate(richness = rarefy(clono_sum, unique(num_clon))) %>%
  mutate(norm_rich = richness/num_clon) %>%
  distinct(subj, norm_div, norm_rich) %>%
  mutate(E2 = ifelse(str_sub(subj,-1,-1) == "p", "pos", "neg")) %>%
  mutate(Neut = ifelse(str_sub(subj, 1, nchar(subj)-1) %in% c("C48", "C117", "C172", "P53", "P54"), "hneut", "lneut")) %>%
  mutate(E2_neut = paste(Neut, E2, sep = "_"))

### Statistical analysis of differences between groups
kruskal.test(norm_div ~ E2_neut, data = HCDR3_div)
kruskal.test(norm_rich ~ E2_neut, data = HCDR3_div)

### Plot
## Specify order of groups
HCDR3_div$E2_neut <- factor(HCDR3_div$E2_neut,levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

## Make plot
# Shannon entropy
HCDR3_shan_p <- ggplot(HCDR3_div, aes(x= E2_neut, y = norm_div, color= E2_neut)) +
  geom_boxplot(outlier.color = NA, show.legend = FALSE, lwd=1)+
  geom_point(size=4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray")) +
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  ggtitle("A")+
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=22),
        axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.185, vjust = -1.3),
        aspect.ratio=1,
        axis.line = element_line(size = 0.75, color = "black"))+
  labs(x= "",
       y= "Normalized Shannon Diversity of IGH")
HCDR3_shan_p

# species richness
HCDR3_rich_p <- ggplot(HCDR3_div, aes(x= E2_neut, y = norm_rich, color= E2_neut)) +
  geom_boxplot(outlier.color = NA, show.legend = FALSE, lwd=1)+
  geom_point(size=4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray")) +
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  ggtitle("C")+
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=22),
        axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.16, vjust = -1.3),
        aspect.ratio=1,
        axis.line = element_line(size = 0.75, color = "black"))+
  labs(x= "",
       y= "Normalized Richness of IGH")
HCDR3_rich_p

#### IGKL
### Quantify Shannon entropy and species richness for each subject (normalized by max possible diversity for each subject)
## IGK
immdataK3_df <- data.frame()
for(i in 1:20){
  subj_data <- immdataK3$data[[i]] %>%
    mutate(V.name = gsub("\\*.*", "", V.name)) %>%
    mutate(J.name = gsub("\\*.*", "", J.name)) %>%
    mutate(clonotype = paste(CDR3.aa, V.name, J.name)) %>%
    group_by(clonotype) %>%
    mutate(clono_sum = n()) 
  subj_name <- sub("\\..*", "", names(immdataK3$data[i]))
  immdataK3_df <- rbind.data.frame(immdataK3_df, mutate(subj_data, subj = subj_name, ), stringsAsFactors = FALSE)
}

LCDR3_divK <- immdataK3_df %>%
  group_by(subj) %>%
  mutate(clono_max =1) %>%
  mutate(max_div = diversity(clono_max)) %>%
  distinct(subj, clonotype, clono_sum, .keep_all = TRUE) %>%
  group_by(subj) %>%
  mutate(shan_div = diversity(clono_sum)) %>%
  mutate(norm_div = shan_div/max_div) %>%
  mutate(num_clon = n()) %>%
  mutate(richness = rarefy(clono_sum, unique(num_clon))) %>%
  mutate(norm_rich = richness/num_clon) %>%
  distinct(subj, norm_div, norm_rich) %>%
  mutate(E2 = ifelse(str_sub(subj,-1,-1) == "p", "pos", "neg")) %>%
  mutate(Neut = ifelse(str_sub(subj, 1, nchar(subj)-1) %in% c("C48", "C117", "C172", "P53", "P54"), "hneut", "lneut")) %>%
  mutate(E2_neut = paste(Neut, E2, sep = "_")) 

# Statistical analysis of differences between groups 
kruskal.test(norm_div ~ E2_neut, data = LCDR3_divK)
kruskal.test(norm_rich ~ E2_neut, data = LCDR3_divK)

## IGL
immdataL3_df <- data.frame()
for(i in 1:20){
  subj_data <- immdataL3$data[[i]] %>%
    mutate(V.name = gsub("\\*.*", "", V.name)) %>%
    mutate(J.name = gsub("\\*.*", "", J.name)) %>%
    mutate(clonotype = paste(CDR3.aa, V.name, J.name)) %>%
    group_by(clonotype) %>%
    mutate(clono_sum = n()) 
  subj_name <- sub("\\..*", "", names(immdataL3$data[i]))
  immdataL3_df <- rbind.data.frame(immdataL3_df, mutate(subj_data, subj = subj_name, ), stringsAsFactors = FALSE)
}

LCDR3_divL <- immdataL3_df %>%
  group_by(subj) %>%
  mutate(clono_max =1) %>%
  mutate(max_div = diversity(clono_max)) %>%
  distinct(subj, clonotype, clono_sum, .keep_all = TRUE) %>%
  group_by(subj) %>%
  mutate(shan_div = diversity(clono_sum)) %>%
  mutate(norm_div = shan_div/max_div) %>%
  mutate(num_clon = n()) %>%
  mutate(richness = rarefy(clono_sum, unique(num_clon))) %>%
  mutate(norm_rich = richness/num_clon) %>%
  distinct(subj, norm_div, norm_rich) %>%
  mutate(E2 = ifelse(str_sub(subj,-1,-1) == "p", "pos", "neg")) %>%
  mutate(Neut = ifelse(str_sub(subj, 1, nchar(subj)-1) %in% c("C48", "C117", "C172", "P53", "P54"), "hneut", "lneut")) %>%
  mutate(E2_neut = paste(Neut, E2, sep = "_")) 

# Statistical analysis of differences between groups
kruskal.test(norm_div ~ E2_neut, data = LCDR3_divL)
kruskal.test(norm_rich ~ E2_neut, data = LCDR3_divL)

### Combine IGK and IGL for plotting
LCDR3_div <- rbind.data.frame(mutate(LCDR3_divK, chain="IGK"),
                              mutate(LCDR3_divL, chain="IGL"), stringsAsFactors = FALSE)

### Plot
## Specify order of groups
LCDR3_div$E2_neut <- factor(LCDR3_div$E2_neut,levels = c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"))

## Make plots
# Shannon entropy
LCDR3_shan_p <- ggplot(LCDR3_div, aes(x= E2_neut, y = norm_div, color= E2_neut)) +
  geom_boxplot(outlier.color = NA, show.legend = FALSE, lwd=1)+
  geom_point(size=4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray")) +
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  facet_wrap(~ chain, ncol=2)+
  ggtitle("B")+
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=22),
        axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.75, color = "black"),
        plot.title = element_text(size = 45, hjust = -0.085, vjust = -1.3),
        strip.text.x = element_text(size=16,face="bold.italic"))+
  labs(x= "",
       y= "Normalized Shannon Diversity")
LCDR3_shan_p

# Species richness
LCDR3_rich_p <- ggplot(LCDR3_div, aes(x= E2_neut, y = norm_rich, color= E2_neut)) +
  geom_boxplot(outlier.color = NA, show.legend = FALSE, lwd=1)+
  geom_point(size=4)+
  scale_color_manual(values=c("#009999", "dark gray", "#b66dff", "dark gray")) +
  scale_x_discrete(breaks=c("hneut_pos", "hneut_neg", "lneut_pos", "lneut_neg"),
                   labels=c("High Neut\nE2+", "High Neut\nE2-", "Low Neut\nE2+", "Low Neut\nE2-"))+
  facet_wrap(~ chain, ncol=2)+
  ggtitle("D")+
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=22),
        axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.75, color = "black"),
        plot.title = element_text(size = 45, hjust = -0.1, vjust = -1.3),
        strip.text.x = element_text(size=16,face="bold.italic"))+
  labs(x= "",
       y= "Normalized Richness")
LCDR3_rich_p


##### FREQUENCY OF E2-REACTIVE B CELLS
#### Compare frequency of E2-reactive B cells between groups
### Input flow cytometry data
E2freq <- data.frame(Subject = c("C48", "C117", "C172", "C176", "C429", "P49", "P157", "P53", "P54", "P111", "LKP1", "LKP2", "LKP3", "LKP4", "LKP5", "LKP6", "LKP7", "LKP8", "LKP9"), 
                     Neut = c("High Neut", "High Neut", "High Neut", "Low Neut", "Low Neut", "Low Neut", "Low Neut", "High Neut", "High Neut", "Low Neut", "Healthy Control", "Healthy Control", "Healthy Control", "Healthy Control", "Healthy Control", "Healthy Control", "Healthy Control", "Healthy Control", "Healthy Control"), 
                     E2_freq = c(0.38, 0.75, 0.71, 0.13, 0.13, 0.11, 0.17, 1.13, 2.76, 0.15, 0.051, 0.056, 0.062, 0.11, 0.029, 0.081, 0, 0.086, 0.1)) %>%
  group_by(Neut) %>%
  mutate(neut_mean = mean(E2_freq)) %>%
  mutate(neut_se = (sd(E2_freq))/sqrt(n())) %>%
  ungroup()

### Statistical analysis of differences between groups
## Test for normality and equal variance
shapiro.test(E2freq[E2freq$Neut=="High Neut",]$E2_freq)
shapiro.test(E2freq[E2freq$Neut=="Low Neut",]$E2_freq)
var.test(E2freq[E2freq$Neut=="High Neut",]$E2_freq, E2freq[E2freq$Neut=="Low Neut",]$E2_freq) #fail

## Compare groups
kruskal.test(E2_freq ~ Neut, data = E2freq)
E2freq_stat <- dunn_test(E2freq, E2_freq ~ Neut, p.adjust.method = "BH")

### Plot
## Specify order of groups
E2freq$Neut <- factor(E2freq$Neut, levels = c("High Neut", "Low Neut", "Healthy Control"))

## Make plot
E2_freq_p <- ggplot(E2freq, aes(x=Neut, y=E2_freq)) +
  geom_boxplot(aes(), color = c("#676767", "#676767", "#676767"), show.legend = FALSE, lwd=1.25)+
  theme(axis.title.x = element_blank(), 
        axis.text=element_text(color="black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"))+
  geom_dotplot(aes(x=Neut,y=E2_freq, fill = Neut, color = Neut), binaxis = "y", stackdir = "center", binwidth = 0.1, dotsize=0.75)+
  scale_color_manual(values = c("#009999", "#b66dff", "light gray"), labels = c("High Neut", "Low Neut", "Healthy Control"))+
  scale_fill_manual(values=c("#009999", "#b66dff", "light gray")) +
  scale_x_discrete(labels = c("Healthy Control" = "HC"))+
  stat_pvalue_manual(E2freq_stat, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(2.9, 3.1), tip.length = 0.005, vjust=0.7, size=11, bracket.size = 0.75) +
  ggtitle("A") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text=element_text(size=20),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.1, vjust = -1),
        aspect.ratio = 1.125)+
  labs(x= "",
       y= "Freq. E2+ B cells (%)")
E2_freq_p

#### Determine association between freqency of E2-reactive B cells and neutralization score
### Input neutralization scores
E2freq_neut <- E2freq[1:10,] %>%
  mutate(Neut_score = c(12, 28,13, 1, 1, 1, 2, 11, 19, 2))

### Plot
E2freq_neut_p <- ggplot(E2freq_neut, aes(x=Neut_score, y=E2_freq)) +
  stat_cor(method = "kendall", size=7, label.x.npc="middle", label.y=3, digits=1)+
  geom_smooth(method = "lm", se = TRUE, color = "black", size=1, fill="light gray")+
  geom_point(size=3, aes(color = Neut), position=position_jitter(h=0.085,w=0.085)) +
  scale_color_manual(values = c("#009999", "#b66dff"))+
  scale_x_continuous(breaks = seq(from = 0, to = 28, by =5), limits = c(0,28.5))+
  ggtitle("B") +
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text=element_text(size=18),
        legend.position="right",
        axis.line = element_line(size = 0.75, color = "black"),
        plot.title = element_text(size = 45, hjust = -0.1, vjust = -1),
        aspect.ratio = 0.85)+
  labs(x= "Neutralization Score",
       y= "% E2+ B cells")
E2freq_neut_p


##### COMPARISON OF E2-REACTIVE PUBLIC CLONOTYPES TO THE BULK POPULATION 
#### CDR3 Length
### IGH
## Quantify CDRH3 length for public clonotypes and bulk population
# Public clonotypes
pubG_len <- as.data.frame(rbind.data.frame(mutate(hneut_pos100VG[,1:3], Group = "hneut_pos"),
                                           mutate(lneut_pos100VG[,1:3], Group = "lneut_pos"), stringsAsFactors = FALSE)) %>%
  mutate(len = nchar(CDR3.aa))

# Bulk population
pubG_all_len <- as.data.frame(rbind.data.frame(mutate(filter(HCDR3_len, hneut_pos>0), Group = "hneut_pos_all"),
                                               mutate(filter(HCDR3_len, lneut_pos>0), Group = "lneut_pos_all"), stringsAsFactors = FALSE)) %>%
  filter(!(paste(V.name, J.name, CDR3.aa) %in% paste(pubG_len$V.name, pubG_len$J.name, pubG_len$CDR3.aa))) %>%
  mutate(len = nchar(CDR3.aa)) %>%
  select(CDR3.aa, V.name, J.name, Group, len)

# Public and bulk populations combined
pubG_comp_len <- as.data.frame(rbind.data.frame(pubG_len, pubG_all_len, stringsAsFactors = FALSE))

## Statistical analysis of differences between groups
pubG_comp_len_stat <- as.data.frame(pubG_comp_len) %>%
  dunn_test(len ~ Group, 
            p.adjust.method = "BH")

## Plot
HCDR3_comp_pub_len_p <- ggplot(pubG_comp_len, aes(x=Group, y = len)) +
  geom_point(aes(color = Group), shape =1, stroke = 1, position = position_jitter(width=0.35), size=1.5)+
  geom_violin(color="#676767", show.legend = FALSE, lwd=1.2, alpha = 0.4) +
  geom_boxplot(outlier.color = NA, color="#676767", show.legend = FALSE, lwd=1, alpha = 0.5, width = 0.25)+
  scale_color_manual(values=c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_x_discrete(breaks=c("hneut_pos", "hneut_pos_all", "lneut_pos", "lneut_pos_all"),
                   labels=c("High Neut\nPublic", "High Neut\nAll", "Low Neut\nPublic", "Low Neut\nAll"))+
  scale_y_continuous(breaks = seq(from = 5, to = 38.95, by = 5), limits = c(5,38.95))+
  ggtitle("A") +
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.12, vjust = -0.25),
        aspect.ratio=1.1)+
  stat_pvalue_manual(pubG_comp_len_stat, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(35, 36.5), tip.length = 0.005,vjust = 0.7, size=10, bracket.size = 0.75) +
  labs(x= "",
       y= "CDRH3 AA Length")
HCDR3_comp_pub_len_p 

### IGKL
## Quantify CDRL3 length for public clonotypes and bulk population
# Public clonotypes
pubKL_len <- as.data.frame(rbind.data.frame(mutate(hneut_pos100VK[,1:3], Group = "hneut_pos"),
                                            mutate(hneut_pos100VL[,1:3], Group = "hneut_pos"),
                                            mutate(lneut_pos100VK[,1:3], Group = "lneut_pos"),
                                            mutate(lneut_pos100VL[,1:3], Group = "lneut_pos"),
                                            stringsAsFactors = FALSE)) %>%
  mutate(len = nchar(CDR3.aa))


# Bulk population
pubKL_all_len <- as.data.frame(rbind.data.frame(mutate(filter(LCDR3_lenK, hneut_pos>0), Group = "hneut_pos_all"),
                                                mutate(filter(LCDR3_lenL, hneut_pos>0), Group = "hneut_pos_all"),
                                                mutate(filter(LCDR3_lenK, lneut_pos>0), Group = "lneut_pos_all"),
                                                mutate(filter(LCDR3_lenL, lneut_pos>0), Group = "lneut_pos_all"),
                                                stringsAsFactors = FALSE)) %>%
  filter(!(paste(V.name, J.name, CDR3.aa) %in% paste(pubKL_len$V.name, pubKL_len$J.name, pubKL_len$CDR3.aa))) %>%
  mutate(len = nchar(CDR3.aa)) %>%
  select(CDR3.aa, V.name, J.name, Group, len)

# Public and bulk populations combined
pubKL_comp_len <- as.data.frame(rbind.data.frame(pubKL_len, pubKL_all_len, stringsAsFactors = FALSE))

## Statistical analysis of differences between groups
pubKL_comp_len_stat <- as.data.frame(pubKL_comp_len) %>%
  dunn_test(len ~ Group, 
            p.adjust.method = "BH")

## Plot
LCDR3_comp_pub_len_p <- ggplot(pubKL_comp_len, aes(x=Group, y = len)) +
  geom_point(aes(color = Group), shape =1, stroke = 1, position = position_jitter(width=0.35), size=1.5)+
  geom_violin(color="#676767", show.legend = FALSE, lwd=1.2, alpha = 0.4) +
  geom_boxplot(outlier.color = NA, color="#676767", show.legend = FALSE, lwd=1, alpha = 0.5, width = 0.25)+
  scale_color_manual(values=c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_x_discrete(breaks=c("hneut_pos", "hneut_pos_all", "lneut_pos", "lneut_pos_all"),
                   labels=c("High Neut\nPublic", "High Neut\nAll", "Low Neut\nPublic", "Low Neut\nAll"))+
  scale_y_continuous(breaks = seq(from = 5, to = 21, by = 5), limits = c(5,21))+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.12, vjust = -0.25),
        aspect.ratio=1.1)+
  labs(x= "",
       y= "CDRL3 AA Length")
LCDR3_comp_pub_len_p 


#### V-gene usage
### IGH
## Quantify IGHV gene usage for public clonotypes and the bulk population
Vpub_comp_statG <- pubG_comp_len %>%
  group_by(V.name, Group) %>%
  mutate(V_count = n()) %>%
  mutate(V.name=str_remove(V.name, "IGHV")) %>%
  select(V.name, Group, V_count) %>%
  distinct() %>%
  ungroup() %>%
  complete(V.name, Group, fill = list(V_count = 0)) %>%
  group_by(Group) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total)

## Statistical analysis of differences between groups
list_pub_compVG <- unique(Vpub_comp_statG$V.name)

Vpub_comp_statsG <- data.frame()
for(i in 1:length(list_pub_compVG)){
  VG <- filter(Vpub_comp_statG, V.name == list_pub_compVG[i])
  hnp_hnpa <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[1], VG$Total[2]-VG$V_count[2], VG$Total[1]-VG$V_count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(VG$V_count[1], VG$V_count[3], VG$Total[1]-VG$V_count[1], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  lnp_lnpa <- fisher.test(matrix(c(VG$V_count[4], VG$V_count[3], VG$Total[4]-VG$V_count[4], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  hnp_lnpa <- fisher.test(matrix(c(VG$V_count[1], VG$V_count[4], VG$Total[1]-VG$V_count[1], VG$Total[4]-VG$V_count[4]), ncol=2))$p.value
  lnp_hnpa <- fisher.test(matrix(c(VG$V_count[3], VG$V_count[2], VG$Total[3]-VG$V_count[3], VG$Total[2]-VG$V_count[2]), ncol=2))$p.value
  hnpa_lnpa <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[4], VG$Total[2]-VG$V_count[2], VG$Total[4]-VG$V_count[4]), ncol=2))$p.value
  VG_df <- cbind.data.frame(V.name = list_pub_compVG[i], hnp_hnpa, hnp_lnp, lnp_lnpa, hnp_lnpa, lnp_hnpa, hnpa_lnpa)
  Vpub_comp_statsG <- rbind.data.frame(Vpub_comp_statsG, VG_df, stringsAsFactors = FALSE)
}

Vpub_comp_sigG <-data.frame(Vpub_comp_statsG$V.name, matrix(p.adjust(as.vector(as.matrix(Vpub_comp_statsG[,-1])), method='BH'),ncol=6))
colnames(Vpub_comp_sigG) <- colnames(Vpub_comp_statsG)

Vpub_comp_filtG <- Vpub_comp_sigG %>%
  filter_all(any_vars(. <0.05)) %>%
  pivot_longer(., c("hnp_hnpa", "hnp_lnp", "lnp_lnpa", "hnp_lnpa", "lnp_hnpa", "hnpa_lnpa"),
               names_to = "Comp",
               values_to = "p.adj") %>%
  mutate(.y = "Prop") %>%
  filter(p.adj < 0.05) %>%
  filter(Comp != "hnpa_lnpa")

## Plot
Vpub_compG <- Vpub_comp_statG %>%
  dplyr::select(V.name, Group, V_count, Prop) %>%
  reshape2::dcast(., Group ~ V.name, value.var = "Prop")
row.names(Vpub_compG) <- c("High Neut Pub.",  "High Neut All", "Low Neut Pub.",  "Low Neut All" )

Vpub_compG_hm <- Vpub_compG[,-1]

# Adjust the order of V-genes
reorder_VG_pub <- data.frame(vname = colnames(Vpub_compG_hm )) %>%
  mutate(numbering = sub("[0-9]-", "", vname)) %>%
  mutate(numbering = sub("-[0-9]", "", numbering)) %>%
  mutate(fam = sub("-.*", "", vname)) %>%
  group_by(fam) %>%
  arrange(as.numeric(numbering), .by_group = TRUE)

VgenesG_pub_arr <- t(Vpub_compG_hm)[reorder_VG_pub$vname,]

pheatmap(VgenesG_pub_arr, scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 13, labels_row = make_italics(rownames(VgenesG_pub_arr)), main = "")
grid.text("B", x=0.325, y=0.94, gp=gpar(fontsize=30))


### IGKL
## Quantify IGKV/IGLV gene usage for public clonotypes and the bulk population
Vpub_comp_statKL <- pubKL_comp_len %>%
  group_by(V.name, Group) %>%
  mutate(V_count = n()) %>%
  select(V.name, Group, V_count) %>%
  distinct() %>%
  ungroup() %>%
  complete(V.name, Group, fill = list(V_count = 0)) %>%
  group_by(Group) %>%
  mutate(Total = sum(V_count)) %>%
  mutate(Prop = V_count/Total)

## Statistical analysis of differences between groups
list_pub_compVKL <- unique(Vpub_comp_statKL$V.name)

Vpub_comp_statsKL <- data.frame()
for(i in 1:length(list_pub_compVKL)){
  VG <- filter(Vpub_comp_statKL, V.name == list_pub_compVKL[i])
  hnp_hnpa <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[1], VG$Total[2]-VG$V_count[2], VG$Total[1]-VG$V_count[1]), ncol=2))$p.value
  hnp_lnp <- fisher.test(matrix(c(VG$V_count[1], VG$V_count[3], VG$Total[1]-VG$V_count[1], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  lnp_lnpa <- fisher.test(matrix(c(VG$V_count[4], VG$V_count[3], VG$Total[4]-VG$V_count[4], VG$Total[3]-VG$V_count[3]), ncol=2))$p.value
  hnp_lnpa <- fisher.test(matrix(c(VG$V_count[1], VG$V_count[4], VG$Total[1]-VG$V_count[1], VG$Total[4]-VG$V_count[4]), ncol=2))$p.value
  lnp_hnpa <- fisher.test(matrix(c(VG$V_count[3], VG$V_count[2], VG$Total[3]-VG$V_count[3], VG$Total[2]-VG$V_count[2]), ncol=2))$p.value
  hnpa_lnpa <- fisher.test(matrix(c(VG$V_count[2], VG$V_count[4], VG$Total[2]-VG$V_count[2], VG$Total[4]-VG$V_count[4]), ncol=2))$p.value
  VG_df <- cbind.data.frame(V.name = list_pub_compVKL[i], hnp_hnpa, hnp_lnp, lnp_lnpa, hnp_lnpa, lnp_hnpa, hnpa_lnpa)
  Vpub_comp_statsKL <- rbind.data.frame(Vpub_comp_statsKL, VG_df, stringsAsFactors = FALSE)
}

Vpub_comp_sigKL <-data.frame(Vpub_comp_statsKL$V.name, matrix(p.adjust(as.vector(as.matrix(Vpub_comp_statsKL[,-1])), method='BH'),ncol=6))
colnames(Vpub_comp_sigKL) <- colnames(Vpub_comp_statsKL)

Vpub_comp_filtKL <- Vpub_comp_sigKL %>%
  filter_all(any_vars(. <0.05)) %>%
  pivot_longer(., c("hnp_hnpa", "hnp_lnp", "lnp_lnpa", "hnp_lnpa", "lnp_hnpa", "hnpa_lnpa"),
               names_to = "Comp",
               values_to = "p.adj") %>%
  mutate(.y = "Prop") %>%
  filter(p.adj < 0.05) %>%
  filter(Comp != "hnpa_lnpa")

## Plot
Vpub_compKL <- Vpub_comp_statKL %>%
  dplyr::select(V.name, Group, V_count, Prop) %>%
  reshape2::dcast(., Group ~ V.name, value.var = "Prop")
row.names(Vpub_compKL) <- c("High Neut Pub.",  "High Neut All", "Low Neut Pub.",  "Low Neut All" )

Vpub_compKL_hm <- Vpub_compKL[,-1]

# Adjust the order of V-genes
reorder_VKL_pub <- data.frame(vname_full = colnames(Vpub_compKL_hm), vname = sub("IG[A-Z]V", "", colnames(Vpub_compKL_hm))) %>%
  mutate(numbering = sub(".{1,2}-", "", vname)) %>%
  mutate(fam = sub("-.*", "", vname)) %>%
  mutate(chain = substr(vname_full, 1,4)) %>%
  group_by(chain, fam) %>%
  arrange(as.numeric(numbering), .by_group = TRUE)

VgenesKL_pub_arr <- t(Vpub_compKL_hm)[reorder_VKL_pub$vname_full,]

# IGK
pheatmap(VgenesKL_pub_arr[1:45,], scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 13, labels_row = make_italics(gsub("IGKV", "", rownames(VgenesKL_pub_arr)[1:45])), main = "")
grid.text("C", x=0.325, y=0.96, gp=gpar(fontsize=30))

# IGL
pheatmap(VgenesKL_pub_arr[46:81,], scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 12, cellwidth = 18, fontsize = 13, fontsize_row = 13, fontsize_col = 13, labels_row = make_italics(gsub("IGLV", "", rownames(VgenesKL_pub_arr)[46:81])), main = "")



#### Somatic hypermutation of the V-gene
### IGH
## Quantify SHM for public clonotypes and bulk population
# Public clonotypes defined in CDR3 length analysis
hneut_pos_pubG_SHM <- IMGT_IGH_all_grp %>%
  filter(paste(V.name, J.name, CDR3.aa) %in% paste(hneut_pos100VG$V.name, hneut_pos100VG$J.name, hneut_pos100VG$CDR3.aa)) %>%
  mutate(group = "hneut_pos")

lneut_pos_pubG_SHM <- IMGT_IGH_all_grp %>%
  filter(paste(V.name, J.name, CDR3.aa) %in% paste(lneut_pos100VG$V.name, lneut_pos100VG$J.name, lneut_pos100VG$CDR3.aa)) %>%
  mutate(group = "lneut_pos")

# Bulk population
pubG_comp_hneut <- IMGT_IGH_all_grp %>%
  filter(!(paste(V.name, J.name, CDR3.aa) %in% paste(hneut_pos_pubG_SHM$V.name, hneut_pos_pubG_SHM$J.name, hneut_pos_pubG_SHM$CDR3.aa))) %>%
  filter(hneut_pos > 0) %>%
  mutate(group = "hneut_pos_all")

pubG_comp_lneut <- IMGT_IGH_all_grp %>%
  filter(!(paste(V.name, J.name, CDR3.aa) %in% paste(lneut_pos_pubG_SHM$V.name, lneut_pos_pubG_SHM$J.name, lneut_pos_pubG_SHM$CDR3.aa))) %>%
  filter(lneut_pos > 0) %>%
  mutate(group = "lneut_pos_all")

# Public and bulk populations combined
pubG_comp_SHM <- rbind.data.frame(pubG_comp_hneut, pubG_comp_lneut, hneut_pos_pubG_SHM, lneut_pos_pubG_SHM, stringsAsFactors = FALSE)


## Statistical analysis of differences between groups
pubG_comp_SHM_stat <- as.data.frame(pubG_comp_SHM) %>%
  dunn_test(V_identity ~ group, 
            p.adjust.method = "BH")

## Plot
pubG_comp_SHM_p <- ggplot(pubG_comp_SHM, aes(x= group, y = V_identity)) +
  geom_point(aes(color = group), shape =1, stroke = 1, position = position_jitter(width=0.35), size=1.5)+
  geom_violin(color="#676767", show.legend = FALSE, lwd=1.2, alpha = 0.4) +
  geom_boxplot(outlier.color = NA, color="#676767", show.legend = FALSE, lwd=1, alpha = 0.5, width = 0.25)+
  scale_color_manual(values=c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_x_discrete(breaks=c("hneut_pos", "hneut_pos_all", "lneut_pos", "lneut_pos_all"),
                   labels=c("High Neut\nPublic", "High Neut\nAll", "Low Neut\nPublic", "Low Neut\nAll"))+
  ggtitle("D") +
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.165, vjust = -0.25),
        aspect.ratio=1.1)+
  stat_pvalue_manual(pubG_comp_SHM_stat, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(100.5), tip.length = 0.005, vjust=0.6, size=12, bracket.size = 1) +
  labs(x= "",
       y= expression("%"~italic(V[H])~" Identity"))
pubG_comp_SHM_p


### IGKL
## Quantify SHM for public clonotypes and bulk population
# Public clonotypes
hneut_pos_pubKL_SHM <- IMGT_IGKL_all_grp %>%
  filter(paste(V.name, J.name, CDR3.aa) %in% paste(hneut_pos100VK$V.name, hneut_pos100VK$J.name, hneut_pos100VK$CDR3.aa) |
           paste(V.name, J.name, CDR3.aa) %in% paste(hneut_pos100VL$V.name, hneut_pos100VL$J.name, hneut_pos100VL$CDR3.aa)) %>%
  mutate(group = "hneut_pos")

lneut_pos_pubKL_SHM <- IMGT_IGKL_all_grp %>%
  filter(paste(V.name, J.name, CDR3.aa) %in% paste(lneut_pos100VK$V.name, lneut_pos100VK$J.name, lneut_pos100VK$CDR3.aa) |
           paste(V.name, J.name, CDR3.aa) %in% paste(lneut_pos100VL$V.name, lneut_pos100VL$J.name, lneut_pos100VL$CDR3.aa)) %>%
  mutate(group = "lneut_pos")

# Bulk population
pubKL_comp_hneut <- IMGT_IGKL_all_grp %>%
  filter(!(paste(V.name, J.name, CDR3.aa) %in% paste(hneut_pos_pubKL_SHM$V.name, hneut_pos_pubKL_SHM$J.name, hneut_pos_pubKL_SHM$CDR3.aa))) %>%
  filter(hneut_pos > 0) %>%
  mutate(group = "hneut_pos_all")

pubKL_comp_lneut <- IMGT_IGKL_all_grp %>%
  filter(!(paste(V.name, J.name, CDR3.aa) %in% paste(lneut_pos_pubKL_SHM$V.name, lneut_pos_pubKL_SHM$J.name, lneut_pos_pubKL_SHM$CDR3.aa))) %>%
  filter(lneut_pos > 0) %>%
  mutate(group = "lneut_pos_all")

# Public and bulk populations combined
pubKL_comp_SHM <- rbind.data.frame(pubKL_comp_hneut, pubKL_comp_lneut, hneut_pos_pubKL_SHM, lneut_pos_pubKL_SHM, stringsAsFactors = FALSE)


## Statistical analysis of differences between groups
pubKL_comp_SHM_stat <- as.data.frame(pubKL_comp_SHM) %>%
  dunn_test(V_identity ~ group, 
            p.adjust.method = "BH")

## Plot
pubKL_comp_SHM_p <- ggplot(pubKL_comp_SHM, aes(x= group, y = V_identity)) +
  geom_point(aes(color = group), shape =1, stroke = 1, position = position_jitter(width=0.35), size=1.5)+
  geom_violin(color="#676767", show.legend = FALSE, lwd=1.2, alpha = 0.4) +
  geom_boxplot(outlier.color = NA, color="#676767", show.legend = FALSE, lwd=1, alpha = 0.5, width = 0.25)+
  scale_color_manual(values=c("#009999", "#009999", "#b66dff", "#b66dff"))+
  scale_x_discrete(breaks=c("hneut_pos", "hneut_pos_all", "lneut_pos", "lneut_pos_all"),
                   labels=c("High Neut\nPublic", "High Neut\nAll", "Low Neut\nPublic", "Low Neut\nAll"))+
  theme(axis.text=element_text(color="black", size=20),
        panel.background = element_rect(fill = "#FDFDFD"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.75, color = "black"),
        axis.title=element_text(size=22),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = -0.15, vjust = -0.25),
        aspect.ratio=1.1)+
  stat_pvalue_manual(pubKL_comp_SHM_stat, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(102.5, 101.5, 100.5), tip.length = 0.005, vjust=0.6, size=12, bracket.size = 1) +
  labs(x= "",
       y= expression("%"~italic(V[K])~"/"~italic(V[L])~" Identity"))
pubKL_comp_SHM_p

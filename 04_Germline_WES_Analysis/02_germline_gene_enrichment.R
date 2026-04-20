library(dplyr)
library(tidyverse)
library(broom)
library(ggrepel)
library(BuenColors)

## Load the preprocessed one-hot encoded sample by pathogenic variant carrier status dataframe. 
germline_df = read.csv("../data/processed/germline/Jun2023_ASC_Case_Control_Sample_Merged_Germline.tsv",sep="\t",check.names=FALSE)

## Define columns to use as covariate
pc_columns = paste0("PCA",1:10)
pc_columns = c("PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6", "PCA7", "PCA8", "PCA9", "PCA10")
pc_formula = "PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10"

## Convert disease status to int and factors
germline_df$is_ASC = recode(as.factor(germline_df$is_ASC),"True"=1,"False"=0)
germline_df$inferred_sex = recode(as.factor(germline_df$inferred_sex),"F"=0,"M"=1)

## Check which genes have at least 2 carriers in cases
case_germline_df = germline_df %>% filter(is_ASC==1)
case_germline_cols = colnames(case_germline_df)[grepl( "Germline_" , colnames(case_germline_df))]
case_germline_cols
## Total number of genes to be tested
length(case_germline_cols)
case_germline_cols_case_enirched = case_germline_cols[colSums(case_germline_df[case_germline_cols])>=2]

## Germline columns to test
#germline_cols = colnames(germline_df)[grepl( "Germline_" , colnames(germline_df))]
## Remove those with colsum of 0 (how did it get through??)
#germline_cols = germline_cols[colSums(germline_df[germline_cols])!=0] #
#germline_cols = germline_cols[colSums(germline_df[germline_cols])>=5] #
#gene_burden_results = gene_burden_testing(germline_cols)


## Perform Gene Burden Test, takes a while (~0.25 seconds per test x ~1000 test ~= 250 seconds ~4-5 minutes)
gene_burden_testing = function(df,burden_cols) {
  gene_burden_results = burden_cols %>%
    map_df(~ {
      gene = .
      print(gene)
      model = glm(as.formula(paste("is_ASC ~", gene, "+",pc_formula, "+ inferred_sex")),
                  data = df, family = "binomial")
      ## Get model coefficients
      model_terms = tidy(model) %>% filter(term == gene)
      conf_int = confint.default(model)
      ## Get 95 CI for effect
      low_ci = conf_int[gene,"2.5 %"]
      high_ci = conf_int[gene,"97.5 %"]
      model_terms$ci_95_lower = low_ci
      model_terms$ci_95_upper = high_ci
      ## Count the number of carriers in contingency table
      value_counts = table(df[,c(gene,"is_ASC")])
      carrier_in_case = value_counts["1","1"]
      carrier_in_control = value_counts["1","0"]
      noncarrier_in_case = value_counts["0","1"]
      noncarrier_in_control = value_counts["0","0"]
      model_terms$carrier_in_case = carrier_in_case
      model_terms$carrier_in_control = carrier_in_control
      model_terms$noncarrier_in_case = noncarrier_in_case
      model_terms$noncarrier_in_control = noncarrier_in_control
      model_terms
    })
  ## Apply BH pvalue adjustment
  gene_burden_results$adjusted_pval = p.adjust(gene_burden_results$p.value,method="BH")
  gene_burden_results$neg_log_pval = -log(gene_burden_results$p.value)
  gene_burden_results$neg_log_adjusted_pval = -log(gene_burden_results$adjusted_pval)
  extract_middle = function(x) {strsplit(x, "_")[[1]][2]}
  gene_burden_results$gene = sapply(gene_burden_results$term,extract_middle)
  return(gene_burden_results)
}

## Visualize the results
case_filtered_gene_burden_results = gene_burden_testing(germline_df,case_germline_cols_case_enirched)

## Alternatively, read in the results
case_filtered_gene_burden_results = read_tsv("./outputs/germline_PV_gene_enrichment_min_2_in_case.tsv")

## Get gene names
case_filtered_gene_burden_results_pos_only = case_filtered_gene_burden_results[case_filtered_gene_burden_results$estimate > 0,]
case_filtered_gene_burden_results_pos_only = case_filtered_gene_burden_results_pos_only %>% mutate(
  case_num_bin = case_when(
    carrier_in_case >= 10 ~ ">= 10",
    carrier_in_case == 8 ~ "= 8",
    carrier_in_case == 7 ~ "= 7",
    carrier_in_case == 6 ~ "= 6",
    carrier_in_case == 5 ~ "= 5",
    carrier_in_case == 4 ~ "= 4",
    carrier_in_case == 3 ~ "= 3",
    carrier_in_case == 2 ~ "= 2",
    TRUE ~ "< 5"
  )
)
show_text_subset = case_filtered_gene_burden_results_pos_only[
  (case_filtered_gene_burden_results_pos_only$p.value < 0.01)|(case_filtered_gene_burden_results_pos_only$carrier_in_case >=10),
  ]


case_filtered_gene_burden_results_pos_only

viz_genes = c("POT1","CHEK2","BRCA1","BRCA2","MUTYH","ERCC3","ERCC4","MSH2","TP53")
case_filtered_gene_burden_results_pos_only$highlight = case_filtered_gene_burden_results_pos_only$gene %in% viz_genes

case_enrichment_plot_clean = ggplot(case_filtered_gene_burden_results_pos_only %>% filter(estimate < 5) %>% arrange(highlight), 
       aes(x = exp(estimate), y = neg_log_adjusted_pval, color = highlight)) +
  geom_point() +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "gray") + 
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "black")) +
  pretty_plot() + 
  L_border() + 
  scale_x_log10() +
  geom_text_repel(data=case_filtered_gene_burden_results_pos_only %>% filter(highlight),aes(label=gene),size=2) +
  theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y = element_blank())
case_enrichment_plot_clean

cowplot::ggsave2("./outputs/plots/case_enrichment_plot_clean.pdf",case_enrichment_plot_clean,dpi=300,height = 1.3,width=3.5)

# case_filtered_enrichment_plot = ggplot(case_filtered_gene_burden_results_pos_only,aes(x=estimate,y=neg_log_adjusted_pval,label=gene,color=case_num_bin)) + 
#   geom_hline(yintercept = -log(0.05),linetype="dashed",color="black") + 
#   geom_point(aes(size=carrier_in_case)) +
#   scale_size_continuous(range = c(2, 10)) +
#   xlim(0,6) +
#   theme_minimal() +
#   labs(size="# of PV Carriers in patients",
#        color="# of PV Carriers in patients",
#        y="-log(adjusted p-values)",
#        x="log odds\n(ancestry, sex-adjusted)"
#        ) +
#   geom_label_repel(data=show_text_subset,max.overlaps =20,box.padding = 1, size=5,fill="white")
# 
# case_filtered_enrichment_plot

## Check POT1 effect size, dont forget to exponentiate
case_filtered_gene_burden_results[case_filtered_gene_burden_results$gene=="POT1",]


ggsave(case_filtered_enrichment_plot,filename = "./outputs/plots/01_germline_PV_gene_enrichment_min_2_in_case.png",dpi=300,width=16,height=5)
ggsave(case_filtered_enrichment_plot,filename = "./outputs/plots/01_germline_PV_gene_enrichment_min_2_in_case.pdf",dpi=300,width=16,height=5)

write.table(case_filtered_gene_burden_results,file="./outputs/germline_PV_gene_enrichment_min_2_in_case.tsv",sep="\t",row.names=FALSE)



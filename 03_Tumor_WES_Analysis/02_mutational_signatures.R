library(dplyr)
library(MutationalPatterns)
library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

maf_df = read.csv(file="../data/processed/tumor_WES/ASC_mutations.maf",sep="\t",check.names=FALSE)
#maf_df = read.csv(file="../data/processed/tumor_WES/ASC_mutations_new.maf",sep="\t",check.names=FALSE)

## Calculate tri-nucleotide frequency 
mut_genomic_ranges = GRanges(
  seqnames = paste0("chr",Rle(maf_df$Chromosome)),
  ranges = IRanges(start = maf_df$Start_position, end = maf_df$End_position),
  strand = Rle(rep("*", nrow(maf_df))), # Mutations are typically strand-agnostic
  ref = maf_df$Reference_Allele,
  alt = maf_df$Tumor_Seq_Allele2
)
genome(mut_genomic_ranges) = "hg19"
genomic_range_list = split(mut_genomic_ranges, f = maf_df$Tumor_Sample_Barcode)

# Extract SNVs
snv_grl = get_mut_type(genomic_range_list, type = "snv")
# Obtain trinucleotide matrix
tri_nuc_mat = MutationalPatterns::mut_matrix(snv_grl,ref_genome = BSgenome.Hsapiens.UCSC.hg19)
tri_nuc_mat_to_write = as.data.frame(t(tri_nuc_mat),check.names=FALSE)
tri_nuc_mat_to_write = cbind(Tumor_Sample_Barcode = rownames(tri_nuc_mat_to_write), tri_nuc_mat_to_write)

if(FALSE) {
  write.table(tri_nuc_mat_to_write,file="../data/processed/tumor_WES/ASC_trinucleotide_mutation_counts.tsv",sep="\t",quote = FALSE,row.names = FALSE)
  #write.table(tri_nuc_mat_to_write,file="../data/processed/tumor_WES/ASC_trinucleotide_mutation_counts_new.tsv",sep="\t",quote = FALSE,row.names = FALSE)
}

# Fit to existing signatures
cosmic_signatures = get_known_signatures()
strict_decomp = fit_to_signatures_strict(
  tri_nuc_mat,
  cosmic_signatures,
  max_delta = 0.2,
  method = "backwards"
)

# Calculate the relative contributions and visualize
strict_decomp_res = strict_decomp$fit_res$contribution
strict_decomp_res_rel = as.data.frame(sweep(strict_decomp_res,2,colSums(strict_decomp_res),`/`),check.names=FALSE)
strict_decomp_res_rel = cbind(SBS = rownames(strict_decomp_res_rel), strict_decomp_res_rel)
# strict_decomp_res_rel$SBS = rownames(strict_decomp_res_rel)
strict_decomp_res_rel_long = strict_decomp_res_rel %>%
  pivot_longer(
    cols = -SBS, 
    names_to = "sample", 
    values_to = "value")

# Filter out entries with 0 contribution
strict_decomp_res_rel_long_no_zero = strict_decomp_res_rel_long %>% filter(value>0)
strict_decomp_res_rel_long_no_zero$SBS = as.character(strict_decomp_res_rel_long_no_zero$SBS)
strict_decomp_res_rel_long_no_zero$SBS = factor(strict_decomp_res_rel_long_no_zero$SBS, levels=unique(strict_decomp_res_rel_long_no_zero$SBS))
strict_decomp_res_rel_long_no_zero$size_bin = cut(
  strict_decomp_res_rel_long_no_zero$value,
  breaks = c(0,0.2,0.4,0.6,0.8,1.1),
  labels = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8",'0.8-1')
)

# Do the signature dotplot
sig_dotplot = ggplot(strict_decomp_res_rel_long_no_zero,aes(x=SBS,y=sample,size=value,color=size_bin)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45,vjust=0.5)) +
  theme_minimal() +
  labs(x="Signature",y="Sample",size="Signature Proportion",color="Proportion Bin")
sig_dotplot

ggsave(sig_dotplot,file="./outputs/plots/02_mutationa_signature_dot_plot.png",height=16,width=16)

majority_only = strict_decomp_res_rel_long_no_zero %>% arrange(-value) %>% distinct(sample,.keep_all = TRUE)
table(majority_only$SBS)

sample_meta_with_sbs = merge(sample_meta_df,majority_only,
                             by.x="entity:sample_id",by.y="sample",all.x=TRUE)
keep_sigs = c("SBS1", "SBS7a", "SBS7b", "SBS15", "SBS87")
sample_meta_with_sbs = sample_meta_with_sbs %>%
  mutate(SBS_grouped = fct_other(SBS, keep = keep_sigs, other_level = "Other SBS"))
sample_meta_with_sbs_stats = sample_meta_with_sbs %>% group_by(`Primary Site (Recombined)`) %>% summarise() 

signature_proportion_plot <- ggplot(sample_meta_with_sbs, aes(x = `Primary Site (Recombined)`, y = value, fill = SBS_grouped)) +
  # position = "fill" turns the counts/sums into proportions (0 to 1)
  geom_bar(stat = "identity", position = "fill") +
  # Flip coordinates if you have many primary sites to keep labels readable
  coord_flip() +
  # Use a color scale that can handle many signatures (like viridis or expanded palette)
  scale_fill_viridis_d(option = "mako", name = "Signature") +
  labs(
    title = "Proportion of SBS Signatures by Primary Site",
    x = "Primary Site",
    y = "Proportion of Total Signature Value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

print(signature_proportion_plot)

# Inspect reconstruction quality
reconstruction_plot = plot_original_vs_reconstructed(tri_nuc_mat,strict_decomp$fit_res$reconstructed)

## Proportion of samples by primary site with sig7b proportion > 0.2?
# Load metadata
sample_meta_df = read.csv("../data/processed/sample_clin_data.tsv",check.names = FALSE,sep="\t") %>% filter(
  `entity:sample_id` %in% maf_df$Tumor_Sample_Barcode
)
strict_decomp_res_rel_long_no_zero_sbs7_only = strict_decomp_res_rel_long_no_zero %>% filter(SBS=="SBS1")#filter(SBS=="SBS7b" | SBS=="SBS7a")
#strict_decomp_res_rel_long_no_zero_sbs7_only = strict_decomp_res_rel_long_no_zero %>% filter(SBS=="SBS1")

sample_meta_with_sbs7 = merge(sample_meta_df,strict_decomp_res_rel_long_no_zero_sbs7_only,
                              by.x="entity:sample_id",by.y="sample",all.x=TRUE) %>% replace_na(list(value=0))
sample_meta_with_sbs7$is_over_05 = sample_meta_with_sbs7$value >= 0.5

sample_meta_with_sbs7_stats = sample_meta_with_sbs7 %>% group_by(`Primary Site (Recombined)`) %>% summarise(mean(is_over_05),sum(is_over_05),length(is_over_05))
sample_meta_with_sbs7_stats_long = sample_meta_with_sbs7_stats %>% pivot_longer(`mean(is_over_05)`)

sbs7_prop_plot = ggplot(sample_meta_with_sbs7_stats,aes(y=fct_reorder(`Primary Site (Recombined)`,`mean(is_over_05)`),x=`mean(is_over_05)`)) +
  geom_bar(stat="identity") +
  labs(y="Primary Site",x="% Samples with >= 50% SBS7") +
  theme_minimal() +
  scale_x_continuous(labels = scales::percent)
sbs7_prop_plot
ggsave(sbs7_prop_plot,file="./outputs/plots/02_over_050_SBS7_proportion_by_site_plot.png",dpi=300,height=4,width=4)


if (FALSE) {
  write.table(strict_decomp_res_rel,file="../data/processed/tumor_WES/ASC_mutational_signature_relative.tsv",sep="\t",quote = FALSE,row.names = FALSE)
  write.table(strict_decomp_res_rel,file="../data/processed/tumor_WES/ASC_mutational_signature_relative_new.tsv",sep="\t",quote = FALSE,row.names = FALSE)
}


# Obtain primary sites
sample_meta_df = read.csv("../data/processed/sample_clin_data.tsv",check.names = FALSE,sep="\t")
#primary_sites = sample_meta_df_subset %>% filter(`entity:sample_id` %in% rownames(mut_type_occurrences)) %>% pull(`Primary Site (Recombined)`)
primary_sites = sample_meta_df_subset[sample_meta_df_subset$`entity:sample_id` %in% rownames(mut_type_occurrences),]$`Primary Site (Recombined)`
# Obtain mutation type counts by samples
mut_type_occurrences = mut_type_occurrences(snv_grl, ref_genome = BSgenome.Hsapiens.UCSC.hg19)
# Plot count distributions by primary sites
point_mutation_spectrum_plot = plot_spectrum(
  mut_type_occurrences,
  by=primary_sites,
  CT = TRUE,
  indv_points = TRUE,
  legend = TRUE
) + theme_bw()
point_mutation_spectrum_plot
ggsave(point_mutation_spectrum_plot,file="./outputs/plots/02_tri_nucleotide_spectrum_plot.png",width=12,height = 8)


library(dplyr)
library(MutationalPatterns)
library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

maf_df = read.csv(file="data/processed/tumor_WES/ASC_mutations.maf",sep="\t",check.names=FALSE)

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
  write.table(tri_nuc_mat_to_write,file="data/processed/tumor_WES/mutational_signatures/ASC_trinucleotide_mutation_counts.tsv",sep="\t",quote = FALSE,row.names = FALSE)
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
ggsave(sig_dotplot,file="03_Tumor_WES_Analysis/outputs/plots/02_mutationa_signature_dot_plot.png",height=16,width=16)

# Inspect reconstruction quality
reconstruction_plot = plot_original_vs_reconstructed(tri_nuc_mat,strict_decomp$fit_res$reconstructed)

if (FALSE) {
  write.table(strict_decomp_res_rel,file="data/processed/tumor_WES/mutational_signatures/ASC_mutational_signature_relative.tsv",sep="\t",quote = FALSE,row.names = FALSE)
}


# Obtain primary sites
sample_meta_df = read.csv("data/processed/sample_clin_data.tsv",check.names = FALSE,sep="\t")
primary_sites = sample_meta_df_subset[sample_meta_df_subset$`entity:sample_id` %in% rownames(mut_type_occurrences),]$`PRIMARY SITE (Combined)`
# Obtain mutation type counts by samples
mut_type_occurrences = mut_type_occurrences(snv_grl, ref_genome = BSgenome.Hsapiens.UCSC.hg19)
# Plot count distributions by primary sites
point_mutation_spectrum_plot = plot_spectrum(
  mut_type_occurrences,
  by=primary_sites,
  CT = TRUE,
  indv_points = TRUE,
  legend = TRUE
)
point_mutation_spectrum_plot
ggsave(point_mutation_spectrum_plot,file="03_Tumor_WES_Analysis/outputs/plots/02_tri_nucleotide_spectrum_plot.png",width=24,height = 12)


library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(tableone)
library(forcats)
library(ggpubr)
library(ggalluvial)
library(stringr)
library(scales)
library(RColorBrewer)
library(reshape2)
library(ggsci)
library(ggbeeswarm)
library(colorRamp2)
library(ComplexHeatmap)
library(viridis)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

# Import processed clinical Data
clin_df = read.csv("data/processed/clinical_data_manual_deident.tsv",sep="\t",check.names = FALSE)

# Make a clinical characteristic table for the entire cohort
#"Biopsy Location (Wagner2024)",
all_vars = c("Age (Combined)","SEX (EHR_EXTRACTED)","Primary Site (Recombined)","Primary Site (Wagner2024)","Biopsy Location (Wagner2024)","RAAS_LAAS_Class","OTHER_CANCER (PRD)")
#all_vars = c("Age (Combined)","SEX (EHR_EXTRACTED)","PRIMARY SITE (Combined)","RAAS_LAAS_Class","OTHER_CANCER (PRD)")
categorical_vars = c("Primary Site (Recombined)","Primary Site (Wagner2024)","Biopsy Location (Wagner2024)","SEX (EHR_EXTRACTED)","RAAS_LAAS_Class","OTHER_CANCER (PRD)")
table_overall = CreateTableOne(data = clin_df,vars = all_vars, factorVars = categorical_vars)
table_overall_mat = print(table_overall, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_overall_mat, file = "01_Clinical_Data_Analysis/outputs/tables/clinical_characteristic_table.csv")

table(clin_df$`Primary Site (Wagner2024)`)
sort(table(clin_df$`Primary Site (Recombined)`))

# Make a clinical characteristic table for the entire cohort, stratified by primary site
all_vars = c("Primary Site (Wagner2024)","Biopsy Location (Wagner2024)","Age (Combined)","SEX (EHR_EXTRACTED)","RAAS_LAAS_Class","OTHER_CANCER (PRD)")
categorical_vars = c("Primary Site (Wagner2024)","Biopsy Location (Wagner2024)","SEX (EHR_EXTRACTED)","RAAS_LAAS_Class","OTHER_CANCER (PRD)")
table_with_strata = CreateTableOne(data = clin_df,vars = all_vars, factorVars = categorical_vars,strata = c("Primary Site (Recombined)"))
table_with_strata_mat = print(table_with_strata, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_with_strata_mat, file = "01_Clinical_Data_Analysis/outputs/tables/clinical_characteristic_table_stratified_by_site.csv")

# Barplot of patient reported primary sites
primary_site_counts = clin_df %>% group_by(`Primary Site (Recombined)`) %>% mutate(count_name_occurr = n())
clin_df[clin_df$RAAS_LAAS_Class == "",]$RAAS_LAAS_Class = "Unknown"
clin_df$RAAS_LAAS_Class = factor(clin_df$RAAS_LAAS_Class, levels=c("RAAS","LAAS","RAAS & LAAS","Non-RAAS/LAAS","Unknown"))

# Use the npg palette from ggsci
#show_col(pal_npg("nrc")(10))
RAAS_class_palette = c(
  "RAAS" = "#DC0000FF",
  "LAAS" = "4DBBD5FF",
  "RAAS & LAAS" = "#00A087FF",
  "Non-RAAS/LAAS" = "#3C5488FF",
  "Unknown" = "#B09C85FF"
)

primary_site_count_plot = ggplot(clin_df,aes(x=fct_rev(fct_infreq(`Primary Site (Recombined)`)),fill=RAAS_LAAS_Class)) + 
  geom_bar() +
  #geom_text(stat='count', aes(label=..count..), vjust=-1) + 
  theme_minimal() +
  labs(x="Primary Sites",y="# of Patients", fill="AS Subtype") + 
  #scale_fill_brewer(palette="Accent") +
  coord_flip() +
  theme(axis.text.y = element_text(hjust=0.5)) +
  scale_x_discrete(labels=c("Breast (NOS)"="Breast\n(NOS)",
                            "Breast (Parenchymal)"="Breast\n(Parenchymal)",
                            "Breast (Cutaneous)"="Breast\n(Cutaneous)",
                            "Other Visceral Organs" = "Other\nVisceral Organs"
                            )) +
  scale_fill_manual(values=RAAS_class_palette)


primary_site_count_plot
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_num_patients_by_primary_site.png",primary_site_count_plot,dpi=300, width=6,height=4)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_num_patients_by_primary_site.pdf",primary_site_count_plot,dpi=300, width=6,height=4)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_num_patients_by_primary_site_tall.pdf",primary_site_count_plot,dpi=300, width=6,height=12)

## Check statistics for other cancers
other_cancer_df = read.csv("data/processed/other_cancer_data.tsv",sep="\t",check.names = FALSE)
clin_other_cancer_df_merged = merge(clin_df,other_cancer_df,by="STUDY ID",all.x=TRUE)
clin_other_cancer_df_merged$`OTHER CANCERS REMAPPED`[is.na(clin_other_cancer_df_merged$`OTHER CANCERS REMAPPED`)] = "NOT FOUND"
prior_cancer_only = clin_other_cancer_df_merged %>% filter(
  (`DATE OTHER CANCERS (DAYS FROM DX)` <= 0) | (`OTHER CANCERS REMAPPED` == "NOT FOUND")
)
prior_cancer_only
prior_cancer_only_counts = as.data.frame(table(prior_cancer_only$`OTHER CANCERS REMAPPED`,prior_cancer_only$`Primary Site (Recombined)`),stringsAsFactors = FALSE)
colnames(prior_cancer_only_counts) = c("Prior Cancer", "Primary Site", "n")
# prior_cancer_only_counts = as.data.frame(prior_cancer_only_counts)
prior_cancer_only_counts$`Prior Cancer`[prior_cancer_only_counts$`Prior Cancer`=="NOT FOUND"] = "NO PRIOR CANCER"
# prior_cancer_only_counts = prior_cancer_only_counts %>% mutate(
#   case_when(
#     `Prior Cancer` == "NOT FOUND" ~ "NO PRIOR CANCER"
#   )
# )

## Just for the count plot we adjust it
prior_cancer_only_counts_adjusted = prior_cancer_only_counts
prior_cancer_only_counts_adjusted$`Prior Cancer`[prior_cancer_only_counts_adjusted$`Prior Cancer`=="NO PRIOR CANCER"] = "NOT FOUND"

# Step 1: Reorder Primary Site factor by descending total n
site_order <- prior_cancer_only_counts_adjusted %>%
  group_by(`Primary Site`) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  arrange(desc(total_n)) %>%
  pull(`Primary Site`)

# Step 2: Apply factor level ordering (reversed so most is at top)
prior_cancer_only_counts_adjusted$`Primary Site` <- factor(
  prior_cancer_only_counts_adjusted$`Primary Site`, 
  levels = rev(site_order)
)

prior_cancer_only_counts_adjusted$`Prior Cancer` <- factor(
  prior_cancer_only_counts_adjusted$`Prior Cancer`,
  levels = c("NOT FOUND", setdiff(unique(prior_cancer_only_counts_adjusted$`Prior Cancer`), "NOT FOUND"))
)

# Get unique Prior Cancer categories
prior_cancer_levels <- unique(prior_cancer_only_counts_adjusted$`Prior Cancer`)
prior_count_colors <- pal_npg("nrc")(length(prior_cancer_levels) - 1)
names(prior_count_colors) <- setdiff(prior_cancer_levels, "NOT FOUND")
prior_count_palette <- c(prior_count_colors, "NOT FOUND" = "gray80")

# Step 3: Plot
prior_cancer_only_counts_plot <- ggplot(prior_cancer_only_counts_adjusted, 
                                        aes(x = `Primary Site`, y = n, fill = `Prior Cancer`)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() +
  labs(x = "Primary Sites", y = "# of Patients", fill = "Prior Cancer") +
  theme(axis.text.y = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c(
    "Breast (Cutaneous)" = "Breast\n(Cutaneous)",
    "Breast (Parenchymal)" = "Breast\n(Parenchymal)",
    "Other Visceral Organs" = "Other\nVisceral Organs"
  )) +
  scale_fill_manual(values = prior_count_palette)

# Step 4: Save
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_prior_cancer_count_plot.pdf", 
       prior_cancer_only_counts_plot, dpi = 300, width = 6, height = 4)
  

prior_cancer_only_counts_filtered = prior_cancer_only_counts %>% filter(`Prior Cancer` != "NO PRIOR CANCER")

prior_history_count_plot = ggplot(data=prior_cancer_only_counts_filtered,aes(axis1=`Prior Cancer`,axis2=`Primary Site`,y=n)) +
  geom_alluvium(aes(fill = `Primary Site`)) +
  geom_stratum() +
  geom_text(stat = "stratum",size=5,
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Prior Cancer", "AS Primary Site"),
                   expand = c(0.15, 0.05)) +
  theme_minimal() +
  scale_y_continuous(breaks = pretty_breaks(10)) +
  labs(y="# of Prior Cancers",x="",fill="Primary Sites") +
  scale_fill_npg() +
  theme(text = element_text(size=20))

prior_history_count_plot
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_prior_cancer_to_primary_alluvial_plot.png",prior_history_count_plot,dpi=300,width=20,height=12)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_prior_cancer_to_primary_alluvial_plot_no_no_prior.pdf",prior_history_count_plot,dpi=300,width=20,height=14)

# Total number of patients with other cancers
num_with_other_cancers = sum(clin_df$`OTHER_CANCER (PRD)` == "YES")
num_with_other_cancers_prop = num_with_other_cancers/dim(clin_df)[1]
print(num_with_other_cancers)
print(num_with_other_cancers_prop)

# Among those with prior cancer, what does the distribution look like?
clin_other_cancer_df_merged_with_cancer = clin_other_cancer_df_merged %>% filter(`OTHER CANCERS REMAPPED` != "NOT FOUND")
clin_other_cancer_df_merged_no_dup = clin_other_cancer_df_merged_with_cancer %>% distinct(`STUDY ID`,`OTHER CANCERS REMAPPED`, .keep_all = TRUE)
other_cancer_counts = sort(table(clin_other_cancer_df_merged_no_dup$`OTHER CANCERS REMAPPED`))
other_cancer_counts
other_cancer_counts/num_with_other_cancers

## How many of angiosarcoma were RAAS or non-RAAS associated? Particular for breast angiosarcomas
breast_as_ehr_df = clin_df[clin_df$`Primary Site (Recombined)` %in% c("Breast (Cutaneous)","Breast (Parenchymal)"),]
breast_as_counts = dim(breast_as_ehr_df)[1]
site_by_raas_lass = table(clin_df$RAAS_LAAS_Class,clin_df$`Primary Site (Recombined)`)
site_by_raas_lass

## Do non-cutaneous breast angiosarcoma patients have prior cancers?
clin_df[(clin_df$`Primary Site (Recombined)`=="Breast (Parenchymal)")&clin_df$`OTHER_CANCER (PRD)`!="NOT FOUND",]
no_prior_cancer_breast_as_counts = sum(breast_as_ehr_df[breast_as_ehr_df$RAAS_LAAS_Class == "Non-RAAS/LAAS",]$`OTHER_CANCER (PRD)`!="YES")
no_prior_cancer_breast_as_counts
table(breast_as_ehr_df[breast_as_ehr_df$RAAS_LAAS_Class == "Non-RAAS/LAAS",]$`OTHER_CANCER (PRD)`)

# Plot age at diagnosis by primary site
age_count_plot = ggplot(clin_df,aes(x=fct_reorder(`Primary Site (Recombined)`,`Age (Combined)`),y=`Age (Combined)`)) +
  geom_violin(width=1.2) +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  geom_quasirandom(aes(color=RAAS_LAAS_Class)) +
  labs(x="Primary Site",y="Age at Diagnosis", color="AS Subtype") +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("Breast (Parenchymal)","Breast (Cutaneous)")),na.rm=TRUE) +
  scale_color_manual(values=RAAS_class_palette) +
  scale_x_discrete(labels=c("Breast (NOS)"="Breast\n(NOS)",
                            "Breast (Parenchymal)"="Breast\n(Parenchymal)",
                            "Breast (Cutaneous)"="Breast\n(Cutaneous)",
                            "Other Visceral Organs" = "Other\nVisceral Organs"
  )) +
  scale_color_manual(values=RAAS_class_palette) +
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=0.5,vjust = 0.5),plot.margin = margin(t = 20))
  #scale_fill_npg()
age_count_plot
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_age_at_diagnosis_by_primary_site.png",plot=age_count_plot,dpi=300, width=12,height=5)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_age_at_diagnosis_by_primary_site.pdf",plot=age_count_plot,dpi=300, width=12,height=5)


# Get the actual medians of each group
clin_df %>% group_by(`Primary Site (Recombined)`) %>% summarise(median(`Age (Combined)`,na.rm = TRUE))

# How many patients experienced some form of metastasis?
clin_df_with_met_ehr = clin_df %>% filter(!((`Mets Sites Ever (Recombined)` %in% c("Not Extracted"))| is.na(clin_df$`Mets Sites Ever (Recombined)`)))
clin_df_with_met_ehr = clin_df_with_met_ehr %>%
  mutate(has_met = !(`Mets Sites Ever (Recombined)` %in% c("No Mets")))
num_with_mets = dim(clin_df_with_met_ehr[clin_df_with_met_ehr$has_met == TRUE,])[1]
num_with_mets_prop = num_with_mets/dim(clin_df_with_met_ehr)[1]
print(num_with_mets)
print(dim(clin_df_with_met_ehr)[1])
print(num_with_mets_prop)
## Calculate 95% confidence interval
binom_test_res = binom.test(num_with_mets, dim(clin_df_with_met_ehr)[1])
binom_test_res

# Which primary site had the highest percentage of metastasis?
met_count_data = clin_df_with_met_ehr %>%
  group_by(`Primary Site (Recombined)`, has_met) %>%
  summarise(Count = n(), .groups = 'drop')
met_count_total = met_count_data %>%
  group_by(`Primary Site (Recombined)`) %>%
  summarise(Total = sum(Count))
met_percentage_data = met_count_data %>%
  left_join(met_count_total, by = "Primary Site (Recombined)") %>%
  mutate(Percentage = Count / Total) %>%
  mutate(has_met_Percentage = ifelse(has_met, Count / Total, 0)) %>%
  ungroup()
## Hide Missing or Unkown
#show_col(pal_npg("nrc")(2))
met_percentage_data_filtered = met_percentage_data %>%
  filter(`Primary Site (Recombined)` != "Missing or Unknown") %>%
  mutate(has_met = case_when(
    has_met==FALSE ~ "No Mets",
    has_met==TRUE ~ "Has Met"
  ))
## Calculate binom exact in each category
met_percentage_data_filtered = met_percentage_data_filtered %>% rowwise() %>% mutate(
  has_met_ci_low=binom.test(Count,Total)$conf.int[[1]],
  has_met_ci_high=binom.test(Count,Total)$conf.int[[2]],
)
met_percentage_data_filtered
met_percentage_data_filtered[c("Primary Site (Recombined)","has_met","Count","Total","has_met_ci_low","has_met_ci_high")]

met_class_palette = c("Has Met" = "#E64B35FF","No Mets" = "#4DBBD5FF")
met_percentage_data_plot = ggplot(met_percentage_data_filtered, aes(x = fct_rev(fct_reorder(`Primary Site (Recombined)`,has_met_Percentage)), y = Count, fill = has_met,label=Percentage)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Primary Site", y = "Count", fill = "Has Metastasized") +
  theme_minimal() +
  geom_text(aes(label =scales::percent(Percentage)), 
            position = position_stack(vjust = 0.5), size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +
  scale_fill_manual(values=met_class_palette) +
  labs(x="Primary Site",y="# of Patients")# +
  #geom_errorbar(aes(ymin = has_met_ci_low*100, ymax = has_met_ci_high*100), height = 0.2)
  #scale_fill_npg()
met_percentage_data_plot

ggsave("01_Clinical_Data_Analysis/outputs/plots/00_primary_site_met_percentages.png",met_percentage_data_plot,dpi=300,width=12,height=4)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_primary_site_met_percentages.pdf",met_percentage_data_plot,dpi=300,width=12,height=4)

met_percentage_data_filtered

# Convert to long form to visualize the most common metastatic events
clin_df_long = clin_df_with_met_ehr %>%
  separate_rows(`Mets Sites Ever (Recombined)`, sep = ",") %>%
  mutate(`Mets Sites Ever (Recombined)` = str_trim(`Mets Sites Ever (Recombined)`)) %>%
  mutate(`Mets Sites Ever (Recombined)` = case_when(
    `Mets Sites Ever (Recombined)` == "" ~ "Unknown",
    TRUE ~ `Mets Sites Ever (Recombined)`
  ))

# Filter and count 
met_df = clin_df_long %>%
  filter(!is.na(`Mets Sites Ever (Recombined)`)) %>%
  filter(! `Mets Sites Ever (Recombined)` %in% c("Unknown"))

met_cooccurence = met_df %>% count(`Primary Site (Recombined)`,`Mets Sites Ever (Recombined)`)
met_cooccurence = met_cooccurence %>% mutate(`Primary Site`=`Primary Site (Recombined)`,`Metastatic Site`=`Mets Sites Ever (Recombined)`)

# Save this as a source data
write.csv(met_cooccurence, file = "01_Clinical_Data_Analysis/outputs/tables/metastatic_occurence_count_table.csv",row.names = FALSE)

## What is the most common metastatic site?
## Also calculate the proportion of metastatic site total (after subtracting no mets events)

# Step 1: Summarize the data
met_cooccurence_totals <- met_cooccurence %>%
  group_by(`Metastatic Site`) %>%
  summarise(Total = sum(n), .groups = "drop")

# Step 2: Compute total across all metastatic sites except "No Mets"
all_mets_total <- met_cooccurence_totals %>%
  filter(`Metastatic Site` != "No Mets") %>%
  summarise(AllMetsTotal = sum(Total)) %>%
  pull(AllMetsTotal)

# Step 3: Compute proportions and confidence intervals
met_cooccurence_totals <- met_cooccurence_totals %>%
  rowwise() %>%
  mutate(
    AllMetsTotal = all_mets_total,
    met_site_total_prop = Total / AllMetsTotal,
    ci = list(binom.test(Total, AllMetsTotal)$conf.int),
    met_site_ci_low = ci[[1]],
    met_site_ci_high = ci[[2]]
  ) %>%
  select(-ci) %>%
  ungroup()

met_cooccurence_totals
# 
# met_cooccurence_totals = met_cooccurence %>%
#   group_by(`Metastatic Site`) %>%
#   summarise(Total = sum(n)) %>%
#   mutate(AllMetsTotal = sum(Total)-met_cooccurence_totals[met_cooccurence_totals$`Metastatic Site`=="No Mets",]$Total) %>%
#   rowwise() %>%
#   mutate(
#     met_site_total_prop = Total/AllMetsTotal,
#     met_site_ci_low=binom.test(Total,AllMetsTotal)$conf.int[[1]],
#     met_site_ci_high=binom.test(Total,AllMetsTotal)$conf.int[[2]]
#   )
# met_cooccurence_totals

met_cooccurence_merged = met_cooccurence %>% left_join(met_cooccurence_totals)
common_met_plot = ggplot(met_cooccurence_merged,aes(y=fct_reorder(`Metastatic Site`,Total),x=n,fill=`Primary Site (Recombined)`)) +
  geom_bar(stat="identity") +
  scale_fill_npg() +
  theme_minimal() +
  labs(x="# of events",y="Metastatic Site")
common_met_plot
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_metastatic_site_counts.png",common_met_plot,dpi=300,width=5,height=4)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_metastatic_site_counts.pdf",common_met_plot,dpi=300,width=5,height=4)

met_cooccurence_merged

sum(met_cooccurence$n)
# The primary site to metastatic site alleuvial plot
met_site_plot = ggplot(data = met_cooccurence,
       aes(axis1 = `Primary Site`, axis2 = `Metastatic Site`, y = n)) +
  geom_alluvium(aes(fill = `Primary Site`)) +
  geom_stratum() +
  geom_text(stat = "stratum",size=8,
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Primary Site", "Metastatic Site"),
                   expand = c(0.15, 0.05)) +
  theme_minimal() +
  scale_y_continuous(breaks = pretty_breaks(10)) +
  labs(y="# of Metastatic Events",x="",fill="Primary Sites") +
  scale_fill_npg() +
  theme(text = element_text(size=20))

met_site_plot
met_cooccurence_totals
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_metastatic_event_alluvial_plot.png",met_site_plot,dpi=300,width=20,height=12)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_metastatic_event_alluvial_plot.pdf",met_site_plot,dpi=300,width=20,height=18)

## Metastatic Site co-occurence
make_cooccurence_matrix = function(met_site_list) {
  met_site_names = unique(unlist(strsplit(met_site_list, ", ")))
  met_site_co_occurrence_matrix = matrix(0, nrow = length(met_site_names), ncol = length(met_site_names),
                                         dimnames = list(met_site_names, met_site_names))
  # Function to update co-occurrence counts
  update_co_occurrence = function(site_pair, co_occurrence_matrix) {
    sites <- unlist(strsplit(site_pair, ", "))
    if (length(sites) > 1) {
      for (i in 1:(length(sites) - 1)) {
        for (j in (i + 1):length(sites)) {
          co_occurrence_matrix[sites[i], sites[j]] = co_occurrence_matrix[sites[i], sites[j]] + 1
          co_occurrence_matrix[sites[j], sites[i]] = co_occurrence_matrix[sites[j], sites[i]] + 1
        }
      }
    }
    return(co_occurrence_matrix)
  }
  # Update the co-occurrence matrix based on the list
  for (site in met_site_list) {
    met_site_co_occurrence_matrix = update_co_occurrence(site, met_site_co_occurrence_matrix)
  }
  return(met_site_co_occurrence_matrix)
}

met_site_list = clin_df_with_met_ehr$`Mets Sites Ever (Recombined)`
met_site_co_occurrence_matrix= make_cooccurence_matrix(met_site_list)
get_lower_tri=function(cormat){
  cormat[upper.tri(cormat)] = NA
  return(cormat)
}
met_site_co_occurrence_matrix_lower = get_lower_tri(met_site_co_occurrence_matrix)
met_site_co_occurrence_matrix_lower_melted = melt(met_site_co_occurrence_matrix_lower, na.rm = TRUE)
met_site_co_occurrence_matrix_lower_melted_filtered = met_site_co_occurrence_matrix_lower_melted %>%
  filter(`Var1` != "No Mets") %>%
  filter(`Var2` != "No Mets") %>%
  filter(`Var2` != `Var1`)

comat_breaks = seq(0, max(met_site_co_occurrence_matrix_lower_melted_filtered$value), by = 2)

met_site_co_occurrence_plot = ggplot(data = met_site_co_occurrence_matrix_lower_melted_filtered, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "#E64B35FF", mid = "white", 
                      space = "Lab", breaks=comat_breaks,
                       name="Metastasis Co-occurence") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                    hjust = 1))+
  coord_fixed() +
  theme(
    #axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(1, 0.0),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(x="",y="")
met_site_co_occurrence_plot
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_metastatic_site_cooccurence_matrix.png",met_site_co_occurrence_plot,dpi=300,width=6,height=6)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_metastatic_site_cooccurence_matrix.pdf",met_site_co_occurrence_plot,dpi=300,width=6,height=6)

## Which originating sites had bone and liver co-occurence
clin_df_with_bone_liver_met = subset(
  clin_df_with_met_ehr,
  grepl("Liver", `Mets Sites Ever (Recombined)`) & grepl("Bone", `Mets Sites Ever (Recombined)`)
  )
originating_table = table(clin_df_with_bone_liver_met$`Primary Site (Recombined)`)
originating_table

## Check if different breast AS type has different met tropism towards liver
cutaneous_breast_met_sites =  met_cooccurence %>% filter(`Primary Site` == "Breast (Cutaneous)",`Metastatic Site` != "No Mets")
parenchymal_breast_met_sites =  met_cooccurence %>% filter(`Primary Site` == "Breast (Parenchymal)",`Metastatic Site` != "No Mets")
cutaneous_breast_to_liver = cutaneous_breast_met_sites %>% filter(`Metastatic Site`=="Liver") %>% pull(n)
cutaneous_breast_to_non_liver = cutaneous_breast_met_sites %>% filter(`Metastatic Site`!="Liver") %>% pull(n) %>% sum()
parenchymal_breast_to_liver = parenchymal_breast_met_sites %>% filter(`Metastatic Site`=="Liver") %>% pull(n)
parenchymal_breast_to_non_liver = parenchymal_breast_met_sites %>% filter(`Metastatic Site`!="Liver") %>% pull(n) %>% sum()
breast_compare_contigency_table = matrix(c(parenchymal_breast_to_liver,cutaneous_breast_to_liver,
                                           parenchymal_breast_to_non_liver,cutaneous_breast_to_non_liver),nrow=2)
fisher_test_res = fisher.test(breast_compare_contigency_table)
fisher_test_res

# Visualize the types of treatment retrieved in different settings
treatment_df = read.csv("data/processed/treatment_data.tsv",sep="\t",check.names = FALSE)
# Only look at treatments that occurred after angiosarcoma DX
treatment_df_filtered = treatment_df %>% filter(
  `START DATE (DAYS FROM DX)` >= 0,
  `MODE (FORMATTED)` != "OTHER CANCER",
  `MODE (FORMATTED)` != "UNKNOWN",
  DRUG != "NOT FOUND IN RECORD"
)
# Merge with other clinical information
treatment_df_filtered_merged = merge(treatment_df_filtered,clin_df,by="STUDY ID",all.x=TRUE)

## Calculate the total number of treatment records per patient
treatment_records_per_patient = treatment_df_filtered_merged %>% group_by(`STUDY ID`) %>% summarise(n_records = n(), .groups = "drop")
## Calculate average
average_num_records = mean(treatment_records_per_patient$n_records)
treatment_records_per_patient_plot = ggplot(treatment_records_per_patient,aes(x=n_records)) + 
  geom_histogram(binwidth=1,fill = "gray", color = "black", size = 0.3) + 
  theme_minimal() +
  geom_vline(xintercept = average_num_records,linetype="dashed") +
  labs(x="Number of Treatment Records per Patient",y="Number of Patients")

treatment_records_per_patient_plot
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_treatment_records_per_patient_histo.png",treatment_records_per_patient_plot,dpi=300,width=6,height=4)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_treatment_records_per_patient_histo.pdf",treatment_records_per_patient_plot,dpi=300,width=6,height=4)

treatment_records_per_patient

## Calculate the average number of treatment records per patient per primary site
treatment_records_per_patient_merged = merge(treatment_records_per_patient,clin_df[,c("STUDY ID","Primary Site (Recombined)")],by="STUDY ID")
treatment_records_per_patient_per_primary_site_plot = ggplot(treatment_records_per_patient_merged,aes(y=`Primary Site (Recombined)`,x=n_records)) + 
  geom_boxplot() + 
  geom_jitter() +
  theme_minimal() +
  labs(y="Primary Site",x="Number of Treatment Records per Patient")
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_treatment_records_per_patient_per_primary_site_boxplot.png",treatment_records_per_patient_per_primary_site_plot,dpi=300,width=4,height=8)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_treatment_records_per_patient_per_primary_site_boxplot.pdf",treatment_records_per_patient_per_primary_site_plot,dpi=300,width=4,height=8)

# Treatment by Primary Sites
met_treatment_only = treatment_df_filtered_merged[treatment_df_filtered_merged$MODE=="METS",]
site_by_drug_table = table(treatment_df_filtered_merged$DRUG,treatment_df_filtered_merged$`Primary Site (Recombined)`)
treatment_col_fun = colorRamp2(c(0, 20), c("gray", "orange"))
treatment_col_sum = colSums(site_by_drug_table)
treatment_row_sum = rowSums(site_by_drug_table)
site_by_drug_table_ordered = site_by_drug_table[order(treatment_row_sum, decreasing = TRUE), order(treatment_col_sum, decreasing = TRUE)]
treatment_col_anno = columnAnnotation(`Records` = anno_barplot(colSums(site_by_drug_table_ordered)))
treatment_row_anno = rowAnnotation(`Records` = anno_barplot(rowSums(site_by_drug_table_ordered)))
# treatment_col_fun = colorRamp2(
#   c(min(site_by_drug_table_ordered), max(site_by_drug_table_ordered)), 
#   viridis(2)
# )

treatment_complex_heatmap = Heatmap(site_by_drug_table_ordered, 
        name = "# of Treatment Records",
        col=treatment_col_fun,
        top_annotation = treatment_col_anno,cluster_rows = FALSE,cluster_columns = FALSE,
        right_annotation = treatment_row_anno,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", site_by_drug_table_ordered[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        column_names_rot = 45
        )
pdf("01_Clinical_Data_Analysis/outputs/plots/treatment_by_primary_site.pdf",width=10,height=14)
draw(treatment_complex_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

treatment_complex_heatmap
# Count occurences by treatment setting
drug_by_mode_count = treatment_df_filtered_merged %>% group_by(DRUG,`MODE (FORMATTED)`) %>% summarize(ModeTotal = n()) %>% ungroup()
drug_total_count = treatment_df_filtered_merged %>% group_by(DRUG) %>% summarize(OverallTotal = n()) %>% ungroup()
drug_total_count$AllTreatmentTotal = drug_total_count$OverallTotal %>% sum()
drug_total_count = drug_total_count %>% rowwise() %>% mutate(
  drug_percentage=OverallTotal/AllTreatmentTotal,
  drug_prop_ci_low=binom.test(OverallTotal,AllTreatmentTotal)$conf.int[[1]],
  drug_prop_ci_high=binom.test(OverallTotal,AllTreatmentTotal)$conf.int[[2]]
) %>% arrange(-OverallTotal)
drug_total_count
  

drug_count_df = drug_total_count %>% left_join(drug_by_mode_count, by = "DRUG") %>% arrange(desc(OverallTotal))
top_drugs = drug_total_count %>% slice_max(order_by=OverallTotal,n=15)
drug_count_df_top_only = drug_count_df %>% filter(DRUG %in% top_drugs$DRUG)

# Make the plot
treatment_plot = ggplot(drug_count_df_top_only,aes(x=ModeTotal,y=reorder(DRUG,OverallTotal),fill=`MODE (FORMATTED)`)) +
  geom_bar(stat="identity") +
  labs(fill="Treatment Setting",x="# of Treatment Records",y="Treatment Name") +
  theme_minimal() +
  scale_fill_npg()
treatment_plot

drug_count_df_top_only
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_top_treatment_received_barplot.png",treatment_plot,dpi=300,width=6,height=4)
ggsave("01_Clinical_Data_Analysis/outputs/plots/00_top_treatment_received_barplot.pdf",treatment_plot,dpi=300,width=6,height=4)

## Which drug was the most common prescribed in mets setting?
drug_count_df_met_only = drug_count_df %>%
  filter(`MODE (FORMATTED)`=="METS") %>%
  mutate(
    met_total = sum(ModeTotal),
    met_percentage = ModeTotal/met_total
    ) %>%
  rowwise() %>%
  mutate(
    met_percentage=ModeTotal/met_total,
    met_prop_ci_low=binom.test(ModeTotal,met_total)$conf.int[[1]],
    met_prop_ci_high=binom.test(ModeTotal,met_total)$conf.int[[2]]
  ) %>%
  arrange(-ModeTotal) 
drug_count_df_met_only

## How much does the top 5 drug category account for all treatment prescribed in the met setting
top_5_mode_total = sum(drug_count_df_met_only[order(-drug_count_df_met_only$ModeTotal),][1:5,]$ModeTotal)
met_drug_total = sum(drug_count_df_met_only$ModeTotal)
top_5_prop_test = binom.test(top_5_mode_total,met_drug_total)
top_5_mode_total/met_drug_total
top_5_prop_test

drug_count_df_met_only %>% filter(ModeTotal >= 10) %>%select(-c("MODE (FORMATTED)"))
drug_count_df_met_only


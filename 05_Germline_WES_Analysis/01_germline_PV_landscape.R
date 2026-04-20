library(maftools)

germline_maf_file_path = "./outputs/tables/filtered_germline_variants_cases_subset.maf"
germline_clin_file_path = "./outputs/tables/germline_sample_clin_df.tsv"

germline_case_maf = read.maf(germline_maf_file_path,clinicalData = germline_clin_file_path)

## Full table annotation
germline_case_table = read.csv("./outputs/tables/filtered_germline_variants_cases.tsv",sep="\t")
subset_genes = unique(germline_case_table[germline_case_table$pass_germline_filter!="unselected_clinvar_pathogenic_high_conf",]$SYMBOL)

germline_case_maf_subset = subsetMaf(germline_case_maf,genes = subset_genes)

## Rename column for visalization
germline_case_maf_subset@clinical.data$PrimarySite = as.character(
  germline_case_maf_subset@clinical.data$`Primary_Site_(Recombined)`
)
annotation_color_simple = list(
  PrimarySite = clean_palette
)

annotation_color = list(`Primary_Site_(Recombined)`=primary_site_palette_maftools)
landscape_save_path = "./outputs/plots/02_germline_PV_landscape.pdf"
pdf(file=landscape_save_path,height=4,width=6)
oncoplot(
  maf = germline_case_maf_subset,
  top = 100,
  topBarData = "Age",
  clinicalFeatures = "PrimarySite",  # Use the new simple name
  annotationColor = annotation_color_simple,
  color = variant_col,
  draw_titv = FALSE,
  bgCol = "#FFFFFF",
  borderCol = "#000000",
  cohortSize = 229,
  fontSize=0.7,
  titleFontSize=1,
  legendFontSize = 0.8,
  annotationFontSize=0.8
)
dev.off()

## Binomial exact test of # of patients carrying germline variants over 
binom.test(39,229)
germline_case_maf_subset@data$Tumor_Sample_Barcode

## Plot a dual lollipop focused on POT1
somatic_maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
somatic_maf_metadata_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"
somatic_maf = read.maf(maf=somatic_maf_path,clinicalData=somatic_maf_metadata_path)





## Overriding Maftool internal function to reannotate POT1 domain
## Domain annotation from:
## Impact of Gln94Glu mutation on the structure and function of protection of telomere 1, a cause of cutaneous familial melanoma

get_lp_data_custom = function(maf, geneID = NULL, AACol = NULL, refSeqID = NULL, proteinID = NULL, defaultYaxis = FALSE, verbose = TRUE){
  
  if(is.null(geneID)){
    stop('Please provide a gene name.')
  }
  
  # --- CUSTOM INTERCEPTION LOGIC FOR POT1 ---
  if(!is.null(geneID) && geneID == "POT1") {
    message("!!! Intercepting POT1: Using custom structural domains")
    
    # Define the domains manually (OB-folds and TPP1 binding)
    # These coordinates follow the canonical NP_056265.2
    custom_domains = data.table(
      HGNC = "POT1",
      refseq.ID = "NM_015450",
      protein.ID = "NP_056265",
      aa.length = 634,
      Label = c("OB1", "OB2", "OB3", "HJR", "OB3"),
      Start = c(5, 149, 320,395,540), 
      End = c(140, 299, 393,538,634)
    )
    
    # We "trick" the function by replacing the internal GFF read logic
    # Instead of reading from the system file, we use our custom_domains
    gff = custom_domains
  } else {
    #Protein domain source.
    gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
    gff = readRDS(file = gff)
    data.table::setDT(x = gff)
  }
  mut = subsetMaf(maf = maf, includeSyn = FALSE, genes = geneID, query = "Variant_Type != 'CNV'", mafObj = FALSE)
  
  if(is.null(AACol)){
    pchange = c('HGVSp_Short', 'Protein_Change', 'AAChange')
    if(pchange[pchange %in% colnames(mut)] > 0){
      pchange = suppressWarnings(pchange[pchange %in% colnames(mut)][1])
      if(verbose){
        message(paste0("Assuming protein change information are stored under column ", pchange,". Use argument AACol to override if necessary."))
      }
      colnames(mut)[which(colnames(mut) == pchange)] = 'AAChange'
    }else{
      message('Available fields:')
      print(colnames(mut))
      stop('AAChange field not found in MAF. Use argument AACol to manually specifiy field name containing protein changes.')
    }
  }else{
    colnames(mut)[which(colnames(mut) == AACol)] = 'AAChange'
  }
  
  prot.dat = mut[Hugo_Symbol %in% geneID, .(Variant_Type, Variant_Classification, AAChange)]
  if(nrow(prot.dat) == 0){
    return(NULL)
    #stop(paste(geneID, 'does not seem to have any mutations!', sep=' '))
  }
  
  prot = gff[HGNC %in% geneID]
  
  if(nrow(prot) == 0){
    stop(paste('Structure for protein', geneID, 'not found.', sep=' '))
  }
  
  if(!is.null(refSeqID)){
    prot = gff[refseq.ID == refSeqID]
  } else if(!is.null(proteinID)){
    prot = gff[protein.ID == proteinID]
  } else{
    txs = unique(prot$refseq.ID)
    if(length(txs) > 1){
      if(verbose){
        message(paste(length(txs), ' transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.', sep = ''))
        print(prot[!duplicated(protein.ID),.(HGNC, refseq.ID, protein.ID, aa.length)])
      }
      prot = prot[which(prot$aa.length == max(prot$aa.length)),]
      if(length(unique(prot$refseq.ID)) > 1){
        prot = prot[which(prot$refseq.ID == unique(prot[,refseq.ID])[1]),]
        if(verbose){
          message(paste('Using longer transcript', unique(prot[,refseq.ID])[1], 'for now.', sep=' '))
        }
      } else{
        if(verbose){
          message(paste('Using longer transcript', unique(prot[,refseq.ID])[1], 'for now.', sep=' '))
        }
      }
    }
  }
  
  #Legth of protein
  len = as.numeric(max(prot$aa.length, na.rm = TRUE))
  #Remove NA's
  #prot = prot[!is.na(Label)]
  prot = prot[,domain_lenght := End - Start][order(domain_lenght, decreasing = TRUE)][,domain_lenght := NULL]
  
  #hard coded colors for variant classification if user doesnt provide any
  
  sampleSize = as.numeric(maf@summary[ID %in% 'Samples', summary])
  mutRate = round(getGeneSummary(x = maf)[Hugo_Symbol %in% geneID, MutatedSamples]/sampleSize*100, digits = 2)
  cbioSubTitle = paste0(mutRate, "%")
  
  #prot.dat = prot.dat[Variant_Classification != 'Splice_Site']
  #Remove 'p.'
  prot.spl = strsplit(x = as.character(prot.dat$AAChange), split = '.', fixed = TRUE)
  prot.conv = sapply(sapply(prot.spl, function(x) x[length(x)]), '[', 1)
  
  prot.dat[,conv := prot.conv]
  #If conversions are in HGVSp_long (default HGVSp) format, we will remove strings Ter followed by anything (e.g; p.Asn1986GlnfsTer13)
  pos = gsub(pattern = 'Ter.*', replacement = '',x = prot.dat$conv)
  
  #Following parsing takes care of most of HGVSp_short and HGVSp_long format
  pos = gsub(pattern = '[[:alpha:]]', replacement = '', x = pos)
  pos = gsub(pattern = '\\*$', replacement = '', x = pos) #Remove * if nonsense mutation ends with *
  pos = gsub(pattern = '^\\*', replacement = '', x = pos) #Remove * if nonsense mutation starts with *
  pos = gsub(pattern = '\\*.*', replacement = '', x = pos) #Remove * followed by position e.g, p.C229Lfs*18
  
  
  #pos = as.numeric(sapply(strsplit(x = pos, split = '_', fixed = TRUE), '[[', 1))
  pos = as.numeric(sapply(X = strsplit(x = pos, split = '_', fixed = TRUE), FUN = function(x) x[1]))
  prot.dat[,pos := abs(pos)]
  
  if(nrow( prot.dat[is.na(pos)]) > 0){
    if(verbose){
      message(paste('Removed', nrow( prot.dat[is.na(prot.dat$pos),]), 'mutations for which AA position was not available', sep = ' '))
    }
    
    #print(prot.dat[is.na(pos)])
    prot.dat = prot.dat[!is.na(pos)]
  }
  
  prot.snp.sumamry = prot.dat[,.N, .(Variant_Classification, conv, pos)]
  colnames(prot.snp.sumamry)[ncol(prot.snp.sumamry)] = 'count'
  maxCount = max(prot.snp.sumamry$count, na.rm = TRUE)
  
  prot.snp.sumamry = prot.snp.sumamry[order(prot.snp.sumamry$pos),]
  #prot.snp.sumamry$distance = c(0,diff(prot.snp.sumamry$pos))
  
  if(maxCount <= 5){
    prot.snp.sumamry$count2 = 1+prot.snp.sumamry$count
    lim.pos = 2:6
    lim.lab = 1:5
  }else{
    prot.snp.sumamry$count2 = 1+(prot.snp.sumamry$count * (5/max(prot.snp.sumamry$count)))
    lim.pos = prot.snp.sumamry[!duplicated(count2), count2]
    lim.lab = prot.snp.sumamry[!duplicated(count2), count]
  }
  
  if(length(lim.pos) > 6){
    lim.dat = data.table::data.table(pos = lim.pos, lab = lim.lab)
    lim.dat[,posRounded := round(pos)]
    lim.dat = lim.dat[!duplicated(posRounded)]
    lim.pos = lim.dat[,pos]
    lim.lab = lim.dat[,lab]
  }
  
  if(!defaultYaxis){
    lim.pos = c(min(lim.pos), max(lim.pos))
    lim.lab = c(min(lim.lab), max(lim.lab))
  }
  
  clusterSize = 10 #Change this later as an argument to user.
  repel = FALSE
  if(repel){
    prot.snp.sumamry = repelPoints(dat = prot.snp.sumamry, protLen = len, clustSize = clusterSize)
  }else{
    prot.snp.sumamry$pos2 = prot.snp.sumamry$pos
  }
  
  xlimPos = pretty(0:max(prot$aa.length))
  xlimPos[length(xlimPos)] = max(prot$aa.length)
  
  if(xlimPos[length(xlimPos)] - xlimPos[length(xlimPos)-1] <= 10){
    xlimPos = xlimPos[-(length(xlimPos)-1)]
  }
  
  #Get domain overlap labels for each pos (Issue: #794)
  data.table::setkey(x = prot, "Start", "End")
  if(nrow(prot[!is.na(Start)]) == 0){
    domain_olaps = NA
  }else{
    domain_olaps = data.table::foverlaps(
      x = prot.snp.sumamry,
      y = prot[, .(Start, End, Label)],
      by.x = c("pos", "pos2"),
      by.y = c("Start", "End"),
      mult = "all"
    )[, .(Start, End, Label, Variant_Classification, conv, pos, count)]
    colnames(domain_olaps) = c("Domain_Start", "Domain_End", "Domain_Label", "Variant_Classification", "Conversion", "Position", "n_variants")
  }
  
  return(list(prot.snp.sumamry, xlimPos, lim.pos, lim.lab, cbioSubTitle, prot, domain_olaps))
}

unlockBinding("get_lp_data", asNamespace("maftools"))
assignInNamespace("get_lp_data", get_lp_data_custom, ns = "maftools")
lockBinding("get_lp_data", asNamespace("maftools"))

#pot1_dual_lollipop_save_path = "./outputs/plots/02_POT1_dual_lollipop.pdf"
pot1_dual_lollipop_save_path = "./outputs/plots/02_POT1_dual_lollipop_small.pdf"

pdf(file=pot1_dual_lollipop_save_path,height=6,width=6)
lollipopPlot2(
  m1 = somatic_maf,
  m2 = germline_case_maf_subset, gene = "POT1",
  #AACol1 = "amino_acid_change",
  #AACol2 = "amino_acid_change",
  m1_name = "Somatic",
  m2_name = "Germline",
  m1_label="all",
  m2_label="all",
  refSeqID = "NM_015450"
)
dev.off()

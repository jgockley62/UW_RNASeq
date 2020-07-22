---
title: "UW-RNASeq"
output: html_document
---



1. [Write your own config file](https://kelshmo.github.io/sageseqr/articles/customize-config.html).

2. Specify the active configuration by setting `R_CONFIG_ACTIVE`.


```r
Sys.setenv(R_CONFIG_ACTIVE = "UW")
```

3. Load the `sageseqr` library and login to [Synapse](https://www.synapse.org/). `rnaseq_plan()` calls the arguments from the config file and creates the `drake` plan. Execute `drake::make(plan)` to compute. Run this code from your project root.


```r
#library(devtools)
#devtools::install_git('https://github.com/kelshmo/sageseqr.git')
#devtools::install_github('th1vairam/CovariateAnalysis@dev')
library(sageseqr)
library(edgeR)
library(ggplot2)
library(CovariateAnalysis) #get the package from devtools::install_github('th1vairam/CovariateAnalysis@dev')
library(data.table)
library(plyr)
library(tidyverse)
library(psych)
library(limma)
library(edgeR)
library(biomaRt)
library(RColorBrewer)
library(cqn)
library(glmnet)
library(knitr)
library(doParallel)
library(foreach)
library(githubr)
#BiocManager::install("WGCNA")

# Login to Synapse. Make a Synapse account and use synapser to login: https://r-docs.synapse.org/articles/manageSynapseCredentials.html
synapser::synLogin()

Sys.setenv(R_CONFIG_ACTIVE = "UW")
# TO DO: how to setup cache for end user
# make_cache

# Run the analysis
plan <- sageseqr::rnaseq_plan(metadata_id = config::get("metadata")$synID,
                    metadata_version = config::get("metadata")$version,
                    counts_id = config::get("counts")$synID,
                    counts_version = config::get("counts")$version,
                    gene_id_input = config::get("counts")$`gene id`,
                    sample_id_input='individualID',
                    factor_input = config::get("factors"),
                    continuous_input = config::get("continuous"),
                    gene_id = config::get("biomart")$`gene id`,
                    biomart_id = config::get("biomart")$synID,
                    biomart_version = config::get("biomart")$version,
                    filters = config::get("biomart")$filters,
                    host = config::get("biomart")$host,
                    organism = config::get("biomart")$organism, 
                    conditions = config::get("conditions")
)
import_metadata <- sageseqr::get_data("syn22254834", 1L)
import_metadata$sample <- paste0( 'S', import_metadata$sample)
import_metadata$sampleID <- paste0( 'S', import_metadata$sampleID)
import_metadata$individualID <- paste0( 'I', import_metadata$individualID)

import_counts <- sageseqr::get_data("syn22249653", 1L)

counts <- tibble::column_to_rownames(import_counts, var = 'V1')
colnames( counts ) <- counts[ 'feature', ]
counts <- counts[ row.names(counts) != 'feature', ]
colnames(counts) <- paste0( 'S', colnames(counts) )

import_metadata <- import_metadata[,c("sampleID", "individualID", "disease_cohort", 'sex', "ind_mut_status", "apoeGenotype", "ind_mutation", "ADAD_fam_mut", "CDR", "age", "RIN", "AlignmentSummaryMetrics__PF_READS_ALIGNED", "AlignmentSummaryMetrics__TOTAL_READS",	"RnaSeqMetrics__INTERGENIC_BASES", "RnaSeqMetrics__INTRONIC_BASES", "RnaSeqMetrics__PCT_CODING_BASES", "RnaSeqMetrics__PCT_INTERGENIC_BASES", "RnaSeqMetrics__PCT_INTRONIC_BASES", "RnaSeqMetrics__PCT_RIBOSOMAL_BASES")]



clean_md <- sageseqr::clean_covariates(md = import_metadata, factors = c("sampleID", "individualID", "disease_cohort", 'sex', "ind_mut_status", "apoeGenotype", "ind_mutation", "ADAD_fam_mut", "CDR" ), sample_identifier= "sampleID", continuous = c("age", "RIN",
  "AlignmentSummaryMetrics__PF_READS_ALIGNED",
  "AlignmentSummaryMetrics__TOTAL_READS",	"RnaSeqMetrics__INTERGENIC_BASES",
  "RnaSeqMetrics__INTRONIC_BASES", "RnaSeqMetrics__PCT_CODING_BASES",
  "RnaSeqMetrics__PCT_INTERGENIC_BASES", "RnaSeqMetrics__PCT_INTRONIC_BASES",
  "RnaSeqMetrics__PCT_RIBOSOMAL_BASES"
))

filtered_counts <- sageseqr::filter_genes( clean_metadata = clean_md, count_df = counts, conditions = c("ind_mut_status", "ADAD_fam_mut", "ind_mutation", "sex"), cpm_threshold=1, conditions_threshold=.5 )

biomart_results <- get_biomart(count_df = counts, gene_id = "ensembl_gene_id", synid = "syn22254975", 
    version = 2L, filters = "ensembl_gene_id", host = "ensembl.org", 
    organism = "hsa")
colnames(biomart_results)[1] <- "ensembl_gene_id"

filtered_counts<- filtered_counts[ row.names(filtered_counts) %in% row.names(biomart_results), ]
table(row.names(filtered_counts) %in% row.names(biomart_results))

cqn_counts <- sageseqr::cqn(filtered_counts, biomart_results)
```

### Sex Chromosome-specific Gene Expression Patterns


```r
COUNT <- as.data.frame( counts )
COUNT$ensembl_gene_id <- do.call(rbind, strsplit(row.names(COUNT), '\\.'))[,1]
COUNT <- COUNT[ grepl('ENSG', row.names(COUNT)), ]
COUNT <- COUNT[ COUNT$ensembl_gene_id %in% biomart_results$ensembl_gene_id, ]
  
METADATA <- as.data.frame( import_metadata )
row.names(METADATA) <- METADATA$sampleID

REPORTED.GENDER.COUNTS = biomart_results %>% 
  left_join(COUNT) %>%
  dplyr::select(-one_of("percentage_gene_gc_content")) %>%
  filter(chromosome_name == "X" |chromosome_name == "Y") %>% 
  tidyr::gather(key = item, value = value, -c( ensembl_gene_id, hgnc_symbol, gene_biotype, chromosome_name, gene_length, ensembl_gene_id)) %>%
  dplyr::mutate(value = log(value)) %>%
  dplyr::rename(`counts(log)`= value) %>% 
  dplyr::rename(sampleID = item) %>%
  left_join(METADATA[,c("sampleID", "sex")]) %>% 
  dplyr::rename(`Reported Gender` = sex) 

my.theme <- theme_bw() %+replace% theme(legend.position = 'top', axis.text.x = element_text(angle = 90, hjust = 1), plot.title=element_text(hjust=0.5))
p = list()
p[[1]] = ggplot(filter(REPORTED.GENDER.COUNTS, chromosome_name == "X"), aes(x = `Reported Gender`, y = `counts(log)`)) + geom_boxplot()
p[[1]] = p[[1]] + ggtitle('X') + my.theme
p[[2]] = ggplot(filter(REPORTED.GENDER.COUNTS, chromosome_name == "Y"), aes(x = `Reported Gender`, y = `counts(log)`)) + geom_boxplot()
p[[2]] = p[[2]] + ggtitle('Y') + my.theme
multiplot(plotlist = p, cols = 2)
##XIST and UTY expression 
#ENSG00000229807.11 and ENSG00000183878.15 

#Plot initial data
FILT <- REPORTED.GENDER.COUNTS[ , c('ensembl_gene_id', 'chromosome_name', 'sampleID','counts(log)', 'Reported Gender')] %>% 
  filter( ensembl_gene_id == "ENSG00000229807" | ensembl_gene_id == "ENSG00000183878") %>% 
  dplyr::select(-one_of("chromosome_name")) %>% 
  tidyr::spread(key = ensembl_gene_id, value = `counts(log)`) %>% 
  mutate(XIST = as.numeric(`ENSG00000229807`)) %>% 
  mutate(UTY = as.numeric(`ENSG00000183878`)) %>% 
  mutate(UTY = ifelse(UTY == -Inf, 0, UTY)) %>% 
  mutate(XIST = ifelse(XIST == -Inf, 0, XIST))
p = ggplot(FILT, aes (x= XIST, y = UTY)) 
p = p + geom_point(aes(color=`Reported Gender`)) + 
  ggtitle("Sex Check Inital Sex: UW Cohort") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(colour = "Reported Gender")
p

###Find Sex of missing values
if( T %in%  is.na(REPORTED.GENDER.COUNTS$`Reported Gender`) ){
  FILT <- REPORTED.GENDER.COUNTS[ is.na(REPORTED.GENDER.COUNTS$`Reported Gender`), c('ensembl_gene_id', 'chromosome_name', 'sampleID','counts(log)', 'Reported Gender')] %>% 
    filter( ensembl_gene_id == "ENSG00000229807" | ensembl_gene_id == "ENSG00000183878") %>% 
    dplyr::select(-one_of("chromosome_name")) %>% 
    tidyr::spread(key = ensembl_gene_id, value = `counts(log)`) %>% 
    mutate(XIST = as.numeric(`ENSG00000229807`)) %>% 
    mutate(UTY = as.numeric(`ENSG00000183878`)) %>% 
    mutate(UTY = ifelse(UTY == -Inf, 0, UTY)) %>% 
    mutate(XIST = ifelse(XIST == -Inf, 0, XIST))
  p = ggplot(FILT, aes (x= XIST, y = UTY)) 
  p = p + geom_point() + 
    ggtitle("Sex Check Unknown Sex: UW Cohort") + 
    theme(plot.title = element_text(hjust = 0.5, size = 15)) +
    labs(colour = "Reported Gender")
  p
  
  Males <- FILT[ FILT$UTY > 5,]$sampleID
  Females <- FILT[ FILT$UTY < 3 ,]$sampleID
  
  import_metadata$Infered.Sex <- import_metadata$sex 
  if( length(Males) > 0 ){
    writeLines( paste0( 'Inferred Male Induviduals: ', paste0(Males,collapse =', ') ))
    import_metadata[ import_metadata$sampleID %in% Males, ]$Infered.Sex <- "male"
  }else{ writeLines( 'Inferred Male Induviduals: NA') }
  if( length(Females) > 0 ){
    writeLines( paste0( 'Inferred Female Induviduals: ', paste0(Females,collapse =', ') ))
    import_metadata[ import_metadata$sampleID %in% Females, ]$Infered.Sex <- "female"
  }else{ writeLines( 'Inferred Female Induviduals: NA') }
  
  clean_md$Infered.Sex <- import_metadata$Infered.Sex
}else{
  import_metadata$Infered.Sex <- import_metadata$sex
  clean_md$Infered.Sex <- import_metadata$Infered.Sex
  
}
clean_md$Infered.Sex <- as.factor(clean_md$Infered.Sex)
```

### Sample clustering
PCA based clustering of samples

```r
# Find principal components of expression to plot
cqn_counts$E.no.na <- cqn_counts$E
METADATA <- clean_md[ , colnames(clean_md)[colnames(clean_md) != 'sex'] ]
fill<-colnames(METADATA)
METADATA$sampleID<- row.names(METADATA)
METADATA <- METADATA[, c('sampleID',fill)]
colnames(METADATA)[ colnames(METADATA) == 'Infered.Sex'] <- 'sex'
PC <- prcomp(cqn_counts$E.no.na, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(sampleID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])
plotdata <- left_join(plotdata, rownameToFirstColumn(METADATA, 'sampleID'))
p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(color=disease_cohort, shape=apoeGenotype, size=RIN))
p <- p + theme_bw() + theme(legend.position="right")
p
```

Tree based clustering of samples

```r
# Eucledian tree based analysis
COVARIATES.tmp = data.matrix(METADATA[,c( 'individualID', "sex", "apoeGenotype", "disease_cohort", 'ADAD_fam_mut', "ind_mut_status", "ind_mutation")])
COVARIATES.tmp[is.na(COVARIATES.tmp)] = 0
tree = hclust(as.dist(t(cqn_counts$E.no.na)))
cols = WGCNA::labels2colors(COVARIATES.tmp);
WGCNA::plotDendroAndColors(tree, 
                           colors = cols, 
                           dendroLabels = FALSE, 
                           abHeight = 0.80, 
                           main = "Sample dendrogram",
                           groupLabels = colnames(COVARIATES.tmp))
```

### Distribution of samples (log cpm)

```r
# Plot abberent distribution of logcpm counts
tmp1 = cqn_counts$E %>%
  rownameToFirstColumn('Gene.ID') %>%
  tidyr::gather(sampleID, logCPM, -Gene.ID) %>%
  left_join(METADATA %>%
              rownameToFirstColumn('sampleID'))
p = ggplot(tmp1, aes(x = logCPM, color = sampleID)) + geom_density() 
p = p + theme(legend.position = 'NONE') + facet_grid(.~disease_cohort, scale = 'free')
p
```

Coexpression of genes 

```r
cr = cor(t(cqn_counts$E.no.na))
hist(cr, main = 'Distribution of correlation between genes', xlab = 'Correlation')
```


### Significant Covariates
Correlation between pca of unadjusted mRNA expression and covariates are used to find significant covariates

```r
# Find correlation between PC's of gene expression with covariates
cqn_counts$E = cqn_counts$E[,row.names(clean_md)]
cqn_counts$E.no.na = cqn_counts$E[,row.names(clean_md)]
cqn_counts$E.no.na[is.na(cqn_counts$E.no.na)] = 0

#clean_md
METADATA <- as.data.frame( clean_md)
row.names(METADATA) <- import_metadata$sampleID 
METADATA <- METADATA[, colnames(METADATA) != 'individualID' ]

Meta_HeatMap <- METADATA[, (colnames(METADATA) %in% c('individualID','sampleID')==F)  ]
Meta_HeatMap$sex <- Meta_HeatMap$Infered.Sex
Meta_HeatMap[ , colnames(Meta_HeatMap)[ colnames(Meta_HeatMap) != 'Infered.Sex']]
Iters <- colnames(Meta_HeatMap)[colSums(is.na(Meta_HeatMap)) > 0]
dim(Meta_HeatMap)
dim(Meta_HeatMap[ complete.cases(Meta_HeatMap),])
writeLines("Total Metadata with Missing variables:")
preAdjustedSigCovars = runPCAandPlotCorrelations(cqn_counts$E.no.na, 
                                                 Meta_HeatMap,
                                                 'NULL design(voom-normalized)', 
                                                 isKeyPlot=TRUE, 
                                                 MIN_PVE_PCT_PC = 1)

preAdjustedSigCovars[["PC_res"]][[2]]$plotData

writeLines("Iterate across covariates with missingness to attempt to find association:")
for( focus in 1:length(Iters)){
  Meta_HeatMap_t <- METADATA
  Meta_HeatMap_t$sampleID <- row.names(METADATA)
  Meta_HeatMap_t <- Meta_HeatMap_t[ !is.na(Meta_HeatMap_t[,Iters[focus] ]), ]
  exp <- cqn_counts$E.no.na[,Meta_HeatMap_t$sampleID]
  Meta_HeatMap_t <- Meta_HeatMap_t[, (colnames(Meta_HeatMap_t) %in% c('sex', 'sampleID')) == F]
  
  if( Iters[focus] %in% c( 'ADAD_fam_mut', 'CDR' )){
    Meta_HeatMap_t <- Meta_HeatMap_t[ , 
                                    (colnames(Meta_HeatMap_t) %in% names( which(apply(Meta_HeatMap_t, 2, var) == 0) ))==F,]
  }else{}
  
  preAdjustedSigCovars = runPCAandPlotCorrelations(exp, 
                                                 Meta_HeatMap_t,
                                                 'NULL design(voom-normalized)', 
                                                 isKeyPlot=TRUE, 
                                                 MIN_PVE_PCT_PC = 1)

  preAdjustedSigCovars[["PC_res"]][[2]]$plotData
  if( Iters[focus] %in% preAdjustedSigCovars$significantCovars ){
    writeLines( paste0( Iters[focus], " Is Significantly Associated"))
  }else{
    writeLines( paste0( Iters[focus], " Is NOT Significantly Associated"))
  }
}
## RIN, Sex, disease_cohort, RnaSeqMetrics__PCT_CODING_BASES, RnaSeqMetrics__PCT_INTERGENIC_BASES, RnaSeqMetrics__INTRONIC_BASES, RnaSeqMetrics__PCT_INTRONIC_BASES
# Ind Mut Status: Yes PC11
# APOE: Yes PC5, PC11
# Ind_Mutation: NO
# ADAD_fam_mut: NO
# CDR: NO 
# Age: NO
```


## Normalisation (iterative design)
Since many covariates are correlated, re-normalising and re-adjusting COUNTS with an iterative design matrix
1. Adding Batch and Sex a priori to variable selection
2. Primary variable of interest Diagnosis is excluded from the pool of available covariates for selection


```r
# Primary variable of interest
#RIN, Sex, disease_cohort, RnaSeqMetrics__PCT_CODING_BASES, RnaSeqMetrics__PCT_INTERGENIC_BASES, RnaSeqMetrics__INTRONIC_BASES, RnaSeqMetrics__PCT_INTRONIC_BASES
# Ind Mut Status: Yes PC11
# APOE

source('~/sageseqr/utility_functions/parallelDuplicateCorrelation.R')
primaryVariable <- c( "disease_cohort")
#exclude APOE for now, can condition on it later
FactorCovariates <- c( "sex" )
ContCovariates <- c("RIN", "RnaSeqMetrics__PCT_CODING_BASES", "RnaSeqMetrics__PCT_INTERGENIC_BASES", "RnaSeqMetrics__INTRONIC_BASES", "RnaSeqMetrics__PCT_INTRONIC_BASES")
               
#postAdjustCovars = c( 'sex', 'race', 'yearsEducation', 'bmi' );
postAdjustCovars = 'sex';
# Assign residual covariates
residualCovars = setdiff(preAdjustedSigCovars$significantCovars, c(postAdjustCovars, primaryVariable))
residualCovars <- residualCovars[ residualCovars != 'Infered.Sex']
residualSigCovars = preAdjustedSigCovars
covariatesEffects = preAdjustedSigCovars$Effects.significantCovars[residualCovars]
covariatesEffects <- covariatesEffects[ names(covariatesEffects) != 'Infered.Sex' ]
postAdjustCovars = c(postAdjustCovars, names(which.max(covariatesEffects))) %>% unique()
postAdjustCovars <- postAdjustCovars[ postAdjustCovars != 'Infered.Sex']
loopCount = 0 
cqn_counts$E <- cqn_counts$E[,row.names(METADATA)]
cqn_counts$counts <- cqn_counts$counts[,row.names(METADATA)]
cqn_counts$E.no.na <- cqn_counts$E.no.na[,row.names(METADATA)]
NEW.COUNTS <- cqn_counts$E
notEstimable = c()
# Re-Assign the induvidualIDs
METADATA$individualID <- clean_md[ row.names(METADATA), ]$individualID
METADATA$sex <- METADATA$Infered.Sex
METADATA <- METADATA[ , colnames(METADATA) != 'Infered.Sex']
while(length(residualSigCovars$significantCovars)!=0 && loopCount <= 20){
  writeLines(paste('Using following covariates in the model:',
                   paste(postAdjustCovars, collapse=', '),
                   'as fixed effects and individualID as random effect'))
  
  #writeLines(paste('Using following covariates in the model:',
                   #paste(postAdjustCovars, collapse=', '),
                   #'as fixed effects'))
  
  # Post adjusted design matrix
  DM1 = getDesignMatrix(METADATA[,postAdjustCovars,drop=F],Intercept = F)
  DM1$design = DM1$design[,linColumnFinder(DM1$design)$indepCols]
  
  # Estimate voom weights for dispersion control
  cnts = NEW.COUNTS
  cnts[is.na(cnts)] = 0
  VOOM.GENE_EXPRESSION = voom(cqn_counts$counts, 
                              design=DM1$design,
                              #block= METADATA$individualID,
                              plot=F)#,
                              #na.rm = T)
  
  VOOM.GENE_EXPRESSION$E = cqn_counts$E
  correlation = parallelDuplicateCorrelation(VOOM.GENE_EXPRESSION,
                                             block = METADATA$individualID,
                                             method = 'lmer')
  
  if (!is.nan(correlation$cor)){
    # Re-estimate voom weights
    VOOM.GENE_EXPRESSION = voom(cqn_counts$counts, 
                                design=DM1$design, plot=F, 
                                block = METADATA$individualID, 
                                correlation = correlation$cor)
    
    # Fit linear model using new weights and new design
    VOOM.ADJUSTED.FIT = lmFit(cqn_counts$E,
                              design = DM1$design,
                              weights = VOOM.GENE_EXPRESSION$weights,
                              block = METADATA$individualID, 
                              correlation = correlation$cor)
    
    # Residuals after normalisation
    RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(VOOM.ADJUSTED.FIT,
                                                  cqn_counts$E)
  } else {
    notEstimable = c(notEstimable, postAdjustCovars[length(postAdjustCovars)])
    postAdjustCovars = postAdjustCovars[1:(length(postAdjustCovars)-1)]
  }
  
  # Residual covariates to choose from
  residCovars <- setdiff(c(FactorCovariates,ContCovariates),
                         c(postAdjustCovars, primaryVariable, 'individualID', notEstimable))
  
  # Find PC of residual gene expression and significant covariates that are highly correlated with PCs
  expr = RESIDUAL.GENE_EXPRESSION; expr[is.na(expr)] = 0
  residualSigCovars = runPCAandPlotCorrelations(expr, 
                                                METADATA[, residCovars, drop=F], 
                                                'adjusted design(voom-normalized)',
                                                isKeyPlot=TRUE)
  
  # Add postadjusted covariates (if any)
  residCovars = setdiff(residualSigCovars$significantCovars, c(postAdjustCovars, primaryVariable, notEstimable))
  covariatesEffects = residualSigCovars$Effects.significantCovars[residCovars]
  
  postAdjustCovars = c(postAdjustCovars, names(which.max(covariatesEffects)))
  loopCount = loopCount + 1
}
modelStr <- paste(paste(gsub('_','\\\\_',postAdjustCovars), collapse=', '),
                  'as fixed effects')
tmp <- paste('Using following covariates in the final model:', modelStr)
```

### Sanity check

```r
# Find PC of residual gene expression and significant covariates that are highly correlated with PCs
residualSigCovars = runPCAandPlotCorrelations(expr, 
                                              METADATA,
                                              'adjusted design(voom-normalized)',
                                              isKeyPlot=TRUE)
residualSigCovars[["PC_res"]][[2]]$plotData
```
Coexpression of genes 

```r
cr = cor(t(expr))
hist(cr, main = 'Distribution of correlation between genes', xlab = 'Correlation')
```
PCA of residual data

```r
# Find principal components of expression to plot
PC <- prcomp(expr, scale.=T, center = T)
# Plot first 4 PCs
plotdata <- data.frame(individualID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])
plotdata <- left_join(plotdata, rownameToFirstColumn(METADATA, 'individualID'))
p <- ggplot(plotdata, aes(x=PC1, y=PC2)) 
p <- p + geom_point(aes(color=disease_cohort, shape=apoeGenotype, size=RIN))
p <- p + theme_bw() + theme(legend.position="right")
p
```
Tree based clustering of residual data

```r
# Eucledian tree based analysis
COVARIATES.tmp = data.matrix(METADATA[,c( 'individualID', "sex", "apoeGenotype", "disease_cohort", 'ADAD_fam_mut', "ind_mut_status", "ind_mutation")])
COVARIATES.tmp[is.na(COVARIATES.tmp)] = 0
tree = hclust(as.dist(t(expr)))
cols = WGCNA::labels2colors(COVARIATES.tmp);
WGCNA::plotDendroAndColors(tree, 
                           colors = cols, 
                           dendroLabels = FALSE, 
                           abHeight = 0.80, 
                           main = "Sample dendrogram",
                           groupLabels = colnames(COVARIATES.tmp))
```

### Adjust data with covariates for Network Analysis
Identified covariates are regressed out from the expression matrix for network analysis

```r
# Get design matrix
DESIGN.NET = getDesignMatrix(METADATA[, postAdjustCovars, drop = F], Intercept = F)
DESIGN.NET = DESIGN.NET$design[,linColumnFinder(DESIGN.NET$design)$indepCols]
# Estimate voom weights for dispersion control
cnts = cqn_counts$counts
cnts[is.na(cnts)] = 0
VOOM.NET.WEIGHTS = voom(cnts, design=DESIGN.NET, plot=F)
# Fit linear model using new weights and new design
VOOM.NET.FIT = lmFit(cqn_counts$E,
                     design = DESIGN.NET,
                     weights = VOOM.NET.WEIGHTS$weights)
# Residuals after normalisation
RESIDUAL.NET.GENE_EXPRESSION = residuals.MArrayLM(VOOM.NET.FIT,
                                                  cqn_counts$E)
```


### Differential expression analysis (with Risk Cohort as primary variable)
Differential expression is performed on the primary variable by controlling for covariates identified above

Interpretation of Diagnosis
0. Control: No AD Risk
1. Familial: Familial AD Risk
2. Sporatic: Sporatic AD Risk

Genes that are differentially expressed at an FDR <= 0.05 and folde change of 1.2 are

```r
# Get design matrix
DESIGN = getDesignMatrix(METADATA[, c(primaryVariable[1], postAdjustCovars), drop = F], Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
# Estimate voom weights for dispersion control
VOOM.WEIGHTS = voom(cqn_counts$counts, design=DESIGN$design, plot=F)
#_#VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# Fit linear model using new weights and new design
FIT = lmFit(cqn_counts$E,
            design = DESIGN$design,
            weights = VOOM.WEIGHTS$weights)
# Fit contrast
contrast = makeContrasts(contrasts=c("disease_cohort0-disease_cohort1", 
                                     "disease_cohort0-disease_cohort2",
                                     "disease_cohort1-disease_cohort2"),
                         levels = colnames(FIT$coefficients))
FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR = eBayes(FIT.CONTR)
# Get differnetial expression
DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = dim(counts)[1], confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c('Control-FamilialAD', 'Control-SporaticAD', 'FamilialAD-SporaticAD' )
DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(biomart_results %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gene_gc_content, gene_length) %>%
              unique()) %>%
  dplyr::mutate(Study = 'UW',
                Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
tmp = DE %>%
  dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
  group_by(Comparison, Direction) %>%
  dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
  spread(Direction, FDR_0_05_FC_1.2) 
kable(tmp)
p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
p = p + scale_color_manual(values = c('green','grey','red'))
p = p + facet_grid(.~Comparison, scales = 'fixed')
p
all.diff.exp=NULL
all.fit=NULL
all.diff.exp = c(all.diff.exp, list(AD_Cohort = DE))
all.fit = c(all.fit, list(AD_Cohort = FIT))
```

### Differential expression analysis (with Mutation Status as primary variable)
Differential expression is performed on the mutation status that is profiled as mutated by controlling for covariates identified above and disease cohort

Interpretation of Mutation Status
0. NoMut: No Mutation in (APP, PSEN1, PSEN2, Other)
1. Mut: Has a Mutation in (APP, PSEN1, PSEN2, Other)
9. NA: No Mutation Status Available (These are mostly controls)

Genes that are differentially expressed at an FDR <= 0.05 and folde change of 1.2 are

```r
#Remove NAs (The Control samples don't have mutation status)
METADATA_Use <- METADATA
METADATA_Use$ind_mut_status <- addNA(METADATA_Use$ind_mut_status)
levels(METADATA_Use$ind_mut_status) <- c(levels(METADATA_Use$ind_mut_status), 9)
METADATA_Use$ind_mut_status[ is.na(METADATA_Use$ind_mut_status)] <- 9
METADATA_Use$ind_mut_status <- as.factor(as.character(METADATA_Use$ind_mut_status))

# Get design matrix
DESIGN = getDesignMatrix(METADATA_Use[, c("ind_mut_status", primaryVariable, postAdjustCovars), drop = F], Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
# Estimate voom weights for dispersion control
VOOM.WEIGHTS = voom(cqn_counts$counts, design=DESIGN$design, plot=F)
#_#VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# Fit linear model using new weights and new design
FIT = lmFit(cqn_counts$E,
            design = DESIGN$design,
            weights = VOOM.WEIGHTS$weights)
# Fit contrast
contrast = makeContrasts(contrasts=c("ind_mut_status0-ind_mut_status1", 
                                     "ind_mut_status0-ind_mut_status9",
                                     "ind_mut_status1-ind_mut_status9"),
                         levels = colnames(FIT$coefficients))
FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR = eBayes(FIT.CONTR)
# Get differnetial expression
DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = dim(counts)[1], confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c('NoMut-HasMut', 'NoMut-MissingData', 'HasMut-Missing' )
DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(biomart_results %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gene_gc_content, gene_length) %>%
              unique()) %>%
  dplyr::mutate(Study = 'UW',
                Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
tmp = DE %>%
  dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
  group_by(Comparison, Direction) %>%
  dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
  spread(Direction, FDR_0_05_FC_1.2) 
kable(tmp)
p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
p = p + scale_color_manual(values = c('grey','green','red'))
p = p + facet_grid(.~Comparison, scales = 'fixed')
p
all.diff.exp = c(all.diff.exp, list(MutationStatus = DE))
all.fit = c(all.fit, list(MutationStatus = FIT))
```

### Differential expression analysis (with Gene Mutation as primary variable)
Differential expression is performed on the specific gene that is profiled as mutated by controlling for covariates identified above and disease cohort

Interpretation of Mutation Status
0. NoMut: No Mutation in (APP, PSEN1, PSEN2, Other)
1. PS1: Has a Mutation in PSEN1
2. PS2: Has a Mutation in PSEN1 PSEN2
3. APP: Has a Mutation in PSEN1 APP
4. Other: A mutation in a gene other than APP, PSEN1, or PSEN2
9. NA: No Mutation Status Available (These are mostly controls)

Genes that are differentially expressed at an FDR <= 0.05 and folde change of 1.2 are

```r
#Remove NAs (The Control samples don't have mutation status)
METADATA_Use <- METADATA
METADATA_Use$ind_mutation <- addNA(METADATA_Use$ind_mutation)
levels(METADATA_Use$ind_mutation) <- c(levels(METADATA_Use$ind_mutation), 9)
METADATA_Use$ind_mutation[ is.na(METADATA_Use$ind_mutation)] <- 9
METADATA_Use$ind_mutation <- as.factor(as.character(METADATA_Use$ind_mutation))

# Get design matrix
DESIGN = getDesignMatrix(METADATA_Use[, c("ind_mutation", primaryVariable, postAdjustCovars), drop = F], Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
# Estimate voom weights for dispersion control
VOOM.WEIGHTS = voom(cqn_counts$counts, design=DESIGN$design, plot=F)
#_#VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# Fit linear model using new weights and new design
FIT = lmFit(cqn_counts$E,
            design = DESIGN$design,
            weights = VOOM.WEIGHTS$weights)
# Fit contrast
contrast = makeContrasts(contrasts=c("ind_mutation0-ind_mutation1", 
                                     "ind_mutation0-ind_mutation2",
                                     "ind_mutation0-ind_mutation3",
                                     "ind_mutation0-ind_mutation4",
                                     "ind_mutation0-ind_mutation9",
                                     "ind_mutation1-ind_mutation2",
                                     "ind_mutation1-ind_mutation3",
                                     "ind_mutation1-ind_mutation4",
                                     "ind_mutation1-ind_mutation9",
                                     "ind_mutation2-ind_mutation3",
                                     "ind_mutation2-ind_mutation4",
                                     "ind_mutation2-ind_mutation9",
                                     "ind_mutation3-ind_mutation4",
                                     "ind_mutation3-ind_mutation9",
                                     "ind_mutation4-ind_mutation9"),
                         levels = colnames(FIT$coefficients))
FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR = eBayes(FIT.CONTR)
# Get differnetial expression
DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = dim(counts)[1], confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c("NoMut-PS1", 
              "NoMut-PS2",
              "NoMut-APP",
              "NoMut-Other",
              "NoMut-NA",
              "PS1-PS2",
              "PS1-APP",
              "PS1-Other",
              "PS1-NA",
              "PS2-APP",
              "PS2-Other",
              "PS2-NA",
              "APP-Other",
              "APP-NA",
              "Other-NA")

DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(biomart_results %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gene_gc_content, gene_length) %>%
              unique()) %>%
  dplyr::mutate(Study = 'UW',
                Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
tmp = DE %>%
  dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
  group_by(Comparison, Direction) %>%
  dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
  spread(Direction, FDR_0_05_FC_1.2) 
kable(tmp)
p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
p = p + scale_color_manual(values = c('green','grey','red'))
p = p + facet_grid(.~Comparison, scales = 'fixed')
p
all.diff.exp = c(all.diff.exp, list(MutationType = DE))
all.fit = c(all.fit, list(MutationType = FIT))
```
### Differential expression analysis (with Sex as primary variable)
Differential expression is performed on the mutation status that is profiled as mutated by controlling for covariates identified above and disease cohort

Interpretation of Mutation Status


Genes that are differentially expressed at an FDR <= 0.05 and folde change of 1.2 are

```r
#Remove NAs (The Control samples don't have mutation status)
METADATA_Use <- METADATA

# Get design matrix
DESIGN = getDesignMatrix(METADATA_Use[, c("sex", primaryVariable, postAdjustCovars), drop = F], Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
# Estimate voom weights for dispersion control
VOOM.WEIGHTS = voom(cqn_counts$counts, design=DESIGN$design, plot=F)
#_#VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# Fit linear model using new weights and new design
FIT = lmFit(cqn_counts$E,
            design = DESIGN$design,
            weights = VOOM.WEIGHTS$weights)
# Fit contrast
contrast = makeContrasts(contrasts=c("sexmale-sexfemale"),
                         levels = colnames(FIT$coefficients))
FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR = eBayes(FIT.CONTR)
# Get differnetial expression
DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = dim(counts)[1], confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c( 'Male-Female' )
DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(biomart_results %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gene_gc_content, gene_length) %>%
              unique()) %>%
  dplyr::mutate(Study = 'UW',
                Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
tmp = DE %>%
  dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
  group_by(Comparison, Direction) %>%
  dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
  spread(Direction, FDR_0_05_FC_1.2) 
kable(tmp)
p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
p = p + scale_color_manual(values = c('green','grey','red'))
p = p + facet_grid(.~Comparison, scales = 'fixed')
p
all.diff.exp = c(all.diff.exp, list(Sex = DE))
all.fit = c(all.fit, list(Sex = FIT))
```


### Differential expression analysis (with Disease Cohort and Sex)
Differential expression is performed on the Disease Cohort x Sex variable by controlling for covariates identified above

Genes that are differentially expressed at an FDR <= 0.05 and folde change of 1.2 are

Interpretation of Diagnosis
0. Control: No AD Risk
1. Familial: Familial AD Risk
2. Sporatic: Sporatic AD Risk

Genes that are differentially expressed at an FDR <= 0.05 and folde change of 1.2 are

```r
#Remove NAs (The Control samples don't have mutation status)
METADATA_Use <- METADATA
METADATA_Use$disease_cohort <- as.character(METADATA_Use$disease_cohort)
METADATA_Use$disease_cohort[ METADATA_Use$disease_cohort == '0' ] <- 'Controls'
METADATA_Use$disease_cohort[ METADATA_Use$disease_cohort == '1' ] <- 'Famelial'
METADATA_Use$disease_cohort[ METADATA_Use$disease_cohort == '2' ] <- 'Sporatic'
METADATA_Use$disease_cohort <- as.factor(METADATA_Use$disease_cohort)

METADATA_Use$disease_cohort.Sex = paste(METADATA_Use$disease_cohort,METADATA_Use$sex, sep = '.') %>%
  factor
# Get design matrix
DESIGN = getDesignMatrix(METADATA_Use[, c('disease_cohort.Sex', setdiff(postAdjustCovars, 'Sex')), drop = F], Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
colnames(DESIGN$design) = gsub('disease_cohort.Sex','',colnames(DESIGN$design))
# Estimate voom weights for dispersion control
VOOM.WEIGHTS = voom(cnts, design=DESIGN$design, plot=F)
# Fit linear model using new weights and new design
FIT = lmFit(cqn_counts$E,
            design = DESIGN$design,
            weights = VOOM.WEIGHTS$weights)
# Fit contrast
contrast = makeContrasts( 
  contrasts=c(#"(Controls.male+Famelial.male+Sporatic.male)/3-(Controls.female+Famelial.female+Sporatic.female)/3", 
              "Famelial.male-Controls.male",
              "Famelial.male-Sporatic.male",
              "Famelial.female-Controls.female",
              "Famelial.female-Sporatic.female"),
                         levels = colnames(FIT$coefficients))
FIT.CONTR = contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR = eBayes(FIT.CONTR)
# Get differnetial expression
DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = dim(counts)[1], confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c('Familial-Controls.IN.MALE','Familial-Sporatic.IN.MALE', 'Familial-Controls.IN.FEMALE','Familial-Sporatic.IN.FEMALE')
DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(biomart_results %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol, percentage_gene_gc_content, gene_length) %>%
              unique()) %>%
  dplyr::mutate(Study = 'UW',
                Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < log2(1.2)] = 'NONE'
writeLines('Differentially expressed genes at an FDR 0.05 and logFC 1.2')
tmp = DE %>%
  dplyr::select(ensembl_gene_id, Comparison, Direction) %>%
  group_by(Comparison, Direction) %>%
  dplyr::summarise(FDR_0_05_FC_1.2 = length(unique(ensembl_gene_id))) %>%
  spread(Direction, FDR_0_05_FC_1.2) 
kable(tmp)
p = ggplot(DE, aes(y = -log10(adj.P.Val), x = logFC, color = Direction)) + geom_point() + xlim(c(-1,1))
p = p + scale_color_manual(values = c('green','grey','red'))
p = p + facet_grid(.~Comparison, scales = 'fixed')
p
all.diff.exp = c(all.diff.exp, list(AD_Cohort.Sex = DE))
all.fit = c(all.fit, list(AD_Cohort.Sex = FIT))
```

### Store files in synapse



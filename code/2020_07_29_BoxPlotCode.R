# Plot Boxplots for Katie
#  familial mutation - non-mutation expression (only including PS1 mutations) 
#  for the following specific genes: IL1, IL6, KC, TMEM176, CCL3, miR24-2, LGALS2 showing each individual's data

#User enter ENSGs and Gene Names
Genes <- c( 'ENSG00000115008', 'ENSG00000160712', 'ENSG00000163739', 'ENSG00000002933', 'ENSG00000277632', 'ENSG00000284387', 'ENSG00000100079' )
names( Genes ) <- c( 'IL1A', 'IL6R', 'CXCL1', 'TMEM176A', 'CCL3', 'miR24-2', 'LGALS2' )

library(synapser)
library(ggplot2)
library(dplyr)
library(tidyr)

#Login to Synapse
synLogin()

#Pull Expression
Exp <- read.table(synapser::synGet('syn22269112')$path, sep='\t', header=T)
row.names(Exp) <- Exp$ensembl_gene_id

#User Specifies samples to use (Will try and automate more if needed)
controls <- c( "S4", "S5", "S7", "S9", "S10", "S38", "S41", "S43", "S59", "S62", "S72", "S73" )
PSEN1 <- c( "S1", "S2", "S3", "S6", "S12", "S17", "S19", "S22", "S24", "S25", "S27", "S34", "S44", "S52", "S60", "S63", "S66", "S71", "S74" )

#Build Meta DF
Meta <- as.data.frame( matrix( NA, length(c(controls,PSEN1)), 2 ))
colnames(Meta) <- c('Indv', "Status")

Meta$Indv <- c(controls,PSEN1)
Meta$Status <- c( rep( 'No Mutation',length(controls)), rep('PSEN1',length(PSEN1)) )

#Filter Expression and rename columns to gene names
row.names(Meta) <- Meta$Indv
exp <- as.data.frame(t( Exp[ row.names(Exp)[ row.names(Exp) %in% Genes ],Meta$Indv] ))

genes <- names(Genes)
names(genes) <- Genes

colnames(exp) <- as.character( genes[ colnames(exp) ])

# Make ggplot2 dataframe
Total <- Meta %>%
  left_join( rownameToFirstColumn( exp, 'Indv') ) %>%
  gather(key= 'Gene', value=Gene_Expression, -Indv, -Status) 

#Plot to PDF and Save
pdf("outs/SelectGeneBoxPlot.pdf")
  ggplot( Total, aes(x=Gene, y=Gene_Expression) ) + 
    geom_violin( ) +
    geom_boxplot( width=0.1, outlier.shape = NA ) +
    geom_jitter( shape=16, position=position_jitter(0.2), aes(colour = Status) ) + 
    ggtitle("Non-Mutation Familial Controls Versus \n PSEN1 Familial AD Mutation Carriers") +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Push to Synapse

parentId ="syn22257765"
folderName = 'Plots'
CODE <- synapser::Folder(name = folderName, parentId = parentId)
CODE <- synapser::synStore(CODE)

activityName = 'Gene Plots'
activityDescription = 'Katie\'s Genes of Interest to Plot'

thisFileName <- 'code/2020_07_29_BoxPlotCode.R'
# Github link
thisRepo <- githubr::getRepo(repository = "jgockley62/UW_RNASeq", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0(thisFileName))
#Set Used SynIDs For Provenance
#Need to push covariates to synapse to automate above metadata parts and then add to provenance
used <- c( 'syn22269112' )


# Set annotations
all.annotations = list(
  dataType = 'mRNA',
  dataSubType = 'geneExp',
  summaryLevel = 'gene',
  assay	 = 'RNAseq',
  tissueTypeAbrv	= 'UW', 
  study = 'UW', 
  organism = 'HomoSapiens',
  consortium	= 'UW',
  normalizationStatus	= TRUE,
  normalizationType	= 'CQN',
  rnaquantification = 'RSEM',
  genomeAssemblyID = 'GRCh38'
)


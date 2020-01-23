# FrAT
Frameshift Annotation Tool (FrAT)
setwd("/Users/joselazaro/Documents/Genomes")

#https://cran.r-project.org/web/packages/biomartr/vignettes/Sequence_Retrieval.html#example-ncbi-refseq-2
#install.packages("tinytex")


library("biomartr")
library(GenomicRanges)
library(dplyr)
library(msa)
library(latexpdf)
library(tinytex)
require(Biostrings)
library(stringi)



#is.genome.available(db = "refseq", organism = "Homo sapiens")

setwd("/Users/joselazaro/Documents/Genomes")

#HS.genome.refseq <- getGenome( db       = "refseq",
#                               organism = "Homo sapiens",
#                               path     = file.path("_ncbi_downloads","genomes") )

#HS.genome.refseq ="_ncbi_downloads/genomes/Homo_sapiens_genomic_refseq.fna.gz"

#Human_Genome <- read_genome(file = HS.genome.refseq)

#HS.cds.refseq <- getCDS( db       = "refseq",
#                         organism = "Homo sapiens",
#                         path     = file.path("_ncbi_downloads","CDS"))

HS.cds.refseq <- "_ncbi_downloads/CDS/Homo_sapiens_cds_from_genomic_refseq.fna.gz"

Human_CDS <- read_cds(file     = HS.cds.refseq,
                      obj.type = "Biostrings")

names(Human_CDS)

# download the GFF file of Homo sapiens from refseq
# and store the corresponding file in '_ncbi_downloads/annotation'
#HS.gff.refseq <- getGFF( db       = "refseq", 
#                         organism = "Homo sapiens", 
#                         path = file.path("_ncbi_downloads","annotation"))


HS.gff.refseq <- "_ncbi_downloads/annotation/Homo_sapiens_genomic_refseq.gff.gz"

# import downloaded GFF file
Human_GFF <- read_gff(file = HS.gff.refseq)

length(Human_CDS)

Frameshif_Annotation_tool <- function(Gene,NM_number, Indel_type,Indel_start,Indel_end,Base_change)
{  
####################################################################
#####INPUT VARIABLES FOR FRAMESHIFT ANALYSIS########################
####################################################################
#Gene = "CUBN"
  
#NM_number <- "NM_001081"
NM_number <- paste0(NM_number,".")

#Indel_type = "deletion"

#Indel_start = 5913
#Indel_end = 5916

Variant_Aminoacid_Position = ceiling((Indel_start/3))

####################################################################
#####GET NP NUMBER FROM NM_ACCESION DESIRED PROTEIN BY NP NUMBER####
####################################################################
NM_Accesion_search <- filter(Human_GFF, type == "CDS")
NM_Accesion_search$Grep <- grepl(NM_number,NM_Accesion_search$attribute)
NM_Accesion_search <- filter(NM_Accesion_search,Grep == TRUE)

NM_Data_info <- as.data.frame(NM_Accesion_search$attribute)
colnames(NM_Data_info) <- "attribute"
NM_Data_info$Grep <- grepl(Gene,NM_Data_info$attribute)
NM_Data_info <- filter(NM_Data_info, Grep == TRUE)

NM_Data_info <- as.character(NM_Data_info$attribute[1])

NM_Data_info <- as.vector(unlist(strsplit(NM_Data_info,";")))

#NP_number = "NP_000177" example

NP_number =  gsub(".*-","",NM_Data_info[1])
##########################################
#####GREP DESIRED PROTEIN BY NP NUMBER####
#########################################


Genetic_names <- as.data.frame(sub(" .*","",names(Human_CDS)))

Genetic_names <-as.data.frame(names(Human_CDS))
colnames(Genetic_names) <- "Names"

Genetic_names$Grep <- grepl(NP_number,Genetic_names$Names)
Genetic_names$index <- as.numeric(rownames(Genetic_names))

index_find <- filter(Genetic_names,Grep =="TRUE")

WT_protein <- Human_CDS[index_find$index[1]]
WT_protein <- as.character(WT_protein)

Mutant_protein <- as.character(unlist(strsplit(WT_protein,""), use.names=FALSE))


################################
####START TRANSLATION ANALYSIS####
#################################

####Remove bases from Mutant Proteint
#Mutant_protein <- Mutant_protein[-2799:-2800] example
if(Indel_type=="Deletion")
  {
  Mutant_protein <- Mutant_protein[-((Indel_start:Indel_end))]
  }

if(Indel_type=="Insertion")
  {
  Base_change <- as.character(unlist(strsplit(Base_change,"")))
  Mutant_protein <- append(Mutant_protein,Base_change,after = Indel_start)
  }

################################
####START TRANSLATION ANALYSIS####
#################################


Mutant_protein <- paste(Mutant_protein,collapse="")


######Translate Proteins Wild Type and Mutant
WT_protein_translated <- translate(DNAStringSet(WT_protein))
write.table (WT_protein_translated, "protein_wt.txt")
WT_protein_translated@ranges@NAMES <- "WT"

Mutant_protein_translated <- translate(DNAStringSet(Mutant_protein))
Mutant_protein_translated@ranges@NAMES <- "Mutant"
write.table (Mutant_protein_translated, "protein_mutant.txt")

#####PROCED TO ALIGNMENT
Combined_protein_translated <- c(WT_protein_translated,Mutant_protein_translated)

Alignment <- msa(Combined_protein_translated)

sink("Alignment.txt")
print(Alignment, show="complete")
sink()

#########GET THE STOP CODON########

First_Stop_position <- as.vector(unlist(stri_locate_all(pattern = '*', as.character(unlist(Mutant_protein_translated)), fixed = TRUE)))[1]

Aminoacid_WT_Position <- unlist(strsplit(as.character(unlist(WT_protein_translated)),""))[Variant_Aminoacid_Position]
Aminoacid_MUTANT_Position <- unlist(strsplit(as.character(unlist(Mutant_protein_translated)),""))[Variant_Aminoacid_Position]

Stop_Postion_in_fs <- (First_Stop_position - Variant_Aminoacid_Position + 1)

Frame_Shift_nomenclature <- paste0("p.",Aminoacid_WT_Position,Variant_Aminoacid_Position,Aminoacid_MUTANT_Position,"fs","*",Stop_Postion_in_fs)

####PRINT ACCORDING TO NOMENCLATURE
print(Frame_Shift_nomenclature)
return(Frame_Shift_nomenclature)

}#####END Frameshift Annotation Tool.

Frameshift_list <- read.csv(file = "Frameshift_list_FrAT.csv", header = TRUE, colClasses = c("character","character","numeric","numeric","character"))
Frameshift_list$Result <- "Waiting"

Frameshift_list$Gene[1] <- "CPLANE1"

for (i in 1:nrow(Frameshift_list))
{
  Frameshift_list$Result[i] <- Frameshif_Annotation_tool(Frameshift_list$Gene[i],Frameshift_list$NM_Number[i],Frameshift_list$IndelType[i],Frameshift_list$Start_Del[i],Frameshift_list$End_Del[i],Frameshift_list$Base_change[i])
}



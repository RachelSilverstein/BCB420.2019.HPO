if (!require("readr")) {
  install.packages("readr")
}
library(readr)

header <- read_lines("../data/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt", n_max=1)
cat(header)

HPOdata <- read_tsv("../data/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt", 
                    comment = "#",
                    col_names = c(
                      'entrezID',
                      'entrezSymb',
                      'HPOterm',
                      'HPOtermID'),
                    col_types = cols(
                      entrezID = col_character(),
                      entrezSymb = col_character(),
                      HPOterm = col_character(),
                      HPOtermID = col_character()
                    ))

head(HPOdata)

# A tibble: 6 x 4
#   entrezID entrezSymb HPOterm                         HPOtermID 
#      <dbl> <chr>      <chr>                           <chr>     
# 1     8192 CLPP       Seizures                        HP:0001250
# 2     8192 CLPP       Short stature                   HP:0004322
# 3     8192 CLPP       Primary amenorrhea              HP:0000786
# 4     8192 CLPP       Autosomal recessive inheritance HP:0000007
# 5     8192 CLPP       Microcephaly                    HP:0000252
# 6     8192 CLPP       Hypoplasia of the uterus        HP:0000013

unique_genes <- unique(HPOdata$entrezSymb)
length(unique_genes) # 3924
HPO <- data.frame(matrix(nrow = 3924, ncol = 3))
colnames(HPO) <- c("entrezSymb", "entrezID", "HPOterms")

# collapse the phenotypes into one string... this may take a few minutes to run
for (i in seq_along(unique_genes)) {
  HPO$entrezSymb[i] <- unique_genes[i]
  sel <- HPOdata$entrezSymb == unique_genes[i]
  rows <- HPOdata[sel,]
  phenotypes <- rows$HPOterm
  phenoStr <- paste(phenotypes, collapse = ",")
  HPO$HPOterms[i] <- phenoStr
  ID <- unique(rows$entrezID) # makes sure they are all the same or else you will get an error since you cannot add a multi elemnt vecctor to another vector
  HPO$entrezID[i] <- ID
}

head(HPO)

# get a correct mapping of current (as of Feb. 4 2019) HGNC symbols to entrez IDs
HGNC <- HPOdata <- read_tsv("../data/HGNCdata.txt", 
                            skip = 1,
                            col_names = c(
                              'sym',
                              'prevSym',
                              'synonyms',
                              'entrezID'),
                            col_types = cols(
                              sym = col_character(),
                              prevSym = col_character(),
                              synonyms = col_character(),
                              entrezID = col_character()
                            ))
# We want this data frame to be searchable by entrezID so they should be unique
sel <- duplicated(HGNC$entrezID)
((HGNC$entrezID)[sel]) # all the duplicated ones are just missing values. Remove these.
sel <- is.na(HGNC$entrezID)
# keep the HGNC symbols with missing IDs in a data frame for later
missingIDs <- HGNC[sel,]
# remove them from HGNC
HGNC <- HGNC[!sel,]
# now there should be no duplicate entrezIDs so we can set them as the rownames
rownames(HGNC) <- HGNC$entrezID

# check that all the HGNC symbols in HPO data are correct
HGNC_syms <- HGNC$sym
incorrect_sym <- !(HPO$entrezSymb %in% HGNC_syms)
sum(incorrect_sym) # 38 of the gene symbols in HPO are not current HGNC symbols

# get the rows with incorrect symbols
incorrect <- HPO[incorrect_sym,]
# and remove them from the HPO data frame
nrow(HPO) # 3924
HPO <- HPO[!incorrect_sym,]
nrow(HPO) # 3924 - 38 = 3886... correct

# Map the incorrect symbols to the correct ones using the entrezID
for (i in seq_along(incorrect$entrezID)) {
  entrezID <- incorrect$entrezID[i]
  HGNCrow <- HGNC[entrezID,]
  correctSym <- HGNCrow$sym # if the entrezID is not in the HGNC data frame this will be NA
  if (!is.na(correctSym)) {
    incorrect$entrezSymb[i] <- correctSym
  }
}

# check whether there are still any entrez symbols that are wrong
incorrect_sym <- !(incorrect$entrezSymb %in% HGNC_syms)
sum(incorrect_sym) # 2
# This step corrected 36 of the 38 wrong HGNC symbols. Let's add the corrected ones to the HPO data frame.
corrected <- incorrect[!incorrect_sym,]
incorrect <- incorrect[incorrect_sym,]
HPO <- rbind(HPO, corrected)

# Now let's look at the ones that are still incorrect:
(incorrect)

# entrezSymb  entrezID
# 551     H19-ICR 105259599
# 1481    HBB-LCR 109580095

sel <- grepl('H19-ICR', HGNC$synonyms)
sum(sel) # 0. H19-ICR is not listed as a synonym
sel <- grepl('H19-ICR', HGNC$prevSym)
sum(sel) # 0. H19-ICR is not listed as a previous symbol
('H19-ICR' %in% missingIDs$sym) # FALSE, it is not one of the genes that had a missing entrezID and was removed from the HGNC data frame.

sel <- grepl('HBB-LCR', HGNC$synonyms)
sum(sel) # 0. HBB-LCR is not listed as a synonym
sel <- grepl('HBB-LCR', HGNC$prevSym)
sum(sel) # 0. HBB-LCR is not listed as a previous symbol
('HBB-LCR' %in% missingIDs$sym) # FALSE, it is not one of the genes that had a missing entrezID and was removed from the HGNC data frame.

# We will omit these 2 genes from the final data frame as we were unable to map them to an HGNC symbol.

# Check that all the HGNC symbols in our final mapping are valid:
sum(!(HPO$entrezSymb %in% HGNC$sym)) # 0... correct

# Now we can remove the entrez IDs from the HPO data frame to make the mapping tool
HGNCtoHPO <- HPO[, c("entrezSymb", 'HPOterms')]
head(HGNCtoHPO)

#And name the entrez gene column "HGNCsym" since they are now all HGNC symbols
colnames(HGNCtoHPO) <- c("HGNCsym", "HPOterms")
rownames(HGNCtoHPO) <- HGNCtoHPO$HGNCsym

# save the mapping tool
save(HGNCtoHPO, file = file.path("inst", "extdata", "HGNCtoHPO.RData"))

# ANALYSIS OF ENTIRE HPO GENE SET

# Coverage
nrow(HGNCtoHPO)/nrow(HGNC)*100 #9.441729
# the HPO database only lists phenotypes for about 10% of human genes.

number_of_phenotypes <- 0
for (i in seq_along(HGNCtoHPO$HGNCsym)) {
  phenoVec <- unlist(strsplit(HGNCtoHPO$HPOterms[i], split = ","))
  phenoNum <- length(phenoVec)
  number_of_phenotypes[i] <- phenoNum
}

hist(number_of_phenotypes, breaks = 100, xlab = "Number of Phenotypes", 
     ylab = "Number of Genes", main = "Number of phenotypes associated \n with all HPO annotated genes")

mean(number_of_phenotypes) # 35

# How many phenotypes are shared by more than one gene?
phenoLst <- strsplit(HGNCtoHPO$HPOterms, split = ",") # list of vectors
phenoVec <- unlist(phenoLst, recursive = FALSE)
phenoTab <- table(phenoVec)
phenoTab <- data.frame(phenoTab)
barplot(height = sort(phenoTab$Freq, decreasing = TRUE), xlab = 'Phenotype', ylab = 'Number of Genes',
        main = 'Number of genes associated \n with each phenotype')

mean(phenoTab$Freq) # average phenotype associated with about 20 genes
(20/nrow(HGNCtoHPO) * 100) # average phenotype associated with about 0.5% of the genes in HPO data set

# list the top 10 most common phenotypes
sel <- order(phenoTab$Freq, decreasing = TRUE)
ordered <- phenoTab[sel,]
(ordered[1:10,])

# phenoVec Freq
# 1296 Autosomal recessive inheritance 2187
# 3933         Intellectual disability 1780
# 1294  Autosomal dominant inheritance 1373
# 3155      Global developmental delay 1203
# 6112                        Seizures 1196
# 6239                   Short stature  895
# 4906                       Nystagmus  723
# 4517                    Microcephaly  713
# 3102           Generalized hypotonia  706
# 4694              Muscular hypotonia  682

# ANNOTATION OF EXAMPLE GENE SET
xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")

# Determine which of these genes have HPO annotations
sum(xSet %in% HGNCtoHPO$HGNCsym) # 21
sum(xSet %in% HGNCtoHPO$HGNCsym) / length(xSet) * 100 # 24.70588
# About 25% of the genes in the example gene set have annotated phenotypes
# this is enriched over 2x compared to all genes

# Create a data frame of these annotations
# list phenotypes as the empty string if there are not HPO annotations
xAnnotations <- data.frame(matrix(ncol = 2, nrow = length(xSet)))
colnames(xAnnotations) <- c("HGNCsym", "HPOterms")
for (i in seq_along(xSet)) {
  sym <- xSet[i]
  xAnnotations$HGNCsym[i] <- sym
  phenotypes <- HGNCtoHPO[sym, "HPOterms"]
  if (!is.na(phenotypes)) {
    xAnnotations$HPOterms[i] <- phenotypes
  } else {
    xAnnotations$HPOterms[i] <- ""
  }
}

# save this annotated example set
write_tsv(xAnnotations, "./inst/extdata/xAnnotations.txt")


# ANALYSIS OF EXAMPLE GENE SET

number_of_phenotypes2 <- 0
nonZero <- xAnnotations[xAnnotations$HPOterms != "",]
for (i in seq_along(nonZero$HGNCsym)) {
  phenoVec <- unlist(strsplit(nonZero$HPOterms[i], split = ","))
  phenoNum <- length(phenoVec)
  number_of_phenotypes2[i] <- phenoNum
}

hist(number_of_phenotypes2, breaks = 6, xlab = "Number of Phenotypes", 
     ylab = "Number of Genes", main = "Number of phenotypes associated \n with Example Gene Set")

mean(number_of_phenotypes2) # 41

# How many phenotypes are shared by more than one gene?
phenoLst <- strsplit(nonZero$HPOterms, split = ",") # list of vectors
phenoVec <- unlist(phenoLst, recursive = FALSE)
phenoTab <- table(phenoVec)
phenoTab <- data.frame(phenoTab)

barplot(height = sort(phenoTab$Freq, decreasing = TRUE), space = 0.5, xlab = 'Phenotype', ylab = 'Number of Genes',
        main = 'Number of genes associated \n with each phenotype in Example Set')

mean(phenoTab$Freq) # average phenotype associated with about 1.5 genes
(1.5/85 * 100) # average phenotype associated with 1.8 % of genes in list (enriched from 0.5% for all HPO data)

# list the top 10 most common phenotypes
sel <- order(phenoTab$Freq, decreasing = TRUE)
ordered <- phenoTab[sel,]
(ordered[1:10,])

# phenoVec Freq
# 323         Intellectual disability   14
# 88  Autosomal recessive inheritance   10
# 87   Autosomal dominant inheritance    9
# 492                        Seizures    8
# 192                      Dysarthria    7
# 263      Global developmental delay    7
# 83                           Ataxia    6
# 59                        Agitation    5
# 65    Amyotrophic lateral sclerosis    5
# 72                          Anxiety    5

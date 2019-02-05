# `BCB420.2019.HPO`


###### [Rachel Silverstein](https://orcid.org/0000-0001-5724-2252), rachel.silverstein@mail.utoronto.ca


<img src="https://cmg.broadinstitute.org/sites/default/files/HPO-logo-black.png" width="200" height="40" alt="HPO logo">

<small>This project uses the Human Phenotype Ontology (Dec 21, 2018 11:56:54 AM)</small>

&nbsp;


### About the Data Source

----

TO DO

&nbsp;

### Data Download and Import

----

The version of HPO used for this project is Build #153 (Dec 21, 2018 11:56:54 AM). A 'build history' is available at [Project hpo.annotations.monthly](http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/), where you can select Build #153. Note that using the [downloads tab](https://hpo.jax.org/app/download/annotation) of the [HPO page](https://hpo.jax.org/app/) may not provide the same results as this page is updated monthly with new versions of the annotations file.

Once you have navigated to [Build #153](http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/153/), select 'Expand all' to see a list of files available for download. Before downloading data, note the [license restrictions](https://hpo.jax.org/app/license) that apply to all files provided by the HPO project. 

Download the file entitled: 
```text
ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt
```
Save this file in a sister directory of your working directory entitled `data`. `file.exists("../data/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt")` should return `TRUE`.

Install the `readr` package in order to read the data into a data frame in R.

```R
if (!require("readr")) {
  install.packages("readr")
}
library(readr)
```

Now we will look at the first line of the text file to determine the structure of the data frame that we want to create.

```R
header <- read_lines("../data/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt", n_max=1)
cat(header)

# Format: entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-Name<tab>HPO-Term-ID
```

We can see that the data file includes 4 columns, the entrez gene ID, entrez gene symbol, HPO term, and the HPO term symbol. Now that we know the contents of the columns, read in the data using `read_tsv()` from `readr`.

```R
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
```
&nbsp;

The data set retrieved from HPO contains multiple rows with the same gene since each phenotype is listed in a separate row. Since we want to make a data frame with one row for each HGNC symbol, we must condense all the different phenotypes for each gene into the same row. We will do this by creating a single string to represent all the phenotypes and separate them by the comma `,` character so that they can be easily split later on.

```R
# See how many unique genes there are in the HPO data
unique_genes <- unique(HPOdata$entrezSymb)
length(unique_genes) # 3924
# Initialize a data frame of the correct dimensions
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
```

The HPO data includes entrez ID and entrez gene symbol. NCBI reports that gene symbols are provided by HGNC so most of these symbols are likely already what we need. However, it is possible that some of the symbols are out of date since we do not know how often NCBI updates their symbols and we do not know how long ago HPO created this list of NCBI symbols.

Thus, we will create our own mapping of HGNC symbols to entrez gene IDs to make sure that we are able to map all of the genes in HPO to the most current HGNC symbol.

First, we must retrieve the relevant data from the HGNC website. Navigate to the [https://www.genenames.org/download/custom/](custom downloads section of the HGNC website). Select only "Approved symbol,	Previous symbols, and	Synonyms" under the data provided by HGNC. Select "NCBI Gene ID(supplied by NCBI)" under data downloaded from external sources. These are the entrez IDs which we will need to use to map the data. Then download the resulting data file and save it as `"HGNCdata.txt"` to a sister directory called "data" of you working directory. This should be saved in the same place as `ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt` which was saved previously.

Now that the data is downloaded, proceed to load this data into R.

&nbsp;

```R
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
                            
```


&nbsp;

### Mapping HGNC symbols to HPO terms

----

Now let's clean up the data frame so that all of the entrez IDs are unique and can be used for mapping...

```R
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
```
&nbsp;

Now that all the entrez IDs in the HGNC data frame are unique, we can use them to map the entrez symbols in the HPO data to their respective HGNC symbols. First, though, we will see which of the HPO data are already mapped to appropriate HGNC symbols.

```R
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
```

&nbsp;

### Analysis of the Entire Gene Set Provided by HPO

----

Now we have created a data frame that maps all of the HPO phenotype annotations to current HPNC symbols. First, let's look at what fraction of all HGNC symbols have HPO annotations.

```R
# Coverage
nrow(HGNCtoHPO)/nrow(HGNC)*100 # 9.441729
# the HPO database only lists phenotypes for about 10% of human genes.
```

The HPO has annotations for about 10% of human genes. This makes sense since we do not know how and if most genes contribute to phenotype in humans.

&nbsp;

Now let's see how many phenotypes are associated with each gene in the HPO data set.

```R
number_of_phenotypes <- 0
for (i in seq_along(HGNCtoHPO$HGNCsym)) {
  phenoVec <- unlist(strsplit(HGNCtoHPO$HPOterms[i], split = ","))
  phenoNum <- length(phenoVec)
  number_of_phenotypes[i] <- phenoNum
}

hist(number_of_phenotypes, breaks = 100, xlab = "Number of Phenotypes", 
     ylab = "Number of Genes", main = "Number of phenotypes associated \n with all HPO annotated genes")
```

![](./inst/img/fig_1.png?sanitize=true "Phenotype distribution")

&nbsp;

The number of phenotypes associated with each gene seems to decay exponentially. A few genes are associated with many phenotypes while many genes are associated with few phenotypes. It is not clear whether this reflects the biology of the genes or whether the genes with many phenotype annotations are simply better-studied than other genes.

```R
mean(number_of_phenotypes) # 35

```

The mean number of phenotypes for each gene is 35.

&nbsp;

Now let's see how many genes are associated with each phenotype.
```R
# How many genes associated with each phenotype?
phenoLst <- strsplit(HGNCtoHPO$HPOterms, split = ",") # list of vectors
phenoVec <- unlist(phenoLst, recursive = FALSE)
phenoTab <- table(phenoVec)
phenoTab <- data.frame(phenoTab)
barplot(height = sort(phenoTab$Freq, decreasing = TRUE), xlab = 'Phenotype', ylab = 'Number of Genes',
        main = 'Number of genes associated \n with each phenotype')
```

![](./inst/img/fig_2.png?sanitize=true "gene distribution")

&nbsp;

What is the mean number of genes for which a phenotype is associated?
```R
mean(phenoTab$Freq) # average phenotype associated with about 20 genes
(20/nrow(HGNCtoHPO) * 100) # average phenotype associated with about 0.5% of the genes in HPO data set
```

The average phenotype is associated with abot 20 genes which corresponds to about 0.5% of the total genes that have HPO annotations.

&nbsp;

What phenotypes are the most common?

```R
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
```

&nbsp;

### Annotations of the Example Gene Set

----

The example gene set (xSet) was copy and pasted from the BCB420 resources repository on GitHub.

These genes are functionally related in that they are associated with the phagosome / lysosome fusion system.

```
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
```

Let's determine how many of the genes in the example gene set have HPO annotated phenotypes.

```R
# Determine which of these genes have HPO annotations
sum(xSet %in% HGNCtoHPO$HGNCsym) # 21
sum(xSet %in% HGNCtoHPO$HGNCsym) / length(xSet) * 100 # 24.70588
# About 25% of the genes in the example gene set have annotated phenotypes
# this is enriched over 2x compared to all genes
```

About 25% of the genes in the example set have HPO phenotypes which is enriched over 2 times compared to the total gene set. Perhaps this is because they are important genes since they are involved in an important biological process and thus more likely to have phenotypes. Alternatively, they may be more well-studied genes than the average gene.

&nbsp;

Finally, let's create a data frame of the annotations for the example gene set and save it.

```R
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
```

The annotated gene set is saved as a tsv as `"./inst/extdata/xAnnotations.txt"`.

&nbsp;

### Analysis of the Example Gene Set

----

How many phenotypes are associated with each gene?

```R
number_of_phenotypes2 <- 0
nonZero <- xAnnotations[xAnnotations$HPOterms != "",]
for (i in seq_along(nonZero$HGNCsym)) {
  phenoVec <- unlist(strsplit(nonZero$HPOterms[i], split = ","))
  phenoNum <- length(phenoVec)
  number_of_phenotypes2[i] <- phenoNum
}

hist(number_of_phenotypes2, breaks = 6, xlab = "Number of Phenotypes", 
     ylab = "Number of Genes", main = "Number of phenotypes associated \n with Example Gene Set")
     
```
![](./inst/img/fig_3.png?sanitize=true "pheno distribution example")

The number of phenotypes per gene also decays (exponentially?) as in the larger HPO data set.

```R
mean(number_of_phenotypes2) # 41
```

The mean number of phenotypes per gene is larger than for the entire HPO data set. Once again, perhaps this is because they are well studied or biologically important genes (or both).

```R
# How many genes per phenotype?
phenoLst <- strsplit(nonZero$HPOterms, split = ",") # list of vectors
phenoVec <- unlist(phenoLst, recursive = FALSE)
phenoTab <- table(phenoVec)
phenoTab <- data.frame(phenoTab)

barplot(height = sort(phenoTab$Freq, decreasing = TRUE), space = 0.5, xlab = 'Phenotype', ylab = 'Number of Genes',
        main = 'Number of genes associated \n with each phenotype in Example Set')
```

![](./inst/img/fig_4.png?sanitize=true "gene distribution example")

```R
mean(phenoTab$Freq) # average phenotype associated with about 1.5 genes
(1.5/85 * 100) # average phenotype associated with 1.8 % of genes in list (enriched from 0.5% for all HPO data)
```

The average phenotype is associated with 1.5 genes which corresponds to 1.8 percent of genes in the example data set. This is more than the 0.5 percent of genes that are associated with the average phenotype in the entire HPO dataset. This is as expected since all of the genes in the example set are involved in the same process. Thus, we would expect them to have more phenoypes in common than a group of unrelated genes.

&nbsp;

What are the most common phenotypes associated with the example gene set?

```R
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
```

We see some phenotypes that are the same as the most common phenotypes the entire HPO data set but also some that are unique to this group.
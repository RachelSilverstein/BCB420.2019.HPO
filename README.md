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

unique_genes <- unique(HPOdata$entrezSymb)
length(unique_genes) # 7090
HPO <- data.frame(matrix(nrow = 7090, ncol = 3))
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
```

Now we have a data frame containing the HPO annotation data that is ready to be mapped to HGNC symbols.

&nbsp;

### Mapping Entrez Genes to HGNC Symbols

----

The data from HPO contains "entrez gene ID" and "entrez gene symbol" for each phenotype-gene mapping. Although the gene symbols for NCBI genes are provided by HGNC, it is possible that the symbol is not up to date in the data provided by HPO. In order to check that the HGNC symbols are matched to the correct and most up-to-date entrez genes, we will use the biomaRt package to map HGNC symbols to entrez genes and compare our sesults to the mapping provided in the HPO data set.

&nbsp;

#### 1. Creating a Mapping tool using biomaRt

----

The goal is to create a data frame that can act as a map to translate HGNC symbols to entrez gene IDs.

<ol>
<li>The data frame must be searchable by HGNC symbol (these will be the rownames and must be unique).</li>
<li> If the HGNC symbol maps to multiple entrez genes, we want to be able to retrieve all of the (possibly different) phenotypes associated with each of the entrez gene IDs with a single search of their common HGNC symbol. In order to facilitate this, we will concatenate entrez genes associated with the same HGNC symbol into the same string such that they can be stored in a single row.</li>
</ol>

We will thus create a data frame, `map`, of the following structure:

2 columns:
<ol>
<li> `hgnc_symbol` (character vector): complete list of unique HGNC symbols from ensembl data base, accessed via biomaRt. These will also be the rownames of the data frame.
<li> `entrezgene` (character vector): Each element is a string of Entrez gene IDs separated by `,` character if there is more than one
</ol>

&nbsp;

##### 1.1 Download Ensembl Genes and Attributes

Start by installing the biomaRt package and BiocInstaller. We will use biomaRt to access the Ensembl database of human genes in order to create a map between the HGNC symbols and entrez gene IDs.

```R
if (!require("BiocInstaller")) {
  install.packages("BiocInstaller")
}
library("BiocInstaller")

if (!require("biomaRt")) {
  biocLite('biomaRt')
}
library("biomaRt")
```

Select ensembl as the database and fetch the `hsapiens_gene_ensembl` data set. Then search for attributes containing the words "hgnc" and "entrez" to determine the names of the attributes we wish to retrieve.

```R
genes <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(genes)$name
(grep("hgnc", attributes, ignore.case=TRUE, value=TRUE))
(grep("entrez", attributes, ignore.case=TRUE, value=TRUE))

# [1] "hgnc_id"         "hgnc_symbol"     "hgnc_trans_name"
# [1] "entrezgene_trans_name" "entrezgene"
```

The attributes that we want are called "hgnc_symbol" and "entrezgene". Fetch these attributes.

```R
map <- getBM(attributes=c("hgnc_symbol", "entrezgene"), mart = genes)
head(map)

# hgnc_symbol entrezgene
# 1       MT-TF         NA
# 2     MT-RNR1       4549
# 3       MT-TV         NA
# 4     MT-RNR2       4550
# 5      MT-TL1         NA
# 6      MT-ND1       4535
```

The resulting data frame contains a mapping from an HGNC symbol to an entrez gene in each row. Importantly, HGNC symbols that do not map to an entrez gene contain `NA` in the `entrezgene` column and entrez genes that do not map to an HGNC symbol contain the empty string `""` in the `hgnc_symbol` column.

&nbsp;

##### 1.2 Resolving Duplicated Entrez Genes

Since we need to search by HGNC symbol, there must be no duplicates in the HGNC symbol column so we will concatinate any entrez gene IDs that share the same HGNC symbol into the same row. Since any HGNC symbols that are missing ( `""`) will appear to be duplicated, we must first remove these rows from the data frame.

```R
(sum(map$hgnc_symbol == "")) # there are 1243 entrez genes that are 
# not mapped to an HGNC symbol

sel <- map$hgnc_symbol == ""
map <- map[!sel,] # delete them
(sum(map$hgnc_symbol == "")) # 0

# now see whether there are any duplicated HGNC symbols
(sum(duplicated(map$hgnc_symbol))) # 228

# what are the symbols that are duplicated?
dup <- duplicated(map$hgnc_symbol)
dup_names <- map[dup,]$hgnc_symbol
(dup_names[1:10]) # take a look at the first 10

# "PPP1R18"  "PDXK"     "HLA-DQA2" "TP53TG3C" "TP53TG3C" "TP53TG3C" "TP53TG3C"
# "ERCC6"    "TP53TG3E" "TP53TG3B"

# select all the HGNC symbols that are not unique (cannot select them directly # with duplicated() as this will leave out the first occurence of each one)
sel <- map$hgnc_symbol %in% dup_names
duplicated <- map[sel,]
nrow(duplicated) # 421
head(duplicated)

#  hgnc_symbol entrezgene
# 218     PPP1R18  107987457
# 219     PPP1R18     170954
# 405        PDXK  105372824
# 406        PDXK       8566
# 637    HLA-DQA2       3118
# 638    HLA-DQA2       3117

# remove the duplicated ones from the map
nrow(map) # 37807
map <- map[!sel,]
nrow(map) # 37807 - 421 = 37386... correct

# create a list of the HGNC symbols that are duplicated
dup_names <- unique(dup_names)
(length(dup_names)) # 193 unique names

# initialize a new data frame with a row for each unique symbol in duplicated
dup_summ <- data.frame(matrix(ncol = 3, nrow = 193))
colnames(dup_summ) <- c("hgnc_symbol", "entrezgene", "multi")

for (i in seq_along(dup_names)) {
  # choose the rows that have the same HGNC symbol and the entrezgene is not NA
  sel <- (duplicated$hgnc_symbol == dup_names[i]) & (!is.na(duplicated$entrezgene))
  rows <- duplicated[sel,]
  entrezIDs <- paste(rows$entrezgene, collapse = ",")
  dup_summ$hgnc_symbol[i] <- dup_names[i]
  dup_summ$entrezgene[i] <- entrezIDs
}

head(dup_summ) # take a look

#   hgnc_symbol                                    entrezgene
# 1     PPP1R18                               107987457,170954
# 2        PDXK                                 105372824,8566
# 3    HLA-DQA2                                      3118,3117
# 4    TP53TG3C 102724127,102724101,102723713,102723655,653550
# 5       ERCC6                                    267004,2074
# 6    TP53TG3E                            102724101,102723655

dup_summ$multi <- rep(TRUE, times = nrow(dup_summ))
map$multi <- rep(FALSE, times = nrow(map))
# add these rows back to the map
map <- rbind(map, dup_summ)
nrow(map) # 37386 + 193 = 37579... correct

# assign row names (and confirm that the symbols are unique)
rownames(map) <- map$hgnc_symbol
```

We have now created the `map` data frame according to the specifications listed above.

&nbsp;

### 2. Assign HPO Terms to HGNC symbols

The data frame, `HPO`, should be available in your global environment. If not, refer to the steps in the "Data Download and Import" section above.

First, determine the set of HGNC symbols for which we want to retrieve HPO phenotypes. We will use the example gene set provided by [https://github.com/hyginn/](Boris Steipe on GitHub).

```R
# load the example gene set
myURL <- paste0("https://github.com/hyginn/", "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))
head(HGNC)
row_number <- nrow(HGNC) # 27087


# initialize a new data frame to hold the final mappings
mapping <- data.frame(matrix(ncol = 2, nrow = row_number))
colnames(mapping) <- c("HGNC", "phenotypes")

mapping$HGNC <- HGNC$sym

for (i in seq_along(mapping$HGNC)) {
  HGNCsymb <- mapping$HGNC[i]
  entrezgene <- map[HGNCsymb, "entrezgene"]
  if (map$multi[i] == TRUE) {
    entrezVector <- unlist(strsplit(entrezgene, split = ","))
    sel <- HPO$entrezID %in% entrezVector
    HPOrows <- HPO[sel,]
    phenotypes <- HPOrows$HPOterms
    phenoVector <- unlist(strsplit(phenotypes, split = ","))
    phenotypes <- unique(phenotypes)
    phenoStr <- paste(phenotypes, collapse = ",")
    phenotypes <- phenoStr 
  } else {
    sel <- HPO$entrezID == entrezgene
    HPOrow <- HPO[sel,]
    phenotypes <- HPOrow$HPOterms
  }
  mapping$phenotypes[i] <- phenotypes
  mapping$HGNC[i] <- HGNCsymb
}




> 


```







<!--end-->
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
```

This returns: 

`#Format: entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-Name<tab>HPO-Term-ID`

We can see that the data file includes 4 columns, the entrez gene ID, entrez gene symbol, HPO term, and the HPO term symbol. Now that we know the contents of the columns, read in the data using `read_tsv()` from `readr`.

```R
HPOdata <- read_tsv("../data/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt", comment = "#",
col_names = c(
  'entrezID',
  'entrezSymb',
  'HPOterm',
  'HPOtermID'),
col_types = cols(
  entrezID = col_double(),
  entrezSymb = col_character(),
  HPOterm = col_character(),
  HPOtermID = col_character()
))

head(HPOdata)
```

```text
# A tibble: 6 x 4
  entrezID entrezSymb HPOterm                         HPOtermID 
     <dbl> <chr>      <chr>                           <chr>     
1     8192 CLPP       Seizures                        HP:0001250
2     8192 CLPP       Short stature                   HP:0004322
3     8192 CLPP       Primary amenorrhea              HP:0000786
4     8192 CLPP       Autosomal recessive inheritance HP:0000007
5     8192 CLPP       Microcephaly                    HP:0000252
6     8192 CLPP       Hypoplasia of the uterus        HP:0000013
```

Now we have a tibble containing the HPO annotation data that is ready to be mapped to HGNC symbols.

&nbsp;

### Mapping Entrez Genes to HGNC Symbols

----

The data from HPO contains "entrez gene ID" and "entrez gene symbol" for each phenotype-gene mapping. Although the gene symbols for NCBI genes are provided by HGNC, it is possible that the symbol is not up to date in the data provided by HPO. In order to check that the HGNC symbols are matched to the correct and most up-to-date entrez genes, we will use the biomaRt package to map HGNC symbols to entrez genes and compare our sesults to the mapping provided in the HPO data set.

&nbsp;

#### 1. Creating a Mapping tool using biomaRt

----

##### Download Ensembl Genes and Attributes

The goal is to create a data frame that can act as a "dictionary" to translate HGNC symbols to entrez gene IDs.

<ol>
<li>The data frame must be searchable by HGNC symbol (these will be the rownames and must be unique).</li>
<li> If the HGNC symbol maps to multiple entrez genes, we will concatenate the entrez genes into the same string such that they can be stored in a single row.</li>
</ol>

Start by installing the biomaRt package and BiocInstaller for its installation. We will use biomaRt to access the Ensembl database of human genes in order to create a map between the HGNC symbols and entrez gene IDs.

```R
if (!require("BiocInstaller")) {
  install.packages("BiocInstaller")
}

if (!require("biomaRt")) {
  biocLite('biomaRt')
}
library("biomaRt")
library("BiocInstaller")
```

Select ensembl as the database and fetch the `hsapiens_gene_ensembl` data set. Then search for attributes containing the words "hgnc" and "entrez" to determine the names of the attributes we wish to retrieve.

```R
genes <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(genes)$name
(grep("hgnc", attributes, ignore.case=TRUE, value=TRUE))
(grep("entrez", attributes, ignore.case=TRUE, value=TRUE))
```

```text
[1] "hgnc_id"         "hgnc_symbol"     "hgnc_trans_name"
[1] "entrezgene_trans_name" "entrezgene"
```

The attributes that we want are called "hgnc_symbol" and "entrezgene". Fetch these attributes.

```R
BM <- getBM(attributes=c("hgnc_symbol", "entrezgene"), mart = genes)
head(BM)

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

##### Resolving Duplicated Entrez Genes

Since we need to search by HGNC symbol, there must be no duplicates in the HGNC symbol column so we will concatinate any entrez gene IDs that share the same HGNC symbol into the same row.

Since any HGNC symbols that are missing ( `""`) will appear to be duplicated, we must first remove these rows from the data frame.

```R
(sum(BM$hgnc_symbol == "")) # there are 1243 entrez genes that are 
# not mapped to an HGNC symbol

sel <- BM$hgnc_symbol == ""
BM <- BM[!sel,] # delete them
(sum(BM$hgnc_symbol == "")) # 0

# now see whether there are any duplicated HGNC symbols
(sum(duplicated(BM$hgnc_symbol))) # 228

# what are the symbols that are duplicated?
dup <- duplicated(BM$hgnc_symbol)
dup_names <- BM[dup,]$hgnc_symbol
(dup_names[1:10]) # take a look

# "PPP1R18"  "PDXK"     "HLA-DQA2" "TP53TG3C" "TP53TG3C" "TP53TG3C" "TP53TG3C"
# "ERCC6"    "TP53TG3E" "TP53TG3B"

# select all the HGNC symbols that are not unique (cannot select them directly with duplicated() as this will leave out the first occurence of each one)
sel <- BM$hgnc_symbol %in% dup_names
duplicated <- BM[sel,]
nrow(duplicated) # 421
head(duplicated)

#  hgnc_symbol entrezgene
# 218     PPP1R18  107987457
# 219     PPP1R18     170954
# 405        PDXK  105372824
# 406        PDXK       8566
# 637    HLA-DQA2       3118
# 638    HLA-DQA2       3117

# remove the duplicated ones from BM
nrow(BM) # 37807
(37807 - 421) # 37386
BM <- BM[!sel,]
nrow(BM) # 37386... correct

# create a list of the HGNC symbols that are duplicated
dup_names <- unique(dup_names)
(length(dup_names)) # 193 unique names

# initialize a new data frame with a row for each unique symbol in duplicated
summ <- data.frame(matrix(ncol = 2, nrow = 193))
colnames(summ) <- c("hgnc_symbol", "entrezgenes")

for (i in seq_along(dup_names)) {
  sel <- duplicated[duplicated$hgnc_symbol == dup_names[i],]
  entrezIDs <- paste(sel$entrezgene, collapse = ",")
  summ$hgnc_symbol[i] <- dup_names[i]
  summ$entrezgenes[i] <- entrezIDs
}

head(summ) # take a look

#   hgnc_symbol                                    entrezgenes
# 1     PPP1R18                               107987457,170954
# 2        PDXK                                 105372824,8566
# 3    HLA-DQA2                                      3118,3117
# 4    TP53TG3C 102724127,102724101,102723713,102723655,653550
# 5       ERCC6                                    267004,2074
# 6    TP53TG3E                            102724101,102723655
```


```R
# get the names of the genes that are duplicated
# Note: we cannot select the duplicated genes directly because this would
# leave out the first occurence of each of the duplicated genes
sel <- duplicated(BM$entrezgene)
dup_names <- BM[sel,]$entrezgene
# find all the occurences of the duplicated genes
sel <- BM$entrezgene %in% dup_names
sum(sel) # 385 duplicated genes (including the first occurence)

dup <- BM[sel,] # put the duplicated genes into a new data frame
BM <- BM[!sel,] # and remove them from BM
```

We can decide which HGNC symbol to use by seeing which entrez-HGNC pairs also share the same ensembl ID.

```R
# get the ensembl IDs of the entrez genes in dup
ensembl <- getBM(
  attributes = c("hgnc_symbol", "entrezgene", "ensembl_gene_id"),
  mart=genes,
  filters = c("entrezgene"),
  values = dup$entrezgene)
head(ensembl)

# hgnc_symbol entrezgene ensembl_gene_id
# 1         DDT  100037417 ENSG00000275003
# 2        DDTL  100037417 ENSG00000275758
# 3         DDT  100037417 ENSG00000099977
# 4        DDTL  100037417 ENSG00000099974
# 5  FAM83H-AS1  100128338 ENSG00000282685
# 6      IQANK1  100128338 ENSG00000282519
```














<!--end-->
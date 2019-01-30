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

#### 1. Creating a "dictionary" using biomaRt

----

##### Download Ensembl Genes and Attributes

The goal is to create a data frame that can act as a "dictionary" to translate entrez gene IDs to HGNC symbols. This must meet the following requirements:

<ol>
<li>The data frame must be searchable by entrez gene ID (these will be the rownames).</li>
<li> The entrez gene IDs must map to a unique HGNC symbol. </li>
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
BM <- getBM(attributes=c("hgnc_symbol", "entrezgene"), mart=genes)
head(BM)
```

```text
hgnc_symbol entrezgene
1       MT-TF         NA
2     MT-RNR1       4549
3       MT-TV         NA
4     MT-RNR2       4550
5      MT-TL1         NA
6      MT-ND1       4535
```

The resulting data frame contains a mapping from an HGNC symbol to an entrez gene in each row. Importantly, HGNC symbols that do not map to an entrez gene contain `NA` in the `entrezgene` column and entrez genes that do not map to an HGNC symbol contain the empty string `""` in the `hgnc_symbol` column.

&nbsp;

##### Resolving Ambiguities

In order for the mapping tool to be effective, every entrez gene must map to a unique HGNC symbol. Any entrez gene that maps to multple HGNC symbols will appear in more than one row of our data frame so we can check whether the mapping is unique by looking for duplicated entrez genes. However, we must first remove any rows in which the entrez gene is `NA` (HGNC synbols for which there is no entrez gene) as these will also appear as duplications.

```R
(sum(is.na(BM$entrezgene))) # there are 13211 HGNC symbols not mapped to an entrez gene
sel <- is.na(BM$entrezgene)
BM <- BM[!sel,] # remove them
(sum(is.na(BM$entrezgene))) # 0

# Now check the number of duplicated entrez genes.
(sum(duplicated(BM$entrezgene))) # there are 201 duplicated entrez genes
```

There are 201 duplicated entrez genes, which is too many to resolve by hand. We will put the entrez genes that are not unique into their own data frame so that we can manipulate them and then add them back to the main data frame.

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














<!--end-->
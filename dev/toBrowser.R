# toBrowser.R
#
# Purpose: render a markdown file
#
# Version: 1.0
# Date:    2019-01-20
# Author: Boris Steipe (ORCID: 0000-0002-1134-6758)
# License: (c) Author (2019) + MIT
#
#
# ==============================================================================

# NO SIDE EFFECTS:
# This script can be safely source()'d to define the functions it contains and
# install.packages()/run library() as required.
# All other code will not be executed unless this is done interactively.

# ====  PACKAGES  ==============================================================
# Load all required packages.

if (! requireNamespace("rmarkdown", quietly=TRUE)) {
  install.packages("rmarkdown")
}
# Package information:
#  library(help = rmarkdown)       # basic information


# ====  FUNCTIONS  =============================================================

toBrowser <- function(FN = "README.md") {
  # Purpose:
  #     Render a markdown file to html in tempdir() and display it in the
  #     user's default browser.
  # Parameters:
  #     FN:     char   filename of a markdown file. Defaults to README.md
  # Value:
  #     result: NULL (invisible). The function is used for its side-effect
  #             of opening an html file.
  
  html <- file.path(tempdir(), "tmp.html")
  rmarkdown::render(FN, output_format = "html_document", output_file = html)
  browseURL(html)
  return(invisible(NULL))
}


# ====  TESTS  =================================================================
if (FALSE) {
  
  toBrowser()  # Executing this line should render and open the file
  
}


# [END]
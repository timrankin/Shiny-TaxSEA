# Shiny-TaxSEA
Shiny-TaxSEA is a Shiny frontend to [TaxSEA: Taxon Set Enrichment Analysis](https://github.com/feargalr/TaxSEA)

## Quick Start
Access at <https://shiny.taxsea.app>

Start by analysing test data, the plots update as you select up to 8 taxa in the table.

When you're ready to analyse your own data, supply a .csv or .xlsx file with the following columns (column title/header doesn't matter, but order does!): Taxa (e.g. species/genus), rank (e.g. log2 fold changes), P value, Padj (FDR). You can also download test data to see the expected format.

## Running locally
The project is not packaged nicely yet, so to run it locally clone the repository:
```
git clone https://github.com/timrankin/Shiny-TaxSEA
```

Then open the project `Shiny-TaxSEA.Rproj` in RStudio and install dependencies:
```{r output}
library(devtools)
install_github("feargalr/TaxSEA")
install.packages(c('shiny', 'tidyverse', 'ggrepel', 'bslib', 'bsicons', 'openxlsx2', 'DT'))
```

Finally, open `app.r` and click 'Run App'.
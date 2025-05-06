# BONOBO

<img src="www/bonobo.png" align="right" width="120" height="139" alt="Bonobo Image" />
BONOBO is a Shiny interactive web application that analyzes and visualizes metagenomics classification results from tools like Kraken2 and compares them to multi-BAM file alignments. It facilitates intuitive comparison and overlay of organism-level taxonomic classifications and sequencing alignment data across multiple samples.
```bash
# run bonobo
Rscript -e "shiny::runApp('app.R', port=5000, host='0.0.0.0')"
```

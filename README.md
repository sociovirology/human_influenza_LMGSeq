# Genotype by Barcode Sequencing (GbBSeq) applied to high-throughput quantification of genetic exchange among diverse human Influenza A virus strains 

This github repository includes code (and links to data) from the manuscript:  
Influenza A virus reassortment is strain dependent
Kishana Y. Taylor | Ilechukwu Agu  | Ivy José | Sari Mäntynen | A.J. Campbell | Courtney Mattson |  Tsui-wen Chou | Bin Zhou | David Gresham | Elodie Ghedin |  Samuel L. Díaz Muñoz


If you are reading or using this, let us know how these data were useful for you. If you use these data and code, please cite the repository or the paper. Always open to collaborate! Please contact us!

### Quick Start
1. Make sure packages are installed (see #2 below) or use gbbseq-env.yml to set up Anaconda environment
2. git clone https://github.com/sociovirology/human_influenza_GbBSeq.git
3. chmod +x demultiplexing_infection_conditions.sh
4. ./demultiplexing_infection_conditions.sh (pairwise)
5. ./demultiplexing_infection_conditions.sh (controls)
6. Rscript aiv_detection_environment_analysis.R (or load interactively in R)
7. Rscript aiv_detection_environment_analysis.R (or load interactively in R)

### CONTENTS
1. Project Description
2. Packages and software used to test code
3. Data
4. Code

### 1. Project Description
Influenza A viruses can exchange genes when more than one strain co-infects a cell, causing novel strains that can spark pandemics andimpact public health. Thus, understanding genetic exchange potential between strains is an important goal. To understand this genetic exchange, we experimentally co-infected a diverse set of human influenza A viruses in the lab. To determine the genotype of many 1000's of progeny at high throughput, we devised the Genotype by Barcode Sequencing appraoch, which genotypes ~100 progeny in each of the 8 genome segments in each of 10 experimental coinfections (i.e. preserves genetic linkage). Using this sequence data we conducted a through evaluation of the potential for genetic exchange outcomes in a diverse set of influenza A viruses. This repository includes all data (or links to) and analyses in the paper.

Abstract:
RNA viruses can exchange genetic material during coinfection, an interaction that creates novel strains with implications for viral evolution and public health. Influenza A viral genetic exchange occurs when genome segments from distinct strains reassort in coinfection. Predicting potential reassortment between influenza strains is a longstanding goal. Experimental coinfection studies have shed light on factors that limit or promote reassortment. However, determining the reassortment potential between diverse Influenza A strains has remained an elusive goal. To fill this gap, we developed a high throughput genotyping approach to quantify reassortment among a diverse panel of human influenza virus strains, encompassing 41 years of epidemics, multiple geographic locations, and both circulating human subtypes A/H1N1 and A/H3N2. We found that the reassortment rate (proportion of reassortants) is an emergent property of a pair of strains where strain identity is a predictor of the reassortment rate. We show little evidence that antigenic subtype drives reassortment as intersubtype (H1N1xH3N2) and intrasubtype reassortment rates were, on average, similar. Instead, our data suggest that certain strains bias the reassortment rate up or down, independently of the coinfecting partner. We also observe that viral productivity is an emergent property of coinfections and that it is not correlated to reassortment rate, thus affecting the total number of reassortant progeny produced. Assortment of individual segments among progeny, and pairwise segment combinations within progeny, were not random and generally favored homologous combinations. This outcome was not related to strain similarity or shared subtype. Reassortment rate was closely correlated to both the proportion of unique genotypes and the proportion of progeny with heterologous pairwise segment combinations. We provide experimental evidence that viral genetic exchange is potentially an individual social trait subject to natural selection, which implies the propensity for reassortment is not evenly shared among strains. This study highlights the need for research incorporating diverse strains to discover the traits that shift the reassortment potential as we work towards the goal of predicting influenza virus evolution resulting from segment exchange.

### 2. Packages and software used to test code
This code was tested using the following software packages:

* cutadapt (2.6)
* PEAR (0.9.6)
* usearch (8.1.1861)
* R (3.6.3 (2020-02-29) -- "Holding the Windsock") with packages:
    + dplyr, ggplot2, tidyr, reshape2, readr, gridExtra, ggthemes

Anaconda environment file is available in gbbseq-env.yml

### 3. Data
Data consists of sequencing output from the illumina MiSeq platform, sample information, reference database, exprimental coinfection titers, and sample barcodes

1) Sequencing file is available in the Sequence Read Archive (Accession SRX7014890)

2) Information on experimental coinfections is in data/cross_data_runA.csv 

3) Barcode information is in shared/barcodes3.fasta, shared/barcodes3rev.fasta shared/barcodes5.fasta, and shared/barcodes5rev.fasta

4) Database of reference sequences for the amplicons from influenza viruses used in the experimental coinfections is in shared/reference_database.fasta

5) Titers of experimental coinfection supernatants are in data/supernatant_titers.csv (Not needed to run code)

### 4. Code
Below are descriptions of the code files used to generate the tables, figures, and statistics in the paper.

1) demultiplexing_infection_conditions.sh: This file is shell script that downloads raw sequencing reads, demultiplexes each read; and generates a flat text summary file used in downstream analyses 

2) pairwise_infections.R: This file is an R script that analyzes reassortment in pairwise experimental coinfections among 5 human influenza A virus strains

3) control_infections.R: This file is an R script that analyzes reassortment in control experimental coinfections to test the methods and the GbBSeq appraoch
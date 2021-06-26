# MOMIC: A Multi-Omics Pipeline for data analysis, integration and interpretation

## Table of contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Jupyter general instructions](#jupyter-general-instructions)
4. [Analysis Pipelines](#analysis-pipelines)
    1. [Transcriptomics. Genome Wide Expression Studies (GWES) from Microarray](#transcriptomics.-genome-wide-expression-studies-(gwes)-from-microarray)
    2. [Transcriptomics. Genome Wide Expression Studies (GWES) from RNASeq](#transcriptomics.-genome-wide-expression-studies-(gwes)-from-rnaseq)
    3. [Genome Wide Association Analysis (GWAS)](#genome-wide-association-analysis-(gwas))
    4. [Expression Proteomics](#prot)
    5. [Meta Analysis](#meta)
        1. [Meta GWES](#meta-gwes)
        2. [Meta GWAS](#meta-gwas)
    6. [Integrative Analysis](#integrative-analysis)
    7. [Enrichment](#enrichment)
5. [References](#references)

## Introduction

MOMIC offers a complete analysis environment for analysing and integrating multi-omics data in a single, easy-to-use platform.

MOMIC currently compiles protocols for whole genome SNP data (GWAS), mRNA expression (both from arrays and from RNAseq experiments) and protein data. Along with enrichment analysis and methods for combining heterogeneous data at different molecular levels. The proposed protocols are developed as Jupyter notebooks guiding the user through the tasks of pre-processing and transforming the data and performing the actual analysis, allowing the user to modify any piece of code needed along the process to adequate it to each project.
These protocols are based on accepted articles and best practices available in the literature.

The data analysis workflows starts with the pre-processing and quality control of the individual datasets. Each of the datasets are then independently analyzed following appropriate protocols. To consolidate these independent results for the same omic, a meta-analysis can be performed, generating single lists of SNPs or genes. After this, one can take a step forward and combine diverse omic results, along with functional genomics analysis, and try to elucidate the potential causative changes that lead to disease. With this integrative approach the user can get a single ranked list of candidate genes summarising the heterogeneous data at different molecular levels. Finally, visualization and pathway analysis tools are available to further explore analysis results. This whole procedure has been summarized in the next figure.

![main workflow](figures/general_types2.png)


MOMIC is presented as a collection of Jupyter notebooks using JupyterHub with JupyterLab deployed, written mainly in R language and containerized in Docker. It is distributed as a docker-compose project that contains the instructions needed to automatically create a fully working machine with JupyterHub, convenient extensions enabled, like git and table of content, docker volumes for data persistence, the pipeline source code and all the necessary libraries and third-party software. An alternative to a local installation is to use the pipeline hosted at [momic.us.es](momic.us.es). Note this alternative is intended for light analysis or quick testing. 

## Installation

The only requisite to install MOMIC locally is to have Docker and docker-compose already installed. Docker is a platform used to develop, deploy, and run applications with containers. Follow the instructions on each project website ([docker](https://docs.docker.com/install) and [docker-compose](https://docs.docker.com/compose/install)).

There are two git repositories you need to clone to replicate the project. One containing the docker instructions to build the service, called MOMIC_server, and another one containing the jupyter notebooks, called MOMIC_notebooks (**TODO** if using momic image this repo is not needed. This paragraph can be removed). 

**TODO** minimun requirements (RAM and disk)

### For the impatients:

- clone MOMIC_server in your local directory via `git clone https://github.com/laumadmar/MOMIC_server.git` 
- cd into that dir
- run `docker pull momic:latest` to download the docker image. You can ensure that the image is installed by using `docker images` (**No esta publicada, hasta entonces, desde el directorio donde tengamos la imagen guardada, hacemos docker load < momic_server_notebooks_image.tar.gz**)
- run `docker-compose up` and once the server is up, press `CTRL+c` to stop the console output
- run `docker-compose start` to keep the service running in the background
- access the tool at http://localhost:8000/jupyter and log in with user momic, pass m0m1c

An alternative is to create your container from the original instructions, which can be fully customised. Rename the file `Dockerfile.steps` to `Dockerfile` (**TODO** create Dockerfile.image and Dockerfile, the default will be the image)
- run `docker pull ubuntu:18.04` to download the ubuntu docker image
- run `docker-compose up` and once the server is up, press `CTRL+c` to stop the console output
- run `docker-compose start` to keep the service running in the background
- ssh into the container executing the access script (`./access`) and type `nohup Rscript /tmp/install_specific_libraries.R &` to install the required R packages
- access the tool at http://localhost:8000/jupyter and log in with user momic, pass m0m1c
- clone MOMIC Notebooks repo from jupyter or from the terminal after ssh into the container. For the former, go to the git tab located in the left menu, click on the button ‘Clone a Repository’ and provide the repo url. For the latter, ssh into the container, cd into momic home directory and type  `git clone https://github.com/laumadmar/MOMIC_notebooks.git`

If you run into any issues it is likely that the port 8000 is already in used or you need to configure any of the other parameters. Follow the step by step instrucctions in this case.



### For a step by step installation follow the next 3 major steps: / for the developers / move this to another page - remove from here - just keep the link.


#### 1. Get MOMIC server and go over the configuration

##### 1.1. Clone the repo

First, clone MOMIC_server in your local directory via `git clone https://github.com/laumadmar/MOMIC_server.git` and inspect the content:

- docker-compose.yml: YAML file defining services, networks and volumes

- Dockerfile: text document that contains the instructions to build the service

- jupyterhub_config.py: configuration file for JupyterHub

- README.md: readme file with instructions on how to deploy the multiomics pipeline

- software: directory containing third-party software


##### 1.2. Customise parameters

Modify if required the following parameters in docker-compose and/or Jupyter configuration file:

- docker-compose.yml:

```json
    ports:
        - "8000:8000"
    volumes:
        -./jupyterhub_config.py:/home/jupyterconfig/jupyterhub_config.py
        - jupyterdata:/mnt/data
        - home:/home
```


Container port 8000 is exposed in your local machine at port 8000, specified as “host port: container port”. Change it if this port is in used in your local machine.

Volumes are a mechanism for persisting data generated by and used by Docker containers. Three volumes are suggested here for:

    - jupyterhub config file – so there is no need to rebuild if you change this file
    
    - the directory where to keep the data to be analysed

    - jupyter home directory. This will contain the notebooks

Syntax for defining volumes is host_machine/absolute/path/to/dir:container/path/to/dir.

In order to have the same permissions in the data directory on host machine and the directory mounted as a volume in the docker container, it is advisable to create users in the container with the same name and uid as in the host machine. You can create users in the container using the Dockerfile or once the container is up and running; ssh in and create a linux user as usual.

- jupyterhub_config.py: The default parameters you can change if needed is: `c.JupyterHub.bind_url = 'http://8000/jupyter'` which sets the protocol, ip and base url on which the proxy will bind. By default, the JupyterLab view is loaded as indicated with the parameter: `c.Spawner.default_url = '/lab'`. Comment this line to launch the Classic Notebook. Note you can change from one to another later on.

#### 2. Build the service with compose

Follow these steps to build and run locally the multiomics pipeline with Docker Compose. Notice you need sudo privileges or a special group for running docker commands; read more on the Docker web site. 

1. From your project directory, start the server by running `docker-compose up`. Compose builds an image from the instructions specified in the Dockerfile, and starts the services defined.

2. Enter `http://localhost:8000/jupyter` in a browser to see the application running. Modify this url accordingly if you have changed the port and base url in docker-compose.yml.

3. Once the server is up, press `CTRL+c` to stop the console output. This will also stop the container.

4. Run `docker-compose start` to keep the service started in the background.

5. Note the RUN directive that installs R libraries within the Dockerfile. This takes very long to execute and it is commented out. As an alternative, install them after the build, accessing the container from the terminal and executing `nohup Rscript /tmp/install_specific_libraries.R &`

##### Access and log

There are two bash files on this repo for quickly accessing the container via terminal and checking the logs.

The access files contains: `docker exec -it jupyter_config_web_1 bash`

The logs file contains: `docker logs jupyter_config_web_1`

Change jupyter_config_web_1 by the name of your service. You can get if from `sudo docker-compose ps`.

##### Notes

If you modify at some point the Dockerfile, you need you build the image – follow:

`docker-compose stop`

`docker-compose build`

`docker-compose up`


First time you fire up the container and log in it can take a bit longer. It will show a message saying: “Your server is starting up. You will be redirected automatically when it's ready for you.” Refresh after a while if the home page does not come up.

Note you need sudo privileges or create a special group; read more on Docker web site. Few useful docker commands:

- `docker-compose stop` to stop the running container

- `docker-compose ps` to check the status

- `docker inspect --format='{{.LogPath}}' momic_server_web_1` to get the path to the log file

- `docker logs momic_server_web _1` to print the log in console

- `docker exec –it momic_server_web _1 bash` to access the running container.

#### 3. Get MOMIC Notebooks

Log into MOMIC Jupyterhub and go to the git tab located in the left menu. Click on the button ‘Clone a Repository’ and provide the url https://github.com/laumadmar/MOMIC_notebooks.git

An alternative to use the gitlab extension is to clone the repository from the terminal, either using Jupyter terminal or via ssh into the container. CD into your jupyter home directory and type `git clone https://github.com/laumadmar/MOMIC_notebooks.git`.

You now have in your home directory a copy of all the notebooks necessary to carry out the analysis presented in MOMIC.

## Jupyter general instrucctions

The Jupyter Notebook is an open-source web application that allows you to create and share documents that contain live code, equations, visualizations and narrative text. JupyterHub brings the power of notebooks to groups of users.

You will find all the notebooks needed to complete the different analysis in your home directory (first screen right after login).

Notebooks are provided as read-only, we call them templates. To create your own notebook from a template, select the desired template and click on the Duplicate button. Select the new notebook and right-click to see the Rename button. An alternative if using the classic theme is to click on the name of the new notebook to open it and change the name once opened.

It is possible to create a fresh empty notebook in whichever directory you want from the plus button located at the top left menu (or New button if using the classic notebook). If you want to replicate an analysis, follow the same steps indicated in the corresponding template. Note that a cell can have different types, we are working here with markdown and code cells. Thi can be changed in the top menu.

In the same way, you can create a new folder, a text file and open a terminal window.

To execute a cell, press `CTRL+ENTER` or click on the Run button located in the top menu. Alternatively, `CTRL+ALT` to create a new empty code cell bellow it. Inspect the menu to go through all Jupyter features or visit the Jupyter Project Documentation website to know more (https://jupyter-notebook.readthedocs.io/en/stable/ and https://jupyterlab.readthedocs.io/en/stable/)

## Analysis pipelines


Notebooks are provided as read-only and we refer to these as templates. As explained in the previous [section](#jupyter-general-instrucctions), create your own notebook, a duplicate or an empty one, and modify paths and/or code according to your needs. Execute code in cells with CTRL+ENTER or clicking on the Run button located at the top menu. Alternatively, CTRL+ALT to create a new empty code cell below it.

Read carefully the comments on the Jupyter templates as they explain in detail the code and protocol used to perform the different steps of each analysis.

In order to modify a core pipeline function, duplicate the original script, rename it, do the desired changes and save it. Note that the cell that imports these functions in your template, the one that contains this code `source("scripts/whatever.R")`, has to be run after the changes.

### Transcriptomics. Genome Wide Expression Studies (GWES) from Microarray

RNA microarrays, are tools that allow the identification and quantification of the mRNA transcripts present in the cells. RNA microarrays can simultaneously measure the expression level of thousands of genes within a particular mRNA sample. Such high-throughput expression profiling can be used to compare the level of gene transcription in clinical or biological conditions in order to find differences in expression levels between predefined groups of samples. This is called differential expression (DE) analysis.

This pipeline starts from raw expression data or an expression matrix, and ends with a set of differentially expressed genes annotated using Entrez and Human Genome Nomenclature Committee (HGNC) offical gene symbols.

![GWES workflow](figures/Microarray_workflow.png)

The Jupyter templates to follow for completing a Microarray analysis can be found under the home directory, MOMIC_notebooks/GWES/Microarray. Create your own notebook as explained in [section 3](#jupyter-general-instrucctions) and follow the protocol detailed.

Before running a DE analysis, it is critical to perform a thoroughly data cleaning and pre-processing. It is also essential that expression and clinical datasets contain the same samples and in the same order. Pre-processing steps involve transforming raw data into a matrix of normalized intensities, using the RMA, limma or preprocessCore R packages depending on the data platform. Bath effect removal is also available through the combat function of the sva package.

For the DE analysis, we use R package limma [ref] as part of a customised expression pipeline. Note the use of `source("/home/guess/scripts/diffExpressionPipeline.R")` to import the custom script. This performs a DE analysis over the dataset supplied and can be adjusted to your needs via the input parameters. Specify the contrast or conditions you want to compare, i.e., case versus control or treatment1 versus treatment2. It is possible to perform the DE on a subset of the data, i.e., only in male population. You can also provide the variables you want to include in the model to adjust for these covariates, excluding their effect. The details for the parameters of this function can be found in the template itself.

Two notebook templates are provided, called compact and step_by_step. Both implement the same functionality but the last one does it in a more detailed and comprehensive way but more tedious to programmatically execute the analysis. In addition to the dataset analysed in the template, another one is done in the same way to be able to perform a posterior meta-analysis (explained in section **x**).

Regarding gene annotation, if there is more than one gene matching the same probe, deduplication can be performed during the pre-processing step (as illustrated in the Affymetrix rma protocol) or after DE, keeping only the probe with the lowest p value, as gene representative towards integrative analysis.

![GWES protocol](figures/MicroarrayProtocolTable.png)

To illustrate this analysis, datasets GSE48350 and GSE15222, from the Gene Expression Omnibus (GEO) database have been used.

### Transcriptomics. Genome Wide Expression Studies (GWES) from RNASeq

RNA-Seq is a particular technology-based sequencing technique which uses next generation sequencing (NGS) to reveal the presence and quantity of RNA in a biological sample at a given moment. Gene expression is quantified by counting the number of reads that mapped to each locus in the transcriptome assembly step.

This pipeline starts from raw sequence reads, and ends with a set of differentially expressed genes.

![RNASeq workflow](figures/RNASeq_workflow.png)

The Jupyter templates to follow for completing a RNASeq analysis can be found under the home directory, MOMIC_notebooks/GWES/RNASeq. Create your own notebook as explained in [section 3](#jupyter-general-instrucctions) and follow the protocol detailed.

An important step is to look at the quality of the raw reads before proceeding with the analysis. This can be done using fastqc java program. The sensible step after this is to remove sequences with low quality to get better alignment in the later steps. To determine where on the human genome our reads originated from, these are aligned to the reference genome using STAR. Read quantification is performed at the same time with STAR.

Last step is DE analysis, done with DESeq2, to find differences in expression levels between predefined groups of samples. Two notebook templates are provided, called compact and step_by_step. Both implement the same functionality but the last one does it in a more detailed and comprehensive way but more tedious to programmatically execute the analysis.

![RNASeq protocol](figures/RNASeqProtocolTable.png)

To illustrate this analysis synthetic data have been generated.

### Genome Wide Association Analysis (GWAS)

Genome Wide Association Studies are hypothesis free methods to identify associations between genetic regions (loci) and traits (including diseases). It has long been known that genetic variation between individuals can cause differences in phenotypes. These causal variants, and those which are tightly linked to their region of the chromosome, are therefore present at higher frequency in cases (individuals with the trait) than controls (individuals without the trait).

The GWAS protocol is available for genome builds GRCh37/hg19 and GRCh38/hg38, including the UCSC liftOver tools for updating genetic positions.

This pipeline starts from genomic data in PLINK format, and ends with SNPs annotated for rs identifiers (dbSNP v150 for the GRCh37 pipeline and dbSNP v151 for the GRCh38 pipeline) and HGNC annotated genes.

![GWAS workflow](figures/GWAS_workflow.png)

A series of templates are presented to complete a whole GWAS study following the guidelines on Anderson et al. protocol. These are under the GWAS folder in the home directory and are provided in read-only mode. Create your own notebook as explained in [section 3](#jupyter-general-instrucctions) and follow the protocol detailed.

An initial standard quality control needs to be performed; it consists on:

1. Exclusion of individuals: a) with more than 3% missing genotypes, b) with excess autosomal heterozygosity (>0.35 or more than ± 3 standard deviations from the mean), c) showing a discrepancy between genotypic and reported sex or d) showing non-European ancestry. Duplicated and related individuals are also removed by means of IBS estimates (pi-hat >0.1875) within and across studies.

2. Exclusion of SNPs: a) with missing genotype rate > 5%, b) with significant missingness test between cases and controls (p<10−6), c) not in Hardy-Weinberg equilibrium (p<10-6 in controls) and d) SNPs with minor allele frequency (MAF) < 1%.

After the QC, genotype imputation of build GRCh37/hg19 can be performed at the Michigan Imputation Server, using the minimac 3 algorithm, the HRC reference panel and the SHAPEIT tool for haplotype phasing. For GRCh38/hg38 data, imputation is performed at the TopMed server, which uses the minimac4 algorithm and Eagle v2.4 for phasing, with the TopMed population as reference. After imputation, SNPs with an R2 quality estimate lower than 0.3 are excluded from further analyses according to the software recommendations.

Following the QC and optional imputation, a case control association study can be performed with PLINK in order to identify genetic variants or SNPs that can be associated with a trait. Gene-wise statistics can be then computed using MAGMA software, which takes into account physical distance and linkage disequilibrium (LD) between markers to aggregate this data to the level of whole genes.

Finally, results are inspected using Manhattan and QQ plots.

The two main protocols followed in this pipeline are: Anderson, et. al. Data quality control in genetic case‐control association studies [ref] and Marees AT, et. al. A tutorial on conducting genome-wide association studies. It is very recommended that you read first these two papers if you had little knowledge about this kind of analysis.

![GWAS protocol](figures/GWASProtocol.png)

The data used to illustrate this analysis is taken from the 1000 Genomes Project. We have slightly modified this dataset updating the identifiers to rsID, removing SNPs with minor allele frequency and adding fake case/control phenotypes. **We have also purposely duplicated one individual to check that the algorithm detects it.**

### Proteomics

Expression Proteomics pipeline involves includes the analysis of differentially expressed proteins between conditions, such as diseased vs. healthy tissue.

This pipeline starts from the intensity matrix, as the quantitation and identification of proteins are not cover here, and ends with a set of differentially expressed proteins.

![GWAS workflow](figures/Proteomics_workflow.png)

The Jupyter templates to follow for completing a Differential Analysis can be found under the home directory, MOMIC_notebooks/Proteomics. Create your own notebook as explained in [section 3](#jupyter-general-instrucctions) and follow the protocol detailed.

Before running a differential protein expression analysis, it is very important to do data cleaning and pre-processing. It is also essential that the intensity matrix and metadata contain the same samples and in the same order. This pre-processing consists on removing decoy matches and matches to contaminant, extracting the LFQ intensities columns, filtering on missing values, log transformation, normalisation, unique peptide counting and visualisation.

For the differential analysis, we use DEqMS R package. DEqMS builds on top of limma and improves it with proteomics data specific properties, accounting for variance dependence on the number of quantified peptides or PSMs for statistical testing of differential protein expression.

![Proteomics protocol](figures/ProteomicsProtocolTable.png)

Proteomics data from postmortem brain tissue have been collected from The National Institute on Aging’s Baltimore Longitudinal Study of Aging (BLSA) (Synapse 10.7303/syn3606086).

### Meta-Analysis

Meta-analysis is the statistical procedure for synthesising data across studies. It can be performed when there are multiple studies addressing the same question and the same molecular level. This analysis is to be conducted after the individual analysis have been completed. Two different protocols are provided here, one for combining multiple GWAS results and another for transcriptomics and proteomics results.

### Meta-Analysis of gene expression data

This one-step pipeline starts from the results of the DE analysis of various studies, and ends with a unified list of differentially express genes or proteins across studies – with its correspondent average effect sizes and p values.

The MetaDE R library is used to perform this analysis. The Random Effects Model (REM) algorithm takes as input the p values, observed effect size (logFC values) and observed variance to compute summary logFCs and associated statistics. . The Variance is calculated as Standard Error (SE)^2, and SE is calculated as the difference of confidence intervals divided by 3.92. Additionally, a Fisher exact test on p-values can also be performed. The original MetaDE algorithm have been slightly modified to account for genes that are not present in all input studies.

The templates for this meta-analysis can be found in the MetaAnalysis/GWES folder under MOMIC_notebooks within the home directory. This portocol has been illustrated using microarray data, but RNASeq and proteomics can be done in the same way. Note the cell that imports the core pipeline functions `(source ("/home/guess/scripts/metaDE.R")`, it is very important that you run this cell first.

The example notebook provided combines the two limma tables obtained from the DE analysis performed following [section 4.A GWES](#transcriptomics.-genome-wide-expression-studies-(gwes)-from-microarray). The core pipeline function used in this template, `prepare_matrix_function`, is ready to take as input limma datasets containing logFC, confidence intervals and p value. Any output from a different DE algorithm can be used as long as you prepare the datasets accordingly to run the MetaDE algorithm.

### Meta-analysis of GWAS data

This one-step pipeline starts from SNP level GWAS results and ends with a single list of  SNPs and associated summary statistics.

The analysis is peformed using METAL, which allows two analyses scehemes: SAMPLESIZE (default approach, uses p-value and direction of effect, weighted according to sample size) and STDERR (the classical approach, uses effect size estimates and standard errors) A Cochran’s Q test and I2 statistics are generated to evaluate the potential effect of study heterogeneity on the results.

The template for this meta-analysis can be found in MetaAnalysis/GWAS under MOMIC_notebooks within the home directory. The template is provided as read-only. As explained in [section 3](#jupyter-general-instrucctions), create your own notebook (a duplicate or an empty one) and modify paths and/or code according to your needs.

## Integrative analysis

Integrative Analysis aims to consolidate heterogeneous data at different omics levels to understand their interrelation and combined influence on the disease processes.

This one-step pipeline starts from the results of the meta-analysis for different studies, and ends with a unified ranked list of differentially expressed genes.

The integration is performed using Robust Rank Aggregation (RRA) method (Kolde R et al., 2012). It detects genes that are ranked consistently better than expected under null hypothesis of uncorrelated inputs and assigns a significance score for each gene.

For each item, the algorithm looks at how the item is positioned in the ranked lists and compares this to the baseline case where all the preference lists are randomly shuffled. As a result, it assigns a P-value for all items, showing how much better it is positioned in the ranked lists than expected by chance. This P-value is used both for re-ranking the items and deciding their significance.

Since the number of informative ranks is not known, RRA defines the final score for the rank vector r as the minimum of P-values and orders all rank vectors according to their ρ scores. The obtained gene list is ranked according to this global score, using the Rank library from the R Basic package. The same rank is assigned to those genes with NA score.

The template for this analysis can be found in the folder IntegrativeAnalysis under the home directory. It combines the results from the meta-analysis of microarray expression data obtained following [section 4.A GWES](#transcriptomics.-genome-wide-expression-studies-(gwes)-from-microarray), with the GWAS analysis results obtained following [section 4.C GWAS](#genome-wide-association-analysis-(gwas)). 

## Enrichment analysis

Enrichment analysis, or pathway analysis, can identify terms which are statistically over or under-represented within the list of interest, by systematically mapping genes and proteins to their associated biological annotations, such as gene ontology GO terms or pathway membership, and then comparing the distribution of the terms within a gene set of interest with the background distribution of these terms (e.g., all genes represented on a microarray chip).

This one-step pipeline starts from the list of genes to explore, and ends with a list of categories and its associated statistics along with visualization plots.

The enrichment is performed using WebGestalt in R and the visualization  with GOplot package and pheatmap from R CRAN.

Various templates can be found in the folder Enrichment under the home directory. In this case the enrichment is done from the results of the meta-analysis case-control completed following [section 4.E.a Meta GWES]. Read carefully the templates as these contain detailed explanation of the code implemented.

- enrichment_GOplots_template.ipynb: this template performs the enrichment for all predefined categories and prints the GO plots.

- enrichment_GOplots_loop_template.ipynb: this template does the same as the previous one but includes a loop through more than one dataframe to perform the enrichment on.

- pheatmap_template.ipynb: this template uses the R library pheatmap to represent combined results from meta expression, meta GWAS, integrative analysis and enrichment of meta expression, completed following the different analysis available in this pipeline.

Templates are provided as read-only. As explained in [section 3](#jupyter-general-instrucctions), create your own notebook (a duplicate or an empty one) and modify paths and/or code according to your needs. Execute code in cells with `CTRL+ENTER` or doing click on the Run button located at the top menu. Alternatively, `CTRL+ALT` to create a new empty code cell bellow it.

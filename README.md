# ComplexDiseases_Pipeline

Pipeline developed during my master thesis.

I used bash to filtrate and perform quality controls on the original case dataset.

Python was used to parse and merge the VCF files, process the datasets, select genes of interest and train and test Machine Learning models.

Also, R was used for the integration of the genes in the Protein-Protein Interactions network and to create different plots.

Most of files were construted to be used from the terminal (see Bash scripts).

More information in the pdf of my master thesis.

### Organization of the files:

.\
├── data\
│   ├── datasets\
│   │   ├── cleaned_dataset.csv.gz\
│   │   ├── dataset_with_genes.csv.gz\
│   │   ├── imputed_dataset.csv.gz\
│   │   ├── merged_dataset.csv.gz\
│   │   ├── Original\
│   │   │   ├── ensembl_v37.db\
│   │   │   └── ensembl_v37.gtf\
│   │   ├── reduced_dataset_network.csv.gz\
│   │   ├── reduced_dataset_pval.csv.gz\
│   │   ├── reduced_dataset_risk.csv.gz\
│   │   ├── top_dataset_network.csv.gz\
│   │   └── top_dataset_pval.csv.gz\
│   ├── features\
│   │   ├── top100_network.csv\
│   │   └── top100_pval.csv\
│   ├── figures\
│   │   ├── all_network\
│   │   ├── all_pval\
│   │   ├── network\
│   │   ├── pval\
│   │   ├── risk\
│   ├── genes\
│   │   ├── geneList.csv\
│   │   ├── network.csv\
│   │   ├── pval.csv\
│   │   ├── risk.csv\
│   │   ├── riskGenes.csv\
│   │   ├── top100_network.csv\
│   │   └── top100_pval.csv\
│   ├── proteins\
│   │   ├── genes.csv\
│   │   ├── interactions_simple.csv\
│   │   ├── r_genes.csv\
│   ├── variants\
│   │   ├── notnormal.csv\
│   │   └── sigVars.csv\
│   └── vcf\
│       ├── cases\
│       │   └── vcfQC\
│       │       ├── All_PT_biallelic.vcf.gz\
│       │       ├── All_PT.vcf.gz\
│       │       ├── count_snp_indel.sh\
│       │       ├── final_dataset_hwe.vcf.gz\
│       │       ├── final_dataset_ind.vcf.gz\
│       │       ├── final_dataset_q20.vcf.gz\
│       │       ├── pics.R\
│       │       ├── relatedness.py\
│       ├── controls\
│       ├── igsr_samples.tsv\
└── src\
&emsp;&ensp;├── bash\
&emsp;│   ├── plots.sh\
&emsp;│   ├── runCases.sh\
&emsp;│   ├── runControls.sh\
&emsp;│   ├── runMain.sh\
&emsp;│   └── runMerge.sh\
    ├── jupyter\
    │   └── Roc curves and tables.ipynb\
    ├── python\
    │   ├── biogridParser.py\
    │   ├── featureExtraction.py\
    │   ├── feature_selector.py\
    │   ├── geneSelection.py\
    │   ├── main.py\
    │   ├── models.py\
    │   ├── network.py\
    │   ├── pickle\
    │   ├── plots.py\
    │   ├── variantSelection.py\
    │   └── vcfParser.py\
    └── R\
        ├── all_data.RData\
        ├── dmGWAS_3.0.tar.gz\
        ├── network.R\
        ├── plotManhattan.R\
        └── qualityControl.R\

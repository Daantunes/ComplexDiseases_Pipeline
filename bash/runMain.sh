# Run main
python3 main.py -d ../../data/datasets/merged_dataset.csv.gz -cV 10
python3 main.py -d ../../data/datasets/cleaned_dataset.csv.gz -iV most_frequent
python3 main.py -w pickle/dataImp.p -t normality
python3 main.py -w pickle/dataImp.p -t chi2
python3 main.py -w pickle/dataImp.p -tG ../../data/variants/sigVars.csv
python3 main.py -w pickle/dataImp.p -aG ../../data/genes/geneList.csv
python3 main.py -w pickle/dataGenes.p -nF pval
python3 main.py -w pickle/dataGenes.p -nF risk
python3 main.py -w pickle/dataGenes.p -nF network
python3 main.py -d ../../data/datasets/reduced_dataset_pval.csv.gz -tF 25
python3 main.py -d ../../data/datasets/reduced_dataset_network.csv.gz -tF 25
python3 main.py -n ../../data/proteins/interactions_simple.csv
python3 main.py -n ../../data/proteins/interactions_simple.csv ../R/all_data.RData

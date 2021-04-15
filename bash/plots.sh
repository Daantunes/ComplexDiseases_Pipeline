# For each dataset  run model
python3 main.py -d ../../data/datasets/top_dataset_network.csv.gz -c svm -g C:0.25 degree:1 gamma:25 kernel:linear tol:0.001
python3 main.py -d ../../data/datasets/top_dataset_network.csv.gz -c tree -g criterion:entropy max_leaf_nodes:50 min_samples_leaf:3 min_samples_split:5 n_estimators:50
python3 main.py -d ../../data/datasets/top_dataset_network.csv.gz -c log -g C:0.4 penalty:l1 solver:liblinear tol:0.001

python3 main.py -d ../../data/datasets/reduced_dataset_risk.csv.gz -c svm -g C:0.25 degree:1 gamma:25 kernel:linear tol:0.001
python3 main.py -d ../../data/datasets/reduced_dataset_risk.csv.gz -c tree -g criterion:entropy max_leaf_nodes:50 min_samples_leaf:3 min_samples_split:5 n_estimators:50
python3 main.py -d ../../data/datasets/reduced_dataset_risk.csv.gz -c log -g C:0.4 penalty:l1 solver:liblinear tol:0.001

python3 main.py -d ../../data/datasets/top_dataset_pval.csv.gz -c svm -g C:0.25 degree:1 gamma:25 kernel:linear tol:0.001
python3 main.py -d ../../data/datasets/top_dataset_pval.csv.gz -c tree -g criterion:entropy max_leaf_nodes:50 min_samples_leaf:3 min_samples_split:5 n_estimators:50
python3 main.py -d ../../data/datasets/top_dataset_pval.csv.gz -c log -g C:0.4 penalty:l1 solver:liblinear tol:0.001

python3 main.py -d ../../data/datasets/reduced_dataset_pval.csv.gz -c svm -g C:0.25 degree:1 gamma:25 kernel:linear tol:0.001
python3 main.py -d ../../data/datasets/reduced_dataset_pval.csv.gz -c tree -g criterion:entropy max_leaf_nodes:50 min_samples_leaf:3 min_samples_split:5 n_estimators:50
python3 main.py -d ../../data/datasets/reduced_dataset_pval.csv.gz -c log -g C:0.4 penalty:l1 solver:liblinear tol:0.001

python3 main.py -d ../../data/datasets/reduced_dataset_network.csv.gz -c svm -g C:0.25 degree:1 gamma:25 kernel:linear tol:0.001
python3 main.py -d ../../data/datasets/reduced_dataset_network.csv.gz -c tree -g criterion:entropy max_leaf_nodes:50 min_samples_leaf:3 min_samples_split:5 n_estimators:50
python3 main.py -d ../../data/datasets/reduced_dataset_network.csv.gz -c log -g C:0.4 penalty:l1 solver:liblinear tol:0.001

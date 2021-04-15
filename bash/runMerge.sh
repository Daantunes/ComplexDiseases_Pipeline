# Merge files
CHR=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X')
for i in "${CHR[@]}";
do python3 main.py -m $i ../../data/vcf/cases/outputPandas_$i.csv.gz ../../data/vcf/controls/outputPandas_$i.csv.gz
done

python3 main.py -m all

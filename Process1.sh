# tail -n +K to output lines starting with the Kth 
mkdir Process1
cat OriginData/NodesForSmallNetworks.txt | cut -f1 | tail -n +2 | sort -u > Process1/Genes-of-Networks.txt
head -1 OriginData/Profiles-ANan-DEGs.csv > Process1/Profiles-DEGs-nodes-networks.csv
grep -f Process1/Genes-of-Networks.txt OriginData/Profiles-ANan-DEGs.csv >> Process1/Profiles-DEGs-nodes-networks.csv


head -1 OriginData/kat-rpkm-expression.csv > Process1/Profiles-nodes-networks.csv
grep -f Process1/Genes-of-Networks.txt OriginData/kat-rpkm-expression.csv >> Process1/Profiles-nodes-networks.csv


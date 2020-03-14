```
for i in `ls *vafs`
do
awk '{print $1,$2,$2+1,$10}' OFS="\t" $i > $i.map
bedtools intersect -wao -a $i.map -b ~/hpc/db/hg19/refGene.hg19.bed | awk '{print $1":"$2,$10,$4}' OFS="\t" | sort -u  > $i.vaf
done

perl maf2matrix.pl > vafmatrix.txt
Rscript fviz_cluster.R

BiocManager::install("MutationalPatterns")

```




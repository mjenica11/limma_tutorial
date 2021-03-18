# Filter the gencode annotation file to get just the sex chr genes
# Will use this to subset the GTEx count data in filter_count_matrix.R 
# V26 Gencode annotation downloaded from https://www.gencodegenes.org/human/release_26.html
cat gencode.v26.annotation.gtf | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$3,$5,$7,$10,$16}}' | tr -d '";' > gencodeGenes.txt
grep "chrY" gencodeGenes.txt > gencodeGenes_Ychr.txt
grep "chrX" gencodeGenes.txt > gencodeGenes_Xchr.txt

#Brain: 26 species; 1617 samples; 21782 genes.
#Cerebellum: X species; X samples; X genes.
#Liver: 23 species; 431 samples; 19222 genes.
#Kidney: 21 species; 180 samples; 19180 genes.
#Heart: 17 species; 375 samples; 18471 genes.
#Testis: 15 species; 233 samples; 18975 genes.
#Spleen: 13 species; 153 samples; 17889 genes.
#Muscle tissue: 12 species; 92 samples; 17284 genes.
#Lung: 9 species; 105 samples; 17345 genes.
for tissue in "brain" "cerebellum" "liver" "kidney" "heart" "testis" "spleen" "muscle tissue" "lung"
do
  echo $tissue
  Rscript 1_retrieve_Bgee_data.R "${tissue}" all-sex post-juvenile
done
Rscript 3_get_gene_expression_per_orthogroup.R 10

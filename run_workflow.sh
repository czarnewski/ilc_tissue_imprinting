#! /bin/bash -l

##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
eval "$(conda shell.bash hook)"
# conda activate base
# conda install -c conda-forge mamba
# mamba env create -n ilc_tissue_imprinting -f sauron_environment_20201209.yml
# conda activate ilc_tissue_imprinting




########################
### DEFINE VARIABLES ###
########################
var_to_plot='Dataset,Tissue,Plate,Donor,Celltype,Sex'
var_to_regress='perc_mito,perc_rpl,perc_rps,perc_Chr_Y,perc_snRNA,Sex'
main=$(pwd)
cd $main
! test -d analysis && $(mkdir analysis)
! test -d log && $(mkdir log)



#################################
### DOWNLOAD CODE FROM GITHUB ###
#################################
! test -d 'sauron' && {
  $(git clone https://github.com/NBISweden/sauron.git --branch v.0.1-beta)
}
script_path=$main/'sauron/scripts'



##############################
### DOWNLOAD DATA FROM GEO ###
##############################
! test -d data && {
  mkdir -p $main/'data/RNAseq'
  mkdir -p $main/'data/VDJseq'

  #Download metadata file
  cd $main/'data'
  curl -o GSE150050_metadata.csv.gz 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150050/suppl/GSE150050%5Fmetadata%2Ecsv%2Egz'
  gunzip GSE150050_metadata.csv.gz

  #Download raw data from GEO
  cd $main/'data/RNAseq'
  curl -o GSE150050_RAW.tar 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150050/suppl/GSE150050%5FRAW%2Etar'
  tar -xvf GSE150050_RAW.tar

  mv $main/'data/RNAseq'/*VDJ* $main/'data/VDJseq'/
  for i in $(ls $main/'data/RNAseq'/*.tar.gz); do tar -xzvf $i; done

  cd $main/'data/VDJseq'
  for i in $(ls $main/'data/VDJseq'/*.tar.gz); do tar -xzvf $i; done

  rm $main/'data/RNAseq'/*.gz
  rm $main/'data/RNAseq'/*.tar
  rm $main/'data/VDJseq'/*.gz
  rm $main/'data/VDJseq'/*.tar

  cd ..
}




#####################
### LOAD DATASETS ###
#####################
Rscript $script_path/00_load_data.R \
--input_path $main/'data/RNAseq' \
--dataset_metadata_path $main/'data/GSE150050_metadata.csv' \
--species_use 'hsapiens' \
--estimate_molecules_from_read_count 'True' \
--mapping_threshold '2.5' \
--sum_to_gene_level 'True' \
--assay 'RNA' \
--output_path $main/'analysis/1_qc' \
2>&1 | tee $main/log/'00_load_data_log.txt'



###########################
### RUN QUALITY CONTROL ###
###########################
Rscript $script_path/01_qc_filter.R \
--Seurat_object_path $main/'analysis/1_qc/raw_seurat_object.rds' \
--columns_metadata $var_to_plot \
--species_use 'hsapiens' \
--remove_non_coding 'False' \
--plot_gene_family 'RPS,RPL,MT-,HB,MT[RS]' \
--remove_gene_family 'mito,MT[RN]' \
--min_gene_count '5' \
--min_gene_per_cell '200' \
--normalization 'LogNormalize' \
--assay 'RNA' \
--output_path $main/analysis/1_qc \
2>&1 | tee $main/log/'01_QC_log.txt'



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################
Rscript $script_path/02_integrate.R \
--Seurat_object_path $main/'analysis/1_qc/filt_seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--var_genes 'scran,0.001' \
--integration_method 'cca,Sex' \
--cluster_use 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/2_clustering' \
2>&1 | tee $main/log/'02_integrate_log.txt'



####################################
### RUN DIMENSIONALITY REDUCTION ###
####################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'top,50' \
--var_genes 'scran,0.001' \
--dim_reduct_use 'umap' \
--k_nearest_neighbour '10' \
--pre_dim_reduct 'pca' \
--cluster_use 'Celltype,ILC1,ILC2,ILC3,NK,T' \
--cluster_method 'louvain,hc,leiden' \
--output_path $main/'analysis/2_clustering' \
--assay 'cca' \
2>&1 | tee $main/log/'03_dr_and_cluster_log.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
--clustering_use 'HC_21' \
--metadata_use 'none' \
--exclude_cluster 'NONE' \
--DGE_method "wilcox" \
--only_positive "true" \
--covariates "NULL" \
--max_cells_per_ident "200" \
--pvalue_threshold "0.01" \
--logfc_threshold "0.25" \
--assay 'RNA' \
--output_path $main/'analysis/2_clustering/3_diff_expr' \
2>&1 | tee $main/'log/4_diff_expr_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
--marker_lists $script_path/'../support_files/cell_markers/main_cell_types.csv' \
--assay 'RNA' \
--output_path $main/'analysis/2_clustering/cell_type_prediction' \
2>&1 | tee $main/'log_sauronv1/cell_type_prediction_log.txt'



# #Clusters 12, 16 and 20 were associated to Fibroblasts and APCs and were removed
# Hierachical Clustering was used because the total number of clusters is easy to find


####################################
### RUN DIMENSIONALITY REDUCTION ###
####################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'top,50' \
--var_genes 'scran,0.001' \
--dim_reduct_use 'umap' \
--k_nearest_neighbour '10' \
--pre_dim_reduct 'pca' \
--cluster_use 'HC_21,1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19,21' \
--cluster_method 'louvain,hc,leiden' \
--output_path $main/'analysis/2_clustering_filt' \
--assay 'cca' \
2>&1 | tee $main/log/'03_dr_and_cluster_log.txt'



########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
Rscript $script_path/'05_cluster_correlation.R' \
--Seurat_object_path $main/'analysis/2_clustering_filt/seurat_object.rds' \
--clustering_use 'HC_7' \
--exclude_cluster 'NONE' \
--merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
--output_path $main/'analysis/2_clustering_filt/cluster_correlations' \
2>&1 | tee $main/'log/4_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/2_clustering_filt/seurat_object.rds' \
--clustering_use 'HC_14' \
--metadata_use 'none' \
--exclude_cluster 'NONE' \
--DGE_method "wilcox" \
--only_positive "true" \
--covariates "NULL" \
--max_cells_per_ident "100" \
--pvalue_threshold "0.1" \
--logfc_threshold "0.1" \
--assay 'RNA' \
--output_path $main/'analysis/2_clustering_filt/3_diff_expr' \
2>&1 | tee $main/'log/4_diff_expr_log.txt'



########################
### RUN VDJ ANALYSIS ###
########################
Rscript $script_path/VDJ_analysis.R \
--Seurat_object_path $main/'analysis/2_clustering_filt/seurat_object.rds' \
--VDJ_annotation_path $main/'data/VDJseq' \
--columns_metadata 'Celltype' \
--top_TCRs '10' \
--paired_only 'false' \
--only_coding_cdr3 'true' \
--same_scale 'true' \
--chains_use 'TRA,TRB,TRD,TRG' \
--assay 'RNA' \
--output_path $main/'analysis/2_clustering_filt/5_vdj_analysis_ig' \
2>&1 | tee $main/log_sauronv1/'05_vdj_analysis_ig_log.txt'




conda deactivate

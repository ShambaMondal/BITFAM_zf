#' BITFAM main function. BITFAM will infer the transcription factor activities from single cell RNA-seq data based on the ChIP-seq data
#'
#' @param data A matrix or dataframe, normalized single cell RNA-seq data
#' @param species mouse or human
#' @param interseted_TF Transcription factors of interests
#' @param scATAC_obj A preprocessed Seurat object of scATAC-seq data
#' @param number of CPU cores
#' @param number of max iteration
#' @param convergence tolerance on the relative norm of the objective
#' @return sampling results of TF inferred activities and TF-gene weights
#' @export
#' @import rstan
#' @import Seurat

BITFAM <- function(data, species, interseted_TF = NA, scATAC_obj = NA,ncores, iter = 8000, tol_rel_obj=0.005){
  if(species == "mouse"){
    TF_targets_dir <- "mouse/"
  }else if(species == "human"){
    TF_targets_dir <- "human/"
  }else if(species == "zebrafish"){
    TF_targets_dir <- "zebrafish/"
  }else{
    stop("The species must be either mouse or human or zebrafish.")
  }

  if(dim(data)[1] > 5000){
    variable_genes <- Seurat::FindVariableFeatures(data)
    variable_genes <- variable_genes[which(x = variable_genes[, 1, drop = TRUE] != 0), ]
    variable_genes <- variable_genes[order(variable_genes$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
    variable_genes <- head(x = rownames(x = variable_genes), n = 5000)
    data <- data[variable_genes, ]
  }
  
  All_TFs <-system.file("extdata", paste0(TF_targets_dir, "all_TFs.txt"), package = "BITFAM")
  #print(paste0("1. ALL_TFs: \n", ALL_TFs))
  print("1. ALL_TFs: ")
  print(All_TFs)
  All_TFs <- read.table(All_TFs, stringsAsFactors = F)$V1
  print("2. ALL_TFs: ")
  print(All_TFs)
  #print(paste0("2. ALL_TFs: \n", ALL_TFs))
  TF_used <- rownames(data)[rownames(data) %in% All_TFs]
  print("3. TF_used: ")
  print(TF_used)
  #print(paste0("3. TF_used: \n", TF_used))
  rownames(data) <- toupper(rownames(data))
  print("4. rownames(data): ")
  print(rownames(data))
  #print(paste0("4. rownames(data): \n", rownames(data)))
  if(is.na(interseted_TF)){
  }else{
    TF_used <- unique(c(TF_used, interseted_TF))
  }
  print("5. TF_used: ")
  print(TF_used)
  #print(paste0("5. TF_used: \n", TF_used))
  gene_list <- list()
  for(i in TF_used){
    TF_targets_path <-system.file("extdata", paste0(TF_targets_dir, i), package = "BITFAM")
    print("6. i: ")
    print(i)
    print("7. TF_targets_path: ")
    print(TF_targets_path)
    #print(paste0("6. i: \n", i))
    #print(paste0("7. TF_targets_path: \n", TF_targets_path))
    tmp_gene <- read.table(TF_targets_path, stringsAsFactors = F)
    print("8. tmp_gene: ")
    print(tmp_gene)
    #print(paste0("8. tmp_gene: \n", tmp_gene))
    tmp_gene <- toupper(tmp_gene$V1)
    print("9. tmp_gene: ")
    print(tmp_gene)
    #print(paste0("9. tmp_gene: \n", tmp_gene))
    gene_list[[which(TF_used == i)]] <- rownames(data)[rownames(data) %in% tmp_gene]
    print("10. gene_list: ")
    print(gene_list)
    #print(paste0("10. gene_list: \n", gene_list))
  }

  TF_used <- TF_used[ unlist(lapply(gene_list, length)) > 10]
  print("11. TF_used: ")
  print(TF_used)
  #print(paste0("11. TF_used: \n", TF_used))
  gene_list <- list()
  for(i in TF_used){
    TF_targets_path <-system.file("extdata", paste0(TF_targets_dir, i), package = "BITFAM")
    print("12. i: ")
    print(i)
    print("13. TF_targets_path: ")
    print(TF_targets_path)
    #print(paste0("12. i: \n", i))
    #print(paste0("13. TF_targets_path: \n", TF_targets_path))
    tmp_gene <- read.table(TF_targets_path, stringsAsFactors = F)
    print("14. tmp_gene: ")
    print(tmp_gene)
    #print(paste0("14. tmp_gene: \n", tmp_gene))
    tmp_gene <- toupper(tmp_gene$V1)
    print("15. tmp_gene: ")
    print(tmp_gene)
    #print(paste0("15. tmp_gene: \n", tmp_gene))
    gene_list[[which(TF_used == i)]] <- rownames(data)[rownames(data) %in% tmp_gene]
    print("16. gene_list: ")
    print(gene_list)
    #print(paste0("16. gene_list: \n", gene_list))
  }
  
  if(is.na(scATAC_obj)){
  }else{
    for(i in TF_used){
      gene_list[[which(TF_used == i)]] <- gene_list[[which(TF_used == i)]][gene_list[[which(TF_used == i)]] %in% BITFAM_scATAC(scATAC_obj)]
    }
  }
  print("17. This is after the scATAC code block: ")
  
  X <- t(as.matrix(data))
  print("18. X: ")
  print(X)
  #print(paste0("18. X: \n", X))
  chipseq_weight <- matrix(1, nrow = length(colnames(X)), ncol = length(TF_used))
  print("19. chipseq_weight: ")
  print(chipseq_weight)
  #print(paste0("19. chipseq_weight: \n", chipseq_weight))
  for(i in 1:length(TF_used)){
    chipseq_weight[, i] <- ifelse(colnames(X) %in% gene_list[[i]], 1, 0)
    print("20. i: ")
    print(i)
    print("21. chipseq_weight: ")
    print(chipseq_weight)
    #print(paste0("20. i: \n", i))
    #print(paste0("21. chipseq_weight: \n", chipseq_weight))
  }
  print("22. This is after the chipseq_weight code block: ")

  Mask_matrix <- chipseq_weight
  print("23. Mask_matrix: ")
  print(Mask_matrix)
  #print(paste0("23. Mask_matrix: \n", Mask_matrix))
  X <- t(as.matrix(data))
  print("24. X: ")
  print(X)
  #print(paste0("24. X: \n", X))
  N <- dim(X)[1]
  print("25. N: ")
  print(N)
  #print(paste0("25. N: \n", N))
  D <- dim(X)[2]
  print("26. D: ")
  print(D)
  #print(paste0("26. D: \n", D))
  K <- length(TF_used)
  print("27. K: ")
  print(K)
  #print(paste0("27. K: \n", K))
  data_to_model <- list(N = N, D = D, K = K, X = X, Mask = Mask_matrix)
  print("28. data_to_model: ")
  print(data_to_model)
  #print(paste0("28. data_to_model: \n", data_to_model))

  print("29. This is immediately before rstan run: ")


  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = ncores)

  set.seed(100)
  pca_beta_piror <- "
data {
int<lower=0> N; // Number of samples
int<lower=0> D; // The original dimension
int<lower=0> K; // The latent dimension
matrix[N, D] X; // The data matrix
matrix[D, K] Mask; // The binary mask of prior knowledge indicate the target of TFs
}

parameters {
matrix<lower=0, upper=1>[N, K] Z; // The latent matrix
matrix[D, K] W; // The weight matrix
real<lower=0> tau; // Noise term
vector<lower=0>[K] alpha; // ARD prior
}

transformed parameters{
matrix<lower=0>[D, K] t_alpha;
real<lower=0> t_tau;
for(wmd in 1:D){
for(wmk in 1:K){
t_alpha[wmd, wmk] = Mask[wmd, wmk] == 1 ? inv(sqrt(alpha[wmk])) : 0.01;
}
}
t_tau = inv(sqrt(tau));
}
model {
tau ~ gamma(1,1);
to_vector(Z) ~ beta(0.5, 0.5);
alpha ~ gamma(1e-3,1e-3);
for(d in 1:D){
for(k in 1:K){
W[d,k] ~ normal(0, t_alpha[d, k]);
}
}
to_vector(X) ~ normal(to_vector(Z*W'), t_tau);
} "
  
  print("30. pca_beta_prior has been computed ")
  m_beta_prior <- stan_model(model_code = pca_beta_piror)
  print("31. m_beta_prior has been computed ")
  suppressWarnings(fit.vb <- vb(m_beta_prior, data = data_to_model, algorithm = "meanfield",
                                  iter = iter, output_samples = 300, tol_rel_obj = tol_rel_obj))
  print("32. fit.vb has been computed ")
  BITFAM_list <- list(Model = fit.vb,
                      TF_used = TF_used,
                      Genes = rownames(data),
                      Cell_names = colnames(data))
  print("33. BITFAM_list has been computed ")
  return(BITFAM_list)
}






# Getting the assay list ready function -----------------------------------

file_clean = function(filtered_cr){
  # Enter the project name and read the raw and filtered datasets as 10x matrix
  filter_read = Read10X(filtered_cr)
  
  # Provide a list of assays present in the dataset
  assays_present = readline(prompt = "Enter all the assay names present with space (donor crispr rna condition tetramer protein): ")
  assay_list = unlist(strsplit(assays_present, " "))
  
  # If the assay_list contains protein
  if("protein" %in% assay_list){
    # Read the Antibody capture matrix as dataframe
    filter_read.df = data.frame(filter_read[["Antibody Capture"]], check.names = FALSE)
    
    # Split the dataframe into Donor and Protein dataframes and replace underscores with hypen
    filter_read.donor = filter_read.df[grepl("hashtag", rownames(filter_read.df)), ]
    filter_read.protein = filter_read.df[!grepl("hashtag", rownames(filter_read.df)), ]
    rownames(filter_read.donor) = gsub("_", "-", rownames(filter_read.donor))
    rownames(filter_read.protein) = gsub("_", "-", rownames(filter_read.protein))
    
    # If there any cell that doesn't have any donor then remove those cells 
    if(any(colSums(filter_read.donor) == 0)){
      filter_read.donor <- filter_read.donor[, colSums(filter_read.donor) != 0]
      barcodes <- intersect(colnames(filter_read.donor), colnames(filter_read[["Gene Expression"]]))
      
      # remove those cells from other assays as well
      filter_read.protein <- filter_read.protein[, barcodes]
      filter_read[["Gene Expression"]] <- filter_read[["Gene Expression"]][,barcodes]
      filter_read[["CRISPR Guide Capture"]] <- filter_read[["CRISPR Guide Capture"]][, barcodes]
    }
    # adding back the split assays to filter_read
    filter_read[["Protein"]] <- as(as.matrix(filter_read.protein), "dgCMatrix")
    filter_read[["Donor"]] <- as(as.matrix(filter_read.donor), "dgCMatrix")
  } # if the assay list contains tetramer
  else if ("tetramer" %in% assay_list){
    filter_read.df = data.frame(filter_read[["Antibody Capture"]], check.names = FALSE)
    rownames(filter_read.df) = gsub("_", "-", rownames(filter_read.df))
    if(any(colSums(filter_read.df) == 0)){
      filter_read.df <- filter_read.df[, colSums(filter_read.df) != 0]
      barcodes <- intersect(colnames(filter_read.df), colnames(filter_read[["Gene Expression"]]))
      
      # remove those cells from other assays as well
      filter_read[["Gene Expression"]] <- filter_read[["Gene Expression"]][,barcodes]
      filter_read[["CRISPR Guide Capture"]] <- filter_read[["CRISPR Guide Capture"]][, barcodes]
    }
    # adding back the tetramer assays to filter_read
    filter_read[["tetramer"]] <- as(as.matrix(filter_read.df), "dgCMatrix")
  }
  else if ("condition" %in% assay_list){
    filter_read.condition <- data.frame(filter_read[["Antibody Capture"]], check.names = FALSE)
    rownames(filter_read.condition) = gsub("_", "-", rownames(filter_read.condition))
    filter_read[["condition"]] <- as(as.matrix(filter_read.condition), "dgCMatrix")
  }
  return(filter_read)
}

# Setting the CRISPR assay for Seurat Object ----------------------------------------------------
# Function for naming clean of sgRNA guides

# Function to remove the numbers in the sgRNA names
name_clean = function(x){
  result_vector <- c()
  for (i in x){
    y <- substr(i, 1,nchar(i)-2)
    result_vector <- c(result_vector, y)
  }
  return(result_vector)
}

# Making columns for custom metadata --------------------------------------
# Making the column names of perturb_seq.df and assigning them as rownames of new_df
new_df_fun = function(df){
  cols <- colnames(df)
  new_df <- data.frame(matrix(ncol = 0, nrow = length(cols)))
  rownames(new_df) <- cols
  
  # Taking the last character of barcodes as sample number
  new_df$sample_num <- substr(cols, nchar(cols), nchar(cols))
  
  # Adding columns to new_df
  new_df$sample_name <- readline(prompt = "Enter sample name ('PS0X') : ")
  new_df$`sgRNAs_3+` <- 'none'
  new_df$`targets_3+`	<- 'none'
  new_df$`sgRNAs_2+`	<- 'none'
  new_df$`targets_2+` <- 'none'
  new_df$guide_ids <- 'none'
  new_df$gene_ids <- 'none'
  new_df$crispr_perturbation <- 'none'
  for (i in guide_set){
    new_df[[i]] <- 'none' # Adding gene columns to new_df
  }
  return(new_df)
}


# Making the custom metadata to add crispr information

make_meta = function(origin.df, new_df){
  cols <- colnames(origin.df)
  for (i in 1:length(cols)){
    # For each barcode in perturb_seq.df, selecting only rows that have non zero values and adding them to guide_ids
    result <- rownames(origin.df)[origin.df[, cols[i]] >=1 ]
    new_df[i,'guide_ids'] <- paste(result, collapse = ", ")
    
    # Cleaning the names of guide RNAs and storing them in gene_ids column
    guide_name <- name_clean(strsplit(new_df[i, 'guide_ids'],',')[[1]])
    new_df[i,'gene_ids'] <- paste(guide_name, collapse = ", ")
    
    #If there are more than 1 guide RNA
    if (length(strsplit(new_df[i, 'guide_ids'],',')[[1]])> 1){
      new_df[i,'sgRNAs_3+'] <- 'multi'
      new_df[i,'sgRNAs_2+'] <- 'multi'
      new_df[i,'targets_3+'] <- 'multi'
      new_df[i,'targets_2+'] <- 'multi'
    }# for one guide RNA storing the name of the guide RNA 
    else if (length(strsplit(new_df[i,  'guide_ids'],',')[[1]])== 1){
      new_df[i,'sgRNAs_3+'] <- new_df[i,'guide_ids']
      new_df[i,'sgRNAs_2+'] <- new_df[i,'guide_ids']
      new_df[i,'targets_3+'] <- substr(new_df[i,'guide_ids'], 1, nchar(new_df[i,'guide_ids'])-2)
      new_df[i,'targets_2+'] <- substr(new_df[i,'guide_ids'], 1, nchar(new_df[i,'guide_ids'])-2)
    }# for no guide RNA
    else{
      new_df[i,'sgRNAs_3+'] <- 'none'
      new_df[i,'sgRNAs_2+'] <- 'none'
    }
    #Adding guide RNAs detected for each barcode in crispr_perturbation column
    new_df[i,'crispr_perturbation'] <- length(strsplit(new_df[i, 'guide_ids'],',')[[1]])
  }
  # Assigning boolean values to each guide for each cell
  for (i in 1:length(cols)){
    guide_list <- trimws(strsplit(new_df[i, 'gene_ids'],',')[[1]])
    for (j in guide_set){
      if (j %in% guide_list){
        new_df[i, j] <- 1
      }else{
        new_df[i, j] <- 0 
      }
    }
  }
  # Resting guide column name 
  for (w in guide_set){
    t <- paste0(w,"_bool")
    names(new_df)[names(new_df) == w] <- t
  }
  return(new_df)
}


# making the seurat_object ------------------------------------------------

make_seurat_object = function(filter_assay, raw_assay_path){
  assay_type <- names(filter_assay)
  for(i in assay_type){
    if(i == "Gene Expression"){
      raw_assay <- Read10X(raw_assay_path)
      decont <- decontX(filter_assay[[i]], background = raw_assay[["Gene Expression"]])
      dataset_name <- deparse(substitute(filter_assay))
      #QC Plots
      decontx_plot <- plotDecontXContamination(decont)
      
      # Save plot
      ggsave(filename = paste0(dataset_name, "_decontx_UMAP.pdf"), plot = decontx_plot, width = 16, height = 10)
      
      seurat_obj <- CreateSeuratObject(counts = round(decont[["decontXcounts"]]), project = readline(prompt = "Enter project name ('PS0X'): "))
    }else if(i == "CRISPR Guide Capture"){
      crispr.df <- data.frame(counts = filter_assay[[i]], check.names = FALSE)
      df <- new_df_fun(crispr.df)
      df <- make_meta(crispr.df, df)
      seurat_obj@meta.data <- cbind(seurat_obj@meta.data, df)
      seurat_obj[["CRISPR"]] <- CreateAssayObject(counts = filter_assay[[i]])
    }else if(i == "Protein"){
      seurat_obj[["PROTEIN"]] <- CreateAssayObject(counts = filter_assay[[i]])
    }else if(i == 'condition'){
      seurat_obj[["HASH"]] <- CreateAssayObject(counts = filter_assay[[i]])
    }else if(i == "Donor"){
      seurat_obj[["DONOR"]] <- CreateAssayObject(counts = filter_assay[[i]])
    }else if(i == "tetramer"){
      seurat_obj[["HTO"]] <- CreateAssayObject(counts = filter_assay[[i]])
    }
  }
  return(seurat_obj)
}

remove_bad_cells <- function(seurat_obj){
  # removing mitochondrial, ribosomal, platelet, hemoglobin RNA's
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")
  seurat_obj[["percent.plat"]] <- PercentageFeatureSet(seurat_obj, pattern = "PECAM1|PF4")
  seurat_obj.old <- seurat_obj
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & 
                         percent.mt < 10 & percent.rb > 5 & percent.hb < 1 & percent.plat < 0.1)
  
  # QC plots
  dataset_name <- readline(prompt ('PS0X') = "Enter project name: ")
  
  # Scatter plot
  scatter_plot <- FeatureScatter(seurat_obj.old, "nCount_RNA", "nFeature_RNA", pt.size = 0.5) + 
    ggtitle("Pearson correlation between cell and transcripts count of Unfiltered data")
  
  # Violin plots
  vln_plot_old <- VlnPlot(seurat_obj.old, 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb", "percent.plat"), 
                          pt.size = 0, 
                          ncol = 3) + 
    NoLegend() + plot_annotation ("Unfiltered data", theme = theme(plot.title = element_text(size = 28)))
  
  vln_plot_new <- VlnPlot(seurat_obj, 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb", "percent.plat"), 
                          pt.size = 0, 
                          ncol = 3) + 
    NoLegend() + plot_annotation ("Filtered data", theme = theme(plot.title = element_text(size = 28)))
  
  # Combine plots
  combined_plot <- plot_grid(vln_plot_old, vln_plot_new, ncol = 2)
  
  # Save plot
  ggsave(filename = paste0(dataset_name, "_QC_mtRNAs.pdf"), plot = combined_plot, width = 16, height = 10)
  ggsave(filename = paste0(dataset_name, "correlation.pdf"), plot = scatter_plot, width = 16, height = 10)
  return(seurat_obj)
}

standard_seurat_rna <- function(seurat_obj){
  DefaultAssay(seurat_obj) <- "RNA"
  variable_to_regress = readline(prompt = "Enter whether you need to regress for mt RNAs, types yes or no: ")
  selection_type <- readline(prompt = "Enter selection.method type (vst, mean.var.plot, dispersion) to find variable features:  ")
  
  # For RNA assay first cycle of transcripts normalization, cell cycle scoring, pca
  if(variable_to_regress == "yes"){
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    all.genes <- rownames(seurat_obj)
    seurat_obj <- seurat_obj %>%
      NormalizeData(verbose = F) %>%
      FindVariableFeatures(selection.method = selection_type,  verbose = F)  %>%
      ScaleData(features = all.genes, vars.to.regress = c("nFeature_RNA", "percent.mt"),  verbose = F) %>%
      CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
      RunPCA(reduction.name = "pca_nocycle", reduction.key = "pc_ncyc", verbose = F)
    dataset_name <- readline(prompt = "Enter project name ('PS0X'): ")
    eblow_pca_1 <- ElbowPlot(seurat_obj, reduction = "pca_nocycle")
    ggsave(filename = paste0(dataset_name, "first_rna_pca_elbow.pdf"), plot = eblow_pca_1, width = 16, height = 10)
  }
  # To check for in RNA assay post dimultiplexing and calculation of perturbation signatures
  else if (variable_to_regress == "no"){
    reduction_name <- readline(prompt = "Enter new pca name (pcadimul, etc): ")
    reduction_key <- readline(prompt = "Enter new pca key (PCdimul, etc): ")
    seurat_obj <- seurat_obj %>%
      NormalizeData( verbose = F) %>%
      FindVariableFeatures(selection.method = selection_type,  verbose = F)  %>%
      ScaleData(features = VariableFeatures(seurat_obj),  verbose = F) %>%
      RunPCA(features = VariableFeatures(perturb.seq), reduction.name = reduction_name, reduction.key = reduction_key, approx = F, verbose = F)
  }
  return(seurat_obj)
}

dimultiplex_assay <- function(seurat_obj){
  assay_name <- readline(prompt = "Enter assay Name (HASH for condition, HTO for tetramer, for rest assay names in Caps): ")
  # For condition assay
  if(assay_name == "HASH"){
    seurat_obj<- seurat_obj %>% NormalizeData(assay = assay_name, normalization.method = "CLR", margin = 2) %>%
      HTODemux(assay = assay_name, positive.quantile = 0.99)
  }
  # For tetramer or donor
  else if(assay_name == "HTO" || assay_name == "DONOR"){
    seurat_obj<- seurat_obj %>% NormalizeData(assay = assay_name, normalization.method = "CLR", margin = 2) %>%
      HTODemux(assay = assay_name, kfunc = "kmeans", positive.quantile = 0.99)
  }
  # For protein assay
  else{
    seurat_obj<- seurat_obj %>% NormalizeData(assay = assay_name, normalization.method = "CLR", margin = 2)
  }
  return(seurat_obj)
}

embeddings_for_dimutliplex <- function(seurat_obj){
  assay_name <- readline(prompt = "Enter assay Name: ")
  seurat_obj <- seurat_obj %>%
    ScaleData(assay = assay_name, features = rownames(seurat_obj[[assay_name]]), verbose = F) %>%
    RunPCA(assay = assay_name, features = rownames(seurat_obj[[assay_name]]), reduction.name = "pca_dimulassay", 
           reduction.key = "pc_dimulassay", approx = F) %>% 
    RunTSNE(assay = assay_name, reduction = "pca_dimulassay",reduction.name = "tsne_dimulassay",
            reduction.key = "tsne_dimulassay", dims = 1:10, perplexity = 100, check_duplicates = F)
  return(seurat_obj)
}

# save_plotsets <- function(seurat_obj) {
#   # Parsing seurat object name
#   object <- readline(prompt = "Enter project name: ")
#   # Which plots
#   options_plot <- readline(prompt = "Enter the plot type: ")
#   # Naming for plots
#   dataset_name <- readline(prompt = "Enter filename: ")
#   
#   #Ridgeplot
#   if(options_plot == "ridge"){
#     # Which assay
#     assay_name <- readline(prompt = "Enter name for which assay to use: ")
#     
#     # Group by the cell based on the metadata
#     group_by_opt <- readline(prompt = "Enter group the cells by: ")
#     ridge_plot <- RidgePlot(seurat_obj, assay = assay_name, features = rownames(seurat_obj[[assay_name]]), 
#                             group.by = group_by_opt, ncol = 2)
#     ggsave(filename = paste0(object, dataset_name, ".pdf"), plot = ridge_plot, width = 10, height = 8)
#   }
#   # ViolinPlot
#   else if(options_plot == "vln"){
#     # Which assay
#     assay_name <- readline(prompt = "Enter name for which assay to use: ")
#     
#     # Group by the cell based on the metadata
#     group_by_opt <- readline(prompt = "Enter group the cells by: ")
#     enter_features <- unlist(strsplit(readline(prompt = "Enter column names from metadata: ")))
#     vln_plot <- VlnPlot(seurat_obj, features = enter_features, group.by = group_by_opt, pt.size = 0.1, log = TRUE)
#     ggsave(filename = paste0(object, dataset_name, ".pdf"), plot = vln_plot, width = 10, height = 8)
#   }
#   # Dimplot
#   else if(options_plot == "dim"){
#     # Which assay
#     assay_name <- readline(prompt = "Enter name for which assay to use: ")
#     
#     # Group by the cell based on the metadata
#     group_by_opt <- readline(prompt = "Enter group the cells by: ")
#     reduction_name = readline(prompt = "Enter reduction: ")
#     dim_plot <- DimPlot(object, seurat_obj, group.by = group_by_opt, reduction = reduction_name)
#     ggsave(filename = paste0(object, dataset_name, ".pdf"), plot = dim_plot, width = 10, height = 8)
#   }
#   #Heatmap
#   else if (options_plot == "heat"){
#     # Which assay
#     assay_name <- readline(prompt = "Enter name for which assay to use: ")
#     
#     # Group by the cell based on the metadata
#     group_by_opt <- readline(prompt = "Enter group the cells by: ")
#     heatmap_plot <- HTOHeatmap(seurat_obj, assay = assay_name)
#     ggsave(filename = paste0(object, dataset_name, ".pdf"), plot = heatmap_plot, width = 10, height = 8)
#   }
#   return(seurat_obj)
# }



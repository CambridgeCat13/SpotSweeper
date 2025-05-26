detectEdgeDryspots <- function(
    spe, 
    qc_metric = "sum_gene",
    samples = "sample_id", 
    mad_threshold = 3,
    edge_threshold = 0.75,
    shifted = FALSE,
    batch_var = "both",
    name = "edge_dryspot") {
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input data must be a SpatialExperiment or inherit from SpatialExperiment.")
  }
  if (!qc_metric %in% colnames(colData(spe))) {
    stop("qc_metric must be present in colData.")
  }
  if (!samples %in% colnames(colData(spe))) {
    stop("samples column must be present in colData.")
  }
  required_cols <- c("in_tissue", "array_row", "array_col")
  missing_cols <- required_cols[!required_cols %in% colnames(colData(spe))]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  if (!is.numeric(mad_threshold) || mad_threshold <= 0) {
    stop("'mad_threshold' must be a positive numeric value.")
  }
  # ===== Convert Jacqui's script to function =====
  lg10_metric <- paste0("lg10_", qc_metric)
  colData(spe)[[lg10_metric]] <- log10(colData(spe)[[qc_metric]])
  if (batch_var %in% c("slide", "both")) {
    if ("slide" %in% colnames(colData(spe))) {
      outlier_slide_col <- paste0(qc_metric, "_3MAD_outlier_slide")
      colData(spe)[[outlier_slide_col]] <- isOutlier(
        colData(spe)[[lg10_metric]], 
        subset = colData(spe)$in_tissue, 
        batch = colData(spe)$slide, 
        type = "lower", 
        nmads = mad_threshold
      )
    } else {
      warning("'slide' column not found, skipping slide-level outlier detection")
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_slide")]] <- FALSE
    }
  }

  if (batch_var %in% c("sample_id", "both")) {
    outlier_sample_col <- paste0(qc_metric, "_3MAD_outlier_sample")
    colData(spe)[[outlier_sample_col]] <- isOutlier(
      colData(spe)[[lg10_metric]], 
      subset = colData(spe)$in_tissue, 
      batch = colData(spe)[[samples]], 
      type = "lower", 
      nmads = mad_threshold
    )
  }
  
  if (batch_var == "both") {
    colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_slide")]] | 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_sample")]]
  } else if (batch_var == "slide") {
    colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_slide")]]
  } else {
    colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
      colData(spe)[[paste0(qc_metric, "_3MAD_outlier_sample")]]
  }
  
  colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]] <- 
    ifelse(colData(spe)$in_tissue == FALSE, FALSE, 
           colData(spe)[[paste0(qc_metric, "_3MAD_outlier_binary")]])
  
  sampleList <- unique(colData(spe)[[samples]])
  names(sampleList) <- sampleList
  message("Detecting edges...")
  genes_edges <- lapply(sampleList, function(x) {
    tmp <- colData(spe)[colData(spe)[[samples]] == x, 
                        c("in_tissue", "array_row", "array_col", 
                          paste0(qc_metric, "_3MAD_outlier_binary"))]
    clumpEdges(tmp[, -1], rownames(tmp)[tmp$in_tissue == FALSE], 
               shifted = shifted, edge_threshold = edge_threshold)
  })
  
  colData(spe)[[paste0(name, "_edge")]] <- FALSE
  edge_spots <- unlist(genes_edges)
  if (length(edge_spots) > 0) {
    colData(spe)[edge_spots, paste0(name, "_edge")] <- TRUE
  }
  message("Number of samples with edges detected: ", 
          sum(sapply(genes_edges, length) > 0))
  samples_with_edges <- names(genes_edges)[sapply(genes_edges, length) > 0]
  if (length(samples_with_edges) > 0) {
    message("Samples with edges detected: ", paste(samples_with_edges, collapse = ", "))
  }
  
  message("Finding problem areas...")
  genes_probs <- lapply(sampleList, function(x) {
    tmp <- colData(spe)[colData(spe)[[samples]] == x, 
                        c("in_tissue", "array_row", "array_col", 
                          paste0(qc_metric, "_3MAD_outlier_binary"))]
    problemAreas(tmp[, -1], 
                 rownames(tmp)[tmp$in_tissue == FALSE], 
                 uniqueIdentifier = x, 
                 shifted = shifted)
  })
  genes_probs <- do.call(rbind, genes_probs)
  colData(spe)[[paste0(name, "_problem_id")]] <- NA
  colData(spe)[[paste0(name, "_problem_size")]] <- 0
  
  if (nrow(genes_probs) > 0) {
    colData(spe)[genes_probs$spotcode, paste0(name, "_problem_id")] <- genes_probs$clumpID
    colData(spe)[genes_probs$spotcode, paste0(name, "_problem_size")] <- genes_probs$clumpSize
  }
  
  message("Edge dryspot detection completed!")
  return(spe)
}


classifyEdgeDryspots <- function(
    spe,
    samples = "sample_id",
    min_spots = 20,
    low_umi_threshold = 100,
    removal_threshold = 0.5,
    name = "edge_dryspot",
    exclude_slides = NULL) {

  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input data must be a SpatialExperiment or inherit from SpatialExperiment.")
  }
  
  required_cols <- c("sum_umi", paste0(name, "_edge"), paste0(name, "_problem_id"))
  missing_cols <- required_cols[!required_cols %in% colnames(colData(spe))]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "),
               "\nRun detectEdgeDryspots() first."))
  }
  
  # ===== Convert Jacqui's classification script =====
  colData(spe)[[paste0(name, "_true_edges")]] <- colData(spe)[[paste0(name, "_edge")]]
  
  if (!is.null(exclude_slides) && "slide" %in% colnames(colData(spe))) {
    colData(spe)[[paste0(name, "_true_edges")]] <- 
      ifelse(colData(spe)$slide %in% exclude_slides, FALSE, 
             colData(spe)[[paste0(name, "_edge")]])
    message("Excluding edges from slides: ", paste(exclude_slides, collapse = ", "))
  }

  cdata <- as.data.frame(colData(spe))
  if (requireNamespace("dplyr", quietly = TRUE)) {
    tmp <- cdata %>%
      dplyr::filter(in_tissue == TRUE, 
                    !!sym(paste0(name, "_true_edges")) == FALSE) %>%
      dplyr::mutate(lowumi = sum_umi <= low_umi_threshold) %>%
      dplyr::group_by(!!sym(paste0(name, "_problem_id"))) %>%
      dplyr::summarise(n_lowumi = sum(lowumi), 
                       n_spots = dplyr::n(), 
                       prop_lowumi = n_lowumi/n_spots) %>%
      dplyr::filter(n_spots > min_spots, !is.na(!!sym(paste0(name, "_problem_id"))))
    
    remove.areas <- unique(dplyr::filter(tmp, prop_lowumi >= removal_threshold)[[paste0(name, "_problem_id")]])
  } else {
    cdata_filtered <- cdata[cdata$in_tissue == TRUE & 
                              cdata[[paste0(name, "_true_edges")]] == FALSE, ]
    cdata_filtered$lowumi <- cdata_filtered$sum_umi <= low_umi_threshold
    problem_ids <- unique(cdata_filtered[[paste0(name, "_problem_id")]])
    problem_ids <- problem_ids[!is.na(problem_ids)]
    
    remove.areas <- c()
    for (pid in problem_ids) {
      subset_data <- cdata_filtered[cdata_filtered[[paste0(name, "_problem_id")]] == pid, ]
      n_spots <- nrow(subset_data)
      n_lowumi <- sum(subset_data$lowumi)
      prop_lowumi <- n_lowumi / n_spots
      
      if (n_spots > min_spots && prop_lowumi >= removal_threshold) {
        remove.areas <- c(remove.areas, pid)
      }
    }
  }
  
  message("Number of problem areas marked for removal: ", length(remove.areas))
  colData(spe)[[paste0(name, "_binary")]] <- 
    ifelse(!is.na(colData(spe)[[paste0(name, "_problem_id")]]), "problem area", "none")
  colData(spe)[colData(spe)[[paste0(name, "_edge")]], paste0(name, "_binary")] <- "edge"

  problem_sizes <- colData(spe)[[paste0(name, "_problem_size")]]
  colData(spe)[[paste0(name, "_classification")]] <- 
    ifelse(problem_sizes <= min_spots, "small", "flag")
  colData(spe)[[paste0(name, "_classification")]] <- 
    ifelse(problem_sizes == 0, "none", 
           colData(spe)[[paste0(name, "_classification")]])

  if (length(remove.areas) > 0) {
    colData(spe)[colData(spe)[[paste0(name, "_problem_id")]] %in% remove.areas, 
                 paste0(name, "_classification")] <- "remove"
  }

  colData(spe)[colData(spe)[[paste0(name, "_true_edges")]], 
               paste0(name, "_classification")] <- "edge"
  message("Edge dryspot classification completed!")
  return(spe)
}
optimized_problemAreas <- function(.xyz, offTissue, uniqueIdentifier=NA, shifted=FALSE) {
  
  #  TIMING START 
  prob_start <- Sys.time()
  cat("  → Starting OPTIMIZED problemAreas...\n")
  
  if(sum(.xyz[,3])==0) {
    cat("  ← optimized problemAreas: No data to process\n")
    return(data.frame())
  }
  
  # Coordinate adjustment if needed (keep original logic)
  if(shifted==TRUE) {
    coord_start <- Sys.time()
    odds = seq(1,max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
    coord_time <- difftime(Sys.time(), coord_start, units = "secs")
    cat("    - Coordinate adjustment:", round(coord_time, 2), "seconds\n")
  }
  
  #  KEEP ORIGINAL RASTER PROCESSING FOR NOW 
  # Step 1: Create raster (keep original)
  step1_start <- Sys.time()
  t1 = raster::rasterFromXYZ(.xyz)
  step1_time <- difftime(Sys.time(), step1_start, units = "secs")
  cat("    - rasterFromXYZ:", round(step1_time, 2), "seconds\n")
  
  # Step 2: Focal transformations (keep original)
  step2_start <- Sys.time()
  t2 = focal_transformations(t1)
  step2_time <- difftime(Sys.time(), step2_start, units = "secs")
  cat("    - focal_transformations total:", round(step2_time, 2), "seconds\n")
  
  # Step 3: Matrix conversion (keep original)
  step3_start <- Sys.time()
  rast = as.matrix(t2)
  step3_time <- difftime(Sys.time(), step3_start, units = "secs")
  cat("    - Matrix conversion:", round(step3_time, 2), "seconds\n")
  
  # Step 4: Main clump operation (keep original)
  step4_start <- Sys.time()
  c1 = raster::clump(t2, directions=8)
  step4_time <- difftime(Sys.time(), step4_start, units = "secs")
  cat("    - MAIN raster::clump:", round(step4_time, 2), "seconds ⚠️\n")
  
  #  OPTIMIZE THE BOTTLENECK: FOR-LOOP PROCESSING 
  step5_start <- Sys.time()
  clumps = as.matrix(c1)
  tot <- max(clumps, na.rm=TRUE)
  
  cat("    - Found", tot, "clumps to process\n")
  
  if(is.na(uniqueIdentifier)) uniqueIdentifier = "X"
  
  # NEW OPTIMIZED APPROACH: Vectorized processing
  cat("    - Using VECTORIZED approach instead of for-loop...\n")
  
  # Create a data.frame with all clump information at once
  clump_indices <- which(!is.na(clumps), arr.ind = TRUE)
  clump_values <- clumps[!is.na(clumps)]
  
  # Create the result data.frame in one go
  result_df <- data.frame(
    row = clump_indices[, 1],
    col = clump_indices[, 2], 
    clump_id = paste(uniqueIdentifier, clump_values, sep="_"),
    stringsAsFactors = FALSE
  )
  
  # Calculate sizes for each clump using table (vectorized)
  clump_sizes <- table(clump_values)
  result_df$size <- clump_sizes[as.character(clump_values)]
  
  step5_time <- difftime(Sys.time(), step5_start, units = "secs")
  cat("    - OPTIMIZED problem area processing:", round(step5_time, 2), "seconds\n")
  
  # Step 6: Lookup processing (keep original logic but optimize)
  step6_start <- Sys.time()
  pAreas = lookupKeyDF(.xyz[,1:2], result_df)
  result <- pAreas[!pAreas$spotcode %in% offTissue,]
  step6_time <- difftime(Sys.time(), step6_start, units = "secs")
  cat("    - Lookup processing:", round(step6_time, 2), "seconds\n")
  
  total_time <- difftime(Sys.time(), prob_start, units = "secs")
  cat("  ← OPTIMIZED problemAreas completed:", round(total_time, 2), "seconds total\n")
  
  return(result)
}
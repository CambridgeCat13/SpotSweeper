test_optimization <- function(.xyz, offTissue, uniqueIdentifier=NA, shifted=FALSE) {
  
  cat("COMPARING OLD VS OPTIMIZED VERSIONS\n")
  
  # Test original version
  cat("Testing ORIGINAL problemAreas:\n")
  original_start <- Sys.time()
  original_result <- problemAreas(.xyz, offTissue, uniqueIdentifier, shifted)
  original_time <- difftime(Sys.time(), original_start, units = "secs")
  
  cat("\nTesting OPTIMIZED problemAreas:\n")
  optimized_start <- Sys.time()
  optimized_result <- optimized_problemAreas(.xyz, offTissue, uniqueIdentifier, shifted)
  optimized_time <- difftime(Sys.time(), optimized_start, units = "secs")
  
  # Compare results
  cat("\n PERFORMANCE COMPARISON \n")
  cat("Original time:", round(original_time, 2), "seconds\n")
  cat("Optimized time:", round(optimized_time, 2), "seconds\n")
  
  # Calculate speedup
  speedup_value <- as.numeric(original_time) / as.numeric(optimized_time)
  cat("Speedup:", round(speedup_value, 2), "x\n")
  
  # Check if results are identical
  cat("\n RESULT VALIDATION \n")
  if(identical(original_result, optimized_result)) {
    cat("Results are IDENTICAL - optimization successful!\n")
  } else {
    cat("Results differ - need to check implementation\n")
    cat("Original rows:", nrow(original_result), "Optimized rows:", nrow(optimized_result), "\n")
    
    cat("\n DETAILED COMPARISON \n")
    cat("Original columns:", paste(colnames(original_result), collapse=", "), "\n")
    cat("Optimized columns:", paste(colnames(optimized_result), collapse=", "), "\n")
    
    cat("\nOriginal first few rows:\n")
    print(head(original_result))
    cat("\nOptimized first few rows:\n")
    print(head(optimized_result))
    
    cat("\n CONTENT VERIFICATION \n")
    original_spots <- sort(original_result$spotcode)
    optimized_spots <- sort(optimized_result$spotcode)
    spots_identical <- identical(original_spots, optimized_spots)
    cat("Same spots (ignoring order):", spots_identical, "\n")
    
    original_sizes <- table(original_result$clumpSize)
    optimized_sizes <- table(optimized_result$clumpSize)
    sizes_identical <- identical(original_sizes, optimized_sizes)
    cat("Same clump size distribution:", sizes_identical, "\n")
    
    cat("Total clumps - Original:", length(unique(original_result$clumpID)),
        "Optimized:", length(unique(optimized_result$clumpID)), "\n")
    
    if(spots_identical && sizes_identical) {
      cat("CONTENT IS IDENTICAL - only order differs (this is OK!)\n")
    }
  }
  
  return(list(
    original_time = as.numeric(original_time),
    optimized_time = as.numeric(optimized_time),
    speedup = speedup_value,
    results_identical = identical(original_result, optimized_result),
    original_result = original_result,
    optimized_result = optimized_result
  ))
}
# ==============================================================================
# Edge Detection Helper Functions
# Adapted from Jacqui's raster_edge_functions.r for SpotSweeper integration
# ==============================================================================

my_fill <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  x[is.na(x)] = 1
  if(sum(x[-5])==8) return(1)
  else return(0)
}

my_fill_star <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  if(sum(x[-5], na.rm=T)==4) return(1)
  else return(0)
}


my_outline <- function(x) {
  if(is.na(x[13])) return(NA_integer_)
  if(x[13]==1) return(1)
  x[is.na(x)] = 1
  sum.edges = sum(x[1:5])+sum(x[21:25])+sum(x[c(6,11,16)])+sum(x[c(10,15,20)])
  if(sum.edges==16) return(1)
  else(return(0))
}

focal_transformations <- function(raster_object) {
  r2 <- raster::focal(raster_object, w=matrix(1,3,3), fun=my_fill, pad=T, padValue=1)
  r3 <- raster::focal(r2, w=matrix(1,5,5), fun=my_outline, pad=T, padValue=1)
  star.m = matrix(rep(c(0,1), length.out=9), 3, 3)
  star.m[5] = 1
  r3_s <- raster::focal(r3, star.m, fun=my_fill_star, pad=T, padValue=1)
  rev_r3 = r3_s
  rev_r3[r3_s==0] = 1
  rev_r3[r3_s==1] = 0
  rev_c3 = raster::clump(rev_r3)
  tbl = table(as.matrix(rev_c3))
  flip_clump = as.numeric(names(tbl)[tbl<40])
  r4 = r3_s
  r4[rev_c3 %in% flip_clump] = 1
  return(r4)
}


lookupKey <- function(.xy, results_df) {
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  t1 = raster::rasterFromXYZ(key1)
  rast = as.matrix(t1)
  
  num_key = sapply(1:nrow(results_df), function(x) {
    v1 = as.numeric(results_df[x,])
    rast[v1[1],v1[2]]
  })
  edge_spots = rownames(key1)[key1$numeric_key %in% num_key] 
  return(edge_spots)
}

lookupKeyDF <- function(.xy, results_df) {
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  t1 = raster::rasterFromXYZ(key1)
  rast = as.matrix(t1)
  
  num_key = do.call(rbind, lapply(1:nrow(results_df), function(x) 
    cbind.data.frame("numeric_key"=rast[results_df[x,1],results_df[x,2]], "clump_id"=results_df[x,3], "size"=results_df[x,4])
  ))
  
  matching = match(num_key$numeric_key, key1$numeric_key)
  return(cbind.data.frame("spotcode"=rownames(key1)[matching],"clumpID"=num_key$clump_id, "clumpSize"=num_key$size))
}


clumpEdges <- function(.xyz, offTissue, shifted=FALSE, edge_threshold=0.75) {
  if(sum(.xyz[,3])==0) return(c())
  if(shifted==TRUE) {
    odds = seq(1,max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  t1 = raster::rasterFromXYZ(.xyz)
  t2 <- focal_transformations(t1)
  rast = as.matrix(t2)
  c1 = raster::clump(t2, directions=8)
  clumps = as.matrix(c1)
  edgeClumps = c()
  north = sum(!is.na(clumps[1,]))/sum(!is.na(rast[1,]))
  if(north>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[1,]))
  east = sum(!is.na(clumps[,ncol(clumps)]))/sum(!is.na(rast[,ncol(rast)]))
  if(east>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[,ncol(clumps)]))
  south = sum(!is.na(clumps[nrow(clumps),]))/sum(!is.na(rast[nrow(rast),]))
  if(south>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[nrow(clumps),]))
  west = sum(!is.na(clumps[,1]))/sum(!is.na(rast[,1]))
  if(west>=edge_threshold) edgeClumps = c(edgeClumps, unique(clumps[,1]))
  
  edgeClumps = edgeClumps[!is.na(edgeClumps)]
  if(length(edgeClumps)==0) return(c())
  
  res <- vector("list", length(edgeClumps))
  names(res) <- as.character(edgeClumps)
  for (i in edgeClumps){
    res[[as.character(i)]] <- as.data.frame(which(clumps == i, arr.ind = TRUE))
  }
  res.df = do.call(rbind, res) 
  edgeSpots = lookupKey(.xyz[,1:2], res.df)
  return(setdiff(edgeSpots, offTissue))
}


problemAreas <- function(.xyz, offTissue, uniqueIdentifier=NA, shifted=FALSE) { 
  if(sum(.xyz[,3])==0) return(data.frame())
  if(shifted==TRUE) {
    odds = seq(1,max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = 
      .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  t1 = raster::rasterFromXYZ(.xyz)
  t2 = focal_transformations(t1)
  rast = as.matrix(t2)
  c1 = raster::clump(t2, directions=8)
  clumps = as.matrix(c1)
  tot <- max(clumps, na.rm=TRUE)
  res <- vector("list",tot)
  if(is.na(uniqueIdentifier)) uniqueIdentifier = "X"
  
  for (i in 1:tot){
    res[i] <- list(which(clumps == i, arr.ind = TRUE))
    res[[i]] <- cbind.data.frame(
      res[[i]],
      "clump_id"=paste(uniqueIdentifier, i, sep="_"),
      "size"=nrow(res[[i]])
    )
  }
  res.df = do.call(rbind.data.frame, res) 
  pAreas = lookupKeyDF(.xyz[,1:2], res.df)
  return(pAreas[!pAreas$spotcode %in% offTissue,])
}
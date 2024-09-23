# A simple function to get mean Rho in genomic ranges
get_r = function(query, ldmap_range, weighted = FALSE, ncores = 8){
  # overlap = findOverlaps(query, ldmap_range)
  # overlap = split(overlap, queryHits(overlap))
  # r = unlist(lapply(overlap, function(x){mean(ldmap_range$Mean_rho[subjectHits(x)])}))
  index = seq(1, length(query))
  res = split(query, as.factor(index))
  res = pbmclapply(res, function(x){findOverlaps(x, ldmap_range)}, mc.cores = ncores)
  
  # DEBUG
  # x = res[993]
  
  if (weighted) {
    weigthed.rho = function(x) {
      queryhit = as.integer(names(x))
      y = x[[1]]
      # Get ldmap ranges
      intervals = ldmap_range[subjectHits(y)]
      # Restrict to query range
      start = start(query)[queryhit]
      end = end(query)[queryhit]
      intervals = restrict(intervals, start, end)
      # Get widths of each interval (i.e. weigths)
      weigth = width(intervals)
      # Get Rho of each interval (i.e. x)
      value = intervals$Mean_rho
      # Compute the weighted mean
      weigthed.mean = sum(weigth*value)/sum(weigth)
      return(weigthed.mean)
    }
    
    res2 = pbmclapply(1:length(res), function(x){weigthed.rho(res[x])})
    r = unlist(res)
  } else {
    res = pbmclapply(res, function(x){mean(ldmap_range$Mean_rho[subjectHits(x)])}, mc.cores = ncores)
    r = unlist(res)
  }
  return(r)
  
  # Check that coordinates are within the LD map
  # if (x < max(ldmap_trimmed$end[which(ldmap_trimmed$chromosome == chr)], na.rm = TRUE)) {
  #   r = ldmap_trimmed$Mean_rho[which((ldmap_trimmed$start < x & ldmap_trimmed$end > x) & ldmap_trimmed$chromosome == chr)]
  # } else {
  #   r = NA
  # }
  # if (length(r) == 0) {
  #   r = NA
  # }
}




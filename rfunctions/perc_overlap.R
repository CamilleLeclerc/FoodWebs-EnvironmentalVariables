perc_overlap = function(x.start, x.end, y.start, y.end){
  #  if(x.start == y.start & x.end == y.end){
  #    return(100)
  #  }
  x.len = abs(x.end - x.start)
  # largest start
  max.start = max(c(x.start, y.start))
  min.end = min(c(x.end, y.end))
  overlap = min.end - max.start
  overlap = ifelse(overlap <= 0, 0, overlap)
  perc_overlap = overlap / x.len * 100
  return(perc_overlap)
}
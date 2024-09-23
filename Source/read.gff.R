read.gff = function(gff_file) {
  gff = read.table(gzfile(gff_file), header = T, sep = "\t", quote = "",
                   na.strings = c(NA, "NaN", ""), fill = TRUE)
  colnames(gff)[1] = "seqname"
  # colnames(gff)[colnames(gff) == "seqnames"] = "seqname"
  gff$seqname = as.factor(gff$seqname)
  gff$feature = as.factor(gff$feature)
  gff$start = as.integer(gff$start)
  gff$end = as.integer(gff$end)
  return(gff)
}

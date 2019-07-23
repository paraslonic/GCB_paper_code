
# og.file <- "LF82_B2_selected_fna/faa/Results_Dec21/Orthogroups.csv"
# c.file <- "LF82_B2_graph/prob_window_complexity_contig_NC_011993.1.txt"
# strain <- "GCF_000284495.1_ASM28449v1_genomic"

calc_tab <- function(c.file, og.file, strain, logTranform){
  og. <- fread(og.file, head=TRUE)
  og. <- data.frame(og.$V1, og.[,strain, with=FALSE], stringsAsFactors = F)
  colnames(og.) <- c("og","gene")
  og. <- subset(og., !grepl(", ",og.$gene) & og.$gene!="")
  c. <- read.delim(c.file, head = FALSE)
  tab <- merge(c., og., by.x="V1",by.y="og")
  colnames(tab) <- c("og","complexity","gene")
  q <- lapply(tab$gene, function(x) { strsplit(x, "\\|")[[1]]} )
  q <- do.call(rbind, q)
  tab$chr <- q[,4]
  tab$pos <- (as.integer(q[,5])+as.integer(q[,6]))/2
  tab$start <- as.integer(q[,5])
  tab$end <- as.integer(q[,6])
  tab. <- tab
  tab. <- tab.[order(tab.$pos),]  
  ## transform
  if(logTranform){
    print("qu")
    c. <- tab.$complexity
    c.log <- log(c.+min(c.[c.!=0]))
    c.log <- c.log-min(c.log)
    tab.$complexity = c.log
  }
  return(tab.)
}

plot_pvogs <- function(pvog.file, h., scale.){
  print(pvog.file)
  pvogs <- read.delim(pvog.file, head = F)
  ## white list from https://github.com/bobeobibo/phigaro/blob/master/phigaro/const.py
  white.list <- read.delim("pvogs/White.lst", head = F, sep=" ", as.is=T)
  pvogs <- subset(pvogs, pvogs$V1 %in% white.list$V1) 
  q <- lapply(as.character(pvogs$V3), function(x) { strsplit(x, "\\|")[[1]]} )
  q <- do.call(rbind, q)
  pos <- (as.integer(q[,5])+as.integer(q[,6]))/2
  for(i in 1:length(pos)){
    points(pos[i]*scale., h., col="red", cex=6, lwd = 6, pch="|")
  }
}




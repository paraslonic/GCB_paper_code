
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
    print("log scale")
    c. <- tab.$complexity
    #c.log <- log(c.+min(c.[c.!=0]))
    c.log <- log(c.+1)
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

plot_complexity <- function(tab, h, .complexity_mult, scale){

  for(i in 1:nrow(tab)){
    h. <- (tab$complexity[i])*.complexity_mult+h
    polygon(c(tab$start[i],tab$start[i],tab$end[i],tab$end[i])*scale,c(h,h.,h.,h), lwd=0.5)
  }
  x. <- seq(0, max(tab$pos), 1e6)
  y. <- rep(h-15, length(x.))
  text(x.*scale, y., x./1e6, cex = 0.8, srt = 0, col = "gray20")
  for(i in 1:length(x.)){lines(c(x.[i], x.[i])*scale, c(h-5, h), col = "gray20")}
  y.tics <- seq(0,floor(max(tab$complexity)/2),length.out = 2)
  axis(side=2, at=c(h+y.tics*.complexity_mult), lwd = 2, col = "gray70", labels = y.tics, cex.axis=0.8)
}

plot_synteny <- function(coords.file, scale1, scale2, h1, h2){
  col. <- rgb(0.1,0.1,0.1,0.18)
  s <- read.delim(coords.file, head = F)
  colnames(s) <- c("S1","E1","S2","E2","L1","L2","id1","id2")
  s <- subset(s, s$L1 > 5000)
  for(i in 1:nrow(s)){
    polygon(c(s$S1[i]*scale2,s$S2[i]*scale1, s$E2[i]*scale1, s$E1[i]*scale2), c(h2, h1, h1, h2),
            col = col., border=F, lwd = 0.01)
  }
}
























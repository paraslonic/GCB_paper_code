library("data.table")
setwd("/data4/bio/operonTravel_d4/PAPER/code/make_figure1B")
source("functions.r")

complexity <- read.delim("LF82_B2_graph/prob_window_complexity_contig_NC_011993.1.txt", head = F)
colnames(complexity) <- c("gene", "complexity")
gene.info <- fread("graph_genes.sif")
gene.info <- subset(gene.info, gene.info$contig == "NC_011993.1")
complexity <- merge(complexity, gene.info)
complexity <- complexity[order(complexity$start),]


# plot --------------------------------------------------------------------
h0 <- 6.2
h.phaster.1 <- 0
h.phaster.2 <- 1.9
h.pai.1 <- 0
h.pai.2 <- 1.9
h.trna.1 <- 2
h.trna.2 <- 3.9
h.rrna.1 <- 2
h.rrna.2 <- 3.9
h.ess.1 <- 4
h.ess.2 <- 5.9
## colors
trna.color <- "#2ecc71"
rrna.color <- "#8e44ad"
phaster.col <- "#f1c40e"
pai.color <- "#e74c3c"
ess.color <- "#0891eb"


svg("figure1_B.svg", height = 4, width = 7)

plot(complexity$start, complexity$complexity+h0, type="n", pch=3, cex = 0.3, col="gray50", lwd = 1,  
     cex.lab=1.5, yaxt="n", axes=F, ylab = "", xlab = "", ylim = c(0,20+h0), xlim = c(0, 4.8e6))
axis(side=2, at=c(h0,h0+10,h0+20), lwd = 2, col = "gray70", labels = c(0,10,20))
mtext("complexity", side = 2, line=2.5, at=15, cex = 1.5, col  ="gray20")
axis(side=1, at=c(seq(0, 4800000, 1000000), 4.8e6, las = 1, cex = 0.5), lwd = 2, col = "gray70",col.ticks = "gray70")
mtext("chromosome position, bp", side = 1, line=2.5, at=2e6, cex = 1.5, las = 1, col = "gray20")
for(i in 1:nrow(complexity)){
  h. <- complexity$complexity[i] + h0
  polygon(c(complexity$start[i],complexity$start[i],complexity$end[i],complexity$end[i]), 
          c(h0,h.,h.,h0), lwd=0.1, col = "black")
}

# grid --------------------------------------------------------------------
max.x <- max(complexity$end)
lines(c(0, max.x), c(h.ess.1, h.ess.1), col = "gray80")
lines(c(0, max.x), c(h.rrna.1, h.rrna.1), col = "gray80")
lines(c(0, max.x), c(h.pai.1, h.pai.1), col = "gray80")
lines(c(0, max.x), c(h.ess.2, h.ess.2), col = "gray80")
for(i in seq(1, max.x, 20000)){
  lines(c(i, i), c(h.pai.1,h.ess.2), col = "gray80")
}

# phaster -----------------------------------------------------------------
phaster <- read.delim("phaster_summary_lf82.txt")
phaster.regions <- strsplit(as.character(phaster$REGION_POSITION), "-")
for(p in phaster.regions){
  rect(p[1],h.phaster.1, p[2], h.phaster.2,col = phaster.col,border = phaster.col)  
}
mtext("Prophage / ", side=2, line=0, at=h.phaster.1+0.5, cex = 0.8, las = 2, col = "#a8880b")

# pai ---------------------------------------------------------------------
pai <- read.delim("pai", head = F)
for(i in 1:nrow(pai)){
  rect(pai[i, 2],h.pai.1, pai[i,3], h.pai.2,col = pai.color,border = pai.color)  
}
mtext("PAI", side=4, line=-2, at=h.pai.1+0.5, cex = 0.8, las = 2, col = "#b62616")

# tRNA --------------------------------------------------------------------
trna <- read.delim("lf82.trna", head = F)
for(i in 1:nrow(trna)){
  rect(trna[i, 1],h.trna.1, trna[i,2], h.trna.2,col = trna.color,border = trna.color)  
}
mtext("tRNA", side=4, line=-2, at=h.trna.1+0.5, cex = 0.8, las = 2, col = "#1c7d44")

# rRNA --------------------------------------------------------------------
rrna <- read.delim("lf82.rrna", head = F)
for(i in 1:nrow(rrna)){
  rect(rrna[i, 1],h.rrna.1, rrna[i,2], h.rrna.2,col = rrna.color,border = rrna.color)  
}
mtext("rRNA /", side=2, line=0, at=h.rrna.1+0.5, cex = 0.8, las = 2, col = rrna.color)

# essential ---------------------------------------------------------------
ess.pec <- read.delim("essential_pec", head = F)
annotation <- read.delim("lf82_ncbi.bed", head = F)
ess <- subset(annotation, annotation$V5 %in% ess.pec$V1)
for(i in 1:nrow(ess)){
  rect(ess[i, 2],h.ess.1, ess[i,3], h.ess.2,col = ess.color, border = ess.color)  
}
mtext("Essential", side=2, line=0, at=h.ess.1+0.5, cex = 0.8, las = 2, col = ess.color)

# PDU operon --------------------------------------------------------------
pdu.chr <- c(2083448, 2101340)
h.pdu <- 6
lines(pdu.chr, c(h0+h.pdu,h0+h.pdu),  lwd = 1, col = "gray50")
text((pdu.chr[1]+pdu.chr[2])/2, h0+h.pdu+1, "Pdu", adj = 0.5, col = "gray50")

# Ori
## data from http://tubic.org/doric/public/index.php/information/bacteria/ori96000365.html
oric <- (3971662+3972039)/2
lines(c(oric, oric), c(h0, h0+20), col = "#4b38ff", lty = 2)
text(oric, h0+18, "oriC", adj = 0)
# ter <- 1536432
# lines(c(ter, ter), c(h0, h0+20), col = "#4b38ff", lty = 2)

dev.off()

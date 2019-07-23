library("ape")
library(phangorn) 
library("data.table")

setwd("/data4/bio/operonTravel_d4/PAPER/code/make_figure3/")
source("functions.r")

tab_o157h7 <- calc_tab("O157H7_E_graph/prob_window_complexity_contig_NC_002695.1.txt",
                       "O157H7_E_graph/Orthogroups.csv",
                       "GCF_000008865.1_ASM886v1_genomic", FALSE)

tab_iai1 <- calc_tab("IAI1_B1_graph/prob_window_complexity_contig_NC_011741.1.txt",
                     "IAI1_B1_graph/Orthogroups.csv",
                     "GCF_000026265.1_ASM2626v1_genomic", FALSE)


tab_k12 <- calc_tab("K12_graph/prob_window_complexity_contig_NC_000913.3.txt",
                    "K12_graph/Orthogroups.csv",
                    "GCF_000005845.2_ASM584v2_genomic", FALSE)

tab_lf82 <- calc_tab("LF82_B2_graph/prob_window_complexity_contig_NC_011993.1.txt",
                     "LF82_B2_graph/Orthogroups.csv",
                     "GCF_000284495.1_ASM28449v1_genomic", FALSE)

tab_umn026 <- calc_tab("UMN026_D_graph/prob_window_complexity_contig_NC_011751.1.txt",
                        "UMN026_D_graph/Orthogroups.csv",
                      "GCF_000026325.1_ASM2632v1_genomic", FALSE)


#### plot ####
complexity_mult <- 12
h_o157 <- 400
h_iai1 <- 300
h_k12 <- 200
h_lf82 <- 0
h_umn026 <- 100

plot(tab_o157h7$pos, (tab_o157h7$complexity)*complexity_mult+h_o157, ylim=c(0,600), 
     type="n", pch=3, cex = 0.3, col="gray50", lwd = 1, xlab="chromosome position, Mbp", 
     ylab="complexity", cex.lab=1.5, yaxt="n", axes=F,) 

#1
plot_complexity(tab_o157h7, h_o157, complexity_mult*1.2, 1)
scale.f_iai1=max(tab_o157h7$pos)/max(tab_iai1$pos)
plot_synteny("alignments/out_1.coords",1, scale.f_iai1, h_o157, h_iai1)

#2
plot_complexity(tab_iai1, h_iai1, complexity_mult*0.5, scale.f_iai1)
scale.f_k12=max(tab_o157h7$pos)/max(tab_k12$pos)
plot_synteny("alignments/out_k12_iai1", scale.f_iai1, scale.f_k12, h_iai1, h_k12)

#3
plot_complexity(tab_k12, h_k12, complexity_mult*0.8, scale.f_k12)
scale.f_lf82=max(tab_o157h7$pos)/max(tab_lf82$pos)
plot_synteny("alignments/out_umn026_k12", scale.f_k12, scale.f_umn026, h_k12, h_umn026)

#4
plot_complexity(tab_umn026, h_umn026, complexity_mult*0.8, scale.f_umn026)
scale.f_umn026=max(tab_o157h7$pos)/max(tab_umn026$pos)
plot_synteny("alignments/out_lf82_umn026", scale.f_umn026, scale.f_lf82, h_umn026, h_lf82)

#5
plot_complexity(tab_lf82, h_lf82, complexity_mult*0.6, scale.f_lf82)


### 
source("read_phaster_results.r")
plot_phaster <- function(phaster, scale, h){
  for(k in 1:ncol(phaster)){
    polygon(c((phaster[1,k]),(phaster[1,k]), phaster[2,k], phaster[2,k])*scale,
            c(h-8, h, h, h-8), col="#ffa54fc8",border="tan1", lwd=0.5)
    # lines(c(phaster[1,k], phaster[2,k])*scale, c(h-5,h-5), col = "gray20", lty = 3, lwd = 2, pch = 17)
  }
}
plot_phaster(phaster.o157, 1, h_o157-5)
plot_phaster(phaster.iai1, scale.f_iai1, h_iai1-5)
plot_phaster(phaster.k12, scale.f_k12, h_k12-5)
plot_phaster(phaster.lf82, scale.f_lf82, h_lf82-5)
plot_phaster(phaster.umn026, scale.f_umn026, h_umn026-5)

#
legend("topright", "Prophage", fill = "tan1", bty = "n", cex = 1, border = "tan1")

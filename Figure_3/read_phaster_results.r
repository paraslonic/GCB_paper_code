phaster.o157 <- read.delim("phaster/summary_o157.txt", head = T, sep = "\t")
phaster.o157 <- strsplit(as.character(phaster.o157$REGION_POSITION), split = '-')
phaster.o157 <- as.data.frame(phaster.o157)
phaster.o157 <- apply(phaster.o157, c(1,2), as.integer)

phaster.iai1 <- read.delim("phaster/summary_iai1.txt", head = T, sep = "\t")
phaster.iai1 <- strsplit(as.character(phaster.iai1$REGION_POSITION), split = '-')
phaster.iai1 <- as.data.frame(phaster.iai1)
phaster.iai1 <- apply(phaster.iai1, c(1,2), as.integer)

phaster.k12 <- read.delim("phaster/summary_k12.txt", head = T, sep = "\t")
phaster.k12 <- strsplit(as.character(phaster.k12$REGION_POSITION), split = '-')
phaster.k12 <- as.data.frame(phaster.k12)
phaster.k12 <- apply(phaster.k12, c(1,2), as.integer)

phaster.lf82 <- read.delim("phaster/summary_lf82.txt", head = T, sep = "\t")
phaster.lf82 <- strsplit(as.character(phaster.lf82$REGION_POSITION), split = '-')
phaster.lf82 <- as.data.frame(phaster.lf82)
phaster.lf82 <- apply(phaster.lf82, c(1,2), as.integer)

phaster.umn026 <- read.delim("phaster/summary_umn026.txt", head = T, sep = "\t")
phaster.umn026 <- strsplit(as.character(phaster.umn026$REGION_POSITION), split = '-')
phaster.umn026 <- as.data.frame(phaster.umn026)
phaster.umn026 <- apply(phaster.umn026, c(1,2), as.integer)


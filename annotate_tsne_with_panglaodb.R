################################################################################
################################################################################
###### How many cell specific genes are used to define cell type clusters? #####
################################################################################
################################################################################
library(ggplot2)

# ---- inputs / outputs --------------------------------------------------------
tsne_file    <- "linseed.cluster.tsneData.k5.All.v2.csv"  # columns: X (gene ID), tSNE1, tSNE2, marker
annot_file   <- "gene_annotation.csv"                      # columns: PROBEID (gene ID), SYMBOL
panglao_file <- "PanglaoDB_markers_27_Mar_2020.tsv"
out_csv      <- "tsne_annotated_panglaodb.csv"
out_pdf      <- "tSNE_celltype_vs_panglaodb.pdf"

all.k5 <- read.csv(tsne_file)
RNAseq.annot <- read.csv(annot_file)
mtcidx <- match(all.k5$X, RNAseq.annot$PROBEID)
all.k5$"SYMBOL" <- RNAseq.annot[mtcidx, "SYMBOL"]

pDB <- read.delim(panglao_file)

RBC <- pDB[pDB$cell.type=="Erythroid-like and erythroid precursor cells" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(RBC) <- c("SYMBOL","Desc","SS","SP")
NEUT <- pDB[pDB$cell.type=="Neutrophils" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(NEUT) <- c("SYMBOL","Desc","SS","SP")
BSO <- pDB[pDB$cell.type=="Basophils" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(BSO) <- c("SYMBOL","Desc","SS","SP")
EOS <- pDB[pDB$cell.type=="Eosinophils" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(EOS) <- c("SYMBOL","Desc","SS","SP")
TCL <- pDB[pDB$cell.type=="T cells" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(TCL) <- c("SYMBOL","Desc","SS","SP")
BCL <- pDB[pDB$cell.type=="B cells" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(BCL) <- c("SYMBOL","Desc","SS","SP")
NKC <- pDB[pDB$cell.type=="NK cells" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(NKC) <- c("SYMBOL","Desc","SS","SP")
MNC <- pDB[pDB$cell.type=="Monocytes" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(MNC) <- c("SYMBOL","Desc","SS","SP")


DEN <- pDB[pDB$cell.type=="Dendritic cells" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(DEN) <- c("SYMBOL","Desc","SS","SP")
MAC <- pDB[pDB$cell.type=="Macrophages" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(MAC) <- c("SYMBOL","Desc","SS","SP")
PLT <- pDB[pDB$cell.type=="Platelets" & grepl("Hs",pDB$species),c("official.gene.symbol","product.description","sensitivity_human","specificity_human")]
colnames(PLT) <- c("SYMBOL","Desc","SS","SP")

# Regions picked by eye on the embedding
target.Lm <- all.k5[(-25 <= all.k5$tSNE1 & all.k5$tSNE1 <= 10 ) & (-25 <= all.k5$tSNE2 & all.k5$tSNE2 <= 17),]   # Lymphocytes
target.Gr <- all.k5[(15 < all.k5$tSNE1 & all.k5$tSNE1 <= 30 ) & (-8 < all.k5$tSNE2 & all.k5$tSNE2 < 10),]    # Granulocytes
target.RB <- all.k5[(5 <= all.k5$tSNE1 & all.k5$tSNE1 <= 15 ) & (20 <= all.k5$tSNE2 & all.k5$tSNE2 <= 30),] # RBCs

############### Impose clusters onto PDB ##############

sort(table(pDB[pDB$official.gene.symbol %in% target.RB$SYMBOL,"cell.type"]), decreasing = TRUE)[1:10]
sort(table(pDB[pDB$official.gene.symbol %in% target.Lm$SYMBOL,"cell.type"]), decreasing = TRUE)[1:10]
sort(table(pDB[pDB$official.gene.symbol %in% target.Gr$SYMBOL,"cell.type"]), decreasing = TRUE)[1:10]

target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="B cells" & grepl("Hs",pDB$species),"official.gene.symbol"], ]
target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="T cells" & grepl("Hs",pDB$species),"official.gene.symbol"], ]
target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="B cells naive" & grepl("Hs",pDB$species),"official.gene.symbol"], ]
target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="Plasma cells" & grepl("Hs",pDB$species),"official.gene.symbol"], ]
target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="T memory cells" & grepl("Hs",pDB$species),"official.gene.symbol"], ]
target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="B memory cells" & grepl("Hs",pDB$species),"official.gene.symbol"], ]
target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="Neurons" & grepl("Hs",pDB$species),"official.gene.symbol"], ]
target.Lm[target.Lm$SYMBOL %in% pDB[pDB$cell.type=="Interneurons" & grepl("Hs",pDB$species),"official.gene.symbol"], ]


# In target.Gr
sort(table(target.Gr$marker), decreasing = TRUE)
sort(table(pDB[pDB$official.gene.symbol %in% target.Gr$SYMBOL,"cell.type"]), decreasing = TRUE)[1:10]
#sort(table(pDB[pDB$official.gene.symbol %in% target.Gr[!is.na(target.Gr$marker) & target.Gr$marker==6,"SYMBOL"],"cell.type"]),decreasing=TRUE)[1:10]
sort(table(pDB[pDB$official.gene.symbol %in% target.Gr[!is.na(target.Gr$marker) & target.Gr$marker==1,"SYMBOL"],"cell.type"]),decreasing=TRUE)[1:10]
sort(table(pDB[pDB$official.gene.symbol %in% target.Gr[!is.na(target.Gr$marker) & target.Gr$marker==2,"SYMBOL"],"cell.type"]),decreasing=TRUE)[1:10]

# In target.Lm
sort(table(target.Lm$marker), decreasing = TRUE)
sort(table(pDB[pDB$official.gene.symbol %in% target.Lm$SYMBOL,"cell.type"]), decreasing = TRUE)[1:10]
sort(table(pDB[pDB$official.gene.symbol %in% target.Lm[!is.na(target.Lm$marker) & target.Lm$marker==3,"SYMBOL"],"cell.type"]), decreasing=TRUE)[1:10]
sort(table(pDB[pDB$official.gene.symbol %in% target.Lm[!is.na(target.Lm$marker) & target.Lm$marker==4,"SYMBOL"],"cell.type"]), decreasing=TRUE)[1:10]
# sort(table(pDB[pDB$official.gene.symbol %in% target.Lm[!is.na(target.Lm$marker) & target.Lm$marker==3,"SYMBOL"],"cell.type"]), decreasing=TRUE)[1:10]
# sort(table(pDB[pDB$official.gene.symbol %in% target.Lm[!is.na(target.Lm$marker) & target.Lm$marker==2,"SYMBOL"],"cell.type"]), decreasing=TRUE)[1:10]


# In target.RB
sort(table(target.RB$marker), decreasing = TRUE)
sort(table(pDB[pDB$official.gene.symbol %in% target.RB$SYMBOL,"cell.type"]), decreasing = TRUE)[1:10]
# sort(table(pDB[pDB$official.gene.symbol %in% target.RB[!is.na(target.RB$marker) & target.RB$marker==8,"SYMBOL"],"cell.type"]), decreasing = TRUE)[1:10]
sort(table(pDB[pDB$official.gene.symbol %in% target.RB[!is.na(target.RB$marker) & target.RB$marker==5,"SYMBOL"],"cell.type"]), decreasing = TRUE)[1:10]


############### Using panglaoDB genes ###############
# How many Lymphocyte related genes (panglaoDB) used to infer cell type clusters?
TCL[TCL$SYMBOL %in% target.Lm[!is.na(target.Lm$marker),"SYMBOL"],]       # in panglaoDB list
TCL[TCL$SYMBOL %in% target.Lm$SYMBOL,]
target.Lm[target.Lm$SYMBOL %in% TCL$SYMBOL & !is.na(target.Lm$marker), ] # in Linseed cluster list
target.Lm[target.Lm$SYMBOL %in% TCL$SYMBOL, ]
target.Lm[target.Lm$SYMBOL %in% TCL$SYMBOL, "pDB"] <- rep("Tcell",nrow(target.Lm[target.Lm$SYMBOL %in% TCL$SYMBOL,]))

BCL[BCL$SYMBOL %in% target.Lm[!is.na(target.Lm$marker),"SYMBOL"],]       # in panglaoDB list
BCL[BCL$SYMBOL %in% target.Lm$SYMBOL,]
target.Lm[target.Lm$SYMBOL %in% BCL$SYMBOL & !is.na(target.Lm$marker), ] # in Linseed cluster list
target.Lm[target.Lm$SYMBOL %in% BCL$SYMBOL, ]
target.Lm[target.Lm$SYMBOL %in% BCL$SYMBOL, "pDB"] <- rep("Bcell",nrow(target.Lm[target.Lm$SYMBOL %in% BCL$SYMBOL,]))

NKC[NKC$SYMBOL %in% target.Lm[!is.na(target.Lm$marker),"SYMBOL"],]       # in panglaoDB list
NKC[NKC$SYMBOL %in% target.Lm$SYMBOL,] 
target.Lm[target.Lm$SYMBOL %in% NKC$SYMBOL & !is.na(target.Lm$marker), ] # in Linseed cluster list
target.Lm[target.Lm$SYMBOL %in% NKC$SYMBOL, ]
target.Lm[target.Lm$SYMBOL %in% NKC$SYMBOL, "pDB"] <- rep("NKcell",nrow(target.Lm[target.Lm$SYMBOL %in% NKC$SYMBOL,]))

target.Lm[target.Lm$SYMBOL %in% MNC$SYMBOL, ]
target.Lm[target.Lm$SYMBOL %in% PLT$SYMBOL,]
target.Lm[target.Lm$SYMBOL %in% DEN$SYMBOL,]
target.Lm[target.Lm$SYMBOL %in% MAC$SYMBOL,]

target.Lm[target.Lm$SYMBOL %in% BSO$SYMBOL, ]
target.Lm[target.Lm$SYMBOL %in% EOS$SYMBOL, ]
target.Lm[target.Lm$SYMBOL %in% NEUT$SYMBOL, ]
target.Lm[target.Lm$SYMBOL %in% RBC$SYMBOL, ]

# How many Granulocyte related genes (panglaoDB) used to infer cell type clusters?
nrow(NEUT[NEUT$SYMBOL %in% target.Gr$SYMBOL,]) #!
NEUT[NEUT$SYMBOL %in% target.Gr[!is.na(target.Gr$marker),"SYMBOL"],] #!
target.Gr[target.Gr$SYMBOL %in% NEUT$SYMBOL & !is.na(target.Gr$marker), ]
target.Gr[target.Gr$SYMBOL %in% NEUT$SYMBOL,]
target.Gr[target.Gr$SYMBOL %in% NEUT$SYMBOL,"pDB"] <- rep("NEUT",nrow(target.Gr[target.Gr$SYMBOL %in% NEUT$SYMBOL,]))

MNC[MNC$SYMBOL %in% target.Gr$SYMBOL,] #!
MNC[MNC$SYMBOL %in% target.Gr[!is.na(target.Gr$marker),"SYMBOL"],]   #!
target.Gr[target.Gr$SYMBOL %in% MNC$SYMBOL & !is.na(target.Gr$marker), ]
target.Gr[target.Gr$SYMBOL %in% MNC$SYMBOL , ]
target.Gr[target.Gr$SYMBOL %in% MNC$SYMBOL,"pDB"] <- rep("MNC",nrow(target.Gr[target.Gr$SYMBOL %in% MNC$SYMBOL,]))

# BSO[BSO$SYMBOL %in% target.Gr$SYMBOL,"SYMBOL"] #!
# BSO[BSO$SYMBOL %in% target.Gr[!is.na(target.Gr$marker),"SYMBOL"],]       # in panglaoDB list
# target.Gr[target.Gr$SYMBOL %in% BSO$SYMBOL & !is.na(target.Gr$marker), ] # in Linseed cluster list
# target.Gr[target.Gr$SYMBOL %in% BSO$SYMBOL, ]
# target.Gr[target.Gr$SYMBOL %in% BSO$SYMBOL, "pDB"] <- rep("BSO",nrow(target.Gr[target.Gr$SYMBOL %in% BSO$SYMBOL,]))
# 
# EOS[EOS$SYMBOL %in% target.Gr[!is.na(target.Gr$marker),"SYMBOL"],]   #!
# target.Gr[target.Gr$SYMBOL %in% EOS$SYMBOL & !is.na(target.Gr$marker), ]
# target.Gr[target.Gr$SYMBOL %in% EOS$SYMBOL , ]
# target.Gr[target.Gr$SYMBOL %in% EOS$SYMBOL , "pDB"] <- rep("EOS",nrow(target.Gr[target.Gr$SYMBOL %in% EOS$SYMBOL,]))
# 
# target.Gr[target.Gr$SYMBOL %in% PLT$SYMBOL,]
# target.Gr[target.Gr$SYMBOL %in% PLT$SYMBOL , "pDB"] <- rep("PLT",nrow(target.Gr[target.Gr$SYMBOL %in% PLT$SYMBOL,]))
# target.Gr[target.Gr$SYMBOL %in% DEN$SYMBOL,]
# target.Gr[target.Gr$SYMBOL %in% DEN$SYMBOL , "pDB"] <- rep("DEN",nrow(target.Gr[target.Gr$SYMBOL %in% DEN$SYMBOL,]))
# target.Gr[target.Gr$SYMBOL %in% MAC$SYMBOL,]
# target.Gr[target.Gr$SYMBOL %in% MAC$SYMBOL , "pDB"] <- rep("MAC",nrow(target.Gr[target.Gr$SYMBOL %in% MAC$SYMBOL,]))

target.Gr[target.Gr$SYMBOL %in% TCL$SYMBOL,]
target.Gr[target.Gr$SYMBOL %in% BCL$SYMBOL,]
target.Gr[target.Gr$SYMBOL %in% NKC$SYMBOL,]
target.Gr[target.Gr$SYMBOL %in% RBC$SYMBOL,]




# How many RBC related genes (panglaoDB) used to infer cell type clusters?
RBC[RBC$SYMBOL %in% target.RB$SYMBOL,]
RBC[RBC$SYMBOL %in% target.RB[!is.na(target.RB$marker),"SYMBOL"],]
target.RB[target.RB$SYMBOL %in% RBC$SYMBOL & !is.na(target.RB$marker),]
target.RB[target.RB$SYMBOL %in% RBC$SYMBOL,]
target.RB[target.RB$SYMBOL %in% RBC$SYMBOL,"pDB"] <- rep("RBC",nrow(target.RB[target.RB$SYMBOL %in% RBC$SYMBOL,]))

target.RB[target.RB$SYMBOL %in% PLT$SYMBOL,]
target.RB[target.RB$SYMBOL %in% DEN$SYMBOL,]
target.RB[target.RB$SYMBOL %in% MAC$SYMBOL,]
target.RB[target.RB$SYMBOL %in% TCL$SYMBOL,]
target.RB[target.RB$SYMBOL %in% BCL$SYMBOL,]
target.RB[target.RB$SYMBOL %in% NKC$SYMBOL,]
target.RB[target.RB$SYMBOL %in% MNC$SYMBOL, ]
target.RB[target.RB$SYMBOL %in% BSO$SYMBOL, ]
target.RB[target.RB$SYMBOL %in% EOS$SYMBOL, ]
target.RB[target.RB$SYMBOL %in% NEUT$SYMBOL, ]

target.Mod <- rbind(target.RB, target.Gr, target.Lm)
nrow(target.Mod[!is.na(target.Mod$marker) & !is.na(target.Mod$pDB),])
target.Mod[!is.na(target.Mod$marker) & !is.na(target.Mod$pDB),"Both.pDB.Linseed"] <- rep(TRUE, nrow(target.Mod[!is.na(target.Mod$marker) & !is.na(target.Mod$pDB),])) 
nrow(target.Mod[is.na(target.Mod$Both.pDB.Linseed),])
nrow(target.Mod) - nrow(target.Mod[!is.na(target.Mod$marker) & !is.na(target.Mod$pDB),])
target.Mod[is.na(target.Mod$Both.pDB.Linseed),"Both.pDB.Linseed"] <- rep(FALSE,nrow(target.Mod[is.na(target.Mod$Both.pDB.Linseed),]))

write.csv(target.Mod, out_csv)

unique(target.Mod$pDB)

target.Mod$pDB <- factor(target.Mod$pDB, levels=c(NA,"RBC","Tcell","Bcell","NKcell","NEUT","MNC"))
target.Mod$marker <- factor(target.Mod$marker, levels=c(NA,1,2,3,4,5))
p<-ggplot(data=target.Mod, aes(x=tSNE1, y=tSNE2, color=marker, size=!is.na(pDB))) +
  geom_point() + theme_bw(base_size=8) +
  #scale_shape_manual(values = c(1, 19)) +
  scale_size_manual(values=c(0.5, 2)) +
  geom_text(data=target.Mod[!is.na(target.Mod$pDB),], aes(x=tSNE1+1.2, y=tSNE2, label=paste(pDB, SYMBOL,sep=".")), size=1.5, color="black") +
  # scale_color_manual(values=c("gray","blue","forestgreen"))+
  # scale_color_manual(values=c("gray","blue","forestgreen","red"))+
  #scale_color_manual(values=c("gray",rgb(248,118,109),rgb(163,165,0),"deepskyblue","forestgreen","red"))+
  #scale_color_manual(values=c("gray","blue","red","green"))+
  theme(legend.key.size = unit(0.1, "in"), aspect.ratio = 1) +
  guides(color=guide_legend(title="CellType"))
ggsave(p, file=out_pdf, height=4.5, width=6)

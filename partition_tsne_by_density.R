#
# Given tSNE plot data, divide clusters at density troughs and make tables for IPA
#
library(ggplot2)

# ---- inputs / outputs --------------------------------------------------------
tsne_file <- "linseed.cluster.tsneData.k5.All.v1.csv"  # columns: X (gene ID), tSNE1, tSNE2, marker
de_file   <- "DE_results_all_genes.csv"                # full DE table; column X holds the same gene IDs
out_dir   <- "tables_for_IPA"
pad       <- 2                                         # padding added beyond the outermost points

k5   <- read.csv(tsne_file)
DEGs <- read.csv(de_file)

# ---- troughs in the density of each axis ------------------------------------
d  <- density(k5$tSNE1)
d2 <- density(k5$tSNE2)

# d$x[which(diff(sign(diff(d$y)))<0)+1]   # Peaks

d$x[which(diff(sign(diff(d$y)))>0)+1]   # Trough tSNE1
d2$x[which(diff(sign(diff(d2$y)))>0)+1] # Trough tSNE2

t1 <- c(min(k5$tSNE1) - pad, d$x[which(diff(sign(diff(d$y)))>0)+1],  max(k5$tSNE1) + pad)
t2 <- c(min(k5$tSNE2) - pad, d2$x[which(diff(sign(diff(d2$y)))>0)+1], max(k5$tSNE2) + pad)

# ---- assign genes to grid regions -------------------------------------------
temp_df <- lapply(c(1:(length(t1)-1)), function(i) {k5[(t1[i] < k5$tSNE1 & k5$tSNE1 <= t1[i+1]) , ]})
k5_list <- lapply(c(1:length(temp_df)), function(x) {
    lapply(c(1:(length(t2)-1)), function(i){temp_df[[x]][(t2[i] < temp_df[[x]]$tSNE2 & temp_df[[x]]$tSNE2 <= t2[i+1]),]})
  }
)

k5_list <- unlist(k5_list, recursive=F)

# ---- one DE table per region ------------------------------------------------
dir.create(out_dir, showWarnings = FALSE)
for (i in c(1:length(k5_list))) {
  write.csv(DEGs[DEGs$X %in% k5_list[[i]]$X,], file=file.path(out_dir, paste0("cluster_", i, ".csv")))
}

# ---- guide plot with break lines and region labels --------------------------
cl <- expand.grid(t1, t2)
names(cl) <- c("X","Y")
cl <- cl[cl$X != t1[length(t1)],]
cl <- cl[cl$Y != t2[length(t2)],]
cl <- cl[order(cl$X),]
cl$text <- paste0("C",c(1:nrow(cl)))

p <- ggplot(k5, aes(x=tSNE1, y=tSNE2)) +
  geom_point() +
  geom_text(data=cl, aes(x=X, y=Y, label=text), hjust=-0.5, vjust=-1, color="green4") +
  geom_vline(xintercept = d$x[which(diff(sign(diff(d$y)))>0)+1], color="blue") +
  geom_hline(yintercept = d2$x[which(diff(sign(diff(d2$y)))>0)+1], color="blue")

ggsave(p, file=file.path(out_dir, "cluster_guidance.pdf"), width=8, height=6.5)

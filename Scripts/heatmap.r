library(pheatmap)
df <- read.csv("Index_hebing.txt", sep='\t', row.names = 1)
pdf("Index_hebing.pdf",width=12,height=12)


breaks <- c(0, seq(0.0000001, 0.9999999, length.out = 98), 1)
# print(breaks)

# gradient color
gradient_colors <- colorRampPalette(c("blue", "red"))(97)

# combine
colors <- c("gray", gradient_colors, "white")

pheatmap(df, color = colors, breaks= breaks, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames=TRUE)
dev.off()

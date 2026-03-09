# Fst heatmap 
library(viridis)
pops <- c("UT", "PT", "ML","MO")
n_pops <- length(pops)
fst_matrix <- matrix(NA, nrow=n_pops, ncol=n_pops, dimnames=list(pops, pops))
fst_matrix["MO", "UT"] <- 0.173
fst_matrix["PT", "UT"] <- 0.212
fst_matrix["ML", "UT"] <- 0.202
fst_matrix["MO","PT"] <- 0.143
fst_matrix["MO", "ML"] <- 0.12
fst_matrix["ML", "PT"] <- 0.149
diag(fst_matrix) <- 0
fst_for_plot <- fst_matrix
fst_for_plot[upper.tri(fst_for_plot, diag = FALSE)] <- NA
layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(4, 1), heights = c(1, 1))
cols <- rocket(100, direction = -1)

par(mar = c(4, 3, 5, 2))

image(
  x = 1:n_pops,
  y = 1:n_pops,
  z = t(fst_for_plot[n_pops:1, ]),
  col = cols,
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "",
  zlim = c(0, 1)
)

axis(1, at = 1:n_pops, labels = colnames(fst_for_plot), las = 2)
axis(2, at = 1:n_pops, labels = rev(rownames(fst_for_plot)), las = 1)
box()

for (i in 1:n_pops) {
  for (j in 1:n_pops) {
    if (!is.na(fst_for_plot[i, j])) {
      text(x = j, y = n_pops - i + 1, labels = round(fst_for_plot[i, j], 3), cex = 1.2)
    }
  }
}

par(mar = c(4, 1, 5, 4))
legend_scale <- seq(0, 1, length.out = 100)

image(y = legend_scale, z = t(legend_scale), col = cols, axes = FALSE, xlab = "", ylab = "")
box()
mtext("Fst", side = 3, line = 0.5)
axis(4, at = seq(0, 1, by = 0.2), las = 1)

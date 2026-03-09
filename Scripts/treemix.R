library(phytools)
tree <- read.tree('../Data/treemix.tree')
plotTree(tree)

source('plotting_funcs.R')
plot_tree('../Data/migration/migration.3')

vioplot(pi_data_list,
        col = "transparent",   # IMPORTANT
        border = "black",
        drawRect = T,
        ylab = "Nucleotide diversity (π)",
        xlab = "Population")


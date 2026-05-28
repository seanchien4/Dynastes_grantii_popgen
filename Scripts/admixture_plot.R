library(viridis)
dat <- read.table('~/Desktop/Projects/Dynastes_grantii_pop/Data/admixture/genetic_prop_2.txt')
order.list <- dat$V1[grepl('PT',dat$V1)]
order.list <- c(order.list,dat$V1[grepl('ML',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('MG',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('MC',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgrC',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgrP',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgrR',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DGS',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgSR',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgUT',dat$V1)])
# Reorder columns
order.dat <- dat[match(order.list, dat$V1), ]
order.dat <- t(order.dat)
gsub('DgrP', 'P', order.dat) -> order.dat
gsub('DgrR', 'R', order.dat) -> order.dat
gsub('DgUT', 'UT', order.dat) -> order.dat
gsub('DgSR', 'UN', order.dat) -> order.dat
gsub('DGS0', 'M', order.dat) -> order.dat
gsub('DgrCO0', 'C', order.dat) -> order.dat

gsub('DgPT', 'PT', order.dat) -> order.dat
gsub('DgML', 'ML', order.dat) -> order.dat
gsub('DgMC', 'MC', order.dat) -> order.dat
gsub('DgMG', 'MG', order.dat) -> order.dat

colnames(order.dat) <- order.dat[1,]
order.dat <- order.dat[-1,]

rownames(order.dat) <- c('Pop1','Pop2')
cols <- viridis(length(rownames(order.dat)))
# par(mfrow = c(2, 2))
# cex.names = 0.4
cex.names = 0.8
par(mar = c(2,2.5,2,1))
barplot(order.dat, 
        beside = FALSE,       
        col = c(cols[1], cols[2], cols[3], cols[4], cols[5]), 
        #main = "K=2", 
        xlab = "",           
        ylab = "",
        yaxt = "n",
        las =2, 
        cex.axis = 1,
        #names.arg = rep("", ncol(order.dat)),
        cex.names = cex.names,las = 2,mgp = c(3, 0, -0.6)
)          
axis(2,at = seq(0, 1.0, by = 0.1),
     labels = seq(0, 100, by = 10))
title(main = paste0("K = ", 2), adj = 0, line = 0.5)
# text(x = 22,
#      y = -0.13,
#      labels = "Pop1",
#      srt = 45,
#      xpd = TRUE)
# text(x = 41,
#      y = -0.13,
#      labels = "Pop2",
#      srt = 45,
#      xpd = TRUE)

# other K 

dat_list <- list()
for (i in c(3,4,5,6)) {
  file_name <- paste0("~/Desktop/Projects/Dynastes_grantii_pop/Data/admixture/genetic_prop_", i, ".txt")
  dat <- read.table(file_name, header = F)
  order.list <- dat$V1[grepl('ML',dat$V1)]
  order.list <- c(order.list,dat$V1[grepl('PT',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('MG',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('MC',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('DgrC',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('DgrP',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('DgrR',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('DGS',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('DgSR',dat$V1)])
  order.list <- c(order.list,dat$V1[grepl('DgUT',dat$V1)])
  # Reorder columns
  order.dat <- dat[match(order.list, dat$V1), ]
  order.dat <- t(order.dat)
  gsub('DgrP', 'P', order.dat) -> order.dat
  gsub('DgrR', 'R', order.dat) -> order.dat
  gsub('DgUT', 'UT', order.dat) -> order.dat
  gsub('DgSR', 'UN', order.dat) -> order.dat
  gsub('DGS0', 'M', order.dat) -> order.dat
  gsub('DgrCO0', 'C', order.dat) -> order.dat
  
  gsub('DgPT', 'PT', order.dat) -> order.dat
  gsub('DgML', 'ML', order.dat) -> order.dat
  gsub('DgMC', 'MC', order.dat) -> order.dat
  gsub('DgMG', 'MG', order.dat) -> order.dat
  
  colnames(order.dat) <- order.dat[1,]
  order.dat <- order.dat[-1,]
  rownames(order.dat) <- paste0("Pop", seq_len(i))
  cols <- viridis(length(rownames(order.dat)))
  barplot(order.dat, 
          beside = FALSE,       
          col = cols[seq_len(i)], 
          main = '', 
          xlab = "",           
          ylab = "",
          yaxt = "n",
          las =2, 
          cex.axis = 1,
          #names.arg = rep("", ncol(order.dat)),
          cex.names = cex.names,las = 2,mgp = c(3, 0, -0.6)
          
  )
  axis(2,at = seq(0, 1.0, by = 0.1),
       labels = seq(0, 100, by = 10))
  title(main = paste0("K = ", i), adj = 0, line = 0.5)
}

# K = 4 
dat <- read.table('../Data/admixture/genetic_prop_4.txt')
order.list <- dat$V1[grepl('PT',dat$V1)]
order.list <- c(order.list,dat$V1[grepl('ML',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('MG',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('MC',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgrC',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgrP',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgrR',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DGS',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgSR',dat$V1)])
order.list <- c(order.list,dat$V1[grepl('DgUT',dat$V1)])
# Reorder columns
order.dat <- dat[match(order.list, dat$V1), ]
order.dat <- t(order.dat)
gsub('DgrP', 'P', order.dat) -> order.dat
gsub('DgrR', 'R', order.dat) -> order.dat
gsub('DgUT', 'UT', order.dat) -> order.dat
gsub('DgSR', 'UN', order.dat) -> order.dat
gsub('DGS0', 'M', order.dat) -> order.dat
gsub('DgrCO0', 'C', order.dat) -> order.dat

gsub('DgPT', 'PT', order.dat) -> order.dat
gsub('DgML', 'ML', order.dat) -> order.dat
gsub('DgMC', 'MC', order.dat) -> order.dat
gsub('DgMG', 'MG', order.dat) -> order.dat

colnames(order.dat) <- order.dat[1,]
order.dat <- order.dat[-1,]

rownames(order.dat) <- c('Pop1','Pop2', 'Pop3','Pop4')
cols <- viridis(length(rownames(order.dat)))
barplot(order.dat, 
        beside = FALSE,       
        col = c(cols[1], cols[3], cols[2], cols[4]), 
        xlab = "",           
        ylab = "",
        yaxt = "n",
        las =2, 
        cex.axis = 1,
        cex.names = 0.6,las = 2,mgp = c(3, 0, -0.6)
)          
axis(2,at = seq(0, 1.0, by = 0.1),
     labels = seq(0, 100, by = 10))
title(main = paste0("K = ", 4), adj = 0, line = 0.5)

# vcf 
par(mfrow = c(2, 2))
# mean_dp <- read.table('../Data/depth_per_individual.tsv', header = T)
hist(mean_dp$MEAN_DEPTH, breaks = 20,
     xlab = 'Mean depth per individual', ylab = 'Frequency', main = '',
     xlim = c(0,30), ylim = c(0,4))
title(main = "(A)", adj = 0, line = 0.5)
# mean_dp_site <- read.table('../Data/depth_per_site.tsv', header = T )
hist(mean_dp_site$MEAN_DEPTH,breaks = 30000,
     xlab = 'Mean depth per site', ylab = 'Frequency', main = '',
     xlim = c(0,40), ylim = c(0,1450000))
mean <- mean(mean_dp_site$MEAN_DEPTH, na.rm = T)
sd <- sd(mean_dp_site$MEAN_DEPTH, na.rm = TRUE)
mean + 2*sd
title(main = "(B)", adj = 0, line = 0.5)
# abline(v=mean + 2*sd, col = 'red')

# missing_ind <- read.table('../Data/missing_per_individual.tsv', header = T)
hist(missing_ind$F_MISS, breaks = 30,
     xlab = 'Missing data per individual(%)', ylab = 'Frequency', main = '',
     xlim = c(0,0.4), ylim = c(0,15),xaxt = "n")
axis(1,at = seq(0, 0.4, by = 0.1),
     labels = seq(0, 40, by = 10))
title(main = "(C)", adj = 0, line = 0.5)

# missing_site <- read.table('../Data/missing_per_site.tsv', header = T)
hist(missing_site$F_MISS,breaks = 40,
     xlab = 'Missing data per site(%)', ylab = 'Frequency', main = '',
     xaxt = "n")
axis(1,at = seq(0, 1.0, by = 0.1),
     labels = seq(0, 100, by = 10))
title(main = "(D)", adj = 0, line = 0.5)



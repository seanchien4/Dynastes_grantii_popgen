##### GONE2 ####
library(viridis)
library(RColorBrewer)
library(scico)
cols <- viridis(4)
cols_fade <- viridis(4, alpha = 0.2)
dat <- read.table('~/Desktop/Projects/Dynastes_grantii_pop/Data/GONE_Ne/Mogollon/Mogollon_1_GONE2_Ne')
colnames(dat) <- c(dat[1,1], dat[1,2])
dat <- dat[-1,]
allNe <- list()
allNe[[1]] <- as.numeric(dat$Ne_diploids)
generations <- dat$Generation
bs <- list.files(path = "~/Desktop/Projects/Dynastes_grantii_pop/Data/GONE_Ne/Mogollon/",
                 pattern = "_GONE2_Ne$",
                 full.names = TRUE)
for(i in 2:length(bs)){
  df <- read.table(bs[i])
  colnames(df) <- c(df[1,1], df[1,2])
  df <- df[-1,]
  allNe[[i]] <- as.numeric(df$Ne_diploids)
}
# convert to matrix: rows = bootstrap, cols = generation
NeMat <- do.call(rbind, allNe)
# calculate median and 95% CI
medNe <- apply(NeMat, 2, median)
lowerNe <- apply(NeMat, 2, quantile, 0.025)
upperNe <- apply(NeMat, 2, quantile, 0.975)
par(mar = c(4,4,1,1))
plot(generations, medNe, type='n',
     xaxt = "n", yaxt = "n",
     ylab = '',
     xlab = "",
     xlim = c(0,150), ylim = c(0,20000)
)
title(ylab = expression(italic(N[e])), line = 2.5)
title(xlab = "Yeasrs ago", line = 2.5)
# add shaded area for 95% CI
polygon(c(generations, rev(generations)),
        c(upperNe, rev(lowerNe)),
        col=cols_fade[1], border=NA)
# overlay median line
lines(generations, medNe, col=cols[1], lwd=2)

g <- seq(from = 0, to = 150, by = 20)
axis(side = 1, at = g, labels = g*2)
ticks_y  <- c(1000,1e+04, 2e+04, 3e+04, 4e+04)
labels_y <- c("1k","10K", "20K", "30K", "40K")
axis(2, at = ticks_y, labels = labels_y, las = 1)

# other populations
subs <- c('Lemmon','Portal','Utah')
for (s in 1:length(subs)) {
  path <- paste0("~/Desktop/Projects/Dynastes_grantii_pop/Data/GONE_Ne/", subs[s], "/")
  
  bs <- list.files(path = path,
                   pattern = "_GONE2_Ne$",
                   full.names = TRUE)
  allNe <- list()
  for(i in 1:length(bs)){
    df <- read.table(bs[i])
    colnames(df) <- c(df[1,1], df[1,2])
    df <- df[-1,]
    allNe[[i]] <- as.numeric(df$Ne_diploids)
  }
  # convert to matrix: rows = bootstrap, cols = generation
  NeMat <- do.call(rbind, allNe)
  # calculate median and 95% CI
  medNe <- apply(NeMat, 2, median)
  lowerNe <- apply(NeMat, 2, quantile, 0.025)
  upperNe <- apply(NeMat, 2, quantile, 0.975)
  # add shaded area for 95% CI
  polygon(c(generations, rev(generations)),
          c(upperNe, rev(lowerNe)),
          col=cols_fade[s+1], border=NA)
  # for(j in 1:nrow(NeMat)){
  #   lines(generations, NeMat[j,], col=cols_fade[s], lwd=0.5)
  # }
  
  # overlay median line
  lines(generations, medNe, col=cols[s+1], lwd=2)
}

# s=800
# labels <- c('Mogollon Rim','Mt. Lemmon','Chiricahua Mt.','Utah')
# for (i in 1:4){
#   x = 117
#   y = 21000-i*s
#   points(x,y, bg = cols[i], pch = 21 )
#   text(x,y, labels[i], pos = 4)
# }


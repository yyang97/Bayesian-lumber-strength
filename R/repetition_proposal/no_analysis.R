library(ggplot2)
library(ggpubr)
non_1 <- readRDS("sim_non_update.rds")
non_2 <- readRDS("sim_non_2.rds")




mean_df <- as.data.frame(rbind(non_1$mean[1:896,], non_2$mean))


# mean 
colMeans(rbind(non_1$mean[1:896,],non_2$mean))


# sd 
colMeans(rbind(non_1$sd[1:896,],non_2$sd))


theta_true <- c(45,5.5,13,1,.7,1,1,1,1,1,1,.7)

upper <- rbind(non_1$upper[1:896,],non_2$upper)
lower <- rbind(non_1$lower[1:896,],non_2$lower)

coverage <- matrix(NA, ncol = 12, nrow = 896)
for (icov in 1:896){
  coverage[icov,] <- (upper[icov,] - theta_true> 0)& (theta_true - lower[icov,] > 0)
}
colMeans(coverage)
colMeans(lower)
colMeans(upper)


# muX
breaks <- pretty(range(mean_df[,1]), n = nclass.FD(mean_df[,1]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_muX <- ggplot(mean_df, aes(x=V1)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 45, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,1], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(mu[X]))

# muY 
breaks <- pretty(range(mean_df[,2]), n = nclass.FD(mean_df[,2]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_muY <- ggplot(mean_df, aes(x=V2)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 5.5, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,2], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(mu[Y]))

# sigmaX
breaks <- pretty(range(mean_df[,3]), n = nclass.FD(mean_df[,3]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_sigmaX <- ggplot(mean_df, aes(x=V3)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 13, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,3], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(sigma[X]))

# sigmaY 
breaks <- pretty(range(mean_df[,4]), n = nclass.FD(mean_df[,4]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_sigmaY <- ggplot(mean_df, aes(x=V4)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 1, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,4], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(sigma[Y]))

# rho
breaks <- pretty(range(mean_df[,5]), n = nclass.FD(mean_df[,5]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_rho <- ggplot(mean_df, aes(x=V5)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.7, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,5], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(rho))


# alpha_R20
breaks <- pretty(range(mean_df[,6]), n = nclass.FD(mean_df[,6]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_alphaR20 <- ggplot(mean_df, aes(x=V6)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.9, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,6], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(alpha[R20]))


# alpha_R40
breaks <- pretty(range(mean_df[,7]), n = nclass.FD(mean_df[,7]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_alphaR40 <- ggplot(mean_df, aes(x=V7)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.8, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,7], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(alpha[R40]))



# alpha_R60
breaks <- pretty(range(mean_df[,8]), n = nclass.FD(mean_df[,8]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_alphaR60 <- ggplot(mean_df, aes(x=V8)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.7, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,8], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(alpha[R60]))




# alpha_T20
breaks <- pretty(range(mean_df[,9]), n = nclass.FD(mean_df[,9]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_alphaT20 <- ggplot(mean_df, aes(x=V9)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.9, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,9], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(alpha[T20]))


# alpha_T40
breaks <- pretty(range(mean_df[,10]), n = nclass.FD(mean_df[,10]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_alphaT40 <- ggplot(mean_df, aes(x=V10)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.8, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,10], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(alpha[T40]))



# alpha_T60
breaks <- pretty(range(mean_df[,11]), n = nclass.FD(mean_df[,11]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_alphaT60 <- ggplot(mean_df, aes(x=V11)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.7, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,11], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(alpha[T60]))




# eta
breaks <- pretty(range(mean_df[,12]), n = nclass.FD(mean_df[,12]), min.n = 2)
bwidth <- breaks[2]-breaks[1]

mean_eta <- ggplot(mean_df, aes(x=V12)) + geom_histogram(binwidth = bwidth)+ 
  geom_vline(xintercept = 0.7, color = "red") + 
  geom_vline(xintercept = quantile(mean_df[,12], probs = c(.025,0.975)), color = "green") +
  labs(x=expression(eta)) + theme(legend.position="bottom")


ggarrange(mean_muX, mean_muY, mean_sigmaX, mean_sigmaY, 
          mean_rho, mean_alphaR20, mean_alphaR40, mean_alphaR60,
          mean_alphaT20, mean_alphaT40, mean_alphaT60, mean_eta,
          ncol = 3, nrow = 4,  common.legend = TRUE)



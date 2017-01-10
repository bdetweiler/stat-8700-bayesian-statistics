# STAT 8700 Final Question 4
# Brian Detweiler
# Thursday, December 15th

library(rjags)
library(knitr)
library(reshape2)
library(ggplot2)
library(dplyr)
source("DBDA2Eprograms/DBDA2E-utilities.R")
set.seed(0)

###############################################################################
# 4 (a)
###############################################################################

# -1 Away win, 0 tie, 1 Home win
epl$Result <- sign(epl$Home.Goals - epl$Away.Goals)

teams <- unique(c(epl$Home.Team, epl$Away.Team))

mean.goals <- mean(c(epl$Home.Goals, epl$Away.Goals))
mean.home.goals <- mean(c(epl$Home.Goals))
sd.goals <- sd(c(epl$Home.Goals, epl$Away.Goals))

epl.M <- matrix(data = NA, nrow = 20, ncol = 20)

for (i in 1:length(epl[,1])) {
  epl.M[epl$Home.ID[i], epl$Away.ID[i]] <- epl$Home.Goals[i]
  epl.M[epl$Away.ID[i], epl$Home.ID[i]] <- epl$Away.Goals[i]
}

a <- c(rep(NA, 20))
d <- c(rep(NA, 20))
for (i in 1:20) {
  epl.mean <- epl %>% filter(Home.ID == i) %>%
    group_by(Home.ID) %>%
    summarise(a_mean = mean(Home.Goals), d_mean = mean(Away.Goals))

  a[i] <- epl.mean$a_mean - mean.goals
  d[i] <- mean.goals - epl.mean$d_mean
}

gam <- sum(a) - sum(d)

fileName <- "Final.4.a"

modelString ="
model {
  for(i in 1:120) {
    H[i] ~ dpois(lambda[HomeTeam[i], AwayTeam[i]])
    A[i] ~ dpois(theta[HomeTeam[i], AwayTeam[i]])
  }
  
  for(i in 1:n_teams) {
    for(j in 1:n_teams) {
      lambda[i, j] <- exp(mu + a[i] - d[j] + gamma) 
      theta[i, j] <- exp(mu + a[j] - d[i])
    }
  }
  
  for(j in 1:(n_teams - 1)) {
    a[j] ~ dnorm(group_attack, group_tau)
    d[j] ~ dnorm(group_defense, group_tau)
  }
  a[n_teams] <- -sum(a[1:(n_teams - 1)])
  d[n_teams] <- -sum(d[1:(n_teams - 1)])

  gamma ~ dgamma(0.01, 0.01)
  
  rank.a <- rank(a)
  rank.d <- rank(d)
  
  for (k in 1:n_teams) {
    for (m in 1:n_teams) {
      #ranks.a[k, m] ~ dbern(pa[k, m])
      #ranks.d[k, m] ~ dbern(pd[k, m])
      #pa[k, m] <- equals(k, rank.a[m])
      #pd[k, m] <- equals(k, rank.d[m])
      ranks.a[k, m] <- equals(k, rank.a[m])
      ranks.d[k, m] <- equals(k, rank.d[m])
    }
  }
  
  mu ~ dnorm(0, 0.0625)

  group_attack ~ dnorm(0, 0.0625)
  group_defense ~ dnorm(0, 0.0625)
  group_tau <- 1 / pow(group_sigma, 2)
  group_sigma ~ dunif(0, 3)
}
"


writeLines(modelString, con=fileName)

data.list <- list(H = epl$Home.Goals, 
                  A = epl$Away.Goals, 
                  HomeTeam = epl$Home.ID,
                  AwayTeam = epl$Away.ID,
                  n_teams = length(teams))


epl.model = jags.model(file=fileName, 
                       data=data.list,
                       n.chains=4)

update(epl.model, n.iter=50000)

epl.samples <- coda.samples(epl.model, 
                            variable.names = c("a",
                                               "d",
                                               "gamma",
                                               "rank.a",
                                               "rank.d",
                                               "ranks.a",
                                               "ranks.d"),
                            n.iter = 50000, 
                            thin = 2)

diagMCMC(epl.samples)

epl.dic <- dic.samples(model = epl.model, n.iter = 20000, thin = 50)

epl.dic

epl.samples.M <- as.matrix(epl.samples)

smry <- summary(epl.samples)


###############################################################################
# 4 (b)
###############################################################################


a.means <- smry$quantiles[1:20, 3]
d.means <- smry$quantiles[21:40, 3]
gamma.means <- smry$quantiles[41, 3]
rank.a.means <- smry$quantiles[42:61, 3]
rank.d.means <- smry$quantiles[62:81, 3]

# sort(a.means)
# teams.with.names$Home.Team[teams.with.names$Home.ID == 9]

best.a <- which(a.means == max(a.means))
worst.a <- which(a.means == min(a.means))

best.d <- which(d.means == max(d.means))
worst.d <- which(d.means == min(d.means))

best.a.team.name <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == best.a])
best.d.team.name <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == best.d])

worst.a.team.name <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == worst.a])
worst.d.team.name <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == worst.d])

###############################################################################
# 4 (c)
###############################################################################

best.rank <- which(rank.d.means == max(rank.d.means))
worst.rank <- which(rank.d.means == min(rank.d.means))

best.team.name.rank1 <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == best.rank[[1]]])
best.team.name.rank2 <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == best.rank[[2]]])
worst.team.name.rank <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == worst.rank[[1]]])

###############################################################################
# 4 (d)
###############################################################################

best.rank <- which(rank.a.means == max(rank.a.means))
worst.rank <- which(rank.a.means == min(rank.a.means))

best.team.name.rank1 <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == best.rank[[1]]])
worst.team.name.rank <- as.character(teams.with.names$Home.Team[teams.with.names$Home.ID == worst.rank[[1]]])

###############################################################################
# 4 (e)
###############################################################################

ranks.df <- data.frame(teams.with.names)
ranks.df$rank1 <- 0
ranks.df$rank2 <- 0
ranks.df$rank3 <- 0
ranks.df$rank4 <- 0
ranks.df$rank5 <- 0
ranks.df$rank6 <- 0
ranks.df$rank7 <- 0
ranks.df$rank8 <- 0
ranks.df$rank9 <- 0
ranks.df$rank10 <- 0
ranks.df$rank11 <- 0
ranks.df$rank12 <- 0
ranks.df$rank13 <- 0
ranks.df$rank14 <- 0
ranks.df$rank15 <- 0
ranks.df$rank16 <- 0
ranks.df$rank17 <- 0
ranks.df$rank18 <- 0
ranks.df$rank19 <- 0
ranks.df$rank20 <- 0

for (i in 1:20) {
  for (j in 1:20) {
    ranks.df[i, paste0("rank", j)] <- mean(epl.samples.M[,paste0("ranks.d[", i, ",", j, "]")])
  }
}

ranks.sub.df <- ranks.df[,-1]

ranks.melt <- melt(ranks.sub.df, id.vars = "Home.Team")
ggplot(ranks.melt, aes(x=variable, y=value)) +
  geom_bar(stat="identity") +
  facet_wrap(~Home.Team) +
  scale_x_discrete("Rank", labels = c("1", "", "", "",
                                      "5", "", "", "", "",
                                      "10", "", "", "", "",
                                      "15", "", "", "", "",
                                      "20"))

ranks.sub.df.rank1 <- ranks.sub.df %>% select(Home.Team, rank1) %>% arrange(-rank1)
ranks.sub.df.rank20 <- ranks.sub.df %>% select(Home.Team, rank20) %>% arrange(-rank20)

worst.d.team.name <- as.character(ranks.sub.df.rank1[1, 1])
best.d.team.name <- as.character(ranks.sub.df.rank20[1, 1])

ggplot(ranks.sub.df.rank1, aes(x=Home.Team, y=rank1)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Simulations with Worst Defense", x="Team", y="Proportion of Top Rank")

ggplot(ranks.sub.df.rank20, aes(x=Home.Team, y=rank20)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Simulations with Best Defense", x="Team", y="Proportion of Top Rank")


###############################################################################
# 4 (f)
###############################################################################

ranks.df <- data.frame(teams.with.names)
ranks.df$rank1 <- 0
ranks.df$rank2 <- 0
ranks.df$rank3 <- 0
ranks.df$rank4 <- 0
ranks.df$rank5 <- 0
ranks.df$rank6 <- 0
ranks.df$rank7 <- 0
ranks.df$rank8 <- 0
ranks.df$rank9 <- 0
ranks.df$rank10 <- 0
ranks.df$rank11 <- 0
ranks.df$rank12 <- 0
ranks.df$rank13 <- 0
ranks.df$rank14 <- 0
ranks.df$rank15 <- 0
ranks.df$rank16 <- 0
ranks.df$rank17 <- 0
ranks.df$rank18 <- 0
ranks.df$rank19 <- 0
ranks.df$rank20 <- 0

for (i in 1:20) {
  for (j in 1:20) {
    ranks.df[i, paste0("rank", j)] <- mean(epl.samples.M[,paste0("ranks.a[", i, ",", j, "]")])
  }
}

ranks.sub.df <- ranks.df[,-1]

ranks.melt <- melt(ranks.sub.df, id.vars = "Home.Team")
ggplot(ranks.melt, aes(x=variable, y=value)) +
  geom_bar(stat="identity") +
  facet_wrap(~Home.Team) +
  scale_x_discrete("Rank", labels = c("1", "", "", "",
                                      "5", "", "", "", "",
                                      "10", "", "", "", "",
                                      "15", "", "", "", "",
                                      "20"))


ranks.sub.df.rank1 <- ranks.sub.df %>% select(Home.Team, rank1) %>% arrange(desc(rank1))
ranks.sub.df.rank20 <- ranks.sub.df %>% select(Home.Team, rank20) %>% arrange(desc(rank20))

best.a.team.name <- as.character(ranks.sub.df.rank20[1, 1])
worst.a.team.name <- as.character(ranks.sub.df.rank1[1, 1])

ggplot(ranks.sub.df.rank20, aes(x=Home.Team, y=rank20)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Simulations with Best Offense", x="Team", y="Proportion of Top Rank")

ggplot(ranks.sub.df.rank1, aes(x=Home.Team, y=rank1)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Simulations With Worst Offense", x="Team", y="Proportion of Bottom Rank")

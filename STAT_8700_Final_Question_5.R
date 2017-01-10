# STAT 8700 Final Question 5
# Brian Detweiler
# Thursday, December 15th

library(rjags)
library(reshape2)
library(ggplot2)
library(dplyr)
source("DBDA2Eprograms/DBDA2E-utilities.R")
set.seed(0)

epl <- read.csv('epl.csv', header=TRUE)
teams.with.names <- epl %>% 
  select(Home.ID, Home.Team) %>%
  arrange(Home.ID) %>%
  unique()

##################################################################
# 5 (a)
##################################################################

# -1 Away win, 0 tie, 1 Home win
epl$Result <- sign(epl$Home.Goals - epl$Away.Goals)

epl.M <- matrix(data = 0,
                nrow = 20,
                ncol = 20)

for (i in 1:20) {
  for (j in 1:20) {
    tmp <- epl %>% filter(Home.ID == i, Away.ID == j)
    epl.M[i, j] <- nrow(tmp)
    if (i == j) {
      epl.M[i, j] <- 1
    }
  }
}



##################################################################
# 5 (b)
##################################################################

epl.M.Home <- matrix(data = 0,
                     nrow = 20,
                     ncol = 20)

epl.M.Away <- matrix(data = 0,
                     nrow = 20,
                     ncol = 20)

epl$Result.Home <- epl$Result + 1
epl$Result.Home[epl$Result.Home == 2] <- 3

# Invert the Home results
epl$Result.Away <- epl$Result * -1
epl$Result.Away <- epl$Result.Away + 1
epl$Result.Away[epl$Result.Away == 2] <- 3

for (i in 1:20) {
  for (j in 1:20) {
    tmp <- epl %>% filter(Home.ID == i, Away.ID == j)
    if (nrow(tmp) > 0) {
      epl.M.Home[i, j] <- tmp$Result.Home
      epl.M.Away[j, i] <- tmp$Result.Away
    }
  }
}

##################################################################
# 5 (c)
##################################################################
chelsea.home.sum <- sum(epl.M.Home[teams.with.names$Home.ID[teams.with.names$Home.Team == "Chelsea"], ])
chelsea.away.sum <- sum(epl.M.Away[teams.with.names$Home.ID[teams.with.names$Home.Team == "Chelsea"], ])

##################################################################
# 5 (d)
##################################################################

teams <- unique(c(epl$Home.Team, epl$Away.Team))

epl.M.Home <- matrix(data = -1,
                     nrow = 20,
                     ncol = 20)

epl.M.Away <- matrix(data = -1,
                     nrow = 20,
                     ncol = 20)

epl$Result.Home <- epl$Result + 1
epl$Result.Home[epl$Result.Home == 2] <- 3

# Invert the Home results
epl$Result.Away <- epl$Result * -1
epl$Result.Away <- epl$Result.Away + 1
epl$Result.Away[epl$Result.Away == 2] <- 3

for (i in 1:20) {
  for (j in 1:20) {
    tmp <- epl %>% filter(Home.ID == i, Away.ID == j)
    if (nrow(tmp) > 0) {
      epl.M.Home[i, j] <- tmp$Result.Home
      epl.M.Away[j, i] <- tmp$Result.Away
    }
  }
}

fileName <- "Final.5.d.jags"

modelString ="
model {
  for(i in 1:120) {
    H[i] ~ dpois(lambda[HomeTeam[i], AwayTeam[i]])
    A[i] ~ dpois(theta[HomeTeam[i], AwayTeam[i]])
  }

  for(i in 1:n_teams) {
    for(j in 1:n_teams) {
      SimGoalsHome[i, j] ~ dpois(lambda[HomeTeam[i], AwayTeam[j]])
      SimGoalsAway[i, j] ~ dpois(theta[HomeTeam[j], AwayTeam[i]])

      # Possible scores are 0, 1, or 3
      SimPointsHome[i, j] <- combinedScoreHome[i, j] - tieHome[i, j]

      # step is 1 if x >= 0, so it'll either be 3 or 0
      combinedScoreHome[i, j] <- step(SimGoalsHome[i, j] - SimGoalsAway[i, j]) * 3 
      # equals is 1 if x == y, so this will either be 2 or 0
      tieHome[i, j] <- equals(SimGoalsHome[i, j], SimGoalsAway[i, j]) * 2

      # Possible scores are 0, 1, or 3
      SimPointsAway[i, j] <- combinedScoreAway[i, j] - tieAway[i, j]

      # step is 1 if x >= 0, so it'll either be 3 or 0
      combinedScoreAway[i, j] <- step(SimGoalsAway[i, j] - SimGoalsHome[i, j]) * 3 
      # equals is 1 if x == y, so this will either be 2 or 0
      tieAway[i, j] <- equals(SimGoalsAway[i, j], SimGoalsHome[i, j]) * 2

      # Data will be -1 if the game has not been played
      FinalPointsHome[i, j] <- step(HomePointsData[i, j]) * HomePointsData[i, j] + equals(HomePointsData[i, j], -1) * SimPointsHome[i, j]

      FinalPointsAway[i, j] <- step(AwayPointsData[i, j]) * AwayPointsData[i, j] + equals(AwayPointsData[i, j], -1) * SimPointsAway[i, j]

    }

    TotalPoints[i] <- sum(FinalPointsHome[i, 1:n_teams]) + sum(FinalPointsAway[i, 1:n_teams]) 
  }

  LeagueRank <- rank(TotalPoints)

  for(k in 1:n_teams) {
    for(m in 1:n_teams) {
      LeagueRanks[k, m] <- equals(k, LeagueRank[m])
    }
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
                  HomePointsData = epl.M.Home,
                  AwayPointsData = epl.M.Away,
                  n_teams = length(teams))


epl.model = jags.model(file=fileName, 
                       data=data.list,
                       n.chains=4)

update(epl.model, n.iter=50000)


epl.samples <- coda.samples(epl.model, 
                            variable.names = c("SimGoalsHome",
                                               "SimGoalsAway",
                                               "SimPointsHome",
                                               "SimPointsAway",
                                               "FinalPointsHome",
                                               "FinalPointsAway",
                                               "TotalPoints",
                                               "LeagueRank",
                                               "LeagueRanks",
                                               "rank.d",
                                               "ranks.a",
                                               "ranks.d"),
                            n.iter = 50000, 
                            thin = 2)

diagMCMC(epl.samples)

epl.DIC <- dic.samples(model = epl.model, n.iter = 200000, thin = 50)
epl.DIC

epl.samples.M <- as.matrix(epl.samples)

smry <- summary(epl.samples)

##################################################################
# 5 (h)
##################################################################

total.points <- rep(0, 20)
for (i in 1:20) {
  total.points[i] <- smry$quantiles[paste0("TotalPoints[", i, "]"), 3]
}

teams.with.names$TotalPoints <- total.points
teams.final <- teams.with.names %>% arrange(desc(TotalPoints))

total.points <- rep(0, 20)
for (i in 1:20) {
  total.points[i] <- smry$quantiles[paste0("TotalPoints[", i, "]"), 3]
}

teams.with.names$TotalPoints <- total.points
teams.final <- teams.with.names %>% arrange(desc(TotalPoints))

##################################################################
# 5 (i)
##################################################################

league.rank <- rep(0, 20)
for (i in 1:20) {
  league.rank[i] <- smry$quantiles[paste0("LeagueRank[", i, "]"), 3]
}

teams.with.names$LeagueRank <- league.rank
league.rank.final <- teams.with.names %>% arrange(desc(LeagueRank))

##################################################################
# 5 (j)
##################################################################

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
    ranks.df[i, paste0("rank", j)] <- mean(epl.samples.M[,paste0("LeagueRanks[", i, ",", j, "]")])
  }
}
ranks.sub.df <- ranks.df[,-1]
ranks.sub.df <- ranks.sub.df[, -2]
ranks.sub.df <- ranks.sub.df[, -2]

ranks.melt <- melt(ranks.sub.df, id.vars = "Home.Team")
ggplot(ranks.melt, aes(x=variable, y=value)) +
  geom_bar(stat="identity") +
  facet_wrap(~Home.Team) +
  scale_x_discrete("Rank", labels = c("1", "", "", "",
                                      "5", "", "", "", "",
                                      "10", "", "", "", "",
                                      "15", "", "", "", "",
                                      "20"))

ranks.sub.df.rank1 <- ranks.sub.df %>% select(Home.Team, rank1) %>% arrange(rank1)
ranks.sub.df.rank20 <- ranks.sub.df %>% select(Home.Team, rank20) %>% arrange(rank20)

ggplot(ranks.sub.df.rank1, aes(x=Home.Team, y=rank1)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Simulations with Lowest Rank", x="Team", y="Proportion of Bottom Rank")

ggplot(ranks.sub.df.rank20, aes(x=Home.Team, y=rank20)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Simulations with Top Rank", x="Team", y="Proportion of Top Rank")

##################################################################
# 5 (k)
##################################################################

team.top.4 <- c()
for (i in 1:20) {
  team <- ranks.df$Home.Team[i]
  percent <- ranks.df$rank20[ranks.df$Home.Team == team] +
    ranks.df$rank19[ranks.df$Home.Team == team] +
    ranks.df$rank18[ranks.df$Home.Team == team] +
    ranks.df$rank17[ranks.df$Home.Team == team]
  
  team.top.4 <- c(team.top.4, percent)
}

ranks.df$PercentTop4 <- team.top.4

ranks.df %>% arrange(desc(PercentTop4)) %>% select(Home.Team, PercentTop4)

##################################################################
# 5 (l)
##################################################################

team.bottom.3 <- c()
for (i in 1:20) {
  team <- ranks.df$Home.Team[i]
  percent <- ranks.df$rank1[ranks.df$Home.Team == team] +
    ranks.df$rank2[ranks.df$Home.Team == team] +
    ranks.df$rank3[ranks.df$Home.Team == team]
  
  team.bottom.3 <- c(team.bottom.3, percent)
}

ranks.df$PercentBottom3 <- team.bottom.3

ranks.df %>% arrange(desc(PercentBottom3)) %>% select(Home.Team, PercentBottom3)

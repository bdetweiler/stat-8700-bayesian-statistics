
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


library(quantmod)
library(PerformanceAnalytics)
library(quadprog)

# pulling adjusted close prices for the 5 largest sector ETFs + SPY as benchmark
# date range: Jan 2010 - Dec 2024 (15 years of monthly data)
tickers <- c("XLK", "XLF", "XLY", "XLV", "XLI", "SPY")
getSymbols(tickers, from = "2010-01-01", to = "2024-12-31", src = "yahoo")

prices <- merge(Ad(XLK), Ad(XLF), Ad(XLY), Ad(XLV), Ad(XLI), Ad(SPY))
colnames(prices) <- tickers

# convert to monthly -- last trading day of each month
# monthly data cuts down on short-term noise, standard in portfolio research
prices_m <- to.monthly(prices, indexAt = "lastof", OHLC = FALSE)
ret_m <- Return.calculate(prices_m, method = "discrete")
ret_m <- na.omit(ret_m)
colnames(ret_m) <- tickers

# annualized returns and covariance matrix
mu_annual <- colMeans(ret_m) * 12
Sigma_annual <- cov(ret_m) * 12
vol_annual <- sqrt(diag(Sigma_annual))

# drop SPY from portfolio construction -- benchmark only
mu_5 <- mu_annual[1:5]
Sigma_5 <- Sigma_annual[1:5, 1:5]

# three portfolios varying only in XLK weight, remainder split evenly
#           XLK    XLF     XLY     XLV     XLI
w_low  <- c(0.10,  0.225,  0.225,  0.225,  0.225)
w_med  <- c(0.25,  0.1875, 0.1875, 0.1875, 0.1875)
w_high <- c(0.50,  0.125,  0.125,  0.125,  0.125)

# portfolio return: w'mu, volatility: sqrt(w'Sigma*w)
ret_low  <- as.numeric(t(w_low)  %*% mu_5)
ret_med  <- as.numeric(t(w_med)  %*% mu_5)
ret_high <- as.numeric(t(w_high) %*% mu_5)

sd_low  <- as.numeric(sqrt(t(w_low)  %*% Sigma_5 %*% w_low))
sd_med  <- as.numeric(sqrt(t(w_med)  %*% Sigma_5 %*% w_med))
sd_high <- as.numeric(sqrt(t(w_high) %*% Sigma_5 %*% w_high))

# Sharpe ratios (rf = 0 for simplicity -- noted as a limitation)
sharpe_low  <- ret_low  / sd_low
sharpe_med  <- ret_med  / sd_med
sharpe_high <- ret_high / sd_high

ret_spy    <- as.numeric(mu_annual["SPY"])
sd_spy     <- as.numeric(vol_annual["SPY"])
sharpe_spy <- ret_spy / sd_spy

summary_table <- data.frame(
  Portfolio  = c("Low Tech (10%)", "Med Tech (25%)", "High Tech (50%)", "SPY Benchmark"),
  XLK_Weight = c("10%", "25%", "50%", "0%"),
  Ann_Return = round(c(ret_low, ret_med, ret_high, ret_spy), 4),
  Ann_Vol    = round(c(sd_low,  sd_med,  sd_high,  sd_spy),  4),
  Sharpe     = round(c(sharpe_low, sharpe_med, sharpe_high, sharpe_spy), 4)
)
print(summary_table)

# Monte Carlo simulation -- 10,000 paths, 5-year horizon
library(MASS)
set.seed(123)

# converting to log returns for compounding
# log returns sum correctly over time; simple returns don't
logret_m_5 <- log(1 + ret_m[, 1:5])
mu_log    <- colMeans(logret_m_5)
Sigma_log <- cov(logret_m_5)

n_sims   <- 10000
T_months <- 60

simulate_5y <- function(w) {
  results <- numeric(n_sims)
  for (s in 1:n_sims) {
    # multivariate normal draws keep historical correlations between ETFs intact
    sim_returns <- mvrnorm(n = T_months, mu = mu_log, Sigma = Sigma_log)
    portfolio_log_ret <- sim_returns %*% w
    results[s] <- exp(sum(portfolio_log_ret)) - 1
  }
  return(results)
}

sim_low  <- simulate_5y(w_low)
sim_med  <- simulate_5y(w_med)
sim_high <- simulate_5y(w_high)

summarize_sim <- function(x, label) {
  data.frame(
    Portfolio  = label,
    Mean_5Y    = round(mean(x), 4),
    Median_5Y  = round(median(x), 4),
    P_Loss     = round(mean(x < 0), 4),
    VaR_5pct   = round(quantile(x, 0.05), 4),
    Best_95pct = round(quantile(x, 0.95), 4)
  )
}

mc_results <- rbind(
  summarize_sim(sim_low,  "Low Tech (10%)"),
  summarize_sim(sim_med,  "Med Tech (25%)"),
  summarize_sim(sim_high, "High Tech (50%)")
)
print(mc_results)

# efficient frontier via quadratic programming
# minimizes portfolio variance for each target return level
n_points    <- 100
target_rets <- seq(min(mu_5) + 0.001, max(mu_5) - 0.001, length.out = n_points)
ef_ret <- numeric(n_points)
ef_sd  <- numeric(n_points)

for (i in 1:n_points) {
  Dmat <- 2 * Sigma_5
  dvec <- rep(0, 5)
  Amat <- cbind(rep(1, 5), mu_5, diag(5))
  bvec <- c(1, target_rets[i], rep(0, 5))
  
  sol <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = 2),
    error = function(e) NULL
  )
  
  if (!is.null(sol)) {
    w_opt     <- sol$solution
    ef_ret[i] <- as.numeric(t(w_opt) %*% mu_5)
    ef_sd[i]  <- as.numeric(sqrt(t(w_opt) %*% Sigma_5 %*% w_opt))
  }
}

valid  <- ef_sd > 0
ef_ret <- ef_ret[valid]
ef_sd  <- ef_sd[valid]

# plot 1: efficient frontier with portfolio overlays
plot(ef_sd, ef_ret,
     type = "l", lwd = 2, col = "black",
     xlim = c(0.12, 0.22), ylim = c(0.12, 0.20),
     xlab = "Annual Volatility (Risk)",
     ylab = "Annual Expected Return",
     main = "Efficient Frontier")

points(sd_low,  ret_low,  pch = 19, col = "blue",      cex = 1.8)
points(sd_med,  ret_med,  pch = 19, col = "darkgreen",  cex = 1.8)
points(sd_high, ret_high, pch = 19, col = "red",        cex = 1.8)
points(sd_spy,  ret_spy,  pch = 17, col = "gray40",     cex = 1.8)

legend("bottomright",
       legend = c("Efficient Frontier", "Low Tech (10%)", "Med Tech (25%)",
                  "High Tech (50%)", "SPY Benchmark"),
       col    = c("black", "blue", "darkgreen", "red", "gray40"),
       lty    = c(1, NA, NA, NA, NA),
       pch    = c(NA, 19, 19, 19, 17),
       lwd    = c(2, NA, NA, NA, NA),
       cex    = 0.8,
       bty    = "n")

# plot 2: monte carlo return distributions
plot(density(sim_low),
     main = "Distribution of 5-Year Portfolio Returns",
     xlab = "5-Year Cumulative Return",
     ylab = "Density",
     col  = "blue", lwd = 2,
     xlim = c(-1, 10))

lines(density(sim_med),  col = "darkgreen", lwd = 2)
lines(density(sim_high), col = "red",       lwd = 2)
abline(v = 0, lty = 2, col = "gray50")

legend("topright",
       legend = c("Low Tech (10%)", "Med Tech (25%)", "High Tech (50%)"),
       col    = c("blue", "darkgreen", "red"),
       lwd = 2, bty = "n")

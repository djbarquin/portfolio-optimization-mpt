# Portfolio Optimization: Tech Sector Concentration Analysis

Tested whether increasing allocation to the technology sector (XLK) produces 
better risk-adjusted returns in a diversified equity portfolio. Built and ran 
this in R using Modern Portfolio Theory and Monte Carlo simulation on 15 years 
of monthly ETF data.

## What I built

- Three portfolios with varying XLK weights (10%, 25%, 50%), with the remaining 
  weight split evenly across XLF, XLY, XLV, and XLI
- MPT framework to calculate annualized returns, volatility, and Sharpe ratios 
  for each portfolio vs. SPY as a benchmark
- 10,000-path Monte Carlo simulation over a 5-year horizon using multivariate 
  normal draws via mvrnorm, which preserves historical correlations between ETFs
- Efficient frontier via quadratic programming to see how each portfolio compares 
  to the theoretical optimum

## Results

| Portfolio | XLK Weight | Ann. Return | Ann. Volatility | Sharpe Ratio |
|-----------|------------|-------------|-----------------|--------------|
| Low Tech  | 10%        | 14.71%      | 15.45%          | 0.952        |
| Med Tech  | 25%        | 15.46%      | 15.42%          | 1.003        |
| High Tech | 50%        | 16.72%      | 15.71%          | 1.064        |
| SPY       | 0%         | 14.37%      | 14.54%          | 0.989        |

Higher tech concentration produced better risk-adjusted returns across the board, 
with volatility staying almost flat across all three portfolios (15.42% to 15.71%). 
The high-tech portfolio also had the lowest simulated probability of loss over 
5 years at 1.51%.

One interesting finding: SPY came closest to the efficient frontier out of all 
four portfolios, meaning a simple index fund was more risk-efficient than any of 
the hand-built allocations. Makes sense in hindsight since manually picking 
weights can't match what formal optimization does.

## Limitations

- The 2010 to 2024 sample period was one of the strongest stretches for tech in 
  market history, so results would likely look very different over 2000 to 2010
- Simulation assumes normally distributed returns; real returns have fat tails 
  and negative skewness, so true downside risk is probably understated
- Risk-free rate was set to zero for simplicity rather than using the actual 
  3-month Treasury yield
- Portfolio weights are held fixed over the full 5-year horizon, which isn't 
  realistic for an actively managed portfolio

## Stack

R: quantmod, PerformanceAnalytics, quadprog, MASS

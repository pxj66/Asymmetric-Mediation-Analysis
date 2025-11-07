# Asymmetric-Mediation-Analysis
R code for Robust Mediation Analysis with Asymmetric Loss
File list:
  1) rlmskew.R: Robust linear regression for skewed data
  2) Basic Functions.R: Basic functions required by rlmskew.R
  3) Example.R: An example for using the rlmskew() function
  4) Simulation:
     i. MSE.R: runing procedure for mean square error for efficiency
     ii. Testing.R: runing procedure for type I error rate and statistical power, Sobel-type and percentile bootstrap confidence intervals
     c) DGPs.R: data generation progress
     d) LocalInfRobustMediation.R: local influence (LI) and robust transformation (RT) approaches from Zu and Yuan (2010, MBR) Local influence and robust procedures for mediation analysis.
     e) Method.R: two functions for bootstrapping, different approach, Sobel-type confidence interval
  5) Real Data Analysis:
     a) LPACOG.rda: Data for light physical activity and cognitive function study
     b) LPACOG-V0.rda: Data for light physical activity and cognitive function study when perceived control is divided into personal mastery and perceived constaints.
     c) PerceivedControlCognitiveDeclines.R: Apply the proposed approach and testing indirect effect
     d) Skewness.R: Compute sample skewness and kurtosis and KS test and symmetric test from Milosevic and Obradociv (2016) compare them with standard normal distribution
     

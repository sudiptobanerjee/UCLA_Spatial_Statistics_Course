##Clear memory
rm(list=ls())

##Load package "surveil"
library(surveil)

##Load a cancer dataset and a standard population data
data(cancer)
cancer_df = subset(cancer, grepl("1999", Year))
data(standard)
pop_df = standard

##The surveil package fits a Bayesian model to estimate disease rates for each age-group and then computes the standardized age-adjusted rate. We will later learn about these Bayesian models. For now, we will demonstrate a direct standardization using a simple empirical estimate of the age-specific cancer incidence rates.  

cancer_rate = cancer_df$Count/cancer_df$Population 
standardized_total_pop = sum(pop_df$standard_pop) ##This sums up to 1 million
pop_rate = pop_df$standard_pop

standardized_incidence_rate = sum(cancer_rate*pop_rate)


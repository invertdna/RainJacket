# RainJacket

R functions for analyzing ddPCR fluorescence observations by training a logistic regression on control samples, using the resulting model to generate a positivity probability for each droplet, and using those probabilities to derive a concentration estimate for each sample.

## Rationale

Often with ddPCR, there is a question about how to set the positivity threshold for a given droplet. That is: was that droplet positive or not? How do I make this decision in a non-arbitrary way?  Here, I train a simple logistic regression on control data, and then apply that model across unknown (experimental) samples.

## Dependencies

RainJacket requires the following R packages:

tidyverse
cluster
rstanarm



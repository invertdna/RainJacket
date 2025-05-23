#core functions for analyzing ddPCR fluorescence observations by training a logistic regression on control samples,
# using the resulting model to generate a positivity probability for each droplet,
# and using those probabilities to derive a concentration estimate for each sample
# RPK 2024 / 2025


#get control data
get_rainJacket_controls <- function(csv_folder, #folder containing only CSV files to be analyzed
                       pos_control_idx, #numeric index (between 1 and the number of files in csv_folder) for samples to use as positive controls
                       neg_control_idx #numeric index (between 1 and the number of files in csv_folder) for samples to use as negative controls
){

  require(tidyverse)
  require(cluster)

  shift <- function(x){(x - quantile(x, 0.05))}
  
  files <- list.files(csv_folder, full.names = T)
  
  POS <- pos_control_idx
  
  NTCs <- neg_control_idx
  
  
  amps <- list()
  amps_shift <- list()
  amps_clust <- list()
  for (i in 1:length(files)){
    amps[[i]] <- read.csv(files[i], skip = 3)[,1] #assumes the BioRad format in which the first three lines of the CSV are comments
    amps_shift[[i]] <- shift(amps[[i]]) #shifts each sample so that 0 level is set at nearly the lowest fluor observed in the sample (this uses the 5th percentile as the zero level)
  }
  
  for (i in 1:length(POS)){
    amps_clust[[i]] <- data.frame(
      f = unlist(amps_shift[[POS[i]]]),
      z = cluster::clara(amps_shift[[POS[i]]], k = 2, stand = T, cluster.only = T) - 1 #classify positive controls using k-medoids clustering algorithm
    )
  }
  pos_controls = do.call(rbind, amps_clust) #combine positive control samples into one df
  
  #define neg controls as not containing target molecule, and combine these into one df
  neg_controls <- data.frame(
    f = unlist(amps_shift[NTCs]),
    z = 0)
  
  #combine all controls and associated pos/neg classifications
  all_controls <- rbind(neg_controls, pos_controls)
  
  return(all_controls)
  
}

# all_controls <- get_rainJacket_controls(csv_folder = here::here("example_data"),
#                                         pos_control_idx = c(9,10),
#                                         neg_control_idx = c(7,8))

  
  
#fit logistic model

fit_rainJacket_model <- function(control_data,
                                 iterations = 1000){
  
  require(tidyverse)
  require(rstanarm)
  options(mc.cores = parallel::detectCores())

  all_controls <- control_data
  
  suppressWarnings(
    m <- stan_glm(z ~ f, family = "binomial", 
                  data = all_controls,
                  chains = 4, iter = iterations)
  )
  
  
  return(m)
  
}

# m <- fit_rainJacket_model(all_controls,
#                           iterations = 500)

#apply model to experimental samples and output

apply_rainJacket <- function(csv_folder, #folder containing only CSV files to be analyzed, each in BioRad format
                             fitted_model, #fitted rstanarm model that is the output of fit_rainJacket_model()
                             template_aliquot_vol = 2, #volume in uL
                             reaction_vol = 22, #volume in uL
                             droplet_vol = 0.00085 #volume in uL
                             ){
  
  require(tidyverse)
  require(rstanarm)
  options(mc.cores = parallel::detectCores())

  m <- fitted_model 
    
  shift <- function(x){(x - quantile(x, 0.05))}
  
  files <- list.files(csv_folder, full.names = T)
  
  amps <- list()
  amps_shift <- list()
  amps_clust <- list()
  for (i in 1:length(files)){
    amps[[i]] <- read.csv(files[i], skip = 3)[,1] #assumes the BioRad format in which the first three lines of the CSV are comments
    amps_shift[[i]] <- shift(amps[[i]]) #shifts each sample so that 0 level is set at nearly the lowest fluor observed in the sample (this uses the 5th percentile as the zero level)
  }
  
  samples_out <- data.frame(matrix(NA, nrow = length(amps), ncol = 9))
  for (i in 1:length(amps)){
    
    print(paste0("Sample ", i))
    
    g <- posterior_epred(m, newdata = data.frame(f = amps_shift[[i]]), type = "response")
    
    samples_out[i,1] <- i
    samples_out[i,2] <- round(mean(rowSums(g)),0) #colSums(g) #expected number positive droplets
    n_hi <- quantile(rowSums(g), .975)
    n_lo <- quantile(rowSums(g), .025)
    
    
    samples_out[i,3] <- length(amps_shift[[i]]) #total droplets
    
    samples_out[i,4] <- samples_out[i,2]/samples_out[i,3] #expected fraction positive
    p_hi <- n_hi/samples_out[i,3]
    p_lo <- n_lo/samples_out[i,3]
    minDetectable <- 1/length(amps_shift[[i]]) #min possible detected concentration = 1 molecule 
    
    
    samples_out[i,5] <- -log(1-samples_out[i,4]) #lambda = copies/droplet
    lambda_hi <- -log(1-p_hi)
    lambda_lo <- -log(1-p_lo)
    
    samples_out[i,6] <- round((samples_out[i,5]/droplet_vol)*(reaction_vol/template_aliquot_vol),2) #expected template copies/uL
    
    
    samples_out[i,7] <- round((lambda_hi/droplet_vol)*(reaction_vol/template_aliquot_vol),2) #hi 95% template copies/uL
    samples_out[i,8] <- round((lambda_lo/droplet_vol)*(reaction_vol/template_aliquot_vol),2) #lo 95% template copies/uL
    
    samples_out[i,9] <- ifelse(samples_out[i,4] > minDetectable, 1, 0) #indicator for mean detection above LOD
  }
  names(samples_out) <- c("sample", "positive", "total", "fraction_positive", "lambda", "template_copies_uL", "hi", "lo", "above_LOD")
  samples_out <- cbind(filename = basename(files), samples_out)
  
  return(samples_out)
  
}

# samples_out <- apply_rainJacket(csv_folder = here::here("example_data"),
#                                 fitted_model = m,
#                                 template_aliquot_vol = 2, #volume in uL
#                                 reaction_vol = 22, #volume in uL
#                                 droplet_vol = 0.00085 #volume in uL
#                                 )


#plot controls w 95CI

plot_rainJacket_controls <- function(modelfit, observations) {
  
  
  m <- modelfit # from rstanarm
  #observations should be a df of two columns, named f and z, containing fluorescence values and a 1/0 indicator for positivity, respectively
  
  params <- summary(m)[1:2]
  inv_logit <- function(x) {
    1 / (1 + exp(-x))
  }
  
  # Define a function
  mod_best_fit <- function(x) {
    inv_logit(params[1] + params[2] * x)
  }
  
  xbar <- seq(
    min(observations$f),
    max(observations$f),
    10
  )
  # ybar <- posterior_epred(m, newdata = data.frame(f = xbar), type = "response")
  # ybar_means <- colMeans(ybar)
  # ybar_hi <- apply(ybar, 2, quantile, 0.975)
  # ybar_lo <- apply(ybar, 2, quantile, 0.025)
  #
  # df <- data.frame(
  #   xbar,
  #   ybar_means,
  #   ybar_hi,
  #   ybar_lo
  # )
  
  base <-
    ggplot() +
    xlim(min(xbar), max(xbar))
  
  return(
    base +
      geom_function(fun = function(x) inv_logit(params[1] + params[2] * x)) +
      geom_point(data = observations[sample(nrow(observations), 1e4), ], aes(x = f, y = z), alpha = 0.2, size = 2) +
      xlab("Fluorescence Relative to Baseline") +
      ylab("Probability of Being True Positive") +
      ggtitle("Control Data") +
      theme_bw()
  )
}

# plot_rainJacket_controls(modelfit = m,
#                          observations = all_controls)


#plot arbitrary samples against model w 95CI

plot_rainJacket_singleSample <- function(modelfit,
                                         csv_folder,
                                         sample_idx){ #within csv folder, give index number of relevant sample to plot
  
  require(tidyverse)
  
  m <- modelfit
  # params <- summary(m)[1:2]
  
  shift <- function(x){(x - quantile(x, 0.05))}
  inv_logit <- function(x) {
    1 / (1 + exp(-x))
  }
  
  # mod_best_fit <- function(x) {
  #   inv_logit(params[1] + params[2] * x)
  # }
  
  files <- list.files(csv_folder, full.names = T)
  
  amps <- list()
  amps_shift <- list()
  
    amps <- read.csv(files[sample_idx], skip = 3)[,1] #assumes the BioRad format in which the first three lines of the CSV are comments
    amps_shift <- shift(amps) #shifts each sample so that 0 level is set at nearly the lowest fluor observed in the sample (this uses the 5th percentile as the zero level)
  
  
    ybar <- posterior_linpred(m, newdata = data.frame(f = amps_shift), type = "response")
    ybar_means <- colMeans(ybar) %>% inv_logit()
    ybar_hi <- apply(ybar, 2, quantile, 0.75) %>% inv_logit()
    ybar_lo <- apply(ybar, 2, quantile, 0.25) %>% inv_logit()

    df <- data.frame(
      x = amps_shift,
      prob = ybar_means,
      prob_hi = ybar_hi,
      prob_lo = ybar_lo
    )
    
    rm(ybar) #for memory management
    
    return(
      df %>% 
        ggplot(aes(x = x, y = prob)) +
        geom_point(alpha = 0.2, size = 2) +
        geom_segment(aes(x = x, xend = x,
                         y = prob_lo, yend = prob_hi)) +
        xlab("Fluorescence Relative to Baseline") +
        ylab("Probability of Being True Positive") +
        ggtitle("Sample Summary") +
        theme_bw()
    )  

}

# plot_rainJacket_singleSample(modelfit = m, 
#                              csv_folder = here::here("example_data"), 
#                              sample_idx = 4)


#plot summary of all samples 

plot_rainJacket_samples <- function(samples_out){ #dataframe resulting from apply_rainJacket() function

  
  require(tidyverse)
  
  p <- samples_out %>% 
    mutate(Sample_Number = 1:nrow(.)) %>% 
    ggplot(aes(x = Sample_Number,
               y = template_copies_uL)) +
    geom_point() +
    geom_segment(aes(x = Sample_Number, xend = Sample_Number,
                     y = lo, yend = hi)) +
    xlab("Sample Number") + ylab("Copies per uL template") +
    theme_bw()
  
  return(p)
  
}

# samples_out %>% 
#   plot_rainJacket_samples()






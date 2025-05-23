#workflow for analyzing ddPCR fluorescence observations via rainJacket
library(here)
library(tidyverse)

#source core functions
source(here("code/rainJacket.R"))


#load data and format controls
all_controls <- get_rainJacket_controls(csv_folder = here::here("example_data"),
                                        pos_control_idx = c(9,10),
                                        neg_control_idx = c(7,8))

#fit logistic regression model
m <- fit_rainJacket_model(all_controls,
                          iterations = 500)

#apply model to each file in data folder
samples_out <- apply_rainJacket(csv_folder = here::here("example_data"),
                                fitted_model = m,
                                template_aliquot_vol = 2, #volume in uL
                                reaction_vol = 22, #volume in uL
                                droplet_vol = 0.00085 #volume in uL
)

#plot control data
plot_rainJacket_controls(modelfit = m,
                         observations = all_controls)

#plot summary of concentrations of all samples
samples_out %>%
  plot_rainJacket_samples()

#if desired, plot individual sample data
plot_rainJacket_singleSample(modelfit = m, 
                             csv_folder = here::here("example_data"), 
                             sample_idx = 4)


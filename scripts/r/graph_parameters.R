#!/usr/bin/env Rscript

#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. April 2020                       ###
#############################################################################
### Graphical parameters for ggplot                                       ###
### Publication quality plots                                             ### 
#############################################################################

# https://github.com/JoseEspinosa/crispr_flow/blob/master/lib/graph_parameters.R

base_size <- 12

dailyInt_theme <- theme_update (
  axis.text.x = element_text (size = base_size * 1.5),
  axis.text.y = element_text (size = base_size * 1.5),
  axis.title.x = element_text (size=base_size * 1.5, face="bold"),
  axis.title.y = element_text (size=base_size * 1.5, angle = 90, face="bold"),
  #strip.text.x = element_text (size=base_size * 1.3, face="bold"),#facet titles size 
  #strip.text.y = element_text (size=base_size * 1.3, face="bold", angle=90),
  plot.title = element_text (size=base_size * 1.5, face="bold", hjust = 0.5), 
  legend.text = element_text (size=base_size * 1.2),             
  #   panel.grid.major = theme_line (colour = "grey90"),
  panel.grid.major = element_blank(),
  #                   panel.grid.minor = element_blank(), 
  panel.grid.minor = element_blank(),
  #panel.border = element_blank(),
  panel.border = element_rect(colour = "black",fill=NA),
  panel.background = element_blank(),
  axis.line = element_line (colour = "black"))
  # axis.ticks = element_blank())
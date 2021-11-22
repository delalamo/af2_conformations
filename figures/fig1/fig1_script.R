#!/usr/bin/env Rscript

#' Author: Diego del Alamo
#' This script reproduces Figure 1B from the manuscript.
#' Note that this figure will have clipped X-axis facet titles.
#' (namely 'MCT1 (without templates)' will be clipped at the edges)
#' It is therefore recommended that you save this as an SVG and manually
#' edit in Inkscape/Adobe Illustrator.
#' This can be achieved by scrolling to the bottom of this file and replacing
#' the destination file ending from ".png" to ".svg"

library(ggplot2)
library(cowplot)
library(viridis)

#' Data set-up
models <- read.csv(
  "fig1_data_models.csv",
  header=TRUE, sep=",")

natives <- read.csv(
  "fig1_data_natives.csv",
  header=TRUE, sep=",")

models <- subset(
  models,
  set_name %in% c( "32 sequences", "128 sequences", "5120 sequences" ) )

models$set_name <- factor(
  models$set_name,
  c( "32 sequences", "128 sequences", "5120 sequences" ) )

bckg_models <- subset(
  models,
  select = -c( set_name ) )

bckg_models$protein <- factor(
  bckg_models$protein,
  c( "CGRPR", "FZD7", "PTH1R", "ASCT2", "Lat1", "STP10", "ZnT8",
    "MCT1 (without templates)", "MCT1 (with templates)" ) )

#' Plot the data
outpng <- ggplot() +
  geom_point( data=bckg_models, aes( x=tm_IF, y=tm_OF ), color="gray",
              size=0.75 ) +
  geom_hline( data=natives, aes( yintercept=tm ), linetype="dotted",
              color="black", size=0.25 ) +
  geom_vline( data=natives, aes( xintercept=tm ), linetype="dotted",
              color="black", size=0.25 ) +
  geom_point( data=models, aes( x=tm_IF, y=tm_OF, fill=set_name ), size=1.0,
              stroke=0.25, pch=21, color="black" ) +
  scale_fill_viridis_d( begin=0.6, option="plasma",
                        name="Number of sequences" ) +
  labs( x="Similarity to inactive/inward-facing structure (TM-score)",
        y="Similarity to active/outward-facing structure (TM-score)" ) +
  theme_minimal() +
  xlim( 0.55, 1.0 ) +
  ylim( 0.55, 1.0 ) +
  facet_grid( set_name ~ protein ) +
  theme( legend.title = element_text( size=5, face="bold" ),
         legend.text = element_text( size=5 ),
         legend.position ="none",
         panel.border = element_rect( colour = "black", fill=NA, size=0.5 ),
         strip.text.x = element_text( size=5, face="bold" ),
         strip.text.y = element_text( size=5, face="bold" ),
         axis.title.y = element_text( size=5, face="bold" ),
         axis.title.x = element_text( size=5, face="bold" ),
         axis.text.y = element_text( size=5 ),
         axis.text.x = element_text( size=5 ),
         aspect.ratio=1 )

# Save the plot
save_plot( "fig1.png",
           plot = outpng, base_width = 6.5, base_height = 2.5, units="in" )

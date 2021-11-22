#!/usr/bin/env Rscript

#' Author: Diego del Alamo
#' This script reproduces Figure 2 from the manuscript.
#' Note that this figure will have clipped facet titles.
#' (namely 'MCT1 (without templates)' will be clipped at the edges)
#' It is therefore recommended that you save this as an SVG and manually
#' edit in Inkscape/Adobe Illustrator.
#' This can be achieved by scrolling to the bottom of this file and replacing
#' the destination file ending from ".png" to ".svg"
#' Finally, note that R-squared values are NOT calculated in this script and 
#' are instead added manually. Nevertheless such values can be manually verified 

library(ggplot2)
library(cowplot)
library(viridis)

#' Data set-up
data <- read.csv(
  "fig2_data.csv",
  header=TRUE, sep="," )

data$protein <- factor(
  data$protein,
  c( "CGRPR", "FZD7", "PTH1R", "ASCT2", "Lat1", "STP10", "ZnT8",
     "MCT1 (without templates)", "MCT1 (with templates)" ) )

#' Plot the data
outpng <- ggplot() +
  geom_point( data=data, aes( x=exp_dist, y=sim_rmsf, fill=avg_plddt ),
              size=1.0, stroke=0.125,pch=21, color="black" ) +
  scale_fill_viridis( name="pLDDT", labels=c( 75, 80, 85, 90, 95, 100 ),
                      limits=c( 75, 100 ), oob = scales::squish ) +
  labs( x="Movement during alternating access (Å)",
        y="RMSF among AlphaFold2 models (Å)" ) +
  theme_minimal() +
  xlim( 0, NA ) +
  ylim( 0, NA ) +
  facet_wrap( . ~ protein, ncol=3, scales="free" ) +
  theme( legend.title = element_text( size=5, face="bold" ),
         legend.text = element_text( size=5 ),
         legend.position ="bottom",
         panel.border = element_rect(colour = "black", fill=NA, size=0.5 ),
         strip.text.x = element_text( size=5, face="bold" ),
         strip.text.y = element_text( size=5, face="bold" ),
         axis.title.y = element_text( size=5, face="bold" ),
         axis.title.x = element_text( size=5, face="bold" ),
         axis.text.y = element_text( size=5 ),
         axis.text.x = element_text( size=5 ),
         aspect.ratio=1 )

save_plot( "fig2.png",
           plot = outpng, base_width = 3.25, base_height = 4.25, units="in" )

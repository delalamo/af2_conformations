#!/usr/bin/env Rscript

#' Author: Diego del Alamo
#' This script reproduces Figure 3B from the manuscript.
#' Note that the figure width is intentionally set to 3.1 in, rather than
#' the more common publication width 3.25 in. This is to provide room to
#' manually introduce another set of Y-axis labels on the right side of the plot
#' To do so, save plot as an SVG and manually edit in Inkscape/Adobe Illustrator
#' This can be achieved by scrolling to the bottom of this file and replacing
#' the destination file ending from ".png" to ".svg". Additionally, unlike the
#' manuscript version of this figure, the top portion of the plot is not flush
#' against the bottom portion (i.e. there is a gap between the two).This can
#' similarly be fixed by manually editing the figure.

library(ggplot2)
library(cowplot)

#' Data set-up 
data <- read.csv(
  "fig3_data.csv",
  header=TRUE, sep=",")

data$protein <- factor(
  data$protein,
  c( "CGRPR", "FZD7", "PTH1R", "ASCT2", "Lat1", "STP10", "ZnT8", "MCT1" ) )

data_models = subset( data, type == "Model" )
data_natives = subset( data, type != "Model" )

outpng1 <- ggplot() +
  geom_density( data=data, aes( x=pc1 ), adjust=3/4, color="grey" ) +
  geom_hline( yintercept=0, color="black", size=0.25 ) +
  labs( x="PC1", y="Frequency" ) +
  facet_wrap( . ~ protein, scales="free", ncol=8 ) +
  theme_minimal() +
  theme( legend.title = element_blank(),
         legend.text = element_text( size=5 ),
         legend.position ="none",
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.text.x = element_text( size=5, face="bold"),
         strip.text.y = element_text( size=5, face="bold"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_blank(),
         axis.text.x = element_blank() )

outpng2 <- ggplot() +
  geom_point( data=data_models, aes( x=pc1, y=TM_OF ), size=0.75, stroke=0.075,
              pch=21, fill="orange", color="black" ) +
  geom_point( data=data_models, aes( x=pc1, y=TM_IF ), size=0.75, stroke=0.075,
              pch=21, fill="aquamarine1", color="black" ) +
  geom_point( data=data_natives, aes( x=pc1, y=max_tm ), color="black",
              size=0.75 ) +
  labs( x="PC1", y="Structural similarity (TM-score)" ) +
  facet_wrap( . ~ protein, scales="free_x", ncol=8 ) +
  theme_minimal() +
  ylim( NA, 1 ) +
  theme( legend.title = element_blank(),
         legend.text = element_text( size=5 ),
         legend.position ="none",
         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.text.x = element_blank(),
         strip.text.y = element_text( size=5 ),
         axis.title.y = element_text( size=5, face="bold"),
         axis.title.x = element_text( size=5, face="bold"),
         axis.text.y = element_text( size=5 ),
         axis.text.x = element_blank() )

outpng <- plot_grid( outpng1, outpng2, rel_heights=c( 0.35, 0.65 ), nrow=2,
                     align="v", axis="lr" )

save_plot( "fig3.png",
           plot = outpng, base_width = 3.10, base_height = 2, units="in" )




#!/usr/bin/env Rscript

library(ggplot2)
library(cowplot)
library(viridis)

data <- read.csv( "fig1_data_models.csv", header=TRUE, sep=",")
data <- subset( data, protein %in% c( "Lat1" ) )
data <- data[ order( data$nseq ), ]

outpng <- ggplot() +
  geom_hline( yintercept=0.86337, linetype="dotted", color="black", size=0.25 ) +
  geom_vline( xintercept=0.86337, linetype="dotted", color="black", size=0.25 ) +
  geom_point( data=data, aes( x=tm_IF, y=tm_OF, fill=nseq ), size=2.0, stroke=0.375,pch=21, color="black" ) +
  #geom_point( data=af2_data, aes( x=tm_IF, y=tm_OF ), fill="black", size=0.75 ) +
  scale_fill_viridis( begin=0.3, end=1, option="cividis", name="Number of sequences", trans="log2" ) +
  scale_x_continuous( limits=c( 0.85,1 ) ) +
  scale_y_continuous( limits=c( 0.85,1 ) ) +
  labs( x="Similarity to inward-facing structure (TM-score)", y="Similarity to outward-facing structure (TM-score)" ) +
  theme_minimal() +
  #facet_grid( set_name ~ protein ) +
  theme(
    legend.title = element_text( size=6, face="bold"),
    legend.text = element_text( size=6 ),
    legend.position ="none",
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    strip.text.x = element_text( size=6, face="bold"),
    strip.text.y = element_text( size=6, face="bold"),
    axis.title.y = element_text( size=6, face="bold"),
    axis.title.x = element_text( size=6, face="bold"),
    axis.text.y = element_text( size=6),
    axis.text.x = element_text( size=6),
    aspect.ratio=1 )

save_plot( "lat1_out.png", plot = outpng, base_width = 2.25, base_height = 2.25, units="in" )
#R is super slow
#cool packages often don't take well to analysing large numbers of WG data
#use external software to generate metrics and then visualise in ggplot2

#load modules
library(dplyr)
library(ggplot2)

#read in PI data as df
pidf <- read.delim("~/Projects/Schisto/GENOMIC/VCF/out.windowed.pi", header = T)

#function to exclude unassembled scaffolds
#doing this manually at the moment
#newpidf <- filter(pidf, CHROM != "")

#calculate cumulative snp positions
don <- pidf %>%
  
  #Compute chromosome size
  group_by(CHROM) %>%
  summarise(chr_len=max(BIN_END)) %>%
  
  #calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  #LOJ this info to initial dataset
  left_join(pidf, .,by=c("CHROM"="CHROM")) %>%
  arrange(CHROM, BIN_END) %>%
  mutate(BIN_ENDcum=BIN_END+tot)

#create x-axis
axisdf = don %>% group_by(CHROM) %>% summarize(center=( max(BIN_ENDcum) + min(BIN_ENDcum) ) / 2 )

#plot!
ggplot(don, aes(x=BIN_ENDcum, y=PI)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

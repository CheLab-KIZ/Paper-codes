library(plotrix)
library(reshape2)
library(phytools)
library(geiger)
library(RColorBrewer)
require(colorspace)
require(scales)
library(devtools)
library(ggtree)
library(ape)

#create colors for 127 states 
color_fn = paste("range_colors.csv")
range_color_list = read.csv(color_fn, header=T, sep=",", colClasses="character")
colors = vector()
for (i in 1:length( sim2$maps ) ) { 
    colors = c(colors, names(sim2$maps[[i]]) )
}
colors = sort(as.numeric(unique(colors)))
col=vector()
for (i in colors){
  col=c(col,range_color_list$color[i])
}
cols = setNames( col, colors)
ranges = vector()
for (i in colors){
  ranges =c(ranges,range_color_list$range[i] )
}
areas= setNames( ranges,colors)

#001_AMPHIBIANS_Hynobiidae_Batrachuperus
character_file = "Batrachuperus.smap.tree"
sim2 = read.simmap(file=character_file, format="phylip")
pdf('001_AMPHIBIANS_Hynobiidae_Batrachuperus.pdf',
    width=6,height=8,paper='special') 
plotSimmap(sim2,cols,
           fsize=0.6,
           ftype="i",lwd=1, split.vertical=T, mar = c(3,0,3,0),
           xlim=c(-5,35),
           ylim = c(0,5.5),
           direction="rightwards")
add.simmap.legend(x=-5,y=5.5,
                  leg = ranges, prompt=F,colors=cols,shape="circle",
                  fsize= 1)
axisPhylo(side=1.5,las = 1,ask = FALSE,
          lwd = .1,
          col = "black",pos=c(0,0)) # plots timescale
title("Batrachuperus", font.main=4)
title(sub='time (Ma)', adj = 0.6,line=1,cex =1.5)
dev.off()

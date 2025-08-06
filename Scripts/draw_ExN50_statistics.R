#!/usr/bin/env Rscript

library(ggplot2)

args<-commandArgs(TRUE)

if (length(args) == 0) {
    stop("\n\n\tusage: plot_ExN50_statistic.Rscript Trinity.ExN50.stats Trinity.ExN50_plot.png \n\n\n")
}

filein = args[1]
fileout = args[2]

message(sprintf("parsing: %s", filein))
data = read.table(filein, header=T, row.names=NULL)

data$inputfile = basename(filein)
E90N50 = paste0("E90N50: ",(data[data$Ex==90, ]$ExN50)*1e-03, "kbp")

p = ggplot(data, aes(x=Ex, y=ExN50, color=inputfile)) + 
    geom_line() + 
    geom_text(x=90, y=max(data$ExN50), label=E90N50, show.legend=FALSE) + 
    theme_bw()
    
ggsave(filename=fileout, plot=p, device="png", dpi=300)

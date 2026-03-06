#The code is sourced from "Directional integration and pathway enrichment analysis for multi-omics data"

library(ggplot2)
library(ggrepel)

data=read.table('DPM_plot.txt',header=T)
data$LogDPM=-log10(data$DPM)
data$LogBrown=-log10(data$Brown)

pdf(file = "scatterplot.pdf", width = 8, height = 7)
ggplot(data) +
  geom_point(
    size = 2.4, shape = 19,
    aes(
      LogBrown, LogDPM,
      color = ifelse(LogBrown <= 1.301, "gray",
                     ifelse(LogDPM > 1.301, "#1F449C", "#F05039"))
    )
  ) +
  geom_text_repel(
    data = subset(data, Label != "-"),
    aes(LogBrown, LogDPM, label = Label),
    size = 4,                     
    max.overlaps = 30,            
    box.padding = 0.5,           
    point.padding = 0.3,          
    segment.color = "grey60",     
    segment.size = 0.3,           
    min.segment.length = 0,       
    show.legend = FALSE
  ) +
  labs(
    title = "",
    x = "Brown P (-log10)",
    y = "DPM P (-log10)"
  ) +
  geom_hline(yintercept = 1.301, linetype = 'dashed', col = 'black', size = 0.5) +
  geom_vline(xintercept = 1.301, linetype = "dashed", col = "black", size = 0.5) +
  geom_abline(size = 0.5, slope = 1, intercept = 0) +
  scale_color_identity() +
  theme(
    plot.title = element_text(size = 23, hjust = 0.5),
    axis.title.x = element_text(size = 18, margin = unit(c(2, 0, 0, 0), "mm")),
    axis.title.y = element_text(size = 18, margin = unit(c(0, 4, 0, 0), "mm")),
    axis.text = element_text(size = 16),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
dev.off()
library(GGally)
library(ggplot2)
require(gridExtra)
# Data set is provided by R natively
# data <- iris

pdf(file="para.pdf", width=12, height=4, onefile=FALSE)
dat<-read.table("link_figure/for_para_4.csv",sep='\t', header=TRUE)
dat_line<-read.table("link_figure/for_para_line_4.csv",sep='\t', header=TRUE)
# Plot
# head(dat)
p1 <- ggparcoord(dat_line,
    showPoints = FALSE, 
    columns = 2:3, groupColumn = 4
    ) 

p2 <- ggplot(dat, aes(x=dat[,2], y=dat[,3], size = Sample_Num, color=Condition)) +
    geom_point(alpha=0.5)+
    scale_size(range = c(1, 10))+
    labs(x = colnames(dat)[2], y = colnames(dat)[3])
    


grid.arrange(p1, p2, ncol=2)
dev.off()
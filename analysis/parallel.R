library(GGally)
library(ggplot2)
# Data set is provided by R natively
# data <- iris
pdf(file="para.pdf", width=6, height=4, onefile=FALSE)
dat<-read.table("for_para.csv",sep='\t', header=TRUE)
# Plot
# head(dat)
# ggparcoord(dat,
#     showPoints = TRUE, 
#     columns = 2:3, groupColumn = 4
#     ) 

p <- ggplot(dat, aes(x=dat[,2], y=dat[,3], size = Sample_Num, color=Condition)) +
    geom_point(alpha=0.5)+
    scale_size(range = c(1, 10))

p + labs(x = colnames(dat)[2], y = colnames(dat)[3])
    


dev.off()
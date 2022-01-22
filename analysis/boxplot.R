# Load ggplot2
library(ggplot2)
require(gridExtra)
 
# The mpg dataset is natively available
dat<-read.table("for_box.csv", header=TRUE)
head(dat)
pdf(file="density_box.pdf", width=8, height=8)

# geom_boxplot proposes several arguments to custom appearance
plot1 <- ggplot(dat, aes(x=Condition, y=Density)) + 
    geom_boxplot(
        
        # custom boxes
        color="blue",
        fill="blue",
        alpha=0.2,
        
        # Notch?
        notch=TRUE,
        notchwidth = 0.8,
        
        # custom outliers
        outlier.colour="red",
        outlier.fill="red",
        outlier.size=3
    
    )

plot2 <- ggplot(dat, aes(x=Condition, y=Transtivity)) + 
    geom_boxplot(
        
        # custom boxes
        color="blue",
        fill="blue",
        alpha=0.2,
        
        # Notch?
        notch=TRUE,
        notchwidth = 0.8,
        
        # custom outliers
        outlier.colour="red",
        outlier.fill="red",
        outlier.size=3
    
    )

plot3 <- ggplot(dat, aes(x=Condition, y=Nodes_Number)) + 
    geom_boxplot(
        
        # custom boxes
        color="blue",
        fill="blue",
        alpha=0.2,
        
        # Notch?
        notch=TRUE,
        notchwidth = 0.8,
        
        # custom outliers
        outlier.colour="red",
        outlier.fill="red",
        outlier.size=3
    
    )

plot4 <- ggplot(dat, aes(x=Condition, y=Pair_Number)) + 
    geom_boxplot(
        
        # custom boxes
        color="blue",
        fill="blue",
        alpha=0.2,
        
        # Notch?
        notch=TRUE,
        notchwidth = 0.8,
        
        # custom outliers
        outlier.colour="red",
        outlier.fill="red",
        outlier.size=3
    
    )

grid.arrange(plot1, plot2, plot3, plot4, ncol=2)
dev.off()

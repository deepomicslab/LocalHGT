# Load ggplot2
library(ggplot2)
require(gridExtra)
 
# The mpg dataset is natively available
dat<-read.table("graph_density.csv", header=TRUE)

pdf(file="density_box.pdf", width=8, height=8)

# geom_boxplot proposes several arguments to custom appearance

new_dat<-dat[dat[,"level"] == 1,]
plot1 <- ggplot(new_dat, aes(x=Status, y=Density)) + 
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
    
    )+
    ggtitle("Phylum")+
    theme(plot.title = element_text(hjust = 0.5))
    # level_dict = {"phylum":1, "class":2, "order":3, "family":4, "genus":5, "species":6}

new_dat<-dat[dat[,"level"] == 2,]
plot2 <- ggplot(new_dat, aes(x=Status, y=Density)) + 
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
    
    )+
    ggtitle("Class")+
    theme(plot.title = element_text(hjust = 0.5))

new_dat<-dat[dat[,"level"] == 3,]
plot3 <- ggplot(new_dat, aes(x=Status, y=Density)) + 
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
    
    )+
    ggtitle("Order")+
    theme(plot.title = element_text(hjust = 0.5))

new_dat<-dat[dat[,"level"] == 4,]
plot4 <- ggplot(new_dat, aes(x=Status, y=Density)) + 
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
    
    )+
    ggtitle("Family")+
    theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot1, plot2, plot3, plot4)
dev.off()

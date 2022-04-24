library(ggplot2)
require(gridExtra)
library(grid)

library("cowplot")


dat<-read.table("cami_comparison.csv",sep=',', header=TRUE)
# dat$Recall
pdf(file="cami.pdf", width=10, height=6, onefile=FALSE)
theme_set(theme_bw())

my_width<-0.005

new_dat <- dat[dat$Complexity == 'high', ]
# new_dat

p1 <- ggplot(new_dat, aes(x=Mutation.Rate , y=CPU.time, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("CPU time") + scale_fill_manual(values=c("slateblue1" , "tomato")) +theme(legend.position="none")+
    ggtitle("High")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )
    
    # theme(axis.text.x = element_text(size = 20))+   
    # theme(axis.text.y = element_text(size = 20))+ 
    

p4 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Peak.RAM, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("Peak RAM")+ scale_fill_manual(values=c("slateblue1" , "tomato")) +theme(legend.position="none")+
    ggtitle("High")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )

p7 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Recall, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("Recall")+ scale_fill_manual(values=c("slateblue1" , "tomato"))+theme(legend.position="none")+
    ggtitle("High")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )

new_dat <- dat[dat$Complexity == 'medium', ]
p2 <- ggplot(new_dat, aes(x=Mutation.Rate , y=CPU.time, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("CPU time")+ scale_fill_manual(values=c("slateblue1" , "tomato"))+theme(legend.position="none")+
    ggtitle("Medium")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )

p5 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Peak.RAM, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("Peak RAM") + scale_fill_manual(values=c("slateblue1" , "tomato"))+theme(legend.position="none")+
    ggtitle("Medium")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )

p8 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Recall, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("Recall") + scale_fill_manual(values=c("slateblue1" , "tomato"))+theme(legend.position="none")+
    ggtitle("Medium")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )


new_dat <- dat[dat$Complexity == 'low', ]
p3 <- ggplot(new_dat, aes(x=Mutation.Rate , y=CPU.time, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("CPU time")+ scale_fill_manual(values=c("slateblue1" , "tomato"))+theme(legend.position="none")+
    ggtitle("Low")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )

p6 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Peak.RAM, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 
    xlab("Mutation rate")+
    ylab("Peak RAM")+ scale_fill_manual(values=c("slateblue1" , "tomato"))+theme(legend.position="none")+
    ggtitle("Low")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )


p9 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Recall, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge", width = my_width)+ 

    xlab("Mutation rate")+
    ylab("Recall")+ scale_fill_manual(values=c("slateblue1" , "tomato"))+theme(legend.position="none")+
    ggtitle("Low")+
    theme(
    plot.title = element_text(size=9, face="italic", hjust = 0.5),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8)
    )

# p10 <- ggplot(dep, aes(x=Depth , y=F1.score, fill = Methods)) + 
#     geom_bar(stat="identity", position = "dodge", width = my_width)+ 
#     xlab("Depth")+
#     ylab("F1 score")+ scale_fill_manual(values=c("slateblue1" , "tomato"))


grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}


# grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=1)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=4)

dev.off()


library(ggplot2)
require(gridExtra)
library(grid)


dat<-read.table("cami_comparison.csv",sep=',', header=TRUE)
dat
dat$Recall
pdf(file="cami.pdf", width=12, height=8, onefile=FALSE)
theme_set(theme_bw())


new_dat <- dat[dat$Complexity == 'high', ]
new_dat

p1 <- ggplot(new_dat, aes(x=Mutation.Rate , y=CPU.time, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )

p4 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Peak.RAM, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )

p7 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Recall, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )

new_dat <- dat[dat$Complexity == 'medium', ]
p2 <- ggplot(new_dat, aes(x=Mutation.Rate , y=CPU.time, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )

p5 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Peak.RAM, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )  

p8 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Recall, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )  


new_dat <- dat[dat$Complexity == 'low', ]
p3 <- ggplot(new_dat, aes(x=Mutation.Rate , y=CPU.time, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )

p6 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Peak.RAM, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )


p9 <- ggplot(new_dat, aes(x=Mutation.Rate , y=Recall, fill = Methods)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    scale_fill_grey(
        start = 0.2,
        end = 0.8,
        na.value = "red",
        aesthetics = "fill"
    )






grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}


grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3)

# grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3)

dev.off()
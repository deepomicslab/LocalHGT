# Load ggplot2
library(ggplot2)

# create dummy data
data <- data.frame(
  Condition=c("CRC", "control"),
  Density=c(0.166, 0.149),
  sd=c(0.031, 0.027)
)

data
pdf(file="density.pdf", width = 10, height = 2)
# Most basic error bar
ggplot(data) +
    geom_bar( aes(x=Condition, y=Density), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=Condition, ymin=Density-sd, ymax=Density+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)


dev.off()
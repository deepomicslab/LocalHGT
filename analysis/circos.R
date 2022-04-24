library(circlize)
set.seed(10)
dat<-read.table("for_circos.csv", header=TRUE)
pdf(file="circles.pdf")




par(cex = 1, mar = c(0, 0, 0, 0))
grid.col = c(Bacteroidales = "pink", Christensenellales = "#6EE2FF", Oscillospirales = "#F7C530",
Lachnospirales = "#95CC5E", Enterobacterales ="#D0DFE6",  other = "lightgrey")
# chordDiagram(as.data.frame(dat),col = c("red", "skyblue", "pink", "yellow"))
# dev.off()


mat<-as.matrix(dat)
mat
cols <- hcl.colors(6, "Temps")


my_order<-c( "Bacteroidales",     "Lachnospirales",  "Oscillospirales", "Christensenellales",  "others",   "Enterobacterales")
chordDiagram(mat, annotationTrack = c("name", "grid"), order=my_order, transparency = 0.2, grid.col = grid.col )
circos.clear()
dev.off()


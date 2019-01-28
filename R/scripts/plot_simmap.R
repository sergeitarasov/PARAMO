library(plotrix)
library(phytools)

character_file = "character.tree"
character_file = "CHAR-1.sstm"

write_pdf = TRUE

if (write_pdf) {
    pdf("simmap.pdf")
}

ls(getNamespace("phytools"), all.names=TRUE)
phytools:::text_to_tree

text <- scan(file="CHAR-14.stm", sep = "\n", what = "character")
regexpr("\\}\t\\(",  text[4])

regexpr("sa", "asa")
describe.simmap

?grep
sim1 = read.simmap(file="character-1.tree", format="phylip")
sim1t = read.simmap(file="CHAR-14-1-t.stm", format="phylip", rev.order=T)
sim2 = read.simmap(file="character-2.tree", format="phylip")

plotSimmap(sim1, setNames(c("blue","red"), c(0,1)), fsize=0.5, lwd=1, split.vertical=TRUE, ftype="i")
plotSimmap(sim2, setNames(c("blue","red"), c(0,1)), fsize=0.5, lwd=1, split.vertical=TRUE, ftype="i")
plotSimmap(sim1t, setNames(c("blue","red"), c(0,1)), fsize=0.5, lwd=1, split.vertical=TRUE, ftype="i")

str(sim2)

sim2$maps

sim2 = read.simmap(file=character_file, format="phylip")
str(sim2)
sim2$maps[[172]] %>% length()
write.simmap(sim2, file = "character-2-rewritte.tree", map.order="right-to-left")

sim2$mapped.edge

#text <- scan(file=character_file, sep = "\n", what = "character")
str(text)
class(text)
write(c(text[1:3], text[5000:5001]), file = "CHAR-14.stm")

colors = vector()
for (i in 1:length( sim2$maps ) ) { 
    colors = c(colors, names(sim2$maps[[i]]) )
}
colors = sort(as.numeric(unique(colors)))

cols = setNames( rainbow(length(colors), start=0.0, end=0.9), colors)

plotSimmap(sim2, cols, fsize=0.5, lwd=1, split.vertical=TRUE, ftype="i")

# add legend
x = 0
y = 15

leg = names(cols)
colors = cols
y = y - 0:(length(leg) - 1)
x = rep(x, length(y))
text(x + 0.0005, y, leg, pos=4, cex=0.75)

mapply(draw.circle, x=x, y=y, col=colors, MoreArgs = list(nv=200, radius=0.001, border="white"))

if (write_pdf) {
    dev.off()
}


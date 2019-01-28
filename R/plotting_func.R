require(graphics)

##########################################
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    #par(mar=c(bottom=3, left=3, top=19, right=35))
  #plot(0:4,  type="n", axes=FALSE, xlab = "", ylab = "", bty="n")
    plot(c(0,1), c(min,max),  type="n",  bty='n', xlab='',  ylab='', main=title, xaxt='n',  yaxt='n')
    #plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}
#########
# # my old functions for color bars
# par(mar=c(bottom=3, left=4, top=20, right=10))
# plot(0:4,  type="n", axes=FALSE, xlab = "", ylab = "", bty="n")
# axis(side=2, at=c(0, 2, 4), labels = c(0, 40, 75), lwd=0, lwd.ticks=1, las=1, padj=0.5,  hadj=1)
# gradient.rect(xleft=0, ybottom=0, xright=1, ytop=4, nslices=max(col$stat), col = rev(bluecols2), gradient="y", border=NA)
# gradient.rect(xleft=0, ybottom=0, xright=1, ytop=4, nslices=max(col$stat), col = rev(color), gradient="y", border=NA)
# ##

# different ways for color bars
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
lut  = c(rev(cool), rev(warm))
par(mar=c(bottom=3, left=3, top=19, right=35))
color.bar(lut, -0.05, 0.05)

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
med=rainbow(50, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('cyan'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
lut  = c(rev(cool), rev(med), rev(warm))
color.bar(lut, -0.05, 0.05)

lut = rev(rainbow(100, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('blue'))[1]))
color.bar(lut, -0.1, 0.1)

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(25, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
lut  = c(rev(cool), rev(warm))
color.bar(lut, -0.05, 0.025, ticks=c(-0.05, -0.025, 0, 0.025))
##########################



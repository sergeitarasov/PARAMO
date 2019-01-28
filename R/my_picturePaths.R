# picture<-petal
# nr = 20
# nc = 5
# label=20
# freeScales = TRUE
# n=1
# str(tg)
# 
# my_picturePaths(petal, nr = 20, nc = 5, freeScales = TRUE, fill = "white")

my_label = function(n) {
  tg <- textGrob(n, x = 0, y = 0, just = c("left", "bottom"),  gp = gpar(fontsize = 12))
  # grid.rect(x = 0, y = 0, height = unit(6, "points"), 
  #           width = grobWidth(tg), just = c("left", "bottom"), 
  #           gp = gpar(fill = "white"))
  
  grid.draw(tg)
}

my_picturePaths<-function(picture, nr, nc, col = "black", fill = "light grey", 
          freeScales = FALSE, xscale = NULL, yscale = NULL, label = function(n) {
            tg <- textGrob(n, x = 0, y = 0, just = c("left", "bottom"), 
                           gp = gpar(fontsize = 6))
            grid.rect(x = 0, y = 0, height = unit(6, "points"), 
                      width = grobWidth(tg), just = c("left", "bottom"), 
                      gp = gpar(fill = "white"))
            grid.draw(tg)
          }, use.gc = TRUE) 
{
  if (missing(nr) || missing(nc)) {
    nrnc <- n2mfrow(picture@summary@numPaths)
    nr <- nrnc[1]
    nc <- nrnc[2]
  }
  if (is.null(xscale) || is.null(yscale)) {
    xscale <- picture@summary@xscale
    yscale <- picture@summary@yscale
  }
  if (freeScales) {
    pushViewport(viewport(layout = grid.layout(nr, nc)))
  }
  else {
    pushViewport(viewport(layout = grid.layout(nr, nc, widths = rep(diff(range(xscale)), 
                                                                    nc), heights = rep(diff(range(yscale)), nr), respect = TRUE)))
  }
  
  ##
  i=1
  for (i in 1:nr) {
    
    #j=1
    for (j in 1:nc) {
      pnum <- (i - 1) * nc + j
      if (pnum <= picture@summary@numPaths) {
        pushViewport(viewport(layout.pos.col = j, layout.pos.row = i))
        if (freeScales) {
          pushViewport(viewport(width = 0.95, height = 0.95))
          grid.rect(gp = gpar(col = col, fill = fill))
          grid.picture(picture[pnum])
          
          #
          my_label(pnum)
          popViewport()
        }
        else {
          pushViewport(viewport(width = 0.95, height = 0.95, 
                                xscale = xscale, yscale = yscale))
          grid.rect(gp = gpar(col = col, fill = fill))
          if (is.function(label)) {
            label(pnum)
          }
          grob <- grobify(picture@paths[[pnum]], use.gc = use.gc)
          if (is.null(grob)) {
            grid.text("EMPTY\nPATH", gp = gpar(fontsize = 6))
          }
          else {
            grid.draw(grob)
          }
          popViewport()
        }
        popViewport()
      }
    }
  }
  popViewport()
}

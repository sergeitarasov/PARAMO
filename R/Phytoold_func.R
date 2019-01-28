z=densityMap(tt1[[1]])
plot(z)
str(z)
z$tree$maps
trees=tt1[[1]]

function (trees=tt1, res = 100, fsize = NULL, ftype = NULL, lwd = 3, 
          check = FALSE, legend = NULL, outline = FALSE, type = "phylogram", 
          direction = "rightwards", plot = TRUE, ...) 
{
  if (hasArg(mar)) 
    mar <- list(...)$mar
  else mar <- rep(0.3, 4)
  if (hasArg(offset)) 
    offset <- list(...)$offset
  else offset <- NULL
  if (hasArg(states)) 
    states <- list(...)$states
  else states <- NULL
  if (hasArg(hold)) 
    hold <- list(...)$hold
  else hold <- TRUE
  if (length(lwd) == 1) 
    lwd <- rep(lwd, 2)
  else if (length(lwd) > 2) 
    lwd <- lwd[1:2]
  tol <- 0.0000000001
  if (!inherits(trees, "multiPhylo") && inherits(trees, "phylo")) 
    stop("trees not \"multiPhylo\" object; just use plotSimmap.")
  if (!inherits(trees, "multiPhylo")) 
    stop("trees should be an object of class \"multiPhylo\".")
  
  h <- sapply(unclass(trees), function(x) max(nodeHeights(x)))
  
  steps <- 0:res/res * max(h)
  trees <- rescaleSimmap(trees, totalDepth = max(h))
  #str(trees)
  #trees[[1]]$maps
  
  if (check) {
    X <- matrix(FALSE, length(trees), length(trees))
    for (i in 1:length(trees)) X[i, ] <- sapply(trees, all.equal.phylo, 
                                                current = trees[[i]])
    if (!all(X)) 
      stop("some of the trees don't match in topology or relative branch lengths")
  }
  tree <- trees[[1]]
  trees <- unclass(trees)
  if (is.null(states)) 
    ss <- sort(unique(c(getStates(tree, "nodes"), getStates(tree, 
                                                            "tips"))))
  else ss <- states
  if (!all(ss == c("0", "1"))) {
    c1 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    c2 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    trees <- lapply(trees, mergeMappedStates, ss[1], c1)
    trees <- lapply(trees, mergeMappedStates, ss[2], c2)
    trees <- lapply(trees, mergeMappedStates, c1, "0")
    trees <- lapply(trees, mergeMappedStates, c2, "1")
  }
  H <- nodeHeights(tree)
  message("sorry - this might take a while; please be patient")
  
  tree$maps <- vector(mode = "list", length = nrow(tree$edge))
  i=1
  for (i in 1:nrow(tree$edge)) {
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > 
                                                   H[i, 1]), which(steps < H[i, 2]))]), c(steps[intersect(which(steps > 
                                                                                                                  H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) - 
      H[i, 1]
    ZZ <- rep(0, nrow(YY))
    
    j=1
    for (j in 1:length(trees)) {
      XX <- matrix(0, length(trees[[j]]$maps[[i]]), 2, 
                   dimnames = list(names(trees$maps[[i]]),  #dimnames = list(names(trees[[j]]$maps[[i]]), 
                                   c("start", "end")))
      XX[1, 2] <- trees[[j]]$maps[[i]][1]
      if (length(trees[[j]]$maps[[i]]) > 1) {
        for (k in 2:length(trees[[j]]$maps[[i]])) {
          XX[k, 1] <- XX[k - 1, 2]
          XX[k, 2] <- XX[k, 1] + trees[[j]]$maps[[i]][k]
        }
      }
      for (k in 1:nrow(YY)) {
        lower <- which(XX[, 1] <= YY[k, 1])
        lower <- lower[length(lower)]
        upper <- which(XX[, 2] >= (YY[k, 2] - tol))[1]
        AA <- 0
        names(lower) <- names(upper) <- NULL
        if (!all(XX == 0)) {
          for (l in lower:upper) AA <- AA + (min(XX[l, 
                                                    2], YY[k, 2]) - max(XX[l, 1], YY[k, 1]))/(YY[k, 
                                                                                                 2] - YY[k, 1]) * as.numeric(rownames(XX)[l])
        }
        else AA <- as.numeric(rownames(XX)[1])
        ZZ[k] <- ZZ[k] + AA/length(trees)
      }
    }
    tree$maps[[i]] <- YY[, 2] - YY[, 1]
    names(tree$maps[[i]]) <- round(ZZ * 1000)
  }
  cols <- rainbow(1001, start = 0.7, end = 0)
  names(cols) <- 0:1000
  tree$mapped.edge <- makeMappedEdge(tree$edge, tree$maps)
  tree$mapped.edge <- tree$mapped.edge[, order(as.numeric(colnames(tree$mapped.edge)))]
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  x <- list(tree = tree, cols = cols, states = ss)
  class(x) <- "densityMap"
  if (plot) 
    plot.densityMap(x, fsize = fsize, ftype = ftype, lwd = lwd, 
                    legend = legend, outline = outline, type = type, 
                    mar = mar, direction = direction, offset = offset, 
                    hold = hold)
  invisible(x)
}

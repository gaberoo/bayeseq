#' Traverse a phylogeny and sum read counts
#' 
#' @param counts Matrix of counts: rows are ASVs and columns are samples.
#' @param tree Phylogeny in \code{ape} format.
#' @param node (optional) Node in phylogeny below which to sum. Defaults to root node.
#' @param prefix (optional) String prefix for internal nodes. Defaults to \code{__node_}.
#' @return Augmented counts matrix with additional rows for internal nodes.

sum.counts <- function(counts,tree,node=tree$edge[1],prefix="__node_") {
  edge.ids <- which(tree$edge[,1] == node)
  # get descendents
  if (length(edge.ids) > 0) {
    children <- tree$edge[edge.ids,2]
    child.cnts <- lapply(children,sum.counts,counts=counts,tree=tree)
    child.cnts <- do.call(rbind,child.cnts)
    node.cnts <- colSums(child.cnts[!grepl("__node_",rownames(child.cnts)),])
    out <- rbind(node.cnts,child.cnts)
    rownames(out)[1] <- sprintf("__node_%d",node)
    out
  } else {
    out <- data.frame(t(counts[tree$tip.label[node],]))
    rownames(out) <- tree$tip.label[node]
    out
  }
}



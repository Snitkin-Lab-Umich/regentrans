#tests for prune_tree
set.seed(1)
test_tr <- ape::rtree(25)
patients <- sample(c('A','B', 'C', 'D', 'E', 'F', 'G'), 25, replace = TRUE)
test_tr <- phytools::midpoint.root(test_tr)
names(patients) <- test_tr$tip.label
plot_pat <- c(patients, rep(NA, ape::Nnode(test_tr)))
ggtree::ggtree(test_tr) + ggtree::geom_tippoint(ggplot2::aes(col=plot_pat), size = 3) + ggtree::geom_tiplab(hjust = -0.5)

### Extract the subtrees
l<-subtrees(test_tr)

### plot all the subtrees
for (i in 1:length(l)) plot(l[[i]], sub=paste("Node", l[[i]]$node.label[1]))
par(mfrow=c(6,4))

test_that("prune_tree works", {
  tr_2 <- prune_tree(test_tr, patients)
  
})
plot_pat <- c(patients[names(patients) %in% tr_2$tip.label], rep(NA, ape::Nnode(tr_2)))
ggtree::ggtree(tr_2) + ggtree::geom_tippoint(ggplot2::aes(col=plot_pat), size = 3) + ggtree::geom_tiplab(hjust = -0.5)

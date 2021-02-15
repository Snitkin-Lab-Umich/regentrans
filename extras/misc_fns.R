# functions to include in the regional genomic transmission analysis package

# some notes:
# we might want to change the function names (should start with a verb, be informative and not too long)
# we might want to change the format/colnames of the output as well
# good to include checks to see if the input is correct (within the function; can have a separate script for all check functions) and unit tests to check the functions (also in a separate script; we'll hopefully find bugs this way!)
# we'll want to make some arguments optional - i've tried to note where this will be, because they're not always optional now
# we'll want to be consistent with what arguments are called across functions (same name = same input)
# you'll have to put the package name in front of functions that come from a different package 

# compare pairwise snv distances returned by ape::dist.dna() using pairwise.deletion = T and F
# example of using dist.dna: dist.dna(dna_aln, model = 'N', pairwise.deletion = T, as.matrix = T)

# get intra- and inter-facility pairwise snv distances (might want to change this name to explain better what we're outputting)
# dists - snv distance matrix returned by dist.dna
# locs - locations of isolates (e.g. facility of isolation), named with the same names as dists
# pt - patient isolate was taken from (optional; will remove pairwise snv distances between the same patient), named with the same name as dists
get_snv_dists <- function(dists, locs, pt = NULL){
  # put checks at beginning
  check_names(rownames(dists), colnames(dists)) # see how many overlap - in this function you'll  put warnings
  # possibly make a subset function (so all are the same)
  loc_sub <- locs[rownames(dists)]
  pt_sub <- pt[rownames(dists)]
  # snps <- dists[lower.tri(dists, diag = FALSE)] <- NA
  snps <- na.omit(data.frame(as.table(as.matrix(dists))))
  snps$loc1 <- loc_sub[snps$Var1]
  snps$loc2 <- loc_sub[snps$Var2]
  snps$pt1 <- pt_sub[snps$Var1]
  snps$pt2 <- pt_sub[snps$Var2]
  snps
  bind_cols(snps %>% filter(pt1 != pt2) %>% mutate(intra=ifelse(loc1==loc2,'Intra-facility pair','Inter-facility pair')))
}

# get fraction of intra-facility pairs for different snv thresholds
# dists - output of get_snv_dists - will want to change either the name of this argument or the other one so the arguments always mean the same input
# threshs - what snv thresholds to use
get_frac_intra <- function(dists, threshs = seq(1,50,1)){
  intra_cts <- t(sapply(threshs, function(i){
    intra <- dists$intra[dists$Freq < i & !is.na(dists$intra)]
    im <- c(i,
            table(factor(intra,levels=c('Intra-facility pair','Inter-facility pair'))),
            round(mean(intra == 'Intra-facility pair'),2),
            round(mean(intra == 'Inter-facility pair'),2))
    names(im) <- c('thresh','inter','intra','frac_intra','frac_inter')
    im
  }))
  data.frame(intra_cts)
}

# get clusters
# tr - tree
# locs - locations - named vector where names are same as tree tip labels # maybe make this more general?
# pureness - how pure the cluster is (<= 0.5)
# maybe add bootstrap value
get_clusters <- function(tr, locs, pureness = 1){ # pureness shouldn't be <= 0.5?
  locs_sub <- locs[tr$tip.label]
  subtrs_sub <- ape::subtrees(tr)
  pure_subtrees <- get_largest_subtree(subtrs_sub, locs_sub, bootstrap = NULL, pureness = pureness) # NOTE: this function is in snitkitr right now, but I think we should migrate it to this package (or at least include it here as well); this _might_ be buggy, so definitely good to add unit tests for it
  pure_subtr_info <- bind_cols(f_id=locs_sub,
                               subtr_size=unlist(pure_subtrees$largest_st),
                               index=unlist(pure_subtrees$largest_st_i))
  # change singletons from 0 to 1
  pure_subtr_info <- pure_subtr_info %>% mutate(subtr_size=ifelse(subtr_size==0 & index == 1, 1, subtr_size))
  # remove duplicates (singletons aren't duplicates)
  pure_subtr_info <- pure_subtr_info[!duplicated(pure_subtr_info) | pure_subtr_info$subtr_size == 1,]
  return(pure_subtr_info=pure_subtr_info) #maybe add pureness of cluster? 
}

## Shortest path functions

# find k shortest paths (yen's algorithm)
# Resources:
#   - pseudocode on wikipedia: https://en.wikipedia.org/wiki/Yen%27s_algorithm
#   - python function: https://gist.github.com/ALenfant/5491853
# input:
#   - graph = igraph graph
#   - from = start vertex
#   - to = target vertex
#   - num_k = number of shortest paths to find
# output:
#   - list of paths and distances for k shortest paths
k_shortest_paths <- function(graph, from, to, num_k, weight=NULL){
  
  # SET GRAPH EDGE WEIGHTS
  if(!is.null(weight)){
    E(graph)$weight = weight
  }
  
  # first shortest path from source to target
  A <- list(as_ids(shortest_paths(graph,from,to,mode='out',output='vpath')$vpath[[1]]))
  A_costs = list(distances(graph,from,to,mode='out',algorithm='dijkstra'))
  
  # initialize heap to store potential kth shortest path
  B = list()
  B_costs = list()
  
  #print(paste('Path 1'))
  #print(as_ids(A[[1]]))
  
  if(num_k == 1){
    return(list(paths = A,dist = A_costs))
  }
  
  for(k in 2:num_k){
    #print(paste('Path',k))
    # spur node ranges from first node to next to last node in shortest path
    for(i in 1:(length((A[[k-1]]))-1)){
      # spur node retrieved from previous k-shortest path, k-1
      spurNode = (A[[k-1]])[i]
      # sequence of nodes from source to spur node of previous k-shortest path
      rootPath = A[[k-1]][1:i]
      
      # store removed edges
      removed_edges = c()
      removed_edge_attr = list()
      
      #c = 0
      for(path in A){
        if(length(path) > i){
          if(all((rootPath) == (path[1:i]))){
            # remove links that are part of previous shortest paths that share the same root path
            edge = paste((path[i:(i+1)]),collapse = '|')
            removed_edges = c(removed_edges, edge)
            #c = c+1
            #removed_edge_attr[[c]] = list(path[i],path[i+1],edge_attr(g,index=edge))
            
          }
        }
      }
      
      V(graph)$names = as_ids(V(graph))
      
      
      newgraph = delete_edges(graph,unique(removed_edges))
      
      for(node in (rootPath)){
        if(node != spurNode){
          newgraph = delete_vertices(newgraph,c(node))
        }
      }
      
      spurNode = (V(graph)$names[as_ids(V(graph)) %in% spurNode])
      spurNode = (V(newgraph)[V(newgraph)$names %in% (spurNode)])
      
      # calculate spur path from spur node to sink
      spurPath = as_ids(shortest_paths(newgraph,names(spurNode),as.character(to),mode='out',output='vpath')$vpath[[1]])
      
      if(length(spurPath) > 1){
        
        totalPath = c((rootPath), (spurPath)[2:length(spurPath)])
        
        totalPathCost = sum(E(graph)$weight[(get.edge.ids(graph,rep((totalPath),each=2)[2:(length((totalPath))*2-1)]))])
        
        
        # add potential k-shortest path to heap
        B[[length(B)+1]] = totalPath
        B_costs[[length(B_costs)+1]] = totalPathCost
      }
      
    }
    
    # sort potential k-shortest paths by cost
    B = B[order(unlist(B_costs))]
    B_costs = B_costs[order(unlist(B_costs))]
    
    # lowest cost path is k-shortest path
    for(i in 1:length(B_costs)){
      cost_ = B_costs[i]
      path_ = B[i]
      if(!(paste((path_[[1]]),collapse = '_') %in% (lapply(A, function(x) paste((unlist(x)),collapse = '_'))))){
        # found new path to add
        A[[k]] = path_[[1]]
        A_costs[[k]] = cost_[[1]]
        break
      }
    }
    
    #print((as_ids(A[[k]])))
  }
  return(list(paths = A,dist = A_costs))
}


# GET DISTANCE MATRICES FOR THE K SHORTEST PATHS 
# inupt: 
#   - k = k shortest paths (integer)
#   - graph = graph (igraph object)
#   - alpha - ignore for now
#   - ltach_names - vertex names of facilities of interest
#   - weight - weight of edges to use to calculate shortest graphs (if NULL, uses E(graph)$weight)
# returns list with each matrix being shortest distance for kth shortest path
get_dist_mat_ksp = function(k, graph, ltach_indices, weight = NULL){
  
  # SET GRAPH EDGE WEIGHTS
  if(!is.null(weight)){
    E(graph)$weight = weight
  }
  
  # GET K SHORTEST PATHS BETWEEN VERTICES OF INTEREST
  ksp = list()
  
  for(i in ltach_indices){
    for(j in ltach_indices){
      if(i != j){
        sp_k = k_shortest_paths(graph,i,j,k,weight)
        ksp[[paste(i,j,sep='_')]] = sp_k
      }
    }
  }
  
  # CALCULATE DISTANCE MATRIX FOR EACH K
  dist_mats = vector('list',k)
  for(x in 1:length(dist_mats)){
    dist_mats[[x]] = matrix(0, nrow=length(ltach_indices),ncol=length(ltach_indices))
    colnames(dist_mats[[x]]) = ltach_indices
    rownames(dist_mats[[x]]) = ltach_indices
  }
  
  paths = list() #vector('list',k)
  
  ksp_count = 0
  for(sps in ksp){
    ksp_count = ksp_count + 1
    i = (strsplit(names(ksp)[ksp_count],'_')[[1]][1])
    j = (strsplit(names(ksp)[ksp_count],'_')[[1]][2])
    c = 0
      paths[[paste(i,j,sep = '_')]] = sps$path
    c = 0
    #distances
    for(d in sps$dist){
      c = c+1
      dist_mats[[c]][i,j] = d
    }
  }
  
  return(list('dist_mats' = dist_mats, 'paths' = paths))
  
}

get_homology_from_phyloseq <- function(phyloseq_object, metric = "bray"){
  return(calculate_homology(as.matrix(distance(taxa_filter(phyloseq_object), metric))))
}

get_persistence_list <- function(phyloseq_object, treatment, metric = "bray", dimension = 2) {
  persistence_list <- list()
  phyloseq_object = taxa_filter(phyloseq_object)
  for(u in unique(phyloseq_object@sam_data[[treatment]])){
    persistence_list[[u]]$phyloseq <- eval(parse(text =
                                                   paste0('subset_samples(phyloseq_object, ', treatment, ' == "', u, '")')
    ))
    # persistence_list[[u]]$phyloseq <- do.call(subset_samples, list(phyloseq_object, treatment == u))
    dist_mat = as.matrix(do.call(distance, list(persistence_list[[u]]$phyloseq, metric)))
    persistence_list[[u]]$hom <-  calculate_homology(dist_mat, dim = dimension, format = "distmat")
    
  }
  return(persistence_list)
}

get_pairwise_permutation_tests <- function(persistence_list, iters = 1000){
  permutation_tests <- list()
  unique_treatments <- names(persistence_list)
  combs <- combn(1:length(unique_treatments), 2)
  for (p in 1:dim(combs)[2]){
    i = combs[1, p]
    j = combs[2, p]
    current <- paste0(names(persistence_list[i]), "-", 
                      names(persistence_list[j]))
    permutation_tests[[current]]$Permutation_Test <- permutation_test(persistence_list[[i]]$hom,
                                                                      persistence_list[[j]]$hom,
                                                                      iterations = iters)
    permutation_tests[[current]]$Permutation_Test[[1]]$permvals <- NULL
    permutation_tests[[current]]$Permutation_Test[[2]]$permvals <- NULL
  }
  return(permutation_tests)
}

get_dimension_heatmap <- function(permutations_object, dimension = 2, triangle = FALSE){
  unique_treatments <- unique(c(str_split(names(permutations_object), "-", simplify = TRUE)))
  dim_perm_values <- matrix(1, nrow = length(unique_treatments), ncol = length(unique_treatments),
                            dimnames = list(unique_treatments, unique_treatments)
  )
  for (i in 1:length(permutations_object)){
    treatments <- str_split(names(permutations_object[i]), '-', simplify = TRUE)
    current <- permutations_object[[i]] 
    dim_perm_values[treatments[1], treatments[2]] <- current$Permutation_Test[[dimension]]$pvalue
    dim_perm_values[treatments[2], treatments[1]] <- current$Permutation_Test[[dimension]]$pvalue
  }
  
  if (triangle){
    dim_perm_values[upper.tri(dim_perm_values)] <- NA
  }
  
  dim_perm_values <- dim_perm_values[ncol(dim_perm_values):1, ]
  
  p <- ggplot(data = na.omit(melt(dim_perm_values)), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "gray") + coord_equal() +
    geom_text(aes(label = value)) +
    scale_fill_gradient2(low = "yellow", high = 'white',
                         midpoint = 0.2, limit = c(0, 1), space = "Lab") + 
    labs(x = element_blank(), y = element_blank(), fill = "p-value", 
         title = "Permutation test p-value between treatments",
         subtitle = paste0("Dimension ", dimension - 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.ticks = element_blank())
  return(p)
}

rotate_triangular <- function(triangle_ggplot)
{
  library(grid)
  triangle_ggplot <- triangle_ggplot + theme(legend.position = "none", 
                                             plot.margin = unit(c(1, 1, 1, 1), "cm"),
                                             axis.text.y = element_text(angle = 180, hjust = 1)) + 
    labs(title = element_blank(), subtitle = element_blank()) 
  q <- ggplot_build(triangle_ggplot)
  q$data[[2]]$angle <- 135
  q <- as.ggplot(ggplot_gtable(q))
  print(q, vp = viewport(angle = -135))
}

get_landscape <- function(homology, max_oe = 20, d = 1, min_x = "auto", max_x = "auto") {
  if (min_x == "auto"){
    t <- homology[, 2:3]
    min_x = min(t[t > 0]) / 2
  }
  if (max_x == "auto"){
    t <- homology[, 2:3]
    max_x = max(t[t > 0])
  }
  x.seq <- seq(min_x, max_x, length = 500)
  lands <- as.data.table(x.seq)
  for (oe in 1:max_oe){
    lands[, deparse(oe) := list(landscape(homology, dimension = d, 
                                          KK = oe, tseq = x.seq))]
    
  }
  colnames(lands) <- c("x", 1:max_oe)
  lands.melted <- melt(lands, "x")
  ggplot(lands.melted, aes(x = x, y = value, color = variable)) + geom_line() +
    labs(x = element_blank(), y = element_blank())  +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_grey()
}



plot_communities <- function(community_matrix) {
  community_matrix %>% 
    t() %>% 
    melt() %>% 
    ggplot(aes(x = Var1, y = fct_rev(Var2), fill = value)) +
    geom_tile(color = "white") +
    labs(x = element_blank(), y = element_blank(), 
         title = paste("Community structure of", deparse(substitute(community_matrix)))
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none") +
    scale_x_discrete(label = function (x) str_to_title(str_replace(x, "_", " ")),
                     expand = c(0,0)) + 
    scale_y_discrete(label = function (x) str_replace(x, "_", " "),
                     expand = c(0,0))
}

otu_mat_from_lists <- function(...){
  otu_mat <- cbind(...)
  row.names(otu_mat) <- paste("OTU", 1:nrow(otu_mat), sep = "_")
  
  return(otu_mat)
}

sim.sample.contig_presence <- function(num_otus, num_present){
  return(c(rep(1, num_present), rep(0, num_otus - num_present)))
}

sim.sample.tax_info <- function(num_otus){
  tax_mat <- matrix(sample(letters, 7 * num_otus, replace = TRUE), 
                    nrow = num_otus, ncol = 7,
                    dimnames = list(paste("OTU", 1:num_otus, sep = "_"), 
                                    c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
  
  return(tax_mat)
}

sim.sample.treatment_data <- function(otu_table){
  sample_matrix <- sample_data(data.frame(
    Treatment = rep("Control", ncol(otu_table)),
    row.names = colnames(sim_A.otu_table),
    stringsAsFactors = FALSE
  ))
  return(sample_matrix)
}

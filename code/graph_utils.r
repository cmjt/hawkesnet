#' A function to create a list of sparse 'time-evolving' adjacency matricies
#' given a vector of dates and a dated known edge list
#' @param x An ordered unique vector of dates/timestamps for each value you want an adjacency matrix
#' @param events A data frame of events, in long format, with the variables \code{date} giving
#' the day of the event and \code{id} with the unique character IDs of the marks/participants
#' involved in that events
#' @param edges A named dataframe of known edges with dates at which the edges were formed:
#' \code(date}, \code{from}, \code{to}.
#' @details Events are connected if they share a mark/id at a previous event and/or are connected historically
#' via edges. The function is set up so thatthe edge data may be from a different source.
#' @return A list of adjacency matricies for each unique event date
create_adjs <- function(events, edges){
    require(Matrix)
    x <- unique(events$date)
    n <- length(x)
    adjs <- list(n)
    adjs[[1]] <-  Matrix(data = 0, ncol = 1, nrow = 1)
    for(i in 2:n){
        adjs[[i]] <- cbind(rbind(adjs[[(i-1)]], rep(0, nrow(adjs[[(i-1)]]))), rep(0, i))
        historic_events <- events[events$date < x[i], ] |> unique()
        historic_relations <- edges[edges$date < x[i],]
        current_event <- events[events$date == x[i], ] |> unique()
        prior <- which(historic_events$id %in% current_event$id)
        if(length(prior) > 0) {
            pidx <- which(x %in% historic_events$date[prior])
            adjs[[i]][i, pidx] <- adjs[[i]][pidx, i] <- 1
        }
        con <- which(historic_relations$id %in% current_event$id)
        if(length(con) > 0) {
            who <- which(historic_events$id %in% historic_relations$to[con])
            if(length(who) > 0) {
                cidx <- which(x %in% historic_events$date[who])
                adjs[[i]][i, cidx] <- adjs[[i]][cidx, i] <- 1
            }
        }
    }
    return(adjs)
}


#' Functon to create a data frame with all edges and the time at which they were formed 
#' @param x A data frame where the first column contains all dated (timestamped) event times,
#' the remaining columns have the unique character IDs of the marks/participants involved in that events.
#' The number of columns equals the maximum of participants involved + 1, an event with fewer participants
#' than the maximum must have NA in the relevant columns.
#' @return A dataframe with all edges/links and the date at which they were formed
find_edges <- function(x){
    idx <- 2:ncol(x)
    dat <- as.matrix(x[, idx])
    elist <- do.call(rbind, apply(dat, 1,
                                  function(x) if (sum(!is.na(x)) > 1)
                                                  matrix(x[!is.na(x)][combn(sum(!is.na(x)), 2)],
                                                         ncol = 2, byrow = TRUE) else NULL))
    times  <- rep(x[,1], times = sapply(apply(dat, 1, function(x) sum(!is.na(x))),
                                        function(y) ifelse(y > 1, combn(y,2) |> ncol(), 0)))
    res <- data.frame(date = times, from = elist[,1], to = elist[,2])
    return(res)
}

#' Make a GIF of fr graphs from a list of adjacency matricies
#' where the layout is 'set' based on the final graphs
#' @param x A list of adjancency matricies or graphs (of class igraph)
#' @param filename File path to output GIF to
#'
make_gif <- function(x, filename = "MY_GIF.gif"){
    require(igraph)
    require(animation)
    n <- length(x)
    if(class(x[[1]]) == "igraph"){
        g <- x
    }else{
        g <- lapply(x,  graph_from_adjacency_matrix, mode = "undirected")
    }
    l <- layout_with_fr (g[[n]])
    saveGIF({
        for (i in 1:n){
            col = rep(c(adjustcolor("black", alpha.f = 1),
                        adjustcolor("black", alpha.f = 0)), c(i,(n-i)))
            plot(g[[i]], layout = l, vertex.color = col, vertex.frame.color = col,
                 vertex.label = NA, edge.size = 2, vertex.size = 3)
        }
    },movie.name = filename, interval = 0.1)
}
#' Function to plot adjacency matrices.
plot_adj <- function(x){
    require(ggplot2)
    require(igraph)
    require(reshape2)
    adj <- as_adjacency_matrix(x, sparse = FALSE)
    colnames(adj) <- rownames(adj) <- 1:nrow(adj)
    long <- melt(adj)
    p <- ggplot(long, aes(x = Var2, y = Var1)) + 
        geom_tile(aes(fill = as.factor(value)), col = "grey") + 
        scale_fill_manual(values= c( "white", "red"))  +
        theme_void() +
        theme(legend.position = "none")
    return(p)
}
#' A fuction to calculate the normalized closeness centrality
#' for each vertex of a graph \source{https://en.wikipedia.org/wiki/Closeness_centrality}
#'
#' Internal use
#' @param graph An object of class \code{igraph}
#' @param min_component the minimum component size to return as the "kitchen sink"
#' default 1
calculate_connectedness <- function(graph, min_component = 1) {
  closeness_vals <- igraph::closeness(graph, normalized = TRUE)
  closeness_vals[is.nan(closeness_vals)] <- 0
  components_info <- igraph::components(graph)
  ## making non-connected their own component
  member_info <- components_info$membership
 member_info[components_info$csize[member_info] <= min_component] <- "single"
  component_sizes <- igraph::sizes(components_info)
  vertex_component_sizes <- component_sizes[components_info$membership]
  closeness_scaled <- closeness_vals * vertex_component_sizes
  return(list("scaled_closeness" = as.numeric(closeness_scaled),
              "components" =  as.factor(member_info)))
}
#' A function that retains the k largest elements of a vector
#' and set others to zero
#'
#' Internal use
#' @param vec A numeric vector
#' @param k The value of the "largest" elements to retain
k_largest_elements <- function(vec, k) {
  if (k == 0) {
    return(rep(0, length(vec)))
  }
  if (k >= length(vec)) {
    return(vec)
  }
  threshold <- sort(vec, decreasing = TRUE)[k]
  vec[vec < threshold] <- 0
  return(vec)
}

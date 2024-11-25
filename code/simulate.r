#' A function that generates dynamic 'marks' accorinding
#' to historic connectedness
#' @inheritParams simulate_marked_hawkes
#' @author Conor Kresin
create_mark_generator <- function(N_min, epsilon, p_max, scale_factor) {
    env <- new.env()
    ## Initialize an empty graph with one vertex
    env$g <- igraph::make_empty_graph(n = 1, directed = FALSE)
    env$vertex_count <- 1
    ## Add the 'number' attribute to the initial vertex
    igraph::V(env$g)$number <- env$vertex_count 
    connectedness <- calculate_connectedness(env$g)
    ## Initialize marks as a list; marks_list[[i]] corresponds to vertex i
    env$marks_list <- list()
    env$marks_list[[env$vertex_count]] <- connectedness[env$vertex_count] * scale_factor
    ## Initialize graphs_list to store incremental graphs
    env$graphs_list <- list()
    env$graphs_list[[1]] <- env$g
    ## Initialize pending_seeds as an empty integer vector
    env$pending_seeds <- integer(0)
    ## Update marks based on current graph
    env$update_marks <- function() {
        connectedness <- calculate_connectedness(env$g)
        for (v in igraph::V(env$g)) {
            env$marks_list[[v]] <<- connectedness[v] * scale_factor
        }
    }
    ## Define the mark generator as a closure that modifies the graph in the environment
    mark_gen <- function() {
        env$vertex_count <- env$vertex_count + 1
        current_vertex <- env$vertex_count
        existing_vertices <- 1:(current_vertex - 1)
        ## Add the new vertex with its 'number' attribute
        env$g <- igraph::add_vertices(env$g, 1, attr = list(number = current_vertex))
        ## Determine if connecting to existing graph or starting a new component
        if (length(env$pending_seeds) > 0) {
            ## Connect to a randomly selected pending seed
            seed_to_connect <- sample(env$pending_seeds, 1)
            env$g <- igraph::add_edges(env$g, c(current_vertex, seed_to_connect))
            ## Remove the seed from pending_seeds
            env$pending_seeds <- env$pending_seeds[env$pending_seeds != seed_to_connect]
        } else {
            ## Decide whether to start a new component
            if (runif(1) < 1 / (current_vertex - 1)) { ## .1% probability of new component
                ## Start a new component: add to pending_seeds
                env$pending_seeds <- c(env$pending_seeds, current_vertex)
                ## Do not connect to any existing vertices
            } else {
                ## Connect to existing graph: ensure at least one connection
                ## Calculate connectedness for existing vertices
                existing_subgraph <- igraph::induced_subgraph(env$g, existing_vertices)
                connectedness <- calculate_connectedness(existing_subgraph)
                ## Determine p_base 
                N_star <- length(existing_vertices) #^2
                if (N_star > 1) {
                    p_base <- 1 / (N_star - 1) + epsilon
                } else {
                    p_base <- 1  ## If only one existing vertex, connect with probability 1
                }
                ## Dynamic probabilities
                max_connectedness <- ifelse(max(connectedness, na.rm = TRUE) == 0, 0,
                                            max(connectedness, na.rm = TRUE))
                if (max_connectedness == 0) {
                    prob_conn <- rep(0, length(existing_vertices))
                } else {
                    prob_conn <- connectedness / max_connectedness
                }
                probabilities <- pmax(prob_conn, p_base)
                ## Cap probabilities at p_max to prevent too many connections
                probabilities <- pmin(probabilities, p_max)
                probabilities <- pmin(pmax(probabilities, 0), 1)
                ## Determine the number of positive edges to add
                numPosEdges <- rpois(1, lambda = 1) + 1
                numPosEdges <- min(numPosEdges, length(existing_vertices))
                ## Select the top-k probabilities
                probabilities_active <- k_largest_elements(vec = probabilities, k = numPosEdges)
                ## Determine connections based on dynamic probabilities
                random_numbers <- runif(length(existing_vertices))
                connections <- existing_vertices[random_numbers < probabilities_active]
                ## Add edges if any connections are determined
                if (length(connections) > 0) {
                    edges_to_add <- as.vector(rbind(rep(current_vertex, length(connections)), connections))
                    env$g <- igraph::add_edges(env$g, edges_to_add)
                }
            }
        }
        ## Update marks based on the new graph
        env$update_marks()
        ## Append the current graph to graphs_list
        env$graphs_list[[current_vertex]] <- env$g
        ## Return the mark for the current vertex
        return(env$marks_list[[current_vertex]])
    }
    return(list(mark_gen = mark_gen, graph_env = env))
}
#' Function to simulate a marked Hawkes process
#' with time-evolving dynamic makes as per \code{create_mark_generator}
#' @param hawkes_params A named vector of marked Hawkes process
#' paramters: \code{lambda0}, the baseline rate,
#' \code{alpha}, the Hawkes jump in intensity,
#' \code{beta1}, temporal decay parameter, and
#' \code{beta2}, temporal decay parameter.
#' @param N_min The minimum number of vertecies to simulate, default \code{100}.
#' @param p_max The maximum connection probability [0,1], default \code{0.05}.
#' @param epsilon A threshold probability for mark connetion probability, default \code{0.01}.
#' @param T_max The maximum value to simulate Hawkes process, default \code{250}.
#' @param scale_factor The scaling factor for marks, default \code{100}.
#' @author Conor Kresin
simulate_marked_hawkes <- function(hawkes_params,
                                   N_min = 100, epsilon = 0.01,
                                   p_max = 0.05, T_max = 250,
                                   scale_factor = 100) {
    lambda0 <- hawkes_params["lambda0"]
    alpha <- hawkes_params["alpha"]
    beta1 <- hawkes_params["beta1"]
    beta2 <- hawkes_params["beta2"]
    ## Initialize the mark generator and graph environment
    mark_gen_obj <- create_mark_generator(N_min, epsilon, p_max, scale_factor)
    mark_gen <- mark_gen_obj$mark_gen
    env <- mark_gen_obj$graph_env  
    ## Initialize event list
    event_times <- c()
    event_vertices <- c()
    ## Initialize the list to store Hawkes process states up to each event
    hawkes_process_list <- list()
    ## Initialize the first event
    current_time <- 0
    event_times <- c(event_times, current_time)
    ## Generate the mark for the first event
    mark <- mark_gen()  ## This adds the first vertex and its mark
    event_vertices <- c(event_vertices, env$vertex_count)
    ## Save the first Hawkes process state
    hawkes_process_list[[1]] <- data.frame(time = event_times[1], 
                                           vertex = event_vertices[1], 
                                           mark = env$marks_list[[event_vertices[1]]])
    ## Event generation loop
    while (current_time < T_max && length(event_times) < N_min) {
        ## Calculate the current intensity
        if (length(event_times) == 0) {
            lambda_t <- lambda0
        } else {
            ## Sum over all past events
            lambda_t <- lambda0
            for (i in 1:length(event_times)) {
                decay_time <- current_time - event_times[i]
                if (decay_time >= 0) {
                    lambda_t <- lambda_t +
                        alpha * exp(-beta1 * decay_time) * exp(-beta2 * env$marks_list[[event_vertices[i]]])
                }
            }
        }
        ## Upper bound for intensity
        lambda_upper <- lambda0 + alpha * sum(exp(-beta2 * unlist(env$marks_list)))
        ## Handle potential overflow in lambda_upper
        if (lambda_upper > 1e10) {  
            lambda_upper <- 1e10
            warning("lambda_upper capped to prevent overflow.")
        }
        ## Generate next event time candidate
        u <- runif(1)
        w <- -log(u) / lambda_upper
        t_candidate <- current_time + w
        if (t_candidate > T_max) {
            break
        }
        ## Calculate intensity at candidate time
        if (length(event_times) == 0) {
            lambda_candidate <- lambda0
        } else {
            lambda_candidate <- lambda0
            for (i in 1:length(event_times)) {
                decay_time <- t_candidate - event_times[i]
                if (decay_time >= 0) {
                    lambda_candidate <- lambda_candidate +
                        alpha * exp(-beta1 * decay_time) * exp(-beta2 * env$marks_list[[event_vertices[i]]])
                }
            }
        }
        ## Accept or reject the candidate
        d <- runif(1)
        if (d <= lambda_candidate / lambda_upper) {
            ## Accept the event
            current_time <- t_candidate
            event_times <- c(event_times, current_time)
            ## Generate the mark using the mark generator
            mark <- mark_gen()
            event_vertices <- c(event_vertices, env$vertex_count)
            ## Save the current Hawkes process state
            hawkes_process_list[[length(event_times)]] <- data.frame(
                time = event_times,
                vertex = event_vertices,
                mark = sapply(event_vertices, function(v) env$marks_list[[v]])
            )
        } else {
            ## Reject the event and continue
            current_time <- t_candidate
        }
    }
    ## Extract the final graph from the environment
    final_graph <- env$g
    ## Put marks into a vector aligned with event order
    event_marks <- sapply(event_vertices, function(v) env$marks_list[[v]])
    ## Retrieve the list of all Hawkes process states
    return(list(
        hawkes_process_list = hawkes_process_list,
        graphs_list = env$graphs_list,
        data = data.frame(time = event_times, vertex = event_vertices, mark = event_marks),
        graph = final_graph
    ))
}


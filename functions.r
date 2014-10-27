# TODO
# - return trace in gm.search

# gm.search(observed: table, graph.init: binary square matrix, forward: bool, backward: bool, score: string)
#   : list(model: list of cliques, score: numeric, call: string)
# Perform a hill climbing search on graphical models to determine the best fitted model
# for the observed data, according to the specified score function.
gm.search = function(observed, graph.init, forward = TRUE, backward = TRUE, score = 'bic')
{
    graph <- graph.init
    modelscore <- graph.assess(observed, graph, score)
    print(paste('Starting with score', modelscore))

    # Repeatedly construct and evaluate neighborhood until no improvement can be made
    repeat
    {
        best_neighbor <- NULL
        best_neighbor_score <- Inf

        # Loop over all edges
        l <- nrow(graph)
        for (v in 1:(l-1))
        {
            for (w in (v+1):l)
            {
                if (graph[v,w] == 0 && forward || graph[v,w] && backward)
                {
                    # Transform graph into a neighbor
                    graph[v,w] <- 1 - graph[v,w]
                    graph[w,v] <- 1 - graph[w,v]

                    # Evaluate the neighbor
                    neighbor_score <- graph.assess(observed, graph, score)
                    if (neighbor_score < modelscore)
                    {
                        best_neighbor <- graph
                        best_neighbor_score <- neighbor_score
                    }

                    # Undo the transformation
                    graph[v,w] <- 1 - graph[v,w]
                    graph[w,v] <- 1 - graph[w,v]
                }
            }
        }

        # Return if no improvement else continue with best neighbors
        if (!is.null(best_neighbor))
        {
            print(paste('Improving to', best_neighbor_score))
            graph <- best_neighbor
            modelscore <- best_neighbor_score
        }
        else
            return(list(model = bronkerbosch(graph), score = modelscore, call = match.call()))
    }
}

# gm.restart(nstart: int, prob: numeric, seed: numeric, observed: table, forward: bool, backward: bool, score: string)
#   : list(model: list of cliques, score: numeric, call: string)
# Repeatedly run gm.search on randomly generated initial graphs.
gm.restart = function(nstart, prob, seed, observed, forward = TRUE, backward = TRUE, score = 'bic')
{
    set.seed(seed)
    graphsize <- log(length(observed), 2)

    best_result <- NULL
    best_score <- Inf
    for(n in 1:nstart)
    {
        print(paste('Run', n))
        graph.init = matrix(0, graphsize, graphsize)

        # Loop over all edges
        l <- nrow(graph.init)
        for (v in 1:(l-1))
        {
            for (w in (v+1):l)
            {
                # Decide randomly if edge should be added
                r <- runif(1, 0.0, 1.0)
                if (r > prob)
                {
                    graph.init[v,w] <- 1
                    graph.init[w,v] <- 1
                }
            }
        }

        # Perform hill climbing search
        result <- gm.search(observed, graph.init, forward, backward, score)
        if (result$score < best_score)
        {
            best_result <- result
            best_score <- result$score
        }
    }

    return(best_result)
}

# graph.assess(observed: table, graph: binary square matrix, score: numeric)
#   : numeric
# Compute the AIC or BIC (whichever specified) for the given graph w.r.t. the observed data.
graph.assess = function(observed, graph, score)
{
    # Detect the cliques from the graph
    cliques <- bronkerbosch(graph)

    # Use loglin to generate fitted model to the cliques
    model <- loglin(observed, cliques, print = FALSE)

    # Compute score and return
    if (score == 'aic')
        return(AIC(model, observed, cliques))
    else if (score == 'bic')
        return(BIC(model, observed, cliques))
    else
        stop(paste('Invalid score parameter', score))
}

# AIC(model: list returned by loglin, observed: table, cliques: list of integer vectors)
#   : numeric
# Compute the AIC score for the given model w.r.t. the given data.
AIC = function(model, observed, cliques)
{
    dimension <- length(observed) - model$df
    return(model$lrt + 2 * dimension)
}

# BIC(model: list returned by loglin, observed: table, cliques: list of integer vectors)
#   : numeric
# Compute the BIC score for the given model w.r.t. the given data.
BIC = function(model, observed, cliques)
{
    dimension <- length(observed) - model$df
    return(model$lrt + log(sum(observed)) * dimension)
}

# bronkerbosch(graph: binary square matrix): list of integer vectors
# Compute all cliques in given graph
bronkerbosch = function(graph)
{
    recursion <- function(include, rest, exclude)
    {
        assert(length(include) > 0 || length(rest) > 0 || length(exclude) > 0)
        assert(length(intersect(include, rest)) == 0)
        assert(length(intersect(include, exclude)) == 0)
        assert(length(intersect(rest, exclude)) == 0)

        if (length(exclude) == 0 && length(rest) == 0)
            return(list(include))

        result <- list()
        for (v in rest)
        {
            neighbors <- which(graph[v,] == 1)
            recursive_result <- recursion(union(include, v),
                                          intersect(rest, neighbors),
                                          intersect(exclude, neighbors))
            rest <- setdiff(rest, v)
            exclude <- union(exclude, v)

            result <- c(result, recursive_result)
        }
        return(result)
    }

    return(recursion(c(), 1:nrow(graph), c()))
}

# Utils
assert = function(bool) if (!bool) stop('Assertion error')


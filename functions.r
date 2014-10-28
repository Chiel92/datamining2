# gm.search(observed: table, graph.init: binary square matrix, forward: bool, backward: bool, score: string)
#   : list(model: list of cliques, score: numeric, call: string)
# Perform a hill climbing search on graphical models to determine the best fitted model
# for the observed data, according to the specified score function.
gm.search = function(observed, graph.init, forward = TRUE, backward = TRUE, score = 'bic')
{
    modelscore <- graph.assess(observed, graph.init, score)
    trace <- paste('Starting with score', modelscore)
    print(trace)

    # Repeatedly construct and evaluate neighborhood until no improvement can be made
    repeat
    {
        graph <- graph.init
        best_neighbor <- NULL
        best_neighbor_score <- NULL
        best_neighbor_msg <- NULL

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
                        best_neighbor_msg <- paste('Flip edge', v, w, 'yielding score', neighbor_score)
                    }

                    # Undo the transformation
                    graph[v,w] <- 1 - graph[v,w]
                    graph[w,v] <- 1 - graph[w,v]
                }
            }
        }

        # Return if no improvement found else continue with best neighbor
        if (is.null(best_neighbor))
        {
            #names <- c("1: cat1", "2: death", "3: swang1", "4: gender")
            names <- c("1: cat1", "2: death", "3: swang1", "4: gender", "5: race", "6: ninsclas", "7: income", "8: ca", "9: age", "10: meanbp1")
            rownames(graph.init) <- names
            colnames(graph.init) <- names

            # Plot the graph
            amgraph <- new("graphAM", adjMat = graph.init, edgemode = "undirected")
            plot(am.graph, attrs = list(node = list(fillcolor = "lightblue", height=2), edge = list()))

            return(list(model = bronkerbosch(graph.init), score = modelscore, trace = trace, call = match.call()))
        }
        else
        {
            print(best_neighbor_msg)
            graph.init <- best_neighbor
            modelscore <- best_neighbor_score
            trace <- c(trace, best_neighbor_msg)
        }
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
        pivot <- resample(union(rest, exclude), 1)
        for (v in setdiff(rest, neighbors(graph, pivot)))
        {
            recursive_result <- recursion(union(include, v),
                                          intersect(rest, neighbors(graph, v)),
                                          intersect(exclude, neighbors(graph, v)))
            rest <- setdiff(rest, v)
            exclude <- union(exclude, v)

            result <- c(result, recursive_result)
        }
        return(result)
    }

    return(recursion(c(), 1:nrow(graph), c()))
}

# neighbors(graph: binary square matrix, vertex: int): vector of integers
# Compute the neighborhood of a vertex.
neighbors = function(graph, vertex) which(graph[vertex,] == 1)

# resample(x: vector, ...): subvector of x
# Same as sample, but does not treat single integers as ranges
resample <- function(x, ...) x[sample.int(length(x), ...)]

# Utils
assert = function(bool) if (!bool) stop('Assertion error')


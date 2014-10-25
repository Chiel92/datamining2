# Utils
assert = function(bool) if (!bool) stop('Assertion error')

graph.empty = function(size)
{
    return(matrix(replicate(size * size, 0), size))
}

# graph.search
graph.search = function(observed, graph)
{
    score <- graph.assess(observed, graph)

    best_neighbor <- NULL
    best_neighbor_score <- 0

    while(TRUE)
    {
        # Construct neighborhood
        for (v in 2:nrow(graph))
        {
            for (w in v:nrow(graph))
            {
                # Transform graph into a neighbor
                graph[v,w] <- 1 - graph[v,w]
                graph[w,v] <- 1 - graph[w,v]

                # Evaluate the neighbor
                neighbor_score <- graph.assess(graph)
                if (neighbor_score > score)
                {
                    best_neighbor <- graph
                    best_neighbor_score <- neighbor_score
                }

                # Undo the transformation
                graph[v,w] <- 1 - graph[v,w]
                graph[w,v] <- 1 - graph[w,v]
            }
        }

        # Return if no improvement else continue with best neighbors
        if (!is.null(best_neighbor))
        {
            graph <- best_neighbor
            score <- best_neighbor_score
        }
        else
            return(graph)
    }

}

graph.assess = function(observed, graph)
{
    # Detect the cliques from the graph
    cliques <- bronkerbosch(graph)

    # Use loglin to generate fitted model to the cliques
    model <- loglin(observed, cliques)

    # Compute score and return
}

bronkerbosch = function(graph)
{
    recursion <- function(include, rest, exclude)
    {
        assert(length(include) > 0 | length(rest) > 0 | length(exclude) > 0)
        assert(length(intersect(include, rest)) == 0)
        assert(length(intersect(include, exclude)) == 0)
        assert(length(intersect(rest, exclude)) == 0)

        if (length(exclude) == 0 & length(rest) == 0)
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

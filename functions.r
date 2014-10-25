# Utils
assert = function(bool) if (!bool) stop('Assertion error')

graph.empty = function(size)
{
    return(matrix(replicate(size * size, 0), size))
}

# graph.search
graph.search = function(observed, graph)
{
    # Detect the cliques from the graph
    cliques <- bronkerbosch(graph)

    # Use loglin to generate fitted model to the cliques

    # Construct neighbourhood and advance
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

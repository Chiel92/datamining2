# gm.search(observed: table, graph.init: binary square matrix, forward: bool, backward: bool, score: string)
#   : list(model: list of cliques, score: numeric, call: string)
# Perform a hill climbing search on graphical models to determine the best fitted model
# for the observed data, according to the specified score function.
gm.search = function(observed, graph.init, forward = TRUE, backward = TRUE, score = 'bic', print = TRUE)
{
    bestgraph <- graph.init
    bestscore <- graph.assess(observed, bestgraph, score)
    trace <- paste('Starting with score', bestscore)
    if (print)
        print(trace)

    # Repeatedly construct and evaluate neighborhood until no improvement can be made
    repeat
    {
        graph <- bestgraph
        best_neighbor <- NULL
        score_threshold <- bestscore

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
                    if (neighbor_score < score_threshold)
                    {
                        score_threshold <- neighbor_score

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
            # Plot the graph
            if (print)
            {
                names <- c("1: cat1", "2: death", "3: swang1", "4: gender", "5: race", "6: ninsclas", "7: income", "8: ca", "9: age", "10: meanbp1")
                rownames(bestgraph) <- names
                colnames(bestgraph) <- names
                my_plot(bestgraph)
            }

            return(list(model = bronkerbosch(bestgraph), score = bestscore, trace = trace, call = match.call(), graph=bestgraph))
        }
        else
        {
            if (print)
                print(best_neighbor_msg)
            bestgraph <- best_neighbor
            bestscore <- best_neighbor_score
            trace <- c(trace, best_neighbor_msg)
        }
    }
}

my_plot = function(graph)
{
    names <- c("1: cat1", "2: death", "3: swang1", "4: gender", "5: race", "6: ninsclas", "7: income", "8: ca", "9: age", "10: meanbp1")
    #names <- c("1: cat1", "2: death", "3: swang1", "4: gender", "5: race")
    rownames(graph) <- names
    colnames(graph) <- names
    amgraph <- new("graphAM", adjMat = graph, edgemode = "undirected")
    plot(amgraph, attrs = list(node = list(fillcolor = "lightblue", height=2), edge = list()))
}


randomgraph = function(prob, seed)
{
    set.seed(seed)
    graph.init = matrix(0, 10, 10)

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

    return(graph.init)
}


get_data = function()
{
    print('---------BIC----------')
    generate_nstart_prob_plot('bic', 'bic_plot.txt')
    print('---------AIC----------')
    generate_nstart_prob_plot('aic', 'aic_plot.txt')
}


generate_nstart_prob_plot = function(score, filename)
{
    csv_output <- ""
    runs <- 1:10
    probvalues <- seq(from = 0.2, to = 0.8, by = 0.2)
    for(prob in probvalues)
    {
        best_result <- NULL
        best_score <- Inf
        for(run in runs)
        {
            seed = run
            graph.init <- randomgraph(prob, seed)
            result <- gm.search(observed, graph.init, forward = TRUE, backward = TRUE, score = score, print = FALSE)

            if (result$score < best_score)
            {
                best_result <- result
                best_score <- result$score
                best_seed <- seed
            }

            msg <- paste(toString(run), ', ', toString(prob), ', ', toString(best_score))
            print(msg)
            csv_output <- cbind(csv_output, msg)
        }

        print(paste('With prob=', toString(prob), 'best score is', toString(best_score), 'with seed=', toString(seed)))
    }

    write(csv_output, file = filename)
}

draw_nstart_prob_plot = function()
{
    plot.dat <- read.csv('aic_plot.txt', header=FALSE, col.names=c('nstart', 'prob', 'score'))
    levelplot(score ~ nstart + prob, data=plot.dat, col.regions=topo.colors)
}

#get_data = function()
#{
    #result <- ""
    #nstartvalues <- c(1, 5, 10)
    #probvalues <- c(0.4, 0.5, 0.6)

    #for(nstart in nstartvalues)
        #for(prob in probvalues)
            #for(score in list('aic', 'bic'))
            #{
                #print(paste(nstart, prob))
                #score_val <- gm.restart(nstart, prob, 37, observed, score=score, print = FALSE)$score
                #msg <- paste(toString(nstart), ', ', toString(prob), ', ', score, '=', toString(score_val))
                #print(msg)
                #result <- cbind(result, msg)
            #}

    #write(result, file = 'output.txt')
#}

# gm.restart(nstart: int, prob: numeric, seed: numeric, observed: table, forward: bool, backward: bool, score: string)
#   : list(model: list of cliques, score: numeric, call: string)
# Repeatedly run gm.search on randomly generated initial graphs.
gm.restart = function(nstart, prob, seed, observed, forward = TRUE, backward = TRUE, score = 'bic', graphsize = 10, print = TRUE)
{
    set.seed(seed)
    #graphsize <- log(length(observed), 2)

    best_result <- NULL
    best_score <- Inf
    for(n in 1:nstart)
    {
        if (print)
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
        result <- gm.search(observed, graph.init, forward, backward, score, print = print)
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
graph.assess = function(observed, graph, score='bic')
{
    # Detect the cliques from the graph
    cliques <- bronkerbosch(graph)

    # Use loglin to generate fitted model to the cliques
    model <- loglin(observed, cliques, print = FALSE)
    dimension <- length(observed) - model$df

    # Compute score and return
    if (score == 'aic')
        return(model$lrt + 2 * dimension)
    else if (score == 'bic')
        return(model$lrt + log(sum(observed)) * dimension)
    else
        stop(paste('Invalid score parameter', score))
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


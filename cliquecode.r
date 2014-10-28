find.cliques <- function (R,P,X,graph,cliques) 
{
    if (length(P)==0 & length(X)==0) {cliques <- list(R)}
    else {
        pivot <- P[sample(length(P),1)]
        for(i in 1:length(P)){
            pivot.nb <- neighbors(graph,pivot)
            if(!is.element(P[i],pivot.nb)){
                P.temp <- setdiff(P,P[i])
                R.new <- union(R,P[i])
                P.new <- intersect(P.temp,neighbors(graph,P[i]))
                X.new <- intersect(X,neighbors(graph,P[i]))
                cliques <- c(cliques,find.cliques(R.new,P.new,X.new,graph,list()))
                X <- union(X,P[i])}
    }}
    cliques
}


neighbors <- function (graph,node) 
{
    nnodes <- dim(graph)[2]
    (1:nnodes)[graph[node,]==1]
}

post.process <- function (cliques) 
{
    unique(lapply(cliques,sort))
}

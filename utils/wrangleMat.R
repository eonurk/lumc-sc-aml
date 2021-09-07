
#' Removes sorts genes, removes duplicated ones and renames the rownames
#' @param matrix genes by (1+samples): +1 is being the gene names
#' @return genes by samples matrix
wrangleMat <- function(matrix){
    # Remove spike-ins for now (may not exists in your data)
    matrix <- matrix [!grepl("^ERCC", matrix[,1]),]
    
    # Eliminate any homologs
    # TODO could be dangerous to do it this way, find a better version...
    matrix[,1] <- sapply(strsplit(matrix[,1], ".", fixed = TRUE), function(x){x[1]})
    
    # Order genes according to their standard deviation in decreasing order
    matrix <- matrix [rev(order(apply(matrix[,-1], 1, stats::sd))),]
    
    # Remove duplicated genes
    matrix <- matrix [!duplicated(matrix[,1]),]
    
    # Make the gene names the row names
    rownames(matrix) <- matrix[,1]
    
    # Filter the genes
    matrix <- matrix[,-1]
    
    # Enforce all counts to be integers
    matrix <- round(matrix, 0)
    
    return(matrix)
}
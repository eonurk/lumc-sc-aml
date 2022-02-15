
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

#' Given the count matrix calculates stemness score for all patients
#' @param counts genes by samples (rownames must be ensembl ids)
#' @references https://sci-hub.se/https://www.nature.com/articles/nature20598
#' @return genes by samples matrix
calculateStemness <- function(counts, gene.symbol = F){
    
    # stemness genes
    if(gene.symbol){
        stem_genes <- c("DNMT3B", "ZBTB46", "NYNRIN", "ARHGAP22", "LAPTM4B" ,
                        "MMRN1", "DPYSL3","KIAA0125", "CDK6", "CPXM1", "SOCS2", 
                        "SMIM24", "EMP1", "NGFRAP1", "CD34", "AKR1C3", "GPR56")
    } else{
        stem_genes <- c("ENSG00000088305", "ENSG00000130584", "ENSG00000205978", 
                        "ENSG00000128805", "ENSG00000104341", "ENSG00000138722", 
                        "ENSG00000113657", "ENSG00000277988", "ENSG00000105810", 
                        "ENSG00000088882", "ENSG00000120833", "ENSG00000095932",
                        "ENSG00000134531", "ENSG00000166681", "ENSG00000174059", 
                        "ENSG00000196139", "ENSG00000205336")
    }
    
    stem_coefficients <- c(0.0874, -0.0347, 0.00865, -0.0138, 0.00582, 0.0258,
                           0.0284, 0.0196, -0.0704, -0.0258, 0.0271, -0.0226, 
                           0.0146, 0.0465, 0.0338, -0.0402, 0.0501)
    
    
    norm.counts <- log2(edgeR::cpm(counts, log = F) + 1)
    return(apply(norm.counts[stem_genes,], 2, function(x){sum(x * stem_coefficients)}))
}


#' Given the count matrix calculates stemness score for all patients
#' @param counts genes by samples (rownames must be ensembl ids or gene symbols)
#' @references https://www.nature.com/articles/s41467-021-22625-y#Sec28
#' @return genes by samples matrix
calculateAPS <- function(counts, gene.symbol = F){
    
    # stemness genes
    if(gene.symbol){
        aps_genes <- c("CD109","CALCRL","CPNE3","MBNL1","ZNF883",	
                        "CCDC6","XYLB","TMEM273","AGRN","SAMD11","MICALL2",	
                        "VSTM1","PPEF1-AS1","LOC100133091","ZSCAN25","CERKL")
    } else{
        aps_genes <- c("ENSG00000156535","ENSG00000064989",
                        "ENSG00000085719","ENSG00000152601",
                        "ENSG00000228623","ENSG00000108091",
                        "ENSG00000093217","ENSG00000204161",
                        "ENSG00000188157","ENSG00000187634",
                        "ENSG00000164877","ENSG00000189068",
                        "ENSG00000237221","ENSG00000205485",
                        "ENSG00000197037","ENSG00000188452")
    }
    
    aps_coefficients <- c(0.103439048,0.100259195,0.098030944,0.082724702,
                           0.076161262,0.069639173,0.043183815,0.039548972,
                           -0.007627419,-0.037055078,-0.057267696,-0.070597079,
                           -0.077714874,-0.104276549,-0.118225984,-0.136155589)
                           
    
    
    norm.counts <- log2(edgeR::cpm(counts, log = F) + 1)
    return(apply(norm.counts[aps_genes,], 2, function(x){sum(x * aps_coefficients)}))
}



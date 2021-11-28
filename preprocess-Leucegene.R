file.list <- list.files("data/Other studies/LEUCEGENE/GSE67040_RAW/", full.names = T)

a <- lapply(file.list, function(x){
    
    y <- fread(x)
    return(y$V4)
}) 

intersected.genes <- Reduce(intersect, a)


cm <- lapply(file.list, function(x){
    
    y <- fread(x)
    res <- as.data.frame(y[match(intersected.genes, y$V4), V6])
    colnames(res) <- strsplit(x, "_", fixed = T)[[1]][3]
    rownames(res) <- y$V4[match(intersected.genes, y$V4)]
    return(res)
}) %>% do.call(cbind, .)

write.csv(as.data.frame(cm),file = "data/Other studies/LEUCEGENE/Leucegene_count_matrix.csv", quote = F )


meta.LEUCE <- as.data.frame(fread("data/Other studies/LEUCEGENE/Leucegene_Annotations.tsv"))
rownames(meta.LEUCE) <- meta.LEUCE$sample_id

meta.LEUCE <- meta.LEUCE[colnames(cm),]

write.csv(meta.LEUCE, file = "data/Other studies/LEUCEGENE/Leucegene_meta_data.csv")



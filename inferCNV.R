library(infercnv)

patients <- list.dirs("data/InferCNV/", full.names = F)[-1]


lapply(patients, function(pt){
    
    ref.list <- c("GMP", "Mono", "Prog", "cDC", "HSC", "ProMono", "earlyEry", 
      "CTL", "T", "Plasma", "ProB", "NK", "B", "lateEry", "pDC")
    
    meta <- read.table(paste0("data/InferCNV/", pt, "/VanGalen_AML_meta.txt"))
    data <- read.table(paste0("data/InferCNV/", pt, "/VanGalen_AML_expr.txt"))
    # Not all samples include all ref cell types
    ref.list <- intersect(meta[,2], ref.list)
    
    if(all(unique(meta[,2]) %in% ref.list)){
        message("\t>> Passing ", pt, " as there is no cancerous cells!")
        return(NULL)
    }
    
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0("data/InferCNV/", pt, "/VanGalen_AML_expr.txt"),
                                        annotations_file=paste0("data/InferCNV/", pt, "/VanGalen_AML_meta.txt"),
                                        delim="\t",
                                        gene_order_file="data/InferCNV/gencode_v19_gene_pos.txt",
                                        ref_group_names=ref.list)
    
    output.dir <- paste0("output/InferCNV/", pt)
    dir.create(output.dir)
    
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff= 1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir=output.dir, 
                                 cluster_by_groups=TRUE, 
                                 denoise=TRUE,
                                 HMM=FALSE)
    
    rm(infercnv_obj, meta)
    gc()
})


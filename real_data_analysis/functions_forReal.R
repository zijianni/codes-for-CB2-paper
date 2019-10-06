barcode_diff <- function(dat_CB2, dat_ED, printing = T, names = c("CB2", "ED")) {
  
  Common_b <- intersect(colnames(dat_ED), colnames(dat_CB2))
  CB2_only <- setdiff(colnames(dat_CB2), colnames(dat_ED))
  ED_only <- setdiff(colnames(dat_ED), colnames(dat_CB2))
  res <- c(length(Common_b), length(CB2_only), length(ED_only))
  names(res) <- c("Common barcodes", paste0("Only in ", names[1]), paste0("Only in ", names[2]))
  if (printing)
    print(res)
  return(list(Common_b = Common_b, CB2_only = CB2_only, ED_only = ED_only, Sum = res))
}

down_sample <- function(x, size) {
  res <- rep(0, length(x))
  names(res) <- names(x)
  tmp <- rep(names(x), x)
  tmp_down <- table(sample(tmp, size))
  res[names(tmp_down)] <- tmp_down
  return(res)
}

k_means <- function(dat, k = 10, Npcs = 100, prefilter = T) {
  if (prefilter) {
    nzero <- which(Matrix::rowSums(dat) != 0)
    dat <- dat[nzero, ]
  }
  
  dat <- t(dat)
  obj_pca <- prcomp_irlba(dat, n = Npcs)
  dat <- obj_pca$x
  return(kmeans(dat, centers = k)$cluster)
}

run_tsne <- function(dat, perplexity = 30, verbose = F, PCA = T, Npcs = 100, prefilter = T, PCA_scale = F) {
  bname <- colnames(dat)
  bc <- Matrix::colSums(dat)
  
  if (prefilter) {
    nzero <- which(Matrix::rowSums(dat) != 0)
    dat <- dat[nzero, ]
  }
  dat <- t(dat)
  if (PCA) {
    obj_pca <- prcomp_irlba(dat, n = Npcs, scale. = PCA_scale)
    dat <- obj_pca$x
  } else{
    dat <- as.matrix(dat)
  }
  if (anyDuplicated(dat)) {
    k <- ncol(dat)
    for (i in which(duplicated(dat))) {
      dat[i, ] <- dat[i, ] + rnorm(k, 0, 10e-5)
    }
  }
  obj <- Rtsne(dat, verbose = verbose, perplexity = perplexity)
  obj <- obj$Y
  rownames(obj) <- bname
  return(obj)
}



exp_heatmap <-
  function(dat,
           main = NULL,
           filter_threshold = 1,
           labCol = F,
           ...) {
    dat <- dat[Matrix::rowSums(dat) >= filter_threshold, ]
    dat <- log(as.matrix((dat)) + 1)
    
    colors = c(seq(min(dat, na.rm = T), max(dat, na.rm = T), length = 200))
    my_palette <- colorRampPalette(c("black", "red"))(n = 199)
    #
    heatmap.2(
      dat,
      col = my_palette,
      breaks = colors,
      dendrogram = 'none',
      symm = F,
      tracecol = NA,
      density.info = "none",
      lwid = c(1, 4),
      lhei = c(1, 6),
      labCol = labCol,
      main = ifelse(
        is.null(main),
        paste0(
          " \nlog(x+1) transformed expression, row column reordered. ",
          dim(dat)[1],
          " genes ",
          dim(dat)[2],
          " cells.\nGene filter threshold=",
          filter_threshold
        ),
        main
      ),
      ...
    )
}

cor_heatmap <- function(dat,
                        main = NULL,
                        dendrogram = 'none',
                        ...) {
  colors = c(seq(-1, 1, length = 200))
  my_palette <-
    colorRampPalette(c("green", "black", "red"))(n = 199)
  cor_temp <- sparse_cor(dat)
  heatmap.2(
    cor_temp,
    col = my_palette,
    breaks = colors,
    dendrogram = dendrogram,
    symm = T,
    tracecol = NA,
    lwid = c(1, 4),
    density.info = "none",
    lhei = c(1, 6),
    main = paste0(main, " \ncorrelation, reordered. ", dim(dat)[2], " cells"),
    ...
  )
}

SIMFUN_bg <- function(raw.mat,
                      threshold = 100,
                      remove_protein = T) {
  if (remove_protein) {
    protein <- grep(pattern = "TotalSeqB", x = rownames(raw.mat))
    if (length(protein) > 0) {
      raw.mat <- raw.mat[-protein, ]
    }
    
  }
  
  stats <- barcodeRanks(raw.mat, lower = threshold)
  
  # Assuming all cells below the inflection point are empty droplets.
  totals <- stats$total
  empty.limit <- stats@metadata$inflection
  ambient.prof <- rowSums(raw.mat[, totals <= empty.limit])
  empty.totals <- totals[totals <= empty.limit]
  
  # Scrambling to generate resampled empty droplets.
  gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  gene.ids <- sample(gene.ids)
  cell.ids <- rep(seq_along(empty.totals), empty.totals)
  resampled <-
    makeCountMatrix(gene.ids, cell.ids, all.genes = rownames(raw.mat))
  colnames(resampled) <- paste0("bg_", 1:ncol(resampled))
  
  # Completed.
  return(resampled)
}


SIM_bg <- function(dat,
                   threshold = 100,
                   remove_protein = T) {
  if (remove_protein) {
    protein <- grep(pattern = "TotalSeqB", x = rownames(dat))
    if (length(protein) > 0) {
      dat <- dat[-protein, ]
    }
    
  }
  brank <- barcodeRanks(dat, lower = threshold)
  inflection <-
    ifelse(is.null(brank$inflection),
           brank@metadata$inflection,
           brank$inflection)
  
  dat_bg1 <- dat[, Matrix::colSums(dat) %in% (1:inflection)]
  g_name <- rownames(dat_bg1)
  bc_dist <- Matrix::colSums(dat_bg1)
  #str(bc_dist)
  bg_dist <- Matrix::rowSums(dat_bg1)
  #str(bg_dist)
  #bg_pool <- rep(names(bg_dist),bg_dist)
  bg_prop <- goodTuringProportions(bg_dist)
  bc_t <- table(bc_dist)
  
  ##dat_bg1 is the new background
  dat_tmp <- dat_bg1[, 0]
  for (i in 1:length(table(bc_dist))) {
    #print(i)
    if (bc_t[i] %/% 5000 > 0) {
      for (j in 1:(bc_t[i] %/% 5000)) {
        #print(j)
        dat_tmp <-
          cbind(dat_tmp, rmultinom(5000, names(bc_t)[i], bg_prop))
      }
      dat_tmp <-
        cbind(dat_tmp, rmultinom(bc_t[i] %% 5000, names(bc_t)[i], bg_prop))
    } else{
      dat_tmp <- cbind(dat_tmp, rmultinom(bc_t[i], names(bc_t)[i], bg_prop))
      
    }
    
    #if(i%%10000==0) print(i)
  }
  
  dat_bg1 <- dat_tmp
  colnames(dat_bg1) <- paste0("bg_", 1:ncol(dat_bg1))
  
  return(dat_bg1)
  
}

reorder_mat <- function(mat, rate = 0.1) {
  N_row <- as.integer(nrow(mat) * rate)
  x <- 1:nrow(mat)
  reorder_row <- sample(x, N_row)
  x[sort(reorder_row)] <- reorder_row
  return(mat[x, ])
}

filter_mr <- function(x) {
  tmp <- grep(pattern = "\\<RP", x = x)
  if (length(tmp) > 0) {
    x <- x[-grep(pattern = "\\<RP", x = x)]
  }
  tmp <- grep(pattern = "\\<MT-", x = x)
  if (length(tmp) > 0) {
    x <- x[-grep(pattern = "\\<MT-", x = x)]
  }
  tmp <- grep(pattern = "\\<Rp", x = x)
  if (length(tmp) > 0) {
    x <- x[-grep(pattern = "\\<Rp", x = x)]
  }
  tmp <- grep(pattern = "\\<mt-", x = x)
  if (length(tmp) > 0) {
    x <- x[-grep(pattern = "\\<mt-", x = x)]
  }
  return(x)
}

Upsample <- function(x, target_count) {
  p <- x / sum(x)
  return(rmultinom(1, target_count, p))
}

SIM2 <-
  function(n = 5,
           dat,
           dat_bg = NULL,
           threshold = 100,
           n_large = 2000,
           n_middle = 2000,
           n_small = 2000,
           new_sim = T,
           remove_protein = T,
           FDR = 0.01,
           reorder_rate = 0.1,
           threegroup = T,
           isEBDrop = F,
           usePackage = T,
           ...) {
    if (remove_protein) {
      protein <- grep(pattern = "TotalSeqB", x = rownames(dat))
      if (length(protein) > 0) {
        dat <- dat[-protein, ]
      }
    }
    brank <- barcodeRanks(dat, lower = threshold)
    inflection <-
      ifelse(is.null(brank$inflection),
             brank@metadata$inflection,
             brank$inflection)
    
    ###background barcodes
    if (is.null(dat_bg)) {
      dat_bg <- SIM_bg(dat, threshold = threshold)
    }
    
    
    all_bg <- sum(Matrix::colSums(dat) %in% 1:inflection)
    th_bg <- sum(Matrix::colSums(dat) %in% (threshold + 1):inflection)
    
    dat_true1 <- dat[, Matrix::colSums(dat) > inflection]
    
    
    TP_CB2 <- matrix(0, 6, n + 1)
    colnames(TP_CB2) <- c(paste0("rep", 1:n), "Ave")
    rownames(TP_CB2) <-
      c("G1", "G1 prop", "G2", "G2 prop", "G3", "G3 prop")
    
    FP_CB2 <- TP_CB2[1:4, ]
    rownames(FP_CB2) <-
      c(
        paste0("all (", all_bg, ")"),
        "all prop",
        paste0(">", threshold, " (", th_bg, ")"),
        paste0(">", threshold, " prop")
      )
    
    TP_ED <- TP_CB2
    FP_ED <- FP_CB2
    
    SEED.rep <- as.integer(runif(n, 1, 100000))
    
    for (sss in 1:n) {
      SEED.tmp <- SEED.rep[sss]
      
      cat(paste0("Run ", sss, "\n"))
      
      
      set.seed(SEED.tmp)
      
      
      
      ###2000 G1 large cells
      G1 <- dat_true1[, sample(ncol(dat_true1), n_large, replace = T)]
      colnames(G1) <- paste0("G1_", 1:ncol(G1))
      G1 <- as(G1, "dgCMatrix")
      
      if (threegroup) {
        ###2000 G2 middle cells (50% downsampling)
        G2 <- dat_true1[, sample(ncol(dat_true1), n_middle, replace = T)]
        G2 <- as(G2, "dgCMatrix")
        G2 <- downsampleMatrix(G2, 0.5)
        colnames(G2) <- paste0("G2_", 1:ncol(G2))
      }
      
      
      ###2000 G3 small cells (10% downsampling)
      G3 <- dat_true1[, sample(ncol(dat_true1), n_small, replace = T)]
      G3 <- as(G3, "dgCMatrix")
      G3 <- downsampleMatrix(G3, 0.1)
      colnames(G3) <- paste0("G3_", 1:ncol(G3))
      
      if (new_sim) {
        G1 <- reorder_mat(G1, rate = reorder_rate)
        if (threegroup)
          G2 <- reorder_mat(G2, rate = reorder_rate)
        G3 <- reorder_mat(G3, rate = reorder_rate)
      }

      if (threegroup) {
        dat_sim1 <- cbind(dat_bg, G1, G2, G3)
      } else{
        dat_sim1 <- cbind(dat_bg, G1, G3)
      }
      
      
      
      set.seed(SEED.tmp)
        if(usePackage){
            res_CB2 <-
                scCB2::CB2FindCell(dat_sim1,
                          pooling_threshold = threshold,
                          FDR_threshold = FDR,
                          ...)
        }else{
            res_CB2 <-
                Find_Cell(dat_sim1,
                          pooling_threshold = threshold,
                          FDR_threshold = FDR,
                          ...)
        }

        res1 <-
          c(colnames(res_CB2$cluster_matrix),
            colnames(res_CB2$cell_matrix))

      
      
      
      TP_CB2[1, sss] <- length(grep(x = res1, pattern = "G1"))
      TP_CB2[2, sss] <- TP_CB2[1, sss] / n_large
      
      TP_CB2[3, sss] <- length(grep(x = res1, pattern = "G2"))
      TP_CB2[4, sss] <- TP_CB2[3, sss] / n_middle
      
      TP_CB2[5, sss] <- length(grep(x = res1, pattern = "G3"))
      TP_CB2[6, sss] <- TP_CB2[5, sss] / n_small
      
      FP_CB2[1, sss] <- length(grep(x = res1, pattern = "bg"))
      FP_CB2[2, sss] <- FP_CB2[1, sss] / all_bg
      
      FP_CB2[3, sss] <- length(grep(x = res1, pattern = "bg"))
      FP_CB2[4, sss] <- FP_CB2[3, sss] / th_bg
      
      
      
      
      eOut <- emptyDrops(dat_sim1, lower = threshold)
      sim_ED <- dat_sim1[, ifelse(is.na(eOut$FDR), FALSE, eOut$FDR <= FDR)]
      res_ED2 <- colnames(sim_ED)
      
      TP_ED[1, sss] <- length(grep(x = res_ED2, pattern = "G1"))
      TP_ED[2, sss] <- TP_ED[1, sss] / n_large
      
      TP_ED[3, sss] <- length(grep(x = res_ED2, pattern = "G2"))
      TP_ED[4, sss] <- TP_ED[3, sss] / n_middle
      
      TP_ED[5, sss] <- length(grep(x = res_ED2, pattern = "G3"))
      TP_ED[6, sss] <- TP_ED[5, sss] / n_small
      
      FP_ED[1, sss] <- length(grep(x = res_ED2, pattern = "bg"))
      FP_ED[2, sss] <- FP_ED[1, sss] / all_bg
      
      FP_ED[3, sss] <- length(grep(x = res_ED2, pattern = "bg"))
      FP_ED[4, sss] <- FP_ED[3, sss] / th_bg
      
    }
    
    
    TP_CB2[, n + 1] <- rowMeans(TP_CB2[, 1:n])
    TP_ED[, n + 1] <- rowMeans(TP_ED[, 1:n])
    FP_CB2[, n + 1] <- rowMeans(FP_CB2[, 1:n])
    FP_ED[, n + 1] <- rowMeans(FP_ED[, 1:n])
    
    Stat <- matrix(NA, 2, 2)
    colnames(Stat) <- c("CB2", "ED")
    rownames(Stat) <- c("Power", "FDR")
    Stat[1, 1] <-
      (TP_CB2[2, n + 1] * n_large + TP_CB2[4, n + 1] * n_middle + TP_CB2[6, n +
                                                                           1] * n_small) / (n_large + n_middle + n_small)
    Stat[1, 2] <-
      (TP_ED[2, n + 1] * n_large + TP_ED[4, n + 1] * n_middle + TP_ED[6, n +
                                                                        1] * n_small) / (n_large + n_middle + n_small)
    Stat[2, 1] <-
      FP_CB2[3, n + 1] / (FP_CB2[3, n + 1] + Stat[1, 1] * (n_large + n_middle +
                                                             n_small))
    Stat[2, 2] <-
      FP_ED[3, n + 1] / (FP_ED[3, n + 1] + Stat[1, 2] * (n_large + n_middle +
                                                           n_small))
    Power_diff <- Stat[1, 1] - Stat[1, 2]
    
    return(
      list(
        TP_CB2 = TP_CB2,
        FP_CB2 = FP_CB2,
        TP_ED = TP_ED,
        FP_ED = FP_ED,
        Stat = Stat,
        Power_diff = Power_diff
      )
    )
  }

Attenuate <- function(x){
    q95 <- quantile(x,0.95)
    x[x>q95] <- q95
    return(x)
}

Htmap <- function(dat,set, Allgene=T, background_threshold=100, Plot="both"){
    dat <- as.matrix(dat)

    if(set=="mbrain"){
        diff_set <- diff_mbrain
        dat_r<- mbrain_r
    }else if(set=="MALT"){
        diff_set <- diff_MALT
        dat_r <- MALT_r
    }else if(set=="PBMC8K"){
        diff_set <- diff_PBMC8K
        dat_r <- PBMC8K_r
    }else if(set=="amy761"){
        diff_set <- diff_amy761
        dat_r<- amy761_raw
    }else if(set=="amy762"){
        diff_set <- diff_amy762
        dat_r<- amy762_raw
    }else if(set=="CionaC"){
        diff_set <- diff_CionaC
        dat_r<- CionaC_r
    }else if(set=="CionaT"){
        diff_set <- diff_CionaT
        dat_r<- CionaT_r
    }else if(set=="mbrain1K"){
        diff_set <- diff_mbrain1K
        dat_r<- mbrain1K_r
    }else if(set=="placenta"){
        diff_set <- diff_placenta
        dat_r<- placenta_r
    }else if(set=="Alz"){
        diff_set <- diff_Alz
        dat_r<- Alz_r
    }
    
    
    Group <- ifelse(colnames(dat)%in%diff_set$CB2_only,"Not captured","Captured")
    
    dat_C <- dat[,Group=="Captured"]
    dat_NC <- dat[,Group=="Not captured"]
    dat_NC <- dat_NC[,order(colSums(dat_NC),decreasing = T)]
    
    EDonly <- colnames(dat_C)%in%diff_set$ED_only
    if(sum(EDonly)>0){
        bset3 <- colnames(dat_C)[EDonly]
    }else{
        bset3 <- NULL
    }

    if(Plot%in%c("both","D")){
        dist_plot(colnames(dat_NC), colnames(dat_C)[!EDonly], bset3,
                  colnames(dat_r)[Matrix::colSums(dat_r)<=background_threshold],rownames(dat),
                  dat=dat_r, Allgene = Allgene)
    }

    
    if(Plot%in%c("both","H")){
        if(!all(Group=="Not captured")){
            
            if(ncol(dat)>400){
                ratio_C <-ncol(dat_C)/ncol(dat)
                ratio_NC <- 1-ratio_C
                
                dat_C <- dat_C[,sample(ncol(dat_C),400*ratio_C)]
                dat_NC <- dat_NC[,sample(ncol(dat_NC),400*ratio_NC)]
            }
            
            dat_heatmap0 <- cbind(dat_C,dat_NC)
            dat_heatmap0 <- t(apply(dat_heatmap0,1,Attenuate))
            dat_heatmap1 <- melt(dat_heatmap0)
            colnames(dat_heatmap1)[1:2] <- c("Gene","Barcode")
            
            Cor <- round(cor(rowSums(dat_C),rowSums(dat_NC)),3)
            
            gp <- ggplot(dat_heatmap1, aes(x=Barcode, y=Gene)) + geom_tile(aes(fill = log(value+1)),colour = "white") + 
                scale_fill_gradient(low = "white",high = "red") + 
                theme(axis.text.y=element_text(size=5), axis.text.x = element_blank()) + 
                geom_vline(xintercept = ncol(dat_C)+0.5)+
                xlab(paste0("Barcode (Left: common cells  Right: CB2 extra cells)\nCorrelation: ",Cor))
            plot(gp)
        }
    }

    # plot(ksmooth(x=1:nrow(dat_C),y=rowMeans(dat_bg[rownames(dat_C),])),type="l",ylab="Ave expr",xlab="Gene index", main="Background")
    # if(ncol(dat_C)!=0) plot(rowMeans(dat_C),type="h",ylab="Ave expr",xlab="Gene index", main="Captured by ED")
    # plot(rowMeans(dat_NC),type="h",ylab="Ave expr",xlab="Gene index", main="Not captured by ED")

}

dist_plot <- function(bset1,bset2=NULL,bset3=NULL,bgset,gene,dat,Allgene = T){
    if(Allgene){
        # dist1 <- (Matrix::rowSums(dat[,bset1])/sum(dat[,bset1]))[gene]
        # if(!is.null(bset2)) dist2 <- (Matrix::rowSums(dat[gene,bset2])/sum(dat[gene,bset2]))[gene]
        # distBG <- (Matrix::rowSums(dat[gene,bgset])/sum(dat[gene,bgset]))[gene]
        dist1 <- edgeR::goodTuringProportions(Matrix::rowSums(dat[,bset1,drop=F]))[gene,1]
        if(!is.null(bset2)){
            dist2 <- edgeR::goodTuringProportions(
                Matrix::rowSums(dat[,bset2,drop=F]))[gene,1] 
        }
        if(!is.null(bset3)){
            dist3 <- edgeR::goodTuringProportions(
                Matrix::rowSums(dat[,bset3,drop=F]))[gene,1] 
        }
        distBG <- edgeR::goodTuringProportions(Matrix::rowSums(dat[,bgset]))[gene,1]

    }
    # else{
    #     dist1 <- Matrix::rowSums(dat[gene,bset1])/sum(dat[gene,bset1])
    #     if(!is.null(bset2)) dist2 <- Matrix::rowSums(dat[gene,bset2])/sum(dat[gene,bset2])
    #     distBG <- Matrix::rowSums(dat[gene,bgset])/sum(dat[gene,bgset])
    #     
    # }
    dist_dat <- data.frame(gene=gene,index=seq_along(gene),dist1=dist1,distBG=distBG)
    if(!is.null(bset2)) dist_dat$dist2 <- dist2
    if(!is.null(bset3)) dist_dat$dist3 <- dist3  
    if(!is.null(bset2)){
        gp2 <- ggplot(data=dist_dat)+theme(legend.position="bottom")+
            geom_bar(stat="identity",aes(x=index, y=dist2,fill="Common"),alpha=0.5) +
            xlab("gene index")+ylab("probability")+ylim(0,max(dist_dat[,-(1:2)]))+ 
            scale_fill_manual(name= "", values=c("Common"="#00c746")) 
        plot(gp2)
    }

    
    if(!is.null(bset3)){
        gp3 <- ggplot(data=dist_dat)+theme(legend.position="bottom")+
            geom_bar(stat="identity",aes(x=index, y=dist3,fill="unique in ED"),alpha=0.5) +
            xlab("gene index")+ylab("probability")+ylim(0,max(dist_dat[,-(1:2)]))+ 
            scale_fill_manual(name= "", values=c("unique in ED"="#56B4E9")) 
        plot(gp3)
    }

    
    gp <- ggplot(data=dist_dat)+theme(legend.position="bottom")+
        geom_bar(stat="identity",aes(x=index, y=dist1,fill="unique in CB2"),alpha=0.5)+
        xlab("gene index")+ylab("probability")+ylim(0,max(dist_dat[,-(1:2)]))+ 
        scale_fill_manual(name= "", values=c("unique in CB2"="#FF9999"))
    plot(gp)
    
    gp_bg <- ggplot(data=dist_dat)+ylim(0,max(dist_dat[,-(1:2)]))+
                geom_bar(stat="identity",aes(x=index,y=distBG,fill="Background"),alpha=0.5)+
        xlab("gene index")+ylab("probability") + scale_fill_manual(name= "", 
            values=c("Background"="#000000"))+ theme(legend.position="bottom")
    plot(gp_bg)
}

enrich_matlist <- function(x, term, topNterm=20, topNgene=200, filter=T){
    if(filter){
        for(i in seq_along(x)){
            RPgene <- grep("\\<RPS|\\<RPL|\\<Rps|\\<Rpl",x[[i]])
            MTgene <- grep("\\<MT-|\\<mt-",x[[i]])
            if(length(c(RPgene,MTgene))!=0){
                x[[i]] <- x[[i]][-c(RPgene,MTgene)]
            }

        }
    }
    enrich_list <- list()
    for(i in seq_along(x)){
        gene_temp <- x[[i]][1:min(topNgene,length(x[[i]]))]
        enrich_list[[i]] <- enrichr(gene_temp, term)
    }
    
    term_list <- list()
    for(j in seq_along(term)){
        enrich_by_term <- lapply(enrich_list,function(x) x[[j]])
        term_all_geneset <- lapply(enrich_by_term, function(x) x$Term[1:topNterm])
        term_unique <- unique(unlist(term_all_geneset))
        term_mat <- matrix(0,length(term_unique),length(x))
        colnames(term_mat) <- names(x)
        rownames(term_mat) <- term_unique
        for(k in seq_along(x)){
            term_mat[term_all_geneset[[k]],k] <- 1
        }
        term_list[[j]] <- term_mat
    }
    names(term_list) <- term
    return(term_list)
}

enrich_heatmap <- function(term_list,...){
    if(is.matrix(term_list)){
        gplots::heatmap.2(term_list, scale = "none",  col = c("white","#FF9999"),
            tracecol = F, key=F,...)
    }else{
        for(i in seq_along(term_list)){
            gplots::heatmap.2(term_list[[i]], scale = "none",  col = c("white","#FF9999"),
                              tracecol = F, main=names(term_list)[i], key=F)
        }
    }
}


mk_plot <- function(mk_gene, prop=F,colorUp=T, diff_dat, dat_tsne_n1, dat_f){
    Ent_mat1 <- dat_f[,diff_dat$CB2_only[diff_dat$CB2_only%in%colnames(dat_f)]]
    
    com_mat1 <- dat_f[,diff_dat$Common_b[diff_dat$Common_b%in%colnames(dat_f)]]
    
    ED_mat1 <- dat_f[,diff_dat$ED_only[diff_dat$ED_only%in%colnames(dat_f)]]
    
    
    
    if(mk_gene=="count"){
        mk <- Matrix::colSums(cbind(com_mat1,Ent_mat1,ED_mat1))
    }else if(mk_gene=="MT"){
        
        mk <- Matrix::colSums(cbind(com_mat1,Ent_mat1,
                                    ED_mat1)[grep("\\<MT-",rownames(com_mat1)),])
    }else if(mk_gene=="mt"){
        
        mk <- Matrix::colSums(cbind(com_mat1,Ent_mat1,
                                    ED_mat1)[grep("\\<mt-",rownames(com_mat1)),])
    }else{
        mk <- cbind(com_mat1,Ent_mat1,ED_mat1)[mk_gene,]
    }
    
    mk[mk>0] <- Attenuate(mk[mk>0])
    tsne_dat <- data.frame(cbind(dat_tsne_n1,mk,Matrix::colSums(cbind(com_mat1,Ent_mat1,ED_mat1))))
    colnames(tsne_dat)[c(1,2,ncol(tsne_dat)-1:0)] <- c("x","y","marker","count")
    
    
    if(prop){
        p2 <- ggplot()+guides(alpha = FALSE)+
            geom_point(data = tsne_dat[mk==0,],aes(x=x,y=y),size=0.2,color="#f2f2f2")
        if(!colorUp){
            p2$layers <- c(geom_point(data = tsne_dat[mk!=0,],aes(x=x,y=y,color=log(marker/count+10e-5),
                                                                   alpha=log(marker/count+10e-5))), p2$layers)
        }else{
            p2$layers <- c(p2$layers, geom_point(data = tsne_dat[mk!=0,],aes(x=x,y=y,color=log(marker/count+10e-5),
                                                                              alpha=log(marker/count+10e-5))))
        }
        
        p2 <- p2+scale_color_gradient(low="blue", high="red",name = paste0("log prop\n ",mk_gene))+ 
            theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
        print(p2)
        
    }else{
        p1 <- ggplot()+
            geom_point(data = tsne_dat[mk==0,],aes(x=x,y=y),size=0.2,color="#f2f2f2")+
            guides(alpha = FALSE)
        if(!colorUp){
            p1$layers <- c(geom_point(data = tsne_dat[mk!=0,],aes(x=x,y=y,color=marker,alpha=marker),size=0.2), p1$layers)
        }else{
            p1$layers <- c(p1$layers, geom_point(data = tsne_dat[mk!=0,],aes(x=x,y=y,color=marker,alpha=marker),size=0.2))
        }
        p1 <- p1+scale_color_gradient(low="blue", high="red",name = mk_gene)+ 
            theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
        print(p1)
    }
}

check_range <- function(x, lower, upper){
    return( (x >= lower) & (x <= upper))
}

CB2Cluster <- function(CB2Out, tsneOut,sizeFilter=0){
    CB2Clust <- list()
    for(i in seq_along(CB2Out$Cluster)){
        if( (CB2Out$ClusterStat$padj[i]<=0.01)&&
            (length(CB2Out$Cluster[[i]])>=sizeFilter)){
            CB2Clust <- c(CB2Clust, list(CB2Out$Cluster[[i]]))
        }
    }
    names(CB2Clust) <- paste0("C",seq_along(CB2Clust))
    
    bcluster <- rep("0",nrow(tsneOut))
    knee_bc <- setdiff(colnames(CB2Out$cell_matrix),
        rownames(CB2Out$BarcodeStat)[CB2Out$BarcodeStat$padj<=0.01])
    
    for(i in 1:nrow(tsneOut)){
        bc <- rownames(tsneOut)[i]
        bc_clust <- unlist(lapply(CB2Clust,function(x) bc%in%x))
        if(any(bc_clust)){
            bcluster[i] <- names(bc_clust)[bc_clust]
        }
        if(bc%in%knee_bc) bcluster[i] <- "knee"
    }

    tsneOut$CB2Cluster <- bcluster
    
    return(tsneOut)
}

Calc_retain <- function(dat,lower){
    dat <- FilterGB(dat)
    totals <- unname(colSums(dat))
    o <- order(totals, decreasing=TRUE)
    
    stuff <- rle(totals[o])
    # Get mid-rank of each run.
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 
    run.totals <- stuff$values
    
    keep <- run.totals > lower
    if (sum(keep)<3) { 
        stop("insufficient unique points for computing knee/inflection points")
    }
    
    x <- log10(run.rank[keep])
    fit <- smooth.spline(log10(run.rank), log10(run.totals), df=20)
    y <- predict(fit)$y[keep]
    #y <- log10(run.totals[keep])
    # Numerical differentiation to identify bounds for spline fitting.
    # The upper/lower bounds are defined at the plateau and inflection, respectively.
    d1n <- diff(y)/diff(x)
    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])
    # We restrict to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    new.keep <- left.edge:right.edge
    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee point.
    d1 <- predict(fit, deriv=1)$y[keep][new.keep]
    d2 <- predict(fit, deriv=2)$y[keep][new.keep]
    curvature <- d2/(1 + d1^2)^1.5
    knee <- 10^(y[new.keep][which.min(curvature)])
    inflection <- 10^(y[right.edge])
    return(list(knee=round(knee),inflection=round(inflection)))
}

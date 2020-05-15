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


Attenuate <- function(x){
    q95 <- quantile(x,0.95)
    x[x>q95] <- q95
    return(x)
}

Htmap <- function(dat,set, Allgene=T, background_threshold=100, 
                  Plot="both", AddBackground=F, equal=F, attenuate=F){
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
    
    # b_BG <-colnames(dat_r)[(!colnames(dat_r)%in%unlist(diff_set)) & (colSums(dat_r)>100)]
    # dat_BG <- dat_r[rownames(dat),sample(b_BG,100)]
    
    Group <- ifelse(colnames(dat)%in%diff_set$CB2_only,"Not captured","Captured")
    
    dat_C <- dat[,Group=="Captured"]
    dat_NC <- dat[,Group=="Not captured"]

    
    EDonly <- colnames(dat_C)%in%diff_set$ED_only
    if(sum(EDonly)>0){
        bset3 <- colnames(dat_C)[EDonly]
    }else{
        bset3 <- NULL
    }

    if(Plot%in%c("both","D")){
        dist_plot(colnames(dat_NC), colnames(dat_C)[!EDonly], bset3,
                  colnames(dat_r)[Matrix::colSums(dat_r)<=background_threshold],rownames(dat),
                  dat=dat_r, Allgene = Allgene, AddBackground=AddBackground)
    }

    
    if(Plot%in%c("both","H")){
        if(!all(Group=="Not captured")){
            
         
            if(ncol(dat)>400){
                ratio_C <-ncol(dat_C)/ncol(dat)
                ratio_NC <- 1-ratio_C

                dat_C <- dat_C[,sample(ncol(dat_C),400*ratio_C)]
                dat_NC <- dat_NC[,sample(ncol(dat_NC),400*ratio_NC)]
            }

            if(equal&ncol(dat_C)>ncol(dat_NC)){
                set.seed(1)
                dat_C <- dat_C[,sample(ncol(dat_C),ncol(dat_NC))]
            }
            
            if(attenuate){
                dat_C <- Attenuate(dat_C)
                dat_NC <- Attenuate(dat_NC)
            }
            
            dat_NC <- dat_NC[,order(Matrix::colSums(dat_NC),decreasing = T)]
            dat_C <- dat_C[,order(Matrix::colSums(dat_C),decreasing = T)]
            
            dat_heatmap0 <- cbind(dat_C,dat_NC)
            dat_heatmap0 <- t(apply(dat_heatmap0,1,Attenuate))
            dat_heatmap1 <- melt(dat_heatmap0)
            colnames(dat_heatmap1)[1:2] <- c("Gene","Barcode")
            
            Cor <- round(cor(Matrix::rowSums(dat_C),Matrix::rowSums(dat_NC)),3)
            #scale_fill_gradient(low = "white",high = "deeppink1")
            gp <- ggplot(dat_heatmap1, aes(x=Barcode, y=Gene)) + geom_tile(aes(fill = log(value+1)),colour = "white", size=10e-5) + 
                scale_fill_gradient(low = "white",high = "green3") + labs(fill="")+
                theme(axis.text.y=element_text(size=5), axis.text.x = element_blank(),
                      axis.title=element_blank()) + 
                geom_vline(xintercept = ncol(dat_C)+0.5, linetype="dashed", size=1, col="black")
            plot(gp)
        }
    }

    # plot(ksmooth(x=1:nrow(dat_C),y=rowMeans(dat_bg[rownames(dat_C),])),type="l",ylab="Ave expr",xlab="Gene index", main="Background")
    # if(ncol(dat_C)!=0) plot(rowMeans(dat_C),type="h",ylab="Ave expr",xlab="Gene index", main="Captured by ED")
    # plot(rowMeans(dat_NC),type="h",ylab="Ave expr",xlab="Gene index", main="Not captured by ED")
}

# Extra plots for manuscript revision
Htmap_revision <- function(barcodeID,dat, background_threshold=100, AddBackground=F){
    
    Top100 <- (names(sort(rowSums(dat[,barcodeID]),decreasing = T)))[1:100]

    dist_plot(barcodeID, NULL, NULL,
            colnames(dat)[Matrix::colSums(dat)<=background_threshold],sort(Top100),
            dat=dat, AddBackground=AddBackground)
}

mk_hist <- function(mk_gene,cluster_cell,dat_r,diff_dat, upper){
    knee_cell <- colnames(dat_r)[colSums(dat_r)>=upper]
    knee_cell <- intersect(knee_cell,cluster_cell)
    common_cell <- setdiff(diff_dat$Common_b,knee_cell)
    common_cell <- intersect(common_cell,cluster_cell)
    CB2_cell <- intersect(c(diff_dat$CB2_only),cluster_cell)
    
    exp_CB2 <- dat_r[,colnames(dat_r)%in%CB2_cell,drop=F] %>%
        rowSums %>%
        goodTuringProportions
    exp_common <- dat_r[,colnames(dat_r)%in%common_cell,drop=F] %>%
        rowSums %>% 
        goodTuringProportions
    exp_knee <- dat_r[,colnames(dat_r)%in%knee_cell,drop=F] %>%
        rowSums %>% 
        goodTuringProportions
    
    dat_bg <- dat_r[,colSums(dat_r)<=100]
    exp_bg <- dat_bg %>% rowSums %>% goodTuringProportions
    
    mk_df <- data.frame(rbind(exp_knee[mk_gene,1],exp_common[mk_gene,1],exp_CB2[mk_gene,1],exp_bg[mk_gene,1]),
                        group=c("common (high count)","common (tested)","unique in CB2","Background"),check.names = F)
    mk_df$group <- factor(mk_df$group,levels = rev(c("common (high count)","common (tested)","unique in CB2","Background")))
    

    mk_df <- melt(mk_df,id="group")
    gp <- ggplot(data=mk_df, aes(x=variable, y=value,fill=group, width=0.8)) + 
        geom_bar(stat="identity", position = "dodge")+coord_flip()+
               scale_fill_manual(values=rev(c( "deeppink3", "deeppink1","gold","#937C37")))+  theme_bw() + 
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), 
        axis.title = element_blank(), legend.position = "none")
    # print(gp)
    return(gp)
    
}

dist_plot <- function(bset1,bset2=NULL,bset3=NULL,bgset,gene,dat,
                      Allgene = T, AddBackground=T){

    distBGall <- edgeR::goodTuringProportions(Matrix::rowSums(dat[,bgset]))[,1]
    if(AddBackground){
        gene <- c(gene, filter_mr(setdiff(names(sort(distBGall,decreasing=T)),gene ))[1:25])
    }
    
    dist1 <- edgeR::goodTuringProportions(Matrix::rowSums(dat[,bset1,drop=F]))[gene,1]
    if(!is.null(bset2)){
        dist2 <- edgeR::goodTuringProportions(
            Matrix::rowSums(dat[,bset2,drop=F]))[gene,1] 
    }
    if(!is.null(bset3)){
        dist3 <- edgeR::goodTuringProportions(
            Matrix::rowSums(dat[,bset3,drop=F]))[gene,1] 
    }
    distBG <- distBGall[gene]

    dist_dat <- data.frame(gene=gene,index=seq_along(gene),dist1=dist1,distBG=distBG)
    if(!is.null(bset2)) dist_dat$dist2 <- dist2
    if(!is.null(bset3)) dist_dat$dist3 <- dist3  
    if(!is.null(bset2)){
        gp2 <- ggplot(data=dist_dat)+theme(legend.position="bottom")+
            geom_bar(stat="identity",aes(x=index, y=dist2,fill="Common")) +
            xlab("Gene Index")+ylab("Probability")+ylim(0,max(dist_dat[,-(1:2)]))+ 
            scale_fill_manual(name= "", values=c("Common"="deeppink1")) + 
            theme_bw() + theme(legend.position = "bottom", 
                               panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
        plot(gp2)
    }

    
    if(!is.null(bset3)){
        gp3 <- ggplot(data=dist_dat)+theme(legend.position="bottom")+
            geom_bar(stat="identity",aes(x=index, y=dist3,fill="unique in ED")) +
            xlab("Gene Index")+ylab("Probability")+ylim(0,max(dist_dat[,-(1:2)]))+ 
            scale_fill_manual(name= "", values=c("unique in ED"="black")) + 
            theme_bw() + theme(legend.position = "bottom", 
                               panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
        plot(gp3)
    }

    
    gp <- ggplot(data=dist_dat)+theme(legend.position="bottom")+
        geom_bar(stat="identity",aes(x=index, y=dist1,fill="unique in CB2"))+
        xlab("Gene Index")+ylab("Probability")+ylim(0,max(dist_dat[,-(1:2)]))+ 
        scale_fill_manual(name= "", values=c("unique in CB2"="yellow"))+ 
        theme_bw() + theme(legend.position = "bottom", 
                           panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    plot(gp)
    
    gp_bg <- ggplot(data=dist_dat)+ylim(0,max(dist_dat[,-(1:2)]))+
                geom_bar(stat="identity",aes(x=index,y=distBG,fill="Background"))+
        xlab("Gene Index")+ylab("Probability") + scale_fill_manual(name= "", 
            values=c("Background"="#937C37"))+ theme(legend.position="bottom")+ 
        theme_bw() + theme(legend.position = "bottom", 
                           panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    plot(gp_bg)
}

mk_plot <- function(mk_gene, prop=F, diff_dat, dat_tsne_n1, dat_f,CB2=F){
    if(!CB2) diff_dat <- NULL
    
    dat_tsne_n1 <- dat_tsne_n1[colnames(dat_f),]
    
    if(mk_gene=="count"){
        mk <- Matrix::colSums(dat_f)
    }else if(mk_gene=="MT"){
        
        mk <- Matrix::colSums(dat_f[grep("\\<MT-",rownames(dat_f)),])
    }else if(mk_gene=="mt"){
        
        mk <- Matrix::colSums(dat_f[grep("\\<mt-",rownames(dat_f)),])
    }else{
        mk <- dat_f[mk_gene,]
    }
    
    mk[mk>0] <- Attenuate(mk[mk>0])
    tsne_dat <- data.frame(dat_tsne_n1,marker=mk,count=Matrix::colSums(dat_f))

    if(CB2){
        tsne_dat <- tsne_dat[rownames(dat_tsne_n1)%in%diff_dat$CB2_only,]
        mk <- tsne_dat$marker
    }
    
    p <- ggplot()+guides(alpha = FALSE)+
        geom_point(data = tsne_dat[mk==0,],aes(x=tSNE1,y=tSNE2),size=0.2,color="#f2f2f2")
    
    if(prop){
        p$layers <- c(p$layers, geom_point(data =tsne_dat[mk!=0,],
                                           aes(x=tSNE1,y=tSNE2,color=log(marker/count+10e-5),
                                               alpha=log(marker/count+10e-5))))
        nm <- paste0("log prop\n ",mk_gene)
    }
    else{
        p$layers <- c(p$layers, geom_point(data = tsne_dat[mk!=0,],
                                           aes(x=tSNE1,y=tSNE2,color=marker,alpha=marker),size=0.2))
        nm <- mk_gene
    }
    
    p <- p+scale_color_gradient(low="orange", high="green",name = nm)+ 
        theme_bw() + 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.8), 
              legend.position = "none",
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"))+ 
        xlim(min(tsne_dat$tSNE1),max(tsne_dat$tSNE1)) +
        ylim(min(tsne_dat$tSNE2),max(tsne_dat$tSNE2))
    
    print(p)
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

DE <- function(celltype="ExN",gene, diff=T,tsne=Alz_tsne2_all, prop=F){
    temp_celltype <- rownames(tsne)[tsne$celltype==celltype]
    
    #temp_celltype <- temp_celltype[colSums(Alz_r[,temp_celltype])>= brange[1] & 
    #                                   colSums(Alz_r[,temp_celltype])<= brange[2]]
    
    if(diff){
        temp_celltype <- setdiff(temp_celltype, temp_common)
    }
    temp_AD <- rownames(tsne)[tsne$Condition=="AD-pathology"]
    temp_C <- rownames(tsne)[tsne$Condition=="Control"]
    
    exp_AD <- Alz_fAll3[gene,intersect(temp_celltype,temp_AD)]
    exp_C <- Alz_fAll3[gene,intersect(temp_celltype,temp_C)]
    
    if(prop){
        exp_AD <- exp_AD/tsne[intersect(temp_celltype,temp_AD),]$size
        exp_C <- exp_C/tsne[intersect(temp_celltype,temp_C),]$size
    }
    # exp_AD <- sample(exp_AD,n)
    # exp_C <- sample(exp_C,n)
    # message(paste0("Sample size: ", n))
    
    tt <- wilcox.test(exp_C,exp_AD)
    #mean(exp_AD)/mean(exp_C)

    return(list(pval=round(-log10(tt$p.value),5), 
                log2FC=log2(mean(exp_AD)/mean(exp_C)),
                size=length(exp_C)+length(exp_AD),
                test=tt))
}



KL <- function(x,bg){ 
    x <- x/sum(x)
    bg <- bg/sum(bg)
    sum(x*log(x/bg),na.rm = T)
}

Ent <- function(x){
    x <- x/sum(x)
    -sum(x*log(x),na.rm = T)
}

reorder_mat <- function(mat, rate = 0.1) {
    N_row <- as.integer(nrow(mat) * rate)
    x <- 1:nrow(mat)
    reorder_row <- sample(x, N_row)
    x[sort(reorder_row)] <- reorder_row
    return(mat[x, ])
}

#This function refers to EmptyDrops' background simulation codes
SIMFUN_bg <- function(raw.mat,  
                      threshold = 100,
                      remove_protein = T) {
    # For Cell Ranger version >=3, proteins will be removed.
    # Reasons are discussed in CB2 paper and scCB2 package.
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
             ED_knee = T,
             run=T,
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
        
        knee <-
            ifelse(is.null(brank$knee),
                   brank@metadata$knee,
                   brank$knee)
        
        ###background barcodes
        if (is.null(dat_bg)) {
            dat_bg <- SIMFUN_bg(dat, threshold = threshold)
        }
        
        all_bg <- sum(Matrix::colSums(dat) %in% 1:inflection)
        th_bg <- sum(Matrix::colSums(dat) %in% (threshold + 1):inflection)
        
        dat_true1 <- dat[, Matrix::colSums(dat) > inflection]
        dat_knee <- dat[, Matrix::colSums(dat) >= knee]
        
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
        
        Size <- matrix(NA,3,n+1)
        rownames(Size) <- c("G1","G2", "G3")
        colnames(Size) <- c(1:n, "Average")
        
        KL_Diver <- Size
        Entropy <- Size
        
        prob_bg <- edgeR::goodTuringProportions(rowSums(dat_bg))[,1]
        Top5000 <- prob_bg>=sort(prob_bg,decreasing = T)[5000]
        prob_knee <- edgeR::goodTuringProportions(rowSums(dat_knee))[,1]
        prob_realbg <- edgeR::goodTuringProportions(
            rowSums(dat[, Matrix::colSums(dat) <= threshold]))[,1]
        Top5000_real <- prob_realbg>=sort(prob_realbg,decreasing = T)[5000]
        KL_knee <- KL(prob_knee[Top5000_real],prob_realbg[Top5000_real])
        Size_knee <- mean(colSums(dat_knee))
        Knee_cells <- Size[,1:n]
        
        
        SEED.rep <- as.integer(runif(n, 1, 100000))
        
        for (sss in 1:n) {
            SEED.tmp <- SEED.rep[sss]
            
            cat(paste0("Run ", sss, "\n"))
            
            set.seed(SEED.tmp)
            
            ###2000 G1 large cells
            G1 <- dat_true1[, sample(ncol(dat_true1), n_large, replace = T)]
            colnames(G1) <- paste0("G1_", 1:ncol(G1))
            G1 <- as(G1, "dgCMatrix")
            
            
            ###2000 G2 middle cells (50% downsampling)
            G2 <- dat_true1[, sample(ncol(dat_true1), n_middle, replace = T)]
            G2 <- as(G2, "dgCMatrix")
            G2 <- downsampleMatrix(G2, 0.5)
            colnames(G2) <- paste0("G2_", 1:ncol(G2))
            
            
            ###2000 G3 small cells (10% downsampling)
            G3 <- dat_true1[, sample(ncol(dat_true1), n_small, replace = T)]
            G3 <- as(G3, "dgCMatrix")
            G3 <- downsampleMatrix(G3, 0.1)
            colnames(G3) <- paste0("G3_", 1:ncol(G3))
            
            if (new_sim) {
                G1 <- reorder_mat(G1, rate = reorder_rate)
                G2 <- reorder_mat(G2, rate = reorder_rate)
                G3 <- reorder_mat(G3, rate = reorder_rate)
            }
            
            Size[1,sss] <- mean(colSums(G1))
            Size[2,sss] <- mean(colSums(G2))
            Size[3,sss] <- mean(colSums(G3))
            
            # KL_Diver[1,sss] <- median(apply(G1[Top5000,],2,KL,bg=prob_bg[Top5000]))
            # KL_Diver[2,sss] <- median(apply(G2[Top5000,],2,KL,bg=prob_bg[Top5000]))
            # KL_Diver[3,sss] <- median(apply(G3[Top5000,],2,KL,bg=prob_bg[Top5000]))
            
            KL_Diver[1,sss] <- KL(edgeR::goodTuringProportions(rowSums(G1))[Top5000,1],prob_bg[Top5000])
            KL_Diver[2,sss] <- KL(edgeR::goodTuringProportions(rowSums(G2))[Top5000,1],prob_bg[Top5000])
            KL_Diver[3,sss] <- KL(edgeR::goodTuringProportions(rowSums(G3))[Top5000,1],prob_bg[Top5000])
            
            # Entropy[1,sss] <- median(apply(G1[Top5000,],2,Ent))
            # Entropy[2,sss] <- median(apply(G2[Top5000,],2,Ent))
            # Entropy[3,sss] <- median(apply(G3[Top5000,],2,Ent))
            
            Entropy[1,sss] <- Ent(edgeR::goodTuringProportions(rowSums(G1))[Top5000,1])
            Entropy[2,sss] <- Ent(edgeR::goodTuringProportions(rowSums(G2))[Top5000,1])
            Entropy[3,sss] <- Ent(edgeR::goodTuringProportions(rowSums(G3))[Top5000,1])
            
            dat_sim1 <- cbind(dat_bg, G1, G2, G3)
            if(run){
                
                set.seed(SEED.tmp)
                
                if(ED_knee){
                    retain_temp <- barcodeRanks(dat_sim1,threshold)@metadata$knee
                }else{
                    Retain <- Calc_retain(dat_sim1, threshold)
                    
                    #check convergence of knee point
                    repeat{
                        retain_temp <- Retain$knee
                        Retain <- Calc_retain(dat_sim1, Retain$inflection + 100)
                        if(Retain$knee==retain_temp) break
                    }
                }
                
                Knee_cell <- colnames(dat_sim1)[colSums(dat_sim1)>=retain_temp]
                Knee_cells[1,sss] <- sum(Knee_cell%in%colnames(G1))
                Knee_cells[2,sss] <- sum(Knee_cell%in%colnames(G2))
                Knee_cells[3,sss] <- sum(Knee_cell%in%colnames(G3))
                
                
                res_CB2 <-
                    CB2FindCell(dat_sim1,
                                lower = threshold,
                                FDR_threshold = FDR, upper = retain_temp,
                                ...)
                
                res1 <- colnames(GetCellMat(res_CB2))
                
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
                
                eOut <- emptyDrops(dat_sim1, lower = threshold,retain = retain_temp)
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
        }
        
        TP_CB2[, n + 1] <- rowMeans(TP_CB2[, 1:n])
        TP_ED[, n + 1] <- rowMeans(TP_ED[, 1:n])
        FP_CB2[, n + 1] <- rowMeans(FP_CB2[, 1:n])
        FP_ED[, n + 1] <- rowMeans(FP_ED[, 1:n])
        
        TP_CB2_adj <- TP_CB2
        colnames(TP_CB2_adj)[n+1] <- "Total"
        TP_ED_adj <- TP_CB2_adj
        
        NonKnee_cells <- c(n_large,n_middle,n_small)-Knee_cells
        TP_CB2_adj[c(1,3,5),1:n] <- TP_CB2[c(1,3,5),1:n]-Knee_cells
        TP_ED_adj[c(1,3,5),1:n] <- TP_ED[c(1,3,5),1:n]-Knee_cells
        TP_CB2_adj[c(2,4,6),1:n] <- TP_CB2_adj[c(1,3,5),1:n]/NonKnee_cells
        TP_ED_adj[c(2,4,6),1:n] <- TP_ED_adj[c(1,3,5),1:n]/NonKnee_cells
        TP_CB2_adj[c(1,3,5),n+1] <- rowSums(TP_CB2_adj[c(1,3,5),1:n])
        TP_ED_adj[c(1,3,5),n+1] <- rowSums(TP_ED_adj[c(1,3,5),1:n])
        TP_CB2_adj[c(2,4,6),n+1] <- TP_CB2_adj[c(1,3,5),n+1]/rowSums(NonKnee_cells)
        TP_ED_adj[c(2,4,6),n+1] <- TP_ED_adj[c(1,3,5),n+1]/rowSums(NonKnee_cells)
        
        Stat <- matrix(NA, 2, 2)
        colnames(Stat) <- c("CB2", "ED")
        rownames(Stat) <- c("Power", "FDR")
        Stat[1, 1] <-
            (TP_CB2[2, n + 1] * n_large + TP_CB2[4, n + 1] * n_middle + TP_CB2[6, n + 1] * n_small) / (n_large + n_middle + n_small)
        Stat[1, 2] <-
            (TP_ED[2, n + 1] * n_large + TP_ED[4, n + 1] * n_middle + TP_ED[6, n + 1] * n_small) / (n_large + n_middle + n_small)
        Stat[2, 1] <-
            FP_CB2[3, n + 1] / (FP_CB2[3, n + 1] + Stat[1, 1] * (n_large + n_middle + n_small))
        Stat[2, 2] <-
            FP_ED[3, n + 1] / (FP_ED[3, n + 1] + Stat[1, 2] * (n_large + n_middle + n_small))
        Power_diff <- Stat[1, 1] - Stat[1, 2]
        
        
        Stat_adj <- Stat
        Stat_adj[1,1] <- sum(TP_CB2_adj[c(1,3,5),n+1])/sum(NonKnee_cells)
        Stat_adj[1,2] <- sum(TP_ED_adj[c(1,3,5),n+1])/sum(NonKnee_cells)
        
        Power_diff_adj <- Stat_adj[1, 1] - Stat_adj[1, 2]
        
        Size[,n+1] <- rowMeans(Size[,1:n])
        KL_Diver[,n+1] <- rowMeans(KL_Diver[,1:n])
        Entropy[,n+1] <- rowMeans(Entropy[,1:n])
        
        
        return(
            list(
                TP_CB2 = TP_CB2,
                FP_CB2 = FP_CB2,
                TP_ED = TP_ED,
                FP_ED = FP_ED,
                Stat = Stat,
                Power_diff = Power_diff,
                TP_CB2_adj = TP_CB2_adj,
                TP_ED_adj = TP_ED_adj,
                Stat_adj = Stat_adj,
                Power_diff_adj = Power_diff_adj,
                NonKnee_cells = NonKnee_cells,
                Size=Size,
                KL_Diver=KL_Diver,
                Entropy=Entropy,
                Size_knee=Size_knee,
                KL_knee=KL_knee,
                session=sessionInfo()
            )
        )
    }


#Calculate knee point of a dataset
#Note: This function is a modified version of barcodeRanks() in 
#package DropletUtils. We fixed a minor bug causing returned knee point
#being larger than actual knee point. We also moved the smooth spline fitting
#to the beginning to avoid unstable knee point estimation when lower threshold
#changes. For its origin, see 
#https://github.com/MarioniLab/DropletUtils/blob/master/R/barcodeRanks.R
Calc_retain <- function(dat,lower){
    dat <- FilterGB(dat)
    totals <- unname(colSums(dat))
    o <- order(totals, decreasing=TRUE)
    
    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values
    
    
    keep <- run.totals > lower
    if (sum(keep)<3) { 
        warning("insufficient unique points for computing knee/inflection points")
        return(NULL)
    }
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])
    fit <- smooth.spline(log10(run.rank), log10(run.totals), df=20)
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
    return(list(knee=knee,inflection=inflection))
}


#####################
# For revision: Add simulation result for simple hard threshold
#####################

# n = 5
# dat
# dat_bg = NULL
# threshold = 100
# n_large = 2000
# n_middle = 2000
# n_small = 2000
# remove_protein = T



SIM_HARD <-
    function(n = 5,
             dat,
             dat_bg = NULL,
             threshold = 100,
             n_large = 2000,
             n_middle = 2000,
             n_small = 2000,
             remove_protein = T,
             inflection_prop=0.5,
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
        
        knee <-
            ifelse(is.null(brank$knee),
                   brank@metadata$knee,
                   brank$knee)
        
        HARD_threshold <- inflection * inflection_prop
        
        ###background barcodes
        if (is.null(dat_bg)) {
            dat_bg <- SIMFUN_bg(dat, threshold = threshold)
        }
        
        all_bg <- sum(Matrix::colSums(dat) %in% 1:inflection)
        th_bg <- sum(Matrix::colSums(dat) %in% (threshold + 1):inflection)
        
        dat_true1 <- dat[, Matrix::colSums(dat) > inflection]
        dat_knee <- dat[, Matrix::colSums(dat) >= knee]
        
        TP_HARD <- matrix(0, 6, n + 1)
        colnames(TP_HARD) <- c(paste0("rep", 1:n), "Ave")
        rownames(TP_HARD) <-
            c("G1", "G1 prop", "G2", "G2 prop", "G3", "G3 prop")
        
        FP_HARD <- TP_HARD[1:4, ]
        rownames(FP_HARD) <-
            c(
                paste0("all (", all_bg, ")"),
                "all prop",
                paste0(">", threshold, " (", th_bg, ")"),
                paste0(">", threshold, " prop")
            )
        
        
        Size <- matrix(NA,3,n+1)
        rownames(Size) <- c("G1","G2", "G3")
        colnames(Size) <- c(1:n, "Average")
        Knee_cells <- Size[,1:n]
        
        SEED.rep <- as.integer(runif(n, 1, 100000))
        
        for (sss in 1:n) {
            SEED.tmp <- SEED.rep[sss]
            
            cat(paste0("Run ", sss, "\n"))
            
            set.seed(SEED.tmp)
            
            ###2000 G1 large cells
            G1 <- dat_true1[, sample(ncol(dat_true1), n_large, replace = T)]
            colnames(G1) <- paste0("G1_", 1:ncol(G1))
            G1 <- as(G1, "dgCMatrix")
            
            
            ###2000 G2 middle cells (50% downsampling)
            G2 <- dat_true1[, sample(ncol(dat_true1), n_middle, replace = T)]
            G2 <- as(G2, "dgCMatrix")
            G2 <- downsampleMatrix(G2, 0.5)
            colnames(G2) <- paste0("G2_", 1:ncol(G2))
            
            
            ###2000 G3 small cells (10% downsampling)
            G3 <- dat_true1[, sample(ncol(dat_true1), n_small, replace = T)]
            G3 <- as(G3, "dgCMatrix")
            G3 <- downsampleMatrix(G3, 0.1)
            colnames(G3) <- paste0("G3_", 1:ncol(G3))
            
            Size[1,sss] <- mean(colSums(G1))
            Size[2,sss] <- mean(colSums(G2))
            Size[3,sss] <- mean(colSums(G3))
            
            dat_sim1 <- cbind(dat_bg, G1, G2, G3)
            
            if(ED_knee){
                retain_temp <- barcodeRanks(dat_sim1,threshold)@metadata$knee
            }else{
                Retain <- Calc_retain(dat_sim1, threshold)
                
                #check convergence of knee point
                repeat{
                    retain_temp <- Retain$knee
                    Retain <- Calc_retain(dat_sim1, Retain$inflection + 100)
                    if(Retain$knee==retain_temp) break
                }
            }
            
            Knee_cell <- colnames(dat_sim1)[colSums(dat_sim1)>=retain_temp]
            Knee_cells[1,sss] <- sum(Knee_cell%in%colnames(G1))
            Knee_cells[2,sss] <- sum(Knee_cell%in%colnames(G2))
            Knee_cells[3,sss] <- sum(Knee_cell%in%colnames(G3))
            
            
            res1 <- colnames(dat_sim1)[colSums(dat_sim1)>=HARD_threshold]
            
            TP_HARD[1, sss] <- length(grep(x = res1, pattern = "G1"))
            TP_HARD[2, sss] <- TP_HARD[1, sss] / n_large
            
            TP_HARD[3, sss] <- length(grep(x = res1, pattern = "G2"))
            TP_HARD[4, sss] <- TP_HARD[3, sss] / n_middle
            
            TP_HARD[5, sss] <- length(grep(x = res1, pattern = "G3"))
            TP_HARD[6, sss] <- TP_HARD[5, sss] / n_small
            
            FP_HARD[1, sss] <- length(grep(x = res1, pattern = "bg"))
            FP_HARD[2, sss] <- FP_HARD[1, sss] / all_bg
            
            FP_HARD[3, sss] <- length(grep(x = res1, pattern = "bg"))
            FP_HARD[4, sss] <- FP_HARD[3, sss] / th_bg
            
        }
        
        TP_HARD[, n + 1] <- rowMeans(TP_HARD[, 1:n])
        FP_HARD[, n + 1] <- rowMeans(FP_HARD[, 1:n])
        
        TP_HARD_adj <- TP_HARD
        colnames(TP_HARD_adj)[n+1] <- "Total"
        
        NonKnee_cells <- c(n_large,n_middle,n_small)-Knee_cells
        TP_HARD_adj[c(1,3,5),1:n] <- TP_HARD[c(1,3,5),1:n]-Knee_cells
        # TP_ED_adj[c(1,3,5),1:n] <- TP_ED[c(1,3,5),1:n]-Knee_cells
        TP_HARD_adj[c(2,4,6),1:n] <- TP_HARD_adj[c(1,3,5),1:n]/NonKnee_cells
        # TP_ED_adj[c(2,4,6),1:n] <- TP_ED_adj[c(1,3,5),1:n]/NonKnee_cells
        TP_HARD_adj[c(1,3,5),n+1] <- rowSums(TP_HARD_adj[c(1,3,5),1:n])
        # TP_ED_adj[c(1,3,5),n+1] <- rowSums(TP_ED_adj[c(1,3,5),1:n])
        TP_HARD_adj[c(2,4,6),n+1] <- TP_HARD_adj[c(1,3,5),n+1]/rowSums(NonKnee_cells)
        # TP_ED_adj[c(2,4,6),n+1] <- TP_ED_adj[c(1,3,5),n+1]/rowSums(NonKnee_cells)
        
        Stat <- matrix(NA, 2, 1)
        colnames(Stat) <- "HARD"
        rownames(Stat) <- c("Power", "FDR")
        Stat[1, 1] <-
            (TP_HARD[2, n + 1] * n_large + TP_HARD[4, n + 1] * n_middle + TP_HARD[6, n + 1] * n_small) / (n_large + n_middle + n_small)
        Stat[2, 1] <-
            FP_HARD[3, n + 1] / (FP_HARD[3, n + 1] + Stat[1, 1] * (n_large + n_middle + n_small))
        
        
        
        Stat_adj <- Stat
        Stat_adj[1,1] <- sum(TP_HARD_adj[c(1,3,5),n+1])/sum(NonKnee_cells)
        
        
        
        return(
            list(
                TP_HARD = TP_HARD,
                FP_HARD = FP_HARD,
                Stat = Stat,
                TP_HARD_adj = TP_HARD_adj,
                Stat_adj = Stat_adj,
                NonKnee_cells = NonKnee_cells
            )
        )
    }
oncotypedx2 <- function (data, annot, do.mapping = FALSE, mapping, verbose = FALSE) 
{
    sig2 <- sig.oncotypedx[sig.oncotypedx[, "group"] != "reference", 
        , drop = FALSE]
    dimnames(sig2)[[1]] <- sig2[, "symbol"]
    gt <- nrow(sig2)
    if (do.mapping) {
        gid1 <- as.numeric(as.character(sig2[, "EntrezGene.ID"]))
        names(gid1) <- dimnames(sig2)[[1]]
        gid2 <- as.numeric(as.character(annot[, "EntrezGene.ID"]))
        names(gid2) <- dimnames(annot)[[1]]
        rm.ix <- is.na(gid1) | duplicated(gid1)
        gid1 <- gid1[!rm.ix]
        rr <- geneid.map(geneid1 = gid2, data1 = data, geneid2 = gid1, 
            verbose = FALSE)
        gm <- length(rr$geneid2)
        mymapping <- c(mapped = gm, total = gt)
        if (length(rr$geneid1) != gt) {
            res <- rep(NA, nrow(data))
            names(res) <- dimnames(data)[[1]]
            warning(sprintf("Probe candidates: %i/%i", gm, gt), 
                "\nIncomplete overlap between the gene signature EntrezGene.IDs", 
                " and the EntrezGene.ID column of annot... Returning all NAs.")
            return(list(score = res, risk = res, mapping = mymapping, 
                probe = NA))
        }
        gid1 <- rr$geneid2
        gid2 <- rr$geneid1
        data <- rr$data1
        myprobe <- cbind(probe = names(gid1), EntrezGene.ID = gid1, 
            new.probe = names(gid2))
        dimnames(data)[[2]] <- names(gid2) <- names(gid1)
    }else {
        myprobe <- NA
        data <- data[, intersect(dimnames(sig2)[[1]], dimnames(data)[[2]])]
        #data <- data[, intersect(sig2$symbol, dimnames(data)[[2]])]
        gm <- ncol(data)
        mymapping <- c(mapped = gm, total = gt)
        if (nrow(sig2) != ncol(data)) {
            res <- rep(NA, nrow(data))
            names(res) <- dimnames(data)[[1]]
            warning(sprintf("Probe candidates: %i/%i", gm, gt), 
                "\nIncomplete overlap between the gene signature EntrezGene.IDs", 
                " and the colnames of data... Returning all NAs.")
            return(list(score = res, risk = res, mapping = mymapping, 
                probe = myprobe))
        }
    }
    dimnames(data)[[2]] <- dimnames(sig2)[[1]] <- sig2[, "symbol"]
    data <- apply(data, 2, function(x) {
        xx <- (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
            min(x, na.rm = TRUE))
        return(xx * 15)
    })
    cc.ix <- complete.cases(data)
    rs <- rs.unscaled <- rsrisk <- NULL
    for (i in 1:nrow(data)) {
        if (cc.ix[i]) {
            grb7.gs <- 0.9 * data[i, "GRB7"] + 0.1 * data[i, 
                "ERBB2"]
            if (grb7.gs < 8) {
                grb7.gs <- 8
            }
            er.gs <- (0.8 * data[i, "ESR1"] + 1.2 * data[i, "PGR"] + 
                data[i, "BCL2"] + data[i, "SCUBE2"])/4
            proliferation.gs <- (data[i, "BIRC5"] + data[i, "MKI67"] + 
                data[i, "MYBL2"] + data[i, "CCNB1"] + data[i, 
                "AURKA"])/5
            if (proliferation.gs < 6.5) {
                proliferation.gs <- 6.5
            }
            invasion.gs <- (data[i, "CTSL2"] + data[i, "MMP11"])/2
            rsu <- 0.47 * (grb7.gs) - 0.34 * (er.gs) + 1.04 * 
                (proliferation.gs) + 0.1 * (invasion.gs) + 0.05 * 
                data[i, "CD68"] - 0.08 * data[i, "GSTM1"] #- 0.07 * data[i, "BAG1"]
            rsu2 <- rsu
            if (rsu >= 0 & rsu <= 100) {
                rsu <- 20 * (rsu - 6.7)
            }
            if (rsu < 0) {
                rsu <- 0
            }
            if (rsu > 100) {
                rsu <- 100
            }
            if (rsu < 18) {
                rsr <- 0
            }
            if (rsu >= 18 & rsu < 31) {
                rsr <- 0.5
            }
            if (rsu >= 31) {
                rsr <- 1
            }
        }
        else {
            rsu <- rsr <- rsu2 <- NA
        }
        rs.unscaled <- c(rs.unscaled, rsu2)
        rs <- c(rs, rsu)
        rsrisk <- c(rsrisk, rsr)
    }
    names(rs) <- names(rs.unscaled) <- names(rsrisk) <- dimnames(data)[[1]]
    return(list(score = rs, risk = rsrisk, mapping = mymapping, 
        probe = myprobe))
}

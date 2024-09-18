#' impute expression
#'
#' @param Exp
#' @param G
#'
#' @return expression profile
#' @export
#'
#' @examples
naImpute <- function(Exp,G=NULL){
        m <- rowMeans(Exp,na.rm=T)
        mm <- mean(m,na.rm=T)
        m[which(is.na(m))] <- mm
        mna <- rowMeans(Exp)
        w <- which(is.na(mna))
        for(i in w){
            try(Exp[i,which(is.na(Exp[i,]))] <- m[i])
        }
        G <- setdiff(G,rownames(Exp))
        if(!is.null(G)) {
            tmp <- as.data.frame(array(mm,dim=c(length(G),ncol(Exp))))
            rownames(tmp) <- G
            names(tmp) <- names(Exp)
            Exp <- rbind(Exp,tmp)
        }
        Exp
}
#
#' change profile
#'
#' @param dataframe
#' @param signaturegene
#'
#' @return converted profile
#' @export
#'
#' @examples
GetRFExp <- function(dataframe, signaturegene){
	datasigexp <- as.data.frame(scale(t(dataframe[signaturegene, ])))
	datasigexp <- as.data.frame(t(naImpute(Exp=as.data.frame(t(datasigexp)),G=signaturegene)))
	comb.mat <- combn(seq_along(signaturegene), 2)   ##???е?����????
	train.pairwise.matrix <- apply(comb.mat, 2, function(pair)
	 datasigexp[,pair[1]] > datasigexp[,pair[2]])
	train.pairwise.vals <- as.data.frame(train.pairwise.matrix)
	return(train.pairwise.vals)
}
#
#' change symbol to entrez
#'
#' @param GeneCount
#'
#' @return gene entrez profile
#' @export
#'
#' @examples
ChangeEntrezFromName <- function(GeneCount){
	GeneNameb1 <- GeneCount[as.character(intersect(rownames(GeneCount), Ensembl2Name2NCBI$Gene.name[which(!is.na(Ensembl2Name2NCBI$NCBI.gene.ID))])), ]
	GeneNameb1Rowname <- as.character(Ensembl2Name2NCBI$NCBI.gene.ID[match(rownames(GeneNameb1), Ensembl2Name2NCBI$Gene.name)])
	Convert <- DeleteDupChoseMax(GeneNameb1, GeneNameb1Rowname)
	return(Convert)
}
#
RowDFtoRawNum <- function(ExpDataframe){
  for(i in 1:dim(ExpDataframe)[2]){
    tmp <- as.numeric(ExpDataframe[, i])
    if(i == 1){
      Matrix <- tmp
    } else {
      Matrix <- cbind(Matrix, tmp)
    }
  }
  colnames(Matrix) <- colnames(ExpDataframe)
  rownames(Matrix) <- rownames(ExpDataframe)
  df <- as.data.frame(Matrix)
  return(df)
}
DeleteDupChoseMax <- function(ExpDataframe, ReplaceName, by="row", ori=FALSE){
	if(by == "row"){
		ReplaceRowname <- ReplaceName
		dupchar <- unique(ReplaceRowname[duplicated(ReplaceRowname)])
		if(length(dupchar) > 0){
			expUnique <- ExpDataframe[!ReplaceRowname %in% dupchar, ]
			rownames(expUnique) <- ReplaceRowname[!ReplaceRowname %in% dupchar]
			expDup <- do.call(rbind, lapply(dupchar, function(x){
				Single <- ExpDataframe[which(ReplaceRowname == x), ]
				tmp <- Single[which.max(rowMeans(Single)), ]
				return(tmp)
			}))
			rownames(expDup) <- dupchar
			expDup <- RowDFtoRawNum(expDup)
			FinalDataframe <- rbind(expUnique, expDup)
		} else {
			FinalDataframe <- ExpDataframe
			rownames(FinalDataframe) <- ReplaceRowname
		}
	}
	if(by == "col"){
		ReplaceColname <- ReplaceName
		dupchar <- unique(ReplaceColname[duplicated(ReplaceColname)])
		expUnique <- ExpDataframe[, !ReplaceColname %in% dupchar]
		colnames(expUnique) <- ReplaceColname[!ReplaceColname %in% dupchar]
		expDup <- do.call(cbind, lapply(dupchar, function(x){
			Single <- ExpDataframe[, which(ReplaceColname == x)]
			tmp <- Single[, which.max(colMeans(Single))]
			return(tmp)
		}))
		colnames(expDup) <- dupchar
		FinalDataframe <- cbind(expUnique, expDup)
		if(ori == TRUE){
			uniquenames <- colnames(ExpDataframe)[!ReplaceColname %in% dupchar]
			remains <- do.call(c, lapply(dupchar, function(x){
				Single <- ExpDataframe[, which(ReplaceColname == x)]
				tmp <- colnames(Single)[which.max(colMeans(Single))]
				return(tmp)
			}))
			colnames(FinalDataframe) <- c(uniquenames, remains)
		}
	}
    return(FinalDataframe)
}
#
#' prediction
#'
#' @param model
#' @param texp
#' @param phenotype
#' @param minPosterior
#'
#' @return probability
#' @export
#'
#' @examples
PredictedInMultiRF <- function(model, texp, phenotype=NULL, minPosterior){
	texp <- as.data.frame(texp)
	colnames(texp) <- gsub("-", "_", colnames(texp))
	p1 <- predict(model, type="prob", newdata = texp)
	p <- predict(model, type="prob", newdata = texp)
	p <- as.data.frame(p)
	p$nearest <- apply(p, 1, function(v) paste0(colnames(p)[which(v == max(v))], collapse=";"))
	p$predicted <- p$nearest
	p$predicted[apply(p, 1, function(x) max(x[1:dim(p1)[2]])) < minPosterior] <- NA
	if(length(phenotype) > 0){
		p$phenotype <- as.character(phenotype)
	}
	return(p)
}
#
#' main function
#'
#' @param Expr
#' @param minPosterior
#' @param scale
#' @param log2transfrom
#'
#' @return probability
#' @export
#'
#' @examples
CMSFFPEClassifier <- function(Expr = NULL, minPosterior = 0.5, log2transfrom=FALSE){
	### check name
	message('The row name should be gene symbol or gene Entrez ID ......')
	if(length(intersect(Ensembl2Name2NCBI$Gene.name, rownames(Expr))) > 10){
		Expr <- ChangeEntrezFromName(GeneCount=Expr)
	}
	if(log2transfrom){
		Expr_log2 <- log2(Expr)
	} else {
		Expr_log2 <- Expr
	}
	if(length(intersect(signatures, rownames(Expr_log2))) < length(signatures)){
		stop('lacking of expression of signature gene ......')
	} else {
		Expr_impute <- Expr_log2
	}
	Expr_val <- GetRFExp(dataframe=Expr_impute, signaturegene=signatures)
	prob <- PredictedInMultiRF(model=finalModel, texp=Expr_val, minPosterior=minPosterior)
	rownames(prob) <- colnames(Expr)
	return(prob)
}

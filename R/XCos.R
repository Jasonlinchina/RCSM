#' XCos
#'
#' The implementation of XCos
#'
#' @details `XCosScore()` returns a data.frame, each row of which contains
#'   score, pValue and adjusted-pValue for one sample in the refMatrix.
#' @references "Cheng J et al. Genome medicine, 2014, 6(12): 95".
#' @param refMatrix A matrix
#' @param query A numeric vector with names.
#' @param permuteNum number of perturbation time for computing pValue
#' @param topN number of genes in top and bottom of reference gene list
#' @param pAdjMethod method to use for computing adjudted-pValue
#' @param mcCore the number of core to use for parallel computing
#' @importFrom parallel mclapply
#' @export
#' @examples
#' set.seed(1234)
#' ref <- matrix(rnorm(1000), nrow = 10,
#'   dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
#' query <- rnorm(4)
#' names(query) <- c("gene1", "gene2", "gene9", "gene10")
#' XCosScore(refMatrix = ref, query = query, topN = 4)

#######################The implementation of XCos#########################
XCosScore <- function(refMatrix, query, topN = 500, permuteNum = 10000,
                            pAdjMethod = "BH", mcCore = 1) {

  if (is.data.frame(refMatrix)) {refMatrix <- as.matrix(refMatrix)}

  if (is.null(colnames(refMatrix)) || is.null(rownames(refMatrix))) {
    stop("Warning: refMatrix must have both rownames and colnames!")
  }

  if (!is.numeric(query)) {stop("Warning: query must be a numeric vector!")}

  if (is.null(names(query))) {stop("Warning: query must have names!")}

  if (topN > nrow(refMatrix) / 2) {stop("Warning: topN is lager than half
                                        the length of gene list!")}

  ## Convert the gene expression matrix to ranked list. topN is the number of
  ## genes in the ranked top, which can be reseted.
  matrixToRankedList <- function(refMatrix, topN) {
    ## Allocate memory for the refList
    refList <- vector("list", ncol(refMatrix))
    for(i in 1:ncol(refMatrix)) {
      ## Sort the reference gene lists based on fold change of gene expression.
      ##  And get head topN genes and tail topN genes.
      refList[[i]] <- c(head(refMatrix[order(refMatrix[, i], decreasing=TRUE), i],
                             n = topN),
                        tail(refMatrix[order(refMatrix[, i], decreasing=TRUE), i],
                             n = topN))
    }
    return(refList)
  }
  ## The core part for computing the XCos score
  XCos <- function(refList, query) {
    reservedRef <- refList[match(intersect(names(refList), names(query)),
                                 names(refList))]
    reservedRef[order(names(reservedRef))]
    reservedQuery <- query[match(intersect(names(refList), names(query)),
                                 names(query))]
    reservedQuery[order(names(reservedQuery))]
    ## Compute the cosine similarity
    if (length(reservedRef) == 0) {
      return(NA)
    }
    else {
    return((crossprod(reservedRef, reservedQuery) /
             sqrt(crossprod(reservedRef) * crossprod(reservedQuery)))[1, 1])
    }
  }

  ## Prepare the ranked reference lists
  refList <- matrixToRankedList(refMatrix, topN = topN)
  ## Compute the scores for each sample in the reference lists. mcCore is the
  ## number of cores to use for parallel computing. Set it based on your computer.
  score <- mclapply(refList, XCos, query = query, mc.cores = mcCore)
  score <- as.vector(do.call(rbind, score))
  ## Allocate memory for the permuteScore that are used to compute the p-value.
  ## The permuteNum can be reseted. Notice large permuteNum means low speed.
  permuteScore <- matrix(0, ncol = permuteNum, nrow = ncol(refMatrix))
  for(n in 1:permuteNum) {
    ## Prepare the random query signatures
    names(query) <- sample(rownames(refMatrix), size = length(query))
    ## Compute the random scores for each sample in the reference lists
    bootScore <- mclapply(refList, XCos, query = query, mc.cores = mcCore)
    permuteScore[, n] <- as.vector(do.call(rbind, bootScore))
  }
  permuteScore[is.na(permuteScore)] <- 0
  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(score)) / permuteNum
  ## Compute the adjusted p-value. The adjusting method can be reseted
  ## (Refer to p.adjust()).
  pAdjust <- p.adjust(pValue, method = pAdjMethod)

  scoreResult <- data.frame(Score = score, pValue = pValue, pAdjValue = pAdjust)
  rownames(scoreResult) <- colnames(refMatrix)
  return(scoreResult)
}
################################################################################

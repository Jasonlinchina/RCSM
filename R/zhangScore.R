#' zhangScore
#'
#' The implementation of zhangScore
#'
#' @details `zhangScore()` returns a data.frame, each row of which contains score,
#' pValue and adjusted-pValue for one sample in the refMatrix.
#' @references "Zhang S D et al. BMC bioinformatics, 2008, 9(1): 258".
#' @param refMatrix A matrix
#' @param queryUp,queryDown character vectors. if there is no up or down signatures,
#' you should set: queryUp = NULL or queryDown = NULL.
#' @param permuteNum number of perturbation time for computing pValue
#' @param pAdjMethod method to use for computing adjudted-pValue
#' @param mcCore the number of core to use for parallel computing
#' @importFrom parallel mclapply
#' @export
#' @examples
#' set.seed(1234)
#' ref <- matrix(rnorm(1000), nrow = 10,
#'   dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
#' Up <- c("gene1", "gene2")
#' Down <- c("gene9", "gene10")
#' zhangScore(refMatrix = ref, queryUp = Up, queryDown = Down)

#######################The implementation of zhangscore#########################
zhangScore <- function(refMatrix, queryUp, queryDown,
                                   permuteNum = 10000, pAdjMethod = "BH",
                                   mcCore = 1) {

  if (is.data.frame(refMatrix)) {refMatrix <- as.matrix(refMatrix)}

  if (is.null(colnames(refMatrix)) || is.null(rownames(refMatrix))) {
    stop("Warning: refMatrix should have both rownames and colnames!")
  }

  if (is.null(queryUp)) {queryUp <- character(0)}
  if (is.null(queryDown)) {queryDown <- character(0)}

  if (!is.character(queryUp)) {queryUp <- as.character(queryUp)}
  if (!is.character(queryDown)) {queryUp <- as.character(queryDown)}

  ## Convert the gene expression matrix to ranked list
  matrixToRankedList <- function(refMatrix) {
    ## Allocate memory for the refList
    refList <- vector("list", ncol(refMatrix))
    for(i in 1:ncol(refMatrix)) {
      ## Sort the reference gene lists based on the absolute value
      refSort <- refMatrix[order(abs(refMatrix[, i]), decreasing=TRUE), i]
      ## Make ranks for the reference gene lists with sign
      ## if two variables are equal, their ranks will be averaged
      refRank <- rank(abs(refSort)) * sign(refSort)
      refList[[i]] <- refRank
    }
    return(refList)
  }

  ## The core part for computing the zhangscore
  computeScore <- function(refRank, queryRank) {
    if(length(intersect(names(refRank), names(queryRank))) > 0) {
      ## Compute the maximal theoretical score
      maxTheoreticalScore <- sum(abs(refRank)[1:length(queryRank)] *
                                   abs(queryRank))
      ## The final score
      score <- sum(queryRank * refRank[names(queryRank)], na.rm=TRUE)/
        maxTheoreticalScore
    }
    else {
      score <- NA
    }
    return(score)
  }

  ## Prepare the ranked reference lists
  refList <- matrixToRankedList(refMatrix)
  ## Prepare the query signatures
  queryVector <- c(rep(1, length(queryUp)), rep(-1, length(queryDown)))
  names(queryVector) <- c(queryUp, queryDown)
  ## Compute the scores for each sample in the reference lists. mcCore is the
  ## number of cores to use for parallel computing. Set it based on your computer.
  score <- mclapply(refList, computeScore, queryRank = queryVector,
                    mc.cores = mcCore)
  score <- as.vector(do.call(rbind, score))
  ## Allocate memory for the permuteScore that are used to compute the p-value.
  ## The permuteNum can be reseted. Notice large permuteNum means low speed.
  permuteScore <- matrix(0, ncol = permuteNum, nrow = ncol(refMatrix))
  for(n in 1:permuteNum) {
    ## Prepare the random query signatures
    bootSample <- sample(c(-1,1), replace = TRUE,
                         size = length(queryUp) + length(queryDown))
    names(bootSample) <- sample(rownames(refMatrix), replace = FALSE,
                                size = length(queryUp) + length(queryDown))
    ## Compute the random scores for each sample in the reference lists
    bootScore <- mclapply(refList, computeScore, queryRank = bootSample,
                          mc.cores = mcCore)
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


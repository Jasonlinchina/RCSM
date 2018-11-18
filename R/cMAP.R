#' cmapScore
#'
#' The implementation of cmapScore
#'
#' @details `cmapScore()` returns a data.frame, each row of which contains score,
#'   pValue and adjusted-pValue for one sample in the refMatrix.
#' @references "Lamb J et al. science, 2006, 313(5795): 1929-1935".
#' @param refMatrix A matrix
#' @param queryUp,queryDown character vectors.
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
#' cmapScore(refMatrix = ref, queryUp = Up, queryDown = Down)

##########################The implementation of cMAP############################
cmapScore <- function(refMatrix, queryUp, queryDown, permuteNum = 10000,
                      pAdjMethod = "BH", mcCore = 1) {

  if (is.data.frame(refMatrix)) {refMatrix <- as.matrix(refMatrix)}

  if (is.null(colnames(refMatrix)) || is.null(rownames(refMatrix))) {
    stop("Warning: refMatrix should have both rownames and colnames!")
  }

  if (!is.character(queryUp)) {queryUp <- as.character(queryUp)}
  if (!is.character(queryDown)) {queryUp <- as.character(queryDown)}

  ## Convert the gene expression matrix to ranked list
  matrixToRankedList <- function(refMatrix) {
    ## Allocate memory for the refList
    refList <- vector("list", ncol(refMatrix))
    for(i in 1:ncol(refMatrix)) {
      ## Sort the reference gene lists based on fold change of gene expression
      refList[[i]] <- names(refMatrix[order(refMatrix[, i], decreasing=TRUE), i])
    }
    return(refList)
  }
  ## Compute the ks value
  ksScore <- function(refList, query) {
    lenRef <- length(refList)
    ## Get the rank of query signatures based on refList
    queryRank <- match(query, refList)
    queryRank <- sort(queryRank[!is.na(queryRank)])
    lenQuery <- length(queryRank)
    ## The core part for computing ks value
    if (lenQuery == 0) {
      return(0)
    } else {
      d <- (1:lenQuery) / lenQuery - queryRank / lenRef
      a <- max(d)
      b <- -min(d) + 1 / lenQuery
      ifelse(a > b, a, -b)
    }
  }
  ## Compute the ks score based on up ks value and down ks value
  ks <- function(refList, queryUp, queryDown) {
    scoreUp <- ksScore(refList, queryUp)
    scoreDown <- ksScore(refList, queryDown)
    ifelse(scoreUp * scoreDown <= 0, scoreUp - scoreDown, 0)
  }

  ## Prepare the ranked reference lists
  refList <- matrixToRankedList(refMatrix)
  ## Prepare the up and down signatures
  queryUp <- intersect(queryUp, rownames(refMatrix))
  queryDown <- intersect(queryDown, rownames(refMatrix))
  ## Compute the scores for each sample in the reference lists. mcCore is the
  ## number of cores to use for parallel computing. Set it based on your computer.
  score <- mclapply(refList, ks, queryUp = queryUp, queryDown = queryDown,
                    mc.cores = mcCore)
  score <- as.vector(do.call(rbind, score))
  ## Allocate memory for the permuteScore that are used to compute the p-value.
  ## The permuteNum can be reseted. Notice large permuteNum means low speed.
  permuteScore <- matrix(0, ncol = permuteNum, nrow = ncol(refMatrix))
  for(n in 1:permuteNum) {
    ## Prepare the random query signatures
    bootUp <- sample(rownames(refMatrix), size = length(queryUp))
    bootDown <- sample(rownames(refMatrix), size = length(queryDown))
    ## Compute the random scores for each sample in the reference lists
    bootScore <- mclapply(refList, ks, queryUp = bootUp, queryDown = bootDown,
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

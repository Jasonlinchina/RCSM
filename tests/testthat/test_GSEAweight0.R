# Test for functions in GSEAweight0.R
set.seed(1234)

# Test for the outputs of GSEAweight0Score
context("GSEAweight0Score")

ref <- matrix(rnorm(1000), nrow = 10,
              dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
Up <- c("gene1", "gene2")
Down <- c("gene9", "gene10")
results <- GSEAweight0Score(refMatrix = ref, queryUp = Up, queryDown = Down)

test_that("Whether GSEAweight0 works well", {
  expect_is(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  expect_equal(results$pAdjValue[48], 0.7270968, tolerance=1e-4)
  expect_equal(results$Score[36], -1.125, tolerance=1e-4)
  expect_equal(results$pValue[99], 0.137, tolerance=1e-4)
})

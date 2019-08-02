# Test for functions in XSum.R
set.seed(1234)

# Test for the outputs of XSumScore
context("XSumScore")

ref <- matrix(rnorm(1000), nrow = 10,
              dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
Up <- c("gene1", "gene2")
Down <- c("gene9", "gene10")
results <- XSumScore(refMatrix = ref, queryUp = Up, queryDown = Down, topN = 4)

test_that("Whether XSum works well", {
  expect_is(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  expect_equal(results$pAdjValue[48], 0.7073684, tolerance=1e-4)
  expect_equal(results$Score[36], -3.702277, tolerance=1e-4)
  expect_equal(results$pValue[99], 0.072, tolerance=1e-4)
})

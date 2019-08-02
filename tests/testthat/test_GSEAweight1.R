# Test for functions in GSEAweight1.R
set.seed(1234)

# Test for the outputs of GSEAweight1Score
context("GSEAweight1Score")

ref <- matrix(rnorm(1000), nrow = 10,
              dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
Up <- c("gene1", "gene2")
Down <- c("gene9", "gene10")
results <- GSEAweight1Score(refMatrix = ref, queryUp = Up, queryDown = Down)

test_that("Whether GSEAweight1 works well", {
  expect_is(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  expect_equal(results$pAdjValue[48], 0.8479592, tolerance=1e-4)
  expect_equal(results$Score[36], -1.74999, tolerance=1e-4)
  expect_equal(results$pValue[99], 1, tolerance=1e-4)
})

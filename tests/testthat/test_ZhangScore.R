# Test for functions in ZhangScore.R
set.seed(1234)

# Test for the outputs of ZhangScore
context("ZhangScore")

ref <- matrix(rnorm(1000), nrow = 10,
              dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
Up <- c("gene1", "gene2")
Down <- c("gene9", "gene10")
results <- ZhangScore(refMatrix = ref, queryUp = Up, queryDown = Down)

test_that("Whether ZhangScore works well", {
  expect_is(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  expect_equal(results$pAdjValue[48], 0.9197701, tolerance=1e-4)
  expect_equal(results$Score[36], -0.4117647, tolerance=1e-4)
  expect_equal(results$pValue[99], 0.3715, tolerance=1e-4)
})

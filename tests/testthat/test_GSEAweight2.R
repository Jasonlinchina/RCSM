# Test for functions in GSEAweight2.R
set.seed(1234)

# Test for the outputs of GSEAweight2Score
context("GSEAweight2Score")

ref <- matrix(rnorm(1000), nrow = 10,
              dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
Up <- c("gene1", "gene2")
Down <- c("gene9", "gene10")
results <- GSEAweight2Score(refMatrix = ref, queryUp = Up, queryDown = Down)

test_that("Whether GSEAweight2 works well", {
  expect_is(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  expect_equal(results$pAdjValue[48], 0.878, tolerance=1e-4)
  expect_equal(results$Score[36], -1.940979, tolerance=1e-4)
  expect_equal(results$pValue[99], 1,, tolerance=1e-4)
})

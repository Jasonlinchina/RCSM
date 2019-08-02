# Test for functions in XCos.R
set.seed(1234)

# Test for the outputs of XCosScore
context("XCosScore")

ref <- matrix(rnorm(1000), nrow = 10,
              dimnames = list(paste0("gene", 1:10), paste0("drug", 1:100)))
query <- rnorm(4)
names(query) <- c("gene1", "gene2", "gene9", "gene10")
results <- XCosScore(refMatrix = ref, query = query, topN = 4)

test_that("Whether XCos works well", {
  expect_is(results, "data.frame")
  expect_equal(nrow(results), 100)
  expect_equal(ncol(results), 3)
  expect_equal(results$pAdjValue[48], 0.8488, tolerance=1e-4)
  expect_equal(results$Score[36], 0.9774505, tolerance=1e-4)
  expect_equal(results$pValue[99], 0.2318, tolerance=1e-4)
})

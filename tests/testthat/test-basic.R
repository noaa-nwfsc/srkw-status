
test_that("test data", {
  skip_on_cran()
  data(orca)
  expect_gt(nrow(orca),200)
  
  required_columns <- c("animal", "pod", "sexF1M2", "birth", "death", "matriline", 
                        "mom", "includeFec", "includeSurv", "population")
  expect_true(all(required_columns %in% names(orca)))
})

test_that("test expanded", {
  skip_on_cran()
  data(orca)
  expanded <- expand(orca)
  expect_gt(nrow(expanded),10800)
  
  required_columns <- c("age", "alive", "gave_birth")
  expect_true(all(required_columns %in% names(expanded)))
})

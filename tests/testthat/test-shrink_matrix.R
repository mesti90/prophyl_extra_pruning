test_that("shrink_matrix() works", {
  mbig <- matrix(1:16, nrow = 4)
  row.names(mbig) <- LETTERS[1:4]
  colnames(mbig) <- LETTERS[1:4]
  msmall <- matrix(rep(1, times = 9), nrow = 3)
  row.names(msmall) <- LETTERS[1:3]
  colnames(msmall) <- LETTERS[1:3]
  shrink <- shrink_matrix(mbig, msmall)
  expect_true(all(row.names(shrink) == row.names(msmall)))
  expect_true(all(colnames(shrink) == colnames(msmall)))
})

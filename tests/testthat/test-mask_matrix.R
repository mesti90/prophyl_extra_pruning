test_that("mask_matrix", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "croatia", "china", "serbia", "thailand")
  )
  # works for a single focus group
  m1 <- mask_matrix(df, focus_by = "continent", focus_on = "europe")
  # some values are NA
  expect_true(all(is.na(m1[c(3,5),c(3,5)])))
  expect_equal(sum(is.na(m1)), 4)
  # the rest are 1
  expect_equal(unname(table(m1)[which(names(table(m1)) == 1)]), 21)
  # works for multiple focus groups
  m2 <- mask_matrix(df, focus_by = "country", focus_on = c("hungary","croatia"))
  # some values are NA
  expect_true(all(is.na(m2[3:5,3:5])))
  expect_equal(sum(is.na(m2)), 9)
  # the rest are 1
  expect_equal(unname(table(m2)[which(names(table(m2)) == 1)]), 16)
  # when all samples are in focus, returns an identity matrix and a warning
  m3 <- suppressWarnings(
    mask_matrix(df, focus_by = "continent", focus_on = c("asia", "europe"))
  )
  # all values are 1
  expect_equal(unname(table(m3)[which(names(table(m3)) == 1)]), 25)
  # there is a warning message
  expect_warning(
    mask_matrix(df, focus_by = "continent", focus_on = c("asia", "europe"))
  )
  # when no samples are in focus, returns a matrix of NAs and a warning
  m4 <- suppressWarnings(
    mask_matrix(df, focus_by = "continent", focus_on = "africa")
  )
  expect_equal(sum(is.na(m4)), 25)
  expect_warning(
    mask_matrix(df, focus_by = "continent", focus_on = "africa")
  )
})

# - TODO function works with missing data in focus by

test_that("varid_matrix() works", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "hungary", "china", "serbia", "china"),
    city = c("budapest", "szeged", "beijing", "belgrade", "beijing")
  )
  # no focus
  m1 <- varid_matrix(df, id_var = "assembly", var = "country")
  expect_equal(class(m1[1,1]), "numeric")
  expect_equal(sum(m1, na.rm = TRUE), 4)
  expect_equal(m1[1,2], 1)
  expect_equal(m1[3,5], 1)
  # focus
  m2 <- varid_matrix(
    df, 
    id_var = "assembly", 
    var = "country", 
    focus_by = "continent", 
    focus_on = "europe"
  )
  expect_equal(sum(m2, na.rm = TRUE), 2)
  expect_equal(sum(is.na(m2)), 7)
  m3 <- varid_matrix(
    df, 
    id_var = "assembly", 
    var = "country", 
    focus_by = "city", 
    focus_on = c("budapest", "szeged")
  )
  expect_equal(sum(is.na(m3)), 11)
})

test_that("varid_matrix() works with missing data", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", NA, "china", "serbia", "china"),
    city = c("budapest", "szeged", "beijing", NA, "beijing")
  )
  m1 <- varid_matrix(df, id_var = "assembly", var = "country")
  expect_equal(sum(is.na(m1)), 13)
  expect_true(all(is.na(m1[2, ])))
  expect_true(all(is.na(m1[,2])))
  m2 <- varid_matrix(df, id_var = "assembly", var = "city")
  expect_equal(sum(is.na(m2)), 13)
  expect_true(all(is.na(m2[4, ])))
  expect_true(all(is.na(m2[,4])))
})

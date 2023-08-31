# - TODO function works with missing data in focus by

test_that("neighbors_matrix standard borders",{
  df <- data.frame(
    assembly = LETTERS[1:6],
    continent = c(
      "europe", "europe", "asia", "europe", "asia", "europe"),
    country = c(
      "hungary", "germany", "china", "serbia", "laos", "albania"),
    country_iso2c = c("HU", "DE", "CN", "RS", "LA", "AL")
  )
  # returns correct results with no focus
  m1 <- neighbors_matrix(df, id_var = "assembly")
  expect_equal(class(m1[1,1]), "numeric")
  expect_equal(sum(m1, na.rm = TRUE), 6)
  expect_equal(m1[1,4], 1)
  expect_equal(m1[3,5], 1)
  expect_equal(m1[4,6], 1)
  # returns correct results with focus
  m2 <- neighbors_matrix(
    df, id_var = "assembly",
    focus_by = "continent",
    focus_on = "europe"
  )
  expect_equal(sum(is.na(m2)), 8)
})

test_that("neighbors_matrix custom borders",{
  df <- data.frame(
    assembly = LETTERS[1:7],
    continent = c(
      "europe", "europe", "asia", "europe", "asia", "europe", "europe"),
    country = c(
      "hungary", "germany", "china", "serbia", "laos", "albania", "kosovo"),
    country_iso2c = c("HU", "DE", "CN", "RS", "LA", "AL", "XK")
  )
  # returns NA for comparison with Kosovo and standard borders
  expect_warning(neighbors_matrix(df, id_var = "assembly"))
  m1 <- suppressWarnings(neighbors_matrix(df, id_var = "assembly"))
  expect_equal(sum(is.na(m1)),  19)
  # returns correct results when using custom borders
  data("custom_country_borders")
  updated_borders <- edit_borders(custom_country_borders)
  m2 <- neighbors_matrix(
    df, id_var = "assembly", country_borders = updated_borders)
  expect_equal(sum(m2, na.rm = TRUE), 8)
  expect_equal(m2[1,4], 1)
  expect_equal(m2[3,5], 1)
  expect_equal(m2[4,7], 1)
  expect_equal(m2[6,7], 1)
})

test_that("neighbors_matrix with NA in country codes", {
  df <- data.frame(
    assembly = LETTERS[1:6],
    continent = c(
      "europe", "europe", "asia", "europe", "asia", "europe"),
    country = c(
      "hungary", "germany", "china", "pirezia", "laos", "albania"),
    country_iso2c = c("HU", "DE", "CN", NA, "LA", "AL")
  )
  m1 <- neighbors_matrix(df, id_var = "assembly")
  expect_equal(sum(is.na(m1)), 16)
})

# - TODO function works with missing data in focus by

test_that("geodist_matrix", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "germany", "china", "hungary", "laos"),
    lat = c(47.16249, 51.16569, 35.86166, 47.16249, 19.85627),
    lon = c(19.50330, 10.451526, 104.19540, 19.50330, 102.4955)
  )
  # returns correct results when there is no focus
  m1 <- geodist_matrix(df, id_var = "assembly")
  expect_equal(class(m1[1,1]), "numeric")
  expect_equal(unname(m1[1,]), c(NA,   795,  6821, 0,    7893))
  expect_equal(unname(m1[2,]), c(795,  NA,   7232, 795,  8450))
  expect_equal(unname(m1[3,]), c(6821, 7232, NA,   6821, 1789))
  expect_equal(unname(m1[4,]), c(0,    795,  6821, NA,   7893))
  expect_equal(unname(m1[5,]), c(7893, 8450, 1789, 7893, NA))
  # some results are masked when there is focus
  m2 <- geodist_matrix(
    df, id_var = "assembly", focus_by = "continent", focus_on = "europe")
  expect_equal(sum(is.na(m2)), 7)
  # stops when there are duplicate rows
  df <- dplyr::bind_rows(
    df,
    data.frame(
      assembly = "A",
      continent = "europe",
      country = "hungary",
      lat = 47.16249,
      lon = 19.50330
    )
  )
  expect_error(
    geodist_matrix(
      df,
      id_var = "assembly"
    )
  )
})

test_that("geodist_matrix with NA in coordinates", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "germany", "china", "hungary", "laos"),
    lat = c(NA, 51.16569, 35.86166, 47.16249, 19.85627),
    lon = c(19.50330, 10.451526, 104.19540, 19.50330, 102.4955)
  )
  m1 <- geodist_matrix(
    df, 
    id_var = "assembly"
  )
  expect_equal(sum(is.na(m1)), 13)
})

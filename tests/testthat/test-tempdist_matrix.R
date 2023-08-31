# - TODO function works with missing data in focus by

test_that("tempdist_matrix", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "germany", "china", "hungary", "laos"),
    collection_date = c(
      "2000-01-01", "2001-01-01", "2005-01-01", "2000-06-15", "2005-06-15" 
    )
  )
  # returns correct results when there is no focus
  m1 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date"
  )
  expect_equal(class(m1[1,1]), "numeric")
  expect_equal(unname(m1[1,]), c(NA,   1,    5,    0.45, 5.45))
  expect_equal(unname(m1[2,]), c(1,    NA,   4,    0.55, 4.45))
  expect_equal(unname(m1[3,]), c(5,    4,    NA,   4.55, 0.45))
  expect_equal(unname(m1[4,]), c(0.45, 0.55, 4.55, NA,   5))
  expect_equal(unname(m1[5,]), c(5.45, 4.45, 0.45, 5,    NA))
  # some results are masked when there is focus
  m2 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date",
    focus_by = "continent",
    focus_on = "europe"
  )
  expect_equal(sum(is.na(m2)), 7)
  expect_true(is.na(m2[3,5]))
  expect_true(is.na(m2[5,3]))
  # stops when there are duplicate rows
  df <- dplyr::bind_rows(
    df,
    data.frame(
      assembly = "A",
      continent = "europe",
      country = "hungary",
      collection_date = "2000-01-01"
    )
  )
  expect_error(
    tempdist_matrix(
      df,
      id_var = "assembly",
      date_var = "collection_date"
    )
  )
})

test_that("tempdist_matrix with uncertain date", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "germany", "china", "hungary", "laos"),
    collection_date = c(
      "2000-01-01", "2001-01-01", "2005-01-01", "2000-06-15",  "2005-06"
    )
  )
  m0 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date"
  )
  expect_equal(sum(is.na(m0)), 13)
  m1 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date",
    estimate_dates = "lower"
  )
  expect_equal(sum(is.na(m1)), 5)
  m2 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date",
    estimate_dates = "upper"
  )
  expect_equal(sum(is.na(m2)), 5)
  m3 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date",
    estimate_dates = "middle"
  )
  expect_equal(sum(is.na(m3)), 5)
  m4 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date",
    estimate_dates = "runif"
  )
  expect_equal(sum(is.na(m4)), 5)
})

test_that("tempdist_matrix with NA in date_var", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "germany", "china", "hungary", "laos"),
    collection_date = c(
      "2000-01-01", "2001-01-01", "2005-01-01", "2000-06-15", NA 
    )
  )
  m1 <- tempdist_matrix(
    df, 
    id_var = "assembly",
    date_var = "collection_date"
  )
  expect_equal(sum(is.na(m1)), 13)
})
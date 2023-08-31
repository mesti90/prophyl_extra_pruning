test_that("phylodist_matrix() works", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "hungary", "china", "serbia", "china"),
    mlst = c("ST1", "ST2", "ST1", "ST2", "ST3")
  )
  set.seed(0)
  tr <- ape::rtree(5, tip.label = df$assembly)
  # works with no focus
  m1 <- phylodist_matrix(tr, df, id_var = "assembly")
  expect_equal(class(m1[1,1]), "numeric")
  # rows and columns are not mixed up
  expect_true(all(row.names(m1) == df$assembly))
  expect_true(all(colnames(m1) == df$assembly))
  # works with focus
  m2 <- phylodist_matrix(
    tr,
    df,
    id_var = "assembly",
    focus_by = "continent",
    focus_on = "europe"
  )
  expect_equal(sum(is.na(m2)), 7)
})

test_that("phylodist_matrix() works with missing data in focus by", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", NA, "asia"),
    country = c("hungary", "hungary", "china", "serbia", "china"),
    mlst = c("ST1", "ST2", "ST1", "ST2", "ST3")
  )
  set.seed(0)
  tr <- ape::rtree(5, tip.label = df$assembly)
  m1 <- phylodist_matrix(
    tr,
    df,
    id_var = "assembly",
    focus_by = "continent",
    focus_on = "europe"
  )
  expect_equal(sum(is.na(m1)), 11)
})
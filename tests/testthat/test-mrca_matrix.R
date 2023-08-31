test_that("mrca_matrix() works", {
  df <- data.frame(
    assembly = LETTERS[1:5],
    continent = c("europe", "europe", "asia", "europe", "asia"),
    country = c("hungary", "germany", "china", "hungary", "laos"),
    collection_date = c(
      "2000-01-01", "2001-01-01", "2005-01-01", "2000-06-15", "2005-06-15" 
    )
  )
  colldist <- tempdist_matrix(
    df, id_var = "assembly", date_var = "collection_date")
  set.seed(2)
  tr <- ape::rtree(5, tip.label = df$assembly, br = function(x) runif(5, 5, 15))
  phylodist <- phylodist_matrix(tree = tr, df = df,id_var = "assembly")
  m1 <- mrca_matrix(phylodist, colldist)
  expect_equal(class(m1[1,1]), "numeric")
  expect_equal(m1[1,1], (phylodist[1,1] - colldist[1,1])/2)
  # negative values can be transformed to 0.
  set.seed(0)
  tr <- ape::rtree(5, tip.label = df$assembly, br = function(x) runif(5, 5, 15))
  phylodist <- phylodist_matrix(tree = tr, df = df,id_var = "assembly")
  expect_warning(
    mrca_matrix(phylodist, colldist, force_nonnegative = TRUE)
  )
  m2 <- suppressWarnings(
    mrca_matrix(phylodist, colldist, force_nonnegative = TRUE)
  )
  expect_equal(m2[2,5], 0)
  expect_equal(m2[5,2], 0)
})

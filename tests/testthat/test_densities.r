context("Density (log lik.) functions")

#Pure multinomial sample 20x3 with 5 draws from p = (0.25, 0.5, 0.25)
X <- structure(c(2L, 1L, 2L, 3L, 0L, 1L, 2L, 1L, 0L, 0L, 3L, 2L, 3L, 
2L, 1L, 1L, 1L, 0L, 3L, 3L, 1L, 2L, 1L, 2L, 3L, 3L, 3L, 2L, 1L, 
3L, 0L, 1L, 2L, 2L, 4L, 3L, 0L, 5L, 1L, 2L, 2L, 2L, 2L, 0L, 2L, 
1L, 0L, 2L, 4L, 2L, 2L, 2L, 0L, 1L, 0L, 1L, 4L, 0L, 1L, 0L), .Dim = c(20L, 
3L))


test_that("single DM density functions work", {
  LL <- ddm(X, phi=0, p=c(0.25, 0.5, 0.25))
  expect_is(LL, "numeric")
  expect_equal(LL, -110.2104, tolerance=1e-3) #value calculated w/ dmultinom and ignoring multinomial coefficient   
  expect_less_than(ddm(X, phi=0.5, p=c(0.25, 0.5, 0.25)), LL)  
})

test_that("mixed DM density functions work", {
  res <- dmdm(X, phi=c(0.01, 0), p= rbind(c(0.25, 0.5, 0.25), c(0.04, 0.32, 0.64)), f=c(0.8, 0.2))
  expect_is(res$w, "matrix")
  expect_is(c(res$w), "numeric")
  expect_is(res$ll, "numeric")
  expect_true(all(sapply(rowSums(res$w), all.equal, 1.0)))
  
  worse <- dmdm(X, phi=c(0.01, 0), p= rbind(c(0.25, 0.5, 0.25), c(0.04, 0.32, 0.64)), f=c(0.2, 0.8))
  expect_less_than(worse$ll, res$ll)
})

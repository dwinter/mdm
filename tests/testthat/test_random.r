#test random generating functions return sensible results


context("random generating fxns")

genotypes <- matrix( 0.01/3, nrow=4, ncol=4)
diag(genotypes) <- 0.99
nfreq <- c(0.4, 0.1, 0.1, 0.4)

test_matrix <- function(m, expected_dims, test_name){
    expect_that(m,is_a("matrix"), label=test_name)
    expect_equal(dim(m), expected_dims, label=test_name)
    expect_true(is.integer(m), test_name)
}

test_that("Random generation functions behave", {
    
    #setup

    m <- rdm(4, 40, phi=0, p=genotypes)
    #test, contents desceibed above)
    test_matrix(m, c(4,4), "documentation")
    test_matrix(rdm(10, 5, p=c(0.9, 0.1),phi=0), c(10,2), "phi paramater")
    test_matrix(rdm(10, 5, scale=c(0.9,0.1)), c(10,2), "scale_paramater")
    
})

test_that("Paramater checking works", {
    expect_error(rdm(4, 40, phi=0), "Must specify 'p'")
    expect_error(rdm(4, 40), "Must specify either")
    expect_error(rdm(4, 40, p=c(0.1,0.7,0.1, 0.1)), "Must specify either")
    
    expect_error(rmdm(4, 40, phi=0), "Must specify 'p'")
    expect_error(rmdm(4, 40), "Must specify either")
    expect_error(rmdm(4, 40, p=c(0.1,0.7,0.1, 0.1)), "Must specify either")
})

test_that("Random mixtures work", {
    test_matrix(rmdm(10, 40, phi=0.01, f=nfreq, p=genotypes), c(10,4), "mixture w/ phi")
    test_matrix(rmdm(10, 40, scale=genotypes, f=nfreq), c(10,4), "mixture w/ scale")
})
    

#test random generating functions return sensible results


context("random")

test_matrix <- function(m, expected_dims, test_name){
    expect_that(m,is_a("matrix"), test_name)
    expet_qual(dim(m), expected_dims, test_name)
    expect_true(is.integer(m), test_name)
}

test_that("Random generation functions behave", {
    
    #setup
    genotypes <- matrix( 0.01/3, nrow=4, ncol=4)
    diag(genotypes) <- 0.99
    m <- rdm(4, 40, phi=0, p=genotypes)
    #test, contents desceibed above)
    test_matrix(m, c(4,4), "documentation")
    test_matrix(rdm(10, 5, p=c(0.9, 0.1),phi=0), c(10,2), "phi paramater")
    test_matrix( rdm(10,5, scale=c(0.9,0.1)), c(10,2), "scale_paramater")
    
})

test_that("Paramater checking works", {
    expect_error(rdm(4, 40, phi=0), "Must specify either")
    expect_error(rdm(4, 40), "Must specify either")
    expect_error(rdm(4, 40, p=c(0.1,0.7,0.1, 0.1), "Must specify either")
})


    

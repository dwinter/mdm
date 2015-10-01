context("Fitting fxns")
set.seed(321)
m <- rdm(4, 40, phi=0, p=c(0.25, 0.5, 0.25))
mod_1 <- fit_dm(m)

mm <- rmdm(100, 40, phi=0, f=c(0.7, 0.3), p=rbind( c(0.25, 0.5, 0.25), c(0.0025, 0.950, 0.9025)) )
mod_2 <- fit_mdm(mm, f=2)



test_that("Single component fitting works", {
    expect_is(mod_1, "mdm_model")
    expect_equal(length(mod_1$params), 4) 
 
})

test_that("Multiple component fitting works", {
    expect_is(mod_2, "mdm_model")
    expect_equal(dim(mod_2$params), c(2,4))
}) 

test_that("mdm model print methods work", {
    expect_output(mod_1, "Dirichelet multinomial  model")
    expect_output(mod_2, "Dirichelet multinomial mixture")

})

test_that("accessors for single mdm models work", {
    expect_is(logLik(mod_1), "numeric")
    expect_is(coef(mod_1), "numeric")
    expect_is(AIC(mod_1), "numeric")
    expect_is(BIC(mod_1), "numeric")
})

test_that("accessors for mixed mdm models work", {
    expect_is(logLik(mod_2), "numeric")
    expect_is(coef(mod_2), "matrix")
    expect_is(AIC(mod_2), "numeric")
    expect_is(BIC(mod_2), "numeric")
})

test_that("we can calculate group-assignments for a model", {
    comps <- get_component_probs(mm, mod_2)
    expect_error(get_component_probs(m, mod_1), "Model has only one component")
    expect_equal(dim(comps), c(dim(mm)[1], length(mod_2$f)))
    expect_true(all(sapply(rowSums(comps), all.equal, target=1.0)))
})

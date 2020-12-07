test_that("ANE is a number", {
  t<-matrix(sample(1:1000, 100), nrow=10, ncol=10)
  colnames(t)<-sample(letters,10)


  expect_length(ANE(t, sample(colnames(t),4)),1)
  expect_type(ENE(t,sample(colnames(t),4)), "double")
})

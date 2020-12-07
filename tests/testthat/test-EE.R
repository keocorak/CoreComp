
test_that("EE is a number", {
  t<-matrix(sample(1:1000, 100), nrow=10, ncol=10)
  colnames(t)<-rownames(t)<-sample(letters,10)
  expect_length(EE(t, sample(colnames(t),4)),1)
  expect_type(EE(t,sample(colnames(t),4)), "double")
})

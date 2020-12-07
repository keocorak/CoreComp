
test_that("noirot contribution value is between 0 and 1", {
  t<-matrix(sample(1:1000, 100), nrow=10, ncol=10)
  colnames(t)<-rownames(t)<-sample(letters,10)

  expect_lt(noirot.contribution(t, sample(colnames(t),4)),1)
  expect_gt(noirot.contribution(t, sample(colnames(t),4)),0)
})

test_that("Basic factor application to formula", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 2, 1, 1, 2, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)

  model <- Cox(a, b, c) ~ loglinear(d, factor(e))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, factor(x = e))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, factor(e, levels = c(1, 2), labels = c("low", "high")))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_high"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, factor(e, levels = c(2, 1), labels = c("high", "low")))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_low"), e$names)
})

test_that("Basic ns and bs application to formula", {
  df <- data.table("a" = 1:100, "b" = 2:101, "c" = c(rep(0, 20), rep(1, 80)), "d" = c(rep(1, 20), rep(2, 50), rep(3, 30)), "e" = 0:99)
  model <- Cox(a, b, c) ~ loglinear(d, ns(e, df = 2))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_ns1", "e_ns2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, ns(e, intercept = TRUE))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_ns1", "e_ns2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, bs(e))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_bs1", "e_bs2", "e_bs3"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, bs(e, Boundary.knots = c(0, 99)))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_bs1", "e_bs2", "e_bs3"), e$names)
})

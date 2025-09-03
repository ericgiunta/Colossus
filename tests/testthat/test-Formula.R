test_that("Coxph time column missing", {
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

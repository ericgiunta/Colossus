## ------------------------------------- ##
## Verify the table codes work
## ------------------------------------- ##
test_that("basic table check", {
  a <- c(0, 1, 2, 3, 4, 5, 6, 2, 2, 3, 4, 2, 1, 5, 6, 4, 2)
  b <- c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 3, 2, 2, 1)
  c <- c(0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1)

  d <- c(1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 4, 4, 2, 1, 1, 2)
  e <- c(1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1)
  f <- c(1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900)
  g <- c(4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 5, 5, 5, 4, 4, 4, 4)
  h <- c(6, 4, 4, 6, 6, 6, 4, 4, 4, 6, 6, 6, 6, 4, 4, 4, 4)
  i <- c(1901, 1902, 1903, 1904, 1905, 1906, 1907, 1903, 1904, 1903, 1904, 1910, 1903, 1904, 1903, 1904, 1910)
  table <- data.table::data.table(
    "a" = a, "b" = b, "c" = c,
    "d" = d, "e" = e, "f" = f,
    "g" = g, "h" = h, "i" = i
  )

  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "mean")
  expect_no_error(Event_Count_Gen(table, categ, summary))
  #
  table$a <- as.factor(table$a)
  categ <- list(
    "a" = "factor",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  time_scale <- list(
    "time AS time_bin" = list(
      "day" = c(1, 1, 1),
      "month" = c(1, 1, 1),
      "year" = c(1899, 1903, 1910)
    )
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  time_scale <- list(
    "time AS time_bin" = list(
      "day" = "d",
      "month" = "e",
      "year" = "f",
      categories = "0/10/20/30/40/50/60/70/80/90"
    )
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  time_scale <- list(
    "age AS agecat" = list(
      "day" = "d",
      "month" = "e",
      "year" = "f",
      "categories" = "0/10/20/30/40/50/60/70/80/90"
    ),
    "time AS time_bin" = list(
      "day" = c(1, 1, 1),
      "month" = c(1, 1, 1),
      "year" = c(1899, 1903, 1910)
    )
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  #
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  time_scale <- list(
    "age AS agecat" = list(
      "day" = "f",
      "month" = "e",
      "year" = "d",
      "categories" = "0/10/20/30/40/50/60/70/80/90"
    )
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  #
  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  time_scale <- list(
    "age AS agecat" = list(
      type = "none",
      "day" = "d",
      "month" = "e",
      "year" = "f",
      "categories" = "0/10/20/30/40/50/60/70/80/90"
    )
  )
  summary <- list("c" = "count", "a" = "mean BY cases", "b" = "weighted_mean", "b" = "sum BY cases AS bcases", "b" = "weighted_sum AS b_weight")
  events <- list("c AS cases")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  summary <- list("c" = "count AS cases")
  new_methods <- c("xsum", "xmean", "weighted_xsum", "weighted_xmean")
  for (i in new_methods) {
    summary <- c(summary, list("b" = i))
  }
  events <- list("c")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
})
test_that("person time, different intervals", {
  a <- c(0, 1, 2, 3, 4, 5, 6, 2, 2, 3, 4, 2, 1, 5, 6, 4, 2)
  b <- c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 3, 2, 2, 1)
  c <- c(0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1)

  d <- c(1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 4, 4, 2, 1, 1, 2)
  e <- c(1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1)
  f <- c(1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900)
  g <- c(4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 5, 5, 5, 4, 4, 4, 4)
  h <- c(6, 4, 4, 6, 6, 6, 4, 4, 4, 6, 6, 6, 6, 4, 4, 4, 4)
  i <- c(1901, 1902, 1903, 1904, 1905, 1906, 1907, 1903, 1904, 1903, 1904, 1910, 1903, 1904, 1903, 1904, 1910)
  table <- data.table::data.table(
    "a" = a, "b" = b, "c" = c,
    "d" = d, "e" = e, "f" = f,
    "g" = g, "h" = h, "i" = i
  )

  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "mean")
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  time_scale <- list("time AS time_bin" = list(
    "day" = c(1, 1, 1),
    "month" = c(1, 1, 1),
    "year" = c(1899, 1903, 1910)
  ))
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  pyr <- list(exit = list(year = "i", month = "h", day = "g"))
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry = list(year = "f", month = "e", day = "d"))
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_no_error(Event_Time_Gen(table, pyr, list(), categ, summary, events, TRUE))
  #
  time_scale <- list("age AS agecat" = list(
    "day" = "d",
    "month" = "e",
    "year" = "f",
    lower = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    upper = c(10, 20, 30, 40, 50, 60, 70, 80, 90)
  ))
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  pyr <- list(exit = list(year = "i", month = "h", day = "g"))
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry = list(year = "f", month = "e", day = "d"))
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_no_error(Event_Time_Gen(table, pyr, list(), categ, summary, events, TRUE))
})
test_that("person time, different intervals and user scale", {
  a <- c(0, 1, 2, 3, 4, 5, 6, 2, 2, 3, 4, 2, 1, 5, 6, 4, 2)
  b <- c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 3, 2, 2, 1)
  c <- c(0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1)
  d <- c(1, 1, 2, 2, 10, 10, 15, 15, 20, 20, 30, 30, 40, 40, 50, 50, 60)
  e <- c(11, 11, 22, 22, 30, 30, 45, 45, 50, 50, 60, 60, 70, 70, 80, 80, 90)
  table <- data.table::data.table(
    "a" = a, "b" = b, "c" = c,
    "d" = d, "e" = e
  )
  pyr <- list(entry = "d", exit = "e")
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  time_scale <- list("time AS time_bin" = "0 / 10 / 50 / 100")
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  pyr <- list(entry = "d")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(exit = "e")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry = "d", exit = "e")
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  #
  time_scale <- list("age AS agecat" = list(
    "day" = "d",
    "month" = "e",
    "year" = "f",
    lower = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    upper = c(10, 20, 30, 40, 50, 60, 70, 80, 90)
  ))
  expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  #
  time_scale <- list("age AS agecat" = list(
    lower = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    upper = c(10, 20, 30, 40, 50, 60, 70, 80, 90)
  ))
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  time_scale <- list("time AS time_bin" = list(categories = "0 / 10 / 50 / 100"))
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
})
test_that("basic table error check", {
  a <- c(0, 1, 2, 3, 4, 5, 6, 2, 2, 3, 4, 2, 1, 5, 6, 4, 2)
  b <- c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 3, 2, 2, 1)
  c <- c(0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1)

  d <- c(1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 4, 4, 2, 1, 1, 2)
  e <- c(1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1)
  f <- c(1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900)
  g <- c(4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 5, 5, 5, 4, 4, 4, 4)
  h <- c(6, 4, 4, 6, 6, 6, 4, 4, 4, 6, 6, 6, 6, 4, 4, 4, 4)
  i <- c(1901, 1902, 1903, 1904, 1905, 1906, 1907, 1903, 1904, 1903, 1904, 1910, 1903, 1904, 1903, 1904, 1910)
  table <- data.table::data.table(
    "a" = a, "b" = b, "c" = c,
    "d" = d, "e" = e, "f" = f,
    "g" = g, "h" = h, "i" = i
  )

  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  categ <- list(
    "a_bad" = "-1/3/5]7",
    "b" = list(
      lower = c(-1, 3, 6), upper = c(3, 6, 10),
      name = c("low", "medium", "high")
    )
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "mean")
  expect_error(Event_Count_Gen(table, categ, summary))

  categ <- list(
    "a_bad" = "-1/3/5]7",
    "b" = list(
      lower = c(-1, 3, 6), upper = c(3, 6, 10),
      name = c("low", "medium", "high")
    )
  )
  time_scale <- list("time AS time_bin" = list(
    "day" = c(1, 1, 1),
    "month" = c(1, 1, 1),
    "year" = c(1899, 1903, 1910)
  ))
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  categ <- list(
    "a" = "-1/3/5]7",
    "b" = list(
      lower = c(-1, 3, 6), upper = c(3, 6, 10),
      name = c("low", "medium", "high")
    )
  )
  expect_error(Event_Time_Gen(data.table(), pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry_not = list(year = "f", month = "e", day = "d"), exit_not = list(year = "i", month = "h", day = "g"))
  expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE, studyid = "bad"))
  pyr <- list(entry = list(duration = "e"), exit = list(year = "i", month = "h", day = "g"))
  expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE, TRUE))
  #
})
test_that("person time, different intervals", {
  a <- c(0, 1, 2, 3, 4, 5, 6, 2, 2, 3, 4, 2, 1, 5, 6, 4, 2)
  b <- c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 3, 2, 2, 1)
  c <- c(0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1)

  d <- c(1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 4, 4, 2, 1, 1, 2)
  e <- c(1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1)
  f <- c(1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900, 1900)
  g <- c(4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 5, 5, 5, 4, 4, 4, 4)
  h <- c(6, 4, 4, 6, 6, 6, 4, 4, 4, 6, 6, 6, 6, 4, 4, 4, 4)
  i <- c(1901, 1902, 1903, 1904, 1905, 1906, 1907, 1903, 1904, 1903, 1904, 1910, 1903, 1904, 1903, 1904, 1910)
  table <- data.table::data.table(
    "a" = a, "b" = b, "c" = c,
    "d" = d, "e" = e, "f" = f,
    "g" = g, "h" = h, "i" = i
  )

  pyr <- list(entry = list(year = "f", month = "e", day = "d"), exit = list(year = "i", month = "h", day = "g"))
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "mean")
  categ <- list(
    "a" = "-1/-1/3/5]7",
    "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
  )
  time_scale <- list("time AS time_bin" = list(
    "day" = c(1, 1, 1),
    "month" = c(1, 1, 1),
    "year" = c(1899, 1903, 1910)
  ))
  summary <- list("c" = "count AS cases", "a" = "mean", "b" = "weighted_mean")
  events <- list("c")
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        pyr_entry <- list()
        pyr_exit <- list()
        time_categ <- list()
        if (i == 1) {
          pyr_entry$day <- "d"
          pyr_exit$day <- "g"
        }
        if (j == 1) {
          pyr_entry$month <- "e"
          pyr_exit$month <- "h"
        }
        if (k == 1) {
          pyr_entry$year <- "f"
          pyr_exit$year <- "i"
        }
        if (i + j + k != 6) {
          pyr <- list(entry = pyr_entry, exit = pyr_exit)
          categ <- list(
            "a" = "-1/-1/3/5]7",
            "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
          )
          time_scale <- list("time AS time_bin" = list(
            "day" = c(1, 1, 1),
            "month" = c(1, 1, 1),
            "year" = c(1899, 1903, 1910)
          ))
          expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
          pyr <- list(entry = pyr_entry)
          expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
          pyr <- list(exit = pyr_exit)
          expect_no_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
          pyr <- list(entry = pyr_entry, exit = pyr_exit)
          categ <- list(
            "a" = "-1/-1/3/5]7",
            "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
          )
          expect_no_error(Event_Time_Gen(table, pyr, list(), categ, summary, events, TRUE))
        } else {
          pyr <- list(entry = pyr_entry, exit = pyr_exit)
          categ <- list(
            "a" = "-1/-1/3/5]7",
            "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
          )
          time_scale <- list("time AS time_bin" = list(
            "day" = c(1, 1, 1),
            "month" = c(1, 1, 1),
            "year" = c(1899, 1903, 1910)
          ))
          expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
          pyr <- list(entry = pyr_entry)
          expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
          pyr <- list(exit = pyr_exit)
          expect_error(Event_Time_Gen(table, pyr, time_scale, categ, summary, events, TRUE))
          pyr <- list(entry = pyr_entry, exit = pyr_exit)
          categ <- list(
            "a" = "-1/-1/3/5]7",
            "b AS b_bin" = list(lower = c(-1, -1, 3, 6), upper = c(-1, 3, 6, "]10"))
          )
          expect_error(Event_Time_Gen(table, pyr, list(), categ, summary, events, TRUE))
        }
      }
    }
  }
})

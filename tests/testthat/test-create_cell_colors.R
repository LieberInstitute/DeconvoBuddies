my_colors <- c("darkorchid4", "deeppink4", "aquamarine3", "darkolivegreen1")

test_that("Error with no color input", {
  expect_error(create_cell_colors(cell_type = c("A", "B", "C", "D"),
                                  palette_name = NULL))
  
})

test_that("Error with bad palette_name put", {
  expect_error(create_cell_colors(cell_type = c("A", "B", "C", "D"),
                                  palette_name = "NOT A palette"))
  
})

test_that("Error with less color than ct", {
  expect_error(create_cell_colors(cell_type = c("A", "B", "C", "D"),
                                  palette = my_colors[1:3]))
})

create_cell_colors(cell_type = c("A"),
                   palette = my_colors[1:3])

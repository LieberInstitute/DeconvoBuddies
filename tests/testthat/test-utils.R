test_that("expected 2-3 colors from gg_color_hue", {
    expect_equal(gg_color_hue(2), c("#F8766D", "#00BFC4"))
    expect_equal(gg_color_hue(3), c("#F8766D", "#00BA38", "#619CFF"))
})

test_that("expected string split from ss", {
    expect_equal(ss("Keep.Drop", "\\."), "Keep")
})

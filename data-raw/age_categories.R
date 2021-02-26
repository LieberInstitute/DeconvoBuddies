## code to prepare `age_categories` dataset goes here
age_categories <- list(Prenatal = list(min = -Inf, max = 0),
                       Infant = list(min = 0, max = 1),
                       Child = list(min = 1, max = 10),
                       Teen = list(min = 10, max = 20),
                       Adult = list(min = 20, max = 50),
                       `50+` = list(min = 50, max = Inf))

usethis::use_data(age_categories, overwrite = TRUE)

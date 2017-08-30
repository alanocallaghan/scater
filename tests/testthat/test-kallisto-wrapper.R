## Testing kallisto wrapper functions

context("test improper input usage")

test_that("failure is as expected for wrong input", {
    directories <- file.path("foo", "bar", 1:10)
    
    expect_error(readKallistoResults(samples = 1:10, directories = directories),
                "Some of the desired directories to import do not exist!")
})

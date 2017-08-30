# Tests for tximport wrapper

context("test improper input usage")

test_that("failure is as expected for wrong input", {
    directories <- file.path("foo", "bar", 1:10)
    files <- file.path(directories, "abundance.tsv")
    expect_that(readTxResults(samples = 1:10, files = files,
                              type = "kallisto"),
                    throws_error("Some of the files provided do not exist!"))
})

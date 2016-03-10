# Tests for tximport wrapper

context("test imporper input usage")

test_that("failure is as expected for wrong input", {
    directories <- file.path("foo", "bar", 1:10)
    expect_that(readTxResults(samples = 1:10, directories = directories, 
                              type = "kallisto"),
                    throws_error("Some of the directories do not contain the expected abundance.tsv file!"))
})

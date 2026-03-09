# Tests for kernel functions

test_that("kernel_bartlett produces correct values", {
  # At x = 0, should be 1
  expect_equal(kernel_bartlett(0), 1)

  # At x = 1, should be 0
  expect_equal(kernel_bartlett(1), 0)

  # At x = 0.5, should be 0.5
  expect_equal(kernel_bartlett(0.5), 0.5)

  # Symmetric in absolute value
  expect_equal(kernel_bartlett(-0.5), kernel_bartlett(0.5))

  # Outside [-1, 1], should be 0
  expect_equal(kernel_bartlett(2), 0)
  expect_equal(kernel_bartlett(-1.5), 0)
})

test_that("kernel_parzen produces correct values", {
  # At x = 0, should be 1
  expect_equal(kernel_parzen(0), 1)

  # At x = 1, should be 0
  expect_equal(kernel_parzen(1), 0)

  # Symmetric in absolute value
  expect_equal(kernel_parzen(-0.3), kernel_parzen(0.3))

  # Non-negative
  x <- seq(0, 1, by = 0.1)
  expect_true(all(kernel_parzen(x) >= 0))
})

test_that("kernel_qs produces correct values", {
  # At x = 0, should be 1
  expect_equal(kernel_qs(0), 1)

  # Symmetric in absolute value
  expect_equal(kernel_qs(-0.5), kernel_qs(0.5))

  # QS kernel can be negative for some values
  # But should be bounded
  x <- seq(0, 5, by = 0.1)
  vals <- kernel_qs(x)
  expect_true(all(abs(vals) <= 1.1))  # Allow small tolerance
})

test_that("get_kernel_function returns correct function", {
  expect_identical(get_kernel_function("bartlett"), kernel_bartlett)
  expect_identical(get_kernel_function("parzen"), kernel_parzen)
  expect_identical(get_kernel_function("qs"), kernel_qs)

  # Case insensitive
  expect_identical(get_kernel_function("BARTLETT"), kernel_bartlett)

  # Invalid kernel

  expect_error(get_kernel_function("invalid"), "Unknown kernel")
})

test_that("get_kernel_exponent returns correct values", {
  expect_equal(get_kernel_exponent("bartlett"), 1)
  expect_equal(get_kernel_exponent("parzen"), 2)
  expect_equal(get_kernel_exponent("qs"), 2)
})

test_that("get_kernel_constant returns correct values", {
  # Andrews (1991) rate constants
  expect_equal(get_kernel_constant("bartlett"), 1.1447)
  expect_equal(get_kernel_constant("parzen"), 2.6614)
  expect_equal(get_kernel_constant("qs"), 1.3221)
})

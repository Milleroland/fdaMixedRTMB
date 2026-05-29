# fdaMixedRTMB

## Installation

To install the package run the following code
```r
# install.packages("remotes")
remotes::install_github("Milleroland/fdaMixedRTMB")
```

## Package introduction
`fdaMixedRTMB` is an R package for operator-based functional mixed models fit
with RTMB. It focuses on the `K.order = 1` operators from `fdaMixed`:

- `operator_k1(lambda = 1)` for `L f = -lambda_d2 * f''`
- `operator_k1(lambda = c(1, 1))` for
  `L f = -lambda_d2 * f'' + lambda_level * f`

Boundary operators are supplied as a numeric length-two vector `c(alpha0, alpha1)`,
where `alpha0` multiplies the function value and `alpha1` multiplies the
derivative at the relevant endpoint. A Dirichlet condition is `c(1, 0)` and a
Neumann condition is `c(0, 1)`.

## Main formula interface

```r
library(fdaMixedRTMB)

source(system.file("examples/boundary_value_example.R", package = "fdaMixedRTMB"))
example_data <- make_boundary_value_example()

fit <- fdaLm_rtmb(
  y | curve ~ 0,
  data = example_data$data,
  boundary_value = TRUE,
  boundary_formula = ~ 1 + boundary_x,
  operator = "laplace",
  left_boundary  = c(1, 0),
  right_boundary = c(1, 0)
)

summary(fit)
coef(fit, component = "boundary")
coef(fit, component = "rho")
plot(fit)
```

The public fitting interface is formula-based: fixed and random effects are
specified in the `y | curve ~ fixed | random` formula, and right-boundary
regressions are specified with `boundary_formula`. The package estimates
boundary values, not boundary conditions. The serial effect is fit in the
homogeneous boundary space, while the nonzero right-boundary value is carried by
a fixed null-space basis and can vary by curve-level covariates. Set
`boundary_value = FALSE` in `fdaLm_rtmb()` to fit the same serial model with a
homogeneous right boundary and no estimated boundary-value fixed component.

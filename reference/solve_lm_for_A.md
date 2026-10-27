# Solve the equation `Y = A B` for `A` or `B`

Solve the equation `Y = A B` for `A` or `B`

## Usage

``` r
solve_lm_for_A(Y, B, w = NULL)

solve_lm_for_B(Y, A, w = NULL)
```

## Arguments

- Y:

  the left side of the equation

- w:

  a vector with weights. If `NULL` it is ignored, otherwise it must be
  of length 1 or have the same length as columns in `Y`. Default: `NULL`

- A, B:

  the known matrix on the right side of the equation

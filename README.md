# sml-itersolve

[![CI](https://github.com/sjqtentacles/sml-itersolve/actions/workflows/ci.yml/badge.svg)](https://github.com/sjqtentacles/sml-itersolve/actions/workflows/ci.yml)

Iterative linear solvers for sparse systems in Standard ML. Built on top of
**sml-sparse** (CSR format) and **sml-matrix** (dense reference solver).

## Implemented solvers

| Solver | Function | Suitable for |
|---|---|---|
| Conjugate Gradient | `cg` | Symmetric positive-definite (SPD) sparse systems |

## API sketch

```sml
(* Build a 3×3 SPD sparse system: A*x = b, with true solution x = [1,2,3] *)
val A : Sparse.csr = Sparse.fromCoo {
  rows = 3, cols = 3,
  entries = [(0,0,4.0),(0,1,1.0),(1,0,1.0),(1,1,3.0),(1,2,1.0),(2,1,1.0),(2,2,2.0)]
}
val b = [6.0, 10.0, 7.0]

(* Solve with CG (max 1000 iterations, starting from zeros) *)
val x0 = [0.0, 0.0, 0.0]
val xCG : real list = IterSolve.cg A b x0 1000
(* ≈ [1.0, 2.0, 3.0] *)

(* Compare CG solution against dense Matrix.solve reference *)
val (xCG, xDense) = IterSolve.compareDense A b
(* |xCG[i] - xDense[i]| < 1e-6 for well-conditioned systems *)
```

## Known limitations

- **CG only applies to SPD matrices**: passing a non-symmetric or indefinite
  matrix produces silently incorrect results or divergence.
- **No preconditioner**: convergence depends on the condition number of `A`;
  ill-conditioned systems may require many iterations or fail to converge.
- **No GMRES / BiCGSTAB**: general non-symmetric systems are not yet supported.
- `compareDense` converts the sparse matrix to dense internally — suitable for
  testing small systems only.

## Installing with smlpkg

```sh
smlpkg add github.com/sjqtentacles/sml-itersolve
smlpkg sync
```

Reference from your `.mlb`:

```
lib/github.com/sjqtentacles/sml-itersolve/itersolve.mlb
```

## Building and testing

```sh
make test        # MLton
make test-poly   # Poly/ML
make all-tests   # both
make clean
```

## Project layout

```
sml.pkg
Makefile
lib/github.com/sjqtentacles/sml-itersolve/
  itersolve.sig     ITERSOLVE signature
  itersolve.sml     CG solver + compareDense
  itersolve.mlb
test/
  test.sml          2×2 and 3×3 SPD system tests vs Matrix.solve
```

## License

MIT. See [LICENSE](LICENSE).

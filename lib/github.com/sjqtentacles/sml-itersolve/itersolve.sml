structure IterSolve :> ITERSOLVE =
struct
  fun dot a b = ListPair.foldl (fn (x, y, s) => s + x * y) 0.0 (a, b)

  fun vecAdd a b = ListPair.map (fn (x, y) => x + y) (a, b)
  fun vecSub a b = ListPair.map (fn (x, y) => x - y) (a, b)
  fun vecScale s v = List.map (fn x => s * x) v

  (* Conjugate Gradient for SPD sparse system Ax = b.
     x0 = initial guess, maxIt = max iterations.
     Convergence: when ||r|| < 1e-10 or maxIt reached. *)
  fun cg (A : Sparse.csr) b x0 maxIt =
    let
      val r0 = vecSub b (Sparse.spmv A x0)
      fun loop x r p rsold iter =
        if iter >= maxIt orelse Math.sqrt rsold < 1E~10 then x
        else
          let
            val ap    = Sparse.spmv A p
            val pap   = dot p ap
            val alpha = rsold / pap
            val x'    = vecAdd x (vecScale alpha p)
            val r'    = vecSub r (vecScale alpha ap)
            val rsnew = dot r' r'
            val beta  = rsnew / rsold
            val p'    = vecAdd r' (vecScale beta p)
          in
            loop x' r' p' rsnew (iter + 1)
          end
    in
      loop x0 r0 r0 (dot r0 r0) 0
    end

  (* Solve Ax = b by both CG and direct dense LU solve (via Matrix.solve).
     Returns (cg_solution, dense_solution) so they can be compared. *)
  fun compareDense (A : Sparse.csr) b =
    let
      val n  = #rows A
      val x0 = List.tabulate (n, fn _ => 0.0)
      val xCG = cg A b x0 1000

      (* Build dense matrix from COO representation *)
      val { rows, cols, entries } = Sparse.toCoo A
      val denseRows = List.tabulate (rows, fn r =>
        List.tabulate (cols, fn c =>
          case List.find (fn (rr, cc, _) => rr = r andalso cc = c) entries of
            SOME (_, _, v) => v | NONE => 0.0))
      val m = Matrix.fromRows denseRows
      val xDense = Matrix.solve m b
    in
      (xCG, xDense)
    end
end

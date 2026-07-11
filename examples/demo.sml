(* demo.sml - conjugate-gradient and dense-LU solves of small sparse SPD
   linear systems, cross-checked against each other. Deterministic: identical
   output on every run and both compilers. *)

structure I = IterSolve

fun fmtR r = Real.fmt (StringCvt.FIX (SOME 4)) (if Real.== (r, 0.0) then 0.0 else r)
fun fmtVec v = "[" ^ String.concatWith ", " (List.map fmtR v) ^ "]"

val () = print "Sparse SPD system A x = b, A = [[4,1,0],[1,3,1],[0,1,2]], b = [1,2,3]:\n"

val coo = { rows = 3, cols = 3
          , entries = [ (0,0,4.0),(0,1,1.0)
                      , (1,0,1.0),(1,1,3.0),(1,2,1.0)
                      , (2,1,1.0),(2,2,2.0) ] }
val A = Sparse.fromCoo coo
val b = [1.0, 2.0, 3.0]
val x0 = [0.0, 0.0, 0.0]

val () = print ("  spmv A [1,1,1]      = " ^ fmtVec (Sparse.spmv A [1.0,1.0,1.0]) ^ "\n")

val xCG1 = I.cg A b x0 1
val () = print ("  CG after 1 iter     = " ^ fmtVec xCG1 ^ "\n")

val (xCG, xDense) = I.compareDense A b
val () = print ("  CG (converged)      = " ^ fmtVec xCG ^ "\n")
val () = print ("  dense LU solve      = " ^ fmtVec xDense ^ "\n")

val () = print "\nSmaller SPD system, A = [[4,1],[1,3]], b = [1,2]:\n"
val coo2 = { rows = 2, cols = 2
           , entries = [(0,0,4.0),(0,1,1.0),(1,0,1.0),(1,1,3.0)] }
val A2 = Sparse.fromCoo coo2
val (xCG2, xDense2) = I.compareDense A2 [1.0, 2.0]
val () = print ("  CG (converged)      = " ^ fmtVec xCG2 ^ "\n")
val () = print ("  dense LU solve      = " ^ fmtVec xDense2 ^ "\n")

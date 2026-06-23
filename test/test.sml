structure Tests = struct open Harness structure I = IterSolve

(* 3x3 SPD system:
   A = [[4,1,0],[1,3,1],[0,1,2]], b = [1,2,3]
   Solution (from Gaussian elimination): approx [0.0435, 0.3043, 1.3478] *)
val coo3 = { rows = 3, cols = 3
           , entries = [ (0,0,4.0),(0,1,1.0)
                       , (1,0,1.0),(1,1,3.0),(1,2,1.0)
                       , (2,1,1.0),(2,2,2.0) ] }

(* 2x2 SPD system: A = [[4,1],[1,3]], b = [1,2]
   Solution: x = [1/11, 7/11] *)
val coo2 = { rows = 2, cols = 2
           , entries = [(0,0,4.0),(0,1,1.0),(1,0,1.0),(1,1,3.0)] }

fun run () = let
  val () = section "CG vs dense on 2x2 SPD"
  val A2  = Sparse.fromCoo coo2
  val b2  = [1.0, 2.0]
  val (xCG2, xD2) = I.compareDense A2 b2
  val () = checkRealTol 1E~3 "x[0] CG vs dense" (List.nth (xD2, 0), List.nth (xCG2, 0))
  val () = checkRealTol 1E~3 "x[1] CG vs dense" (List.nth (xD2, 1), List.nth (xCG2, 1))
  (* Known exact values *)
  val () = checkRealTol 1E~3 "x[0] exact" (1.0/11.0, List.nth (xCG2, 0))

  val () = section "CG vs dense on 3x3 SPD"
  val A3  = Sparse.fromCoo coo3
  val b3  = [1.0, 2.0, 3.0]
  val (xCG3, xD3) = I.compareDense A3 b3
  (* CG should match dense within 1e-3 *)
  val () = checkRealTol 1E~3 "x3[0] CG vs dense" (List.nth (xD3, 0), List.nth (xCG3, 0))
  val () = checkRealTol 1E~3 "x3[1] CG vs dense" (List.nth (xD3, 1), List.nth (xCG3, 1))
  val () = checkRealTol 1E~3 "x3[2] CG vs dense" (List.nth (xD3, 2), List.nth (xCG3, 2))
in Harness.run () end end

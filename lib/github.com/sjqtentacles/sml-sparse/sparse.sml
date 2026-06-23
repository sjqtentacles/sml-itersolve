structure Sparse :> SPARSE =
struct
  type coo = { rows : int, cols : int, entries : (int * int * real) list }
  type csr = { rows : int, cols : int, vals : real array, colIdx : int array, rowPtr : int array }

  fun cmpEntry ((r1,c1,_),(r2,c2,_)) =
    case Int.compare (r1, r2) of EQUAL => Int.compare (c1, c2) | ord => ord

  fun sortEntries entries =
    let fun insert (x, []) = [x]
          | insert (x, y::ys) = if cmpEntry (x,y) = GREATER then y :: insert (x, ys) else x :: y :: ys
    in List.foldr insert [] entries end

  fun fromCoo { rows, cols, entries } =
    let
      val sorted = sortEntries entries
      val vals = Array.fromList (List.map (fn (_,_,v) => v) sorted)
      val colIdx = Array.fromList (List.map (fn (_,c,_) => c) sorted)
      val rowPtr = Array.array (rows + 1, 0)
      val () = List.app (fn (r,_,_) => Array.update (rowPtr, r + 1, Array.sub (rowPtr, r + 1) + 1)) sorted
      val () = List.app (fn i => Array.update (rowPtr, i, Array.sub (rowPtr, i) + Array.sub (rowPtr, i - 1)))
                        (List.tabulate (rows, fn i => i + 1))
    in { rows = rows, cols = cols, vals = vals, colIdx = colIdx, rowPtr = rowPtr } end

  fun toCoo { rows, cols, vals, colIdx, rowPtr } =
    let
      fun rowEntries r =
        let val start = Array.sub (rowPtr, r) val stop = Array.sub (rowPtr, r + 1)
        in List.tabulate (stop - start, fn k => (r, Array.sub (colIdx, start + k), Array.sub (vals, start + k))) end
    in { rows = rows, cols = cols, entries = List.concat (List.tabulate (rows, rowEntries)) } end

  fun transpose csr =
    let val { rows, cols, entries } = toCoo csr
    in fromCoo { rows = cols, cols = rows, entries = List.map (fn (r,c,v) => (c,r,v)) entries } end

  fun spmv { rows, cols, vals, colIdx, rowPtr } x =
    let
      val () = if List.length x <> cols then raise Matrix.Dim "spmv: vector length" else ()
      val xa = Array.fromList x
      fun dotRow r =
        let val start = Array.sub (rowPtr, r) val stop = Array.sub (rowPtr, r + 1)
            fun loop (i, acc) = if i >= stop then acc else loop (i + 1, acc + Array.sub (vals, i) * Array.sub (xa, Array.sub (colIdx, i)))
        in loop (start, 0.0) end
    in List.tabulate (rows, dotRow) end

  fun denseParity csr x =
    let
      val { rows, cols, entries } = toCoo csr
      val denseRows = List.tabulate (rows, fn r =>
        List.tabulate (cols, fn c =>
          let fun find [] = 0.0 | find ((rr,cc,v)::rest) = if rr = r andalso cc = c then v else find rest
          in find entries end))
      val m = Matrix.fromRows denseRows
      val b = Matrix.fromRows [x]
      val prod = Matrix.mul (m, Matrix.transpose b)
      val dense = List.tabulate (rows, fn i => Matrix.sub (prod, i, 0))
    in (spmv csr x, dense) end
end

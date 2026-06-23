(* sparse.sig — CSR/COO sparse matrices. *)

signature SPARSE =
sig
  type coo = { rows : int, cols : int, entries : (int * int * real) list }
  type csr = { rows : int, cols : int, vals : real array, colIdx : int array, rowPtr : int array }

  val fromCoo : coo -> csr
  val toCoo : csr -> coo
  val transpose : csr -> csr
  val spmv : csr -> real list -> real list
  val denseParity : csr -> real list -> real list * real list
end

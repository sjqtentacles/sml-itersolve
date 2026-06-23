signature ITERSOLVE =
sig
  val cg : Sparse.csr -> real list -> real list -> int -> real list
  val compareDense : Sparse.csr -> real list -> real list * real list
end

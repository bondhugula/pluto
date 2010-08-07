[PIP2-like future input] Please enter:
- the context matrix,
0 3
- the bignum column (start at 0, -1 if no bignum),
-1
- the constraint matrix.
4 5
 1 -1 0 1 0
 1 0 -1 1 0
 1 -1 3 -2 12
 1 2 -1 -1 3

(if #[ -1 3]
 (list
  #[ 0 0]
  #[ 0 0]
 )
 (if #[ -1 5]
  (newparm 1 (div #[ 1 1] 2))
  (list
   #[ 1 -1 -1]
   #[ 0 0 0]
  )
  (list
   #[ 1 -4]
   #[ 1 -5]
  )
 )
)

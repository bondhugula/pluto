[PIP2-like future input] Please enter:
- the context matrix,
0 4
- the bignum column (start at 0, -1 if no bignum),
-1
- the constraint matrix.
7 6
 1 2 6 0 0 -9
 1 5 -3 0 0 0
 1 2 -10 0 0 15
 1 -2 6 0 0 -3
 1 -2 -6 0 0 17
 1 0 1 -1 0 0
 1 1 0 0 -1 0

(if #[ -1 0 1]
 (if #[ 0 -1 0]
  ()
  (if #[ 0 -1 5]
   (if #[ -3 -1 5]
    (if #[ 0 -1 1]
     (newparm 2 (div #[ 0 1 1] 3))
     (newparm 3 (div #[ 0 1 5 3] 6))
     (newparm 4 (div #[ 0 2 0 1 0] 3))
     ()
     (if #[ 0 -1 3]
      (newparm 2 (div #[ 0 2 1] 3))
      ()
      ()
     )
    )
    (if #[ 0 -1 3]
     (newparm 2 (div #[ 0 2 1] 3))
     ()
     ()
    )
   )
   ()
  )
 )
 ()
)

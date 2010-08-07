#ifndef _EXT_EHRHART_H_
#define _EXT_EHRHART_H_

#ifndef PARAMS
# if (defined(__STDC__) || defined(__cplusplus))
#  define PARAMS(protos) protos
# else /* no (defined(__STDC__) || defined(__cplusplus)) */
#  define PARAMS(protos) ()
# endif /* no (defined(__STDC__) || defined(__cplusplus)) */
#endif

extern Enumeration *Domain_Enumerate
    PARAMS((Polyhedron *D, Polyhedron *C, unsigned MAXRAYS, char **pn));

extern void new_eadd PARAMS((evalue *e1,evalue *res));

#endif

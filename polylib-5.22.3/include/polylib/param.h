#ifndef _param_H_
#define _param_H_
#if (defined(__STDC__) || defined(__cplusplus))

#if defined(__cplusplus)
extern "C" {
#endif

extern char **Read_ParamNames(FILE *in, int m);

#if defined(__cplusplus)
}
#endif

#else /* (defined(__STDC__) || defined(__cplusplus)) */

extern char **Read_ParamNames(/* FILE *in, int m */);

#endif /* (defined(__STDC__) || defined(__cplusplus)) */
#endif /* _param_H_ */

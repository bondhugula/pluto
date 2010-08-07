#ifndef ISL_STREAM_H
#define ISL_STREAM_H

#include <stdio.h>

#if defined(__cplusplus)
extern "C" {
#endif

enum isl_token_type { ISL_TOKEN_UNKNOWN = 256, ISL_TOKEN_VALUE,
			ISL_TOKEN_IDENT, ISL_TOKEN_GE,
			ISL_TOKEN_LE, ISL_TOKEN_TO, ISL_TOKEN_AND,
			ISL_TOKEN_EXISTS };

struct isl_token {
	enum isl_token_type  type;

	unsigned int on_new_line : 1;
	int line;
	int col;

	union {
		isl_int	v;
		char	*s;
	} u;
};

struct isl_stream {
	struct isl_ctx	*ctx;
	FILE        	*file;
	const char  	*str;
	int	    	line;
	int	    	col;
	int	    	eof;

	char	    	*buffer;
	size_t	    	size;
	size_t	    	len;
	int	    	c;

	struct isl_token	*tokens[5];
	int	    	n_token;
};

struct isl_stream* isl_stream_new_file(struct isl_ctx *ctx, FILE *file);
struct isl_stream* isl_stream_new_str(struct isl_ctx *ctx, const char *str);
void isl_stream_free(struct isl_stream *s);

struct isl_token *isl_stream_next_token(struct isl_stream *s);
void isl_stream_push_token(struct isl_stream *s, struct isl_token *tok);
int isl_stream_eat(struct isl_stream *s, int type);

#if defined(__cplusplus)
}
#endif

#endif

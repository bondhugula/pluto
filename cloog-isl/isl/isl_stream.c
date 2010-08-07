#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <isl_ctx.h>
#include "isl_stream.h"

static struct isl_token *isl_token_new(struct isl_ctx *ctx,
	int line, int col, unsigned on_new_line)
{
	struct isl_token *tok = isl_alloc_type(ctx, struct isl_token);
	if (!tok)
		return NULL;
	tok->line = line;
	tok->col = col;
	tok->on_new_line = on_new_line;
	return tok;
}

void isl_token_free(struct isl_token *tok)
{
	if (!tok)
		return;
	if (tok->type == ISL_TOKEN_VALUE)
		isl_int_clear(tok->u.v);
	else if (tok->type == ISL_TOKEN_IDENT)
		free(tok->u.s);
	free(tok);
}

void isl_stream_error(struct isl_stream *s, struct isl_token *tok, char *msg)
{
	int line = tok ? tok->line : s->line;
	int col = tok ? tok->col : s->col;
	fprintf(stderr, "syntax error (%d, %d): %s\n", line, col, msg);
	if (tok) {
		if (tok->type < 256)
			fprintf(stderr, "got '%c'\n", tok->type);
		else
			fprintf(stderr, "got token type %d\n", tok->type);
	}
}

static struct isl_stream* isl_stream_new(struct isl_ctx *ctx)
{
	int i;
	struct isl_stream *s = isl_alloc_type(ctx, struct isl_stream);
	if (!s)
		return NULL;
	s->ctx = ctx;
	isl_ctx_ref(s->ctx);
	s->size = 256;
	s->file = NULL;
	s->str = NULL;
	s->buffer = isl_alloc_array(ctx, char, s->size);
	if (!s->buffer)
		goto error;
	s->len = 0;
	s->line = 1;
	s->col = 0;
	s->eof = 0;
	s->c = -1;
	for (i = 0; i < 5; ++i)
		s->tokens[i] = NULL;
	s->n_token = 0;
	return s;
error:
	isl_stream_free(s);
	return NULL;
}

struct isl_stream* isl_stream_new_file(struct isl_ctx *ctx, FILE *file)
{
	struct isl_stream *s = isl_stream_new(ctx);
	if (!s)
		return NULL;
	s->file = file;
	return s;
}

struct isl_stream* isl_stream_new_str(struct isl_ctx *ctx, const char *str)
{
    struct isl_stream *s = isl_stream_new(ctx);
    s->str = str;
    return s;
}

static int isl_stream_getc(struct isl_stream *s)
{
	int c;
	if (s->eof)
		return -1;
	if (s->file)
		c = fgetc(s->file);
	else {
		c = *s->str++;
		if (c == '\0')
			c = -1;
	}
	if (c == -1)
		s->eof = 1;
	if (!s->eof) {
		if (s->c == '\n') {
			s->line++;
			s->col = 0;
		} else
			s->col++;
	}
	s->c = c;
	return c;
}

static void isl_stream_ungetc(struct isl_stream *s, int c)
{
	if (s->file)
		ungetc(c, s->file);
	else
		--s->str;
	s->c = -1;
}

static int isl_stream_push_char(struct isl_stream *s, int c)
{
	if (s->len >= s->size) {
		s->size = (3*s->size)/2;
		s->buffer = isl_realloc_array(ctx, s->buffer, char, s->size);
		if (!s->buffer)
			return -1;
	}
	s->buffer[s->len++] = c;
	return 0;
}

void isl_stream_push_token(struct isl_stream *s, struct isl_token *tok)
{
	isl_assert(s->ctx, s->n_token < 5, return);
	s->tokens[s->n_token++] = tok;
}

struct isl_token *isl_stream_next_token(struct isl_stream *s)
{
	int c;
	struct isl_token *tok = NULL;
	int line, col;
	int old_line = s->line;

	if (s->n_token)
		return s->tokens[--s->n_token];

	s->len = 0;

	/* skip spaces */
	while ((c = isl_stream_getc(s)) != -1 && isspace(c))
		/* nothing */
		;

	line = s->line;
	col = s->col;

	if (c == -1)
		return NULL;
	if (c == '(' ||
	    c == ')' ||
	    c == '+' ||
	    c == '/' ||
	    c == '*' ||
	    c == '^' ||
	    c == '=' ||
	    c == ',' ||
	    c == ':' ||
	    c == '[' ||
	    c == ']' ||
	    c == '{' ||
	    c == '}') {
		tok = isl_token_new(s->ctx, line, col, old_line != line);
		if (!tok)
			return NULL;
		tok->type = (enum isl_token_type)c;
		return tok;
	}
	if (c == '-') {
		int c;
		if ((c = isl_stream_getc(s)) == '>') {
			tok = isl_token_new(s->ctx, line, col, old_line != line);
			if (!tok)
				return NULL;
			tok->type = ISL_TOKEN_TO;
			return tok;
		}
		if (c != -1)
			isl_stream_ungetc(s, c);
	}
	if (c == '-' || isdigit(c)) {
		tok = isl_token_new(s->ctx, line, col, old_line != line);
		if (!tok)
			return NULL;
		tok->type = ISL_TOKEN_VALUE;
		isl_int_init(tok->u.v);
		if (isl_stream_push_char(s, c))
			goto error;
		while ((c = isl_stream_getc(s)) != -1 && isdigit(c))
			if (isl_stream_push_char(s, c))
				goto error;
		if (c != -1)
			isl_stream_ungetc(s, c);
		if (s->len == 1 && s->buffer[0] == '-')
			isl_int_set_si(tok->u.v, -1);
		else {
			isl_stream_push_char(s, '\0');
			isl_int_read(tok->u.v, s->buffer);
		}
		return tok;
	}
	if (isalpha(c)) {
		tok = isl_token_new(s->ctx, line, col, old_line != line);
		if (!tok)
			return NULL;
		isl_stream_push_char(s, c);
		while ((c = isl_stream_getc(s)) != -1 && isalnum(c))
			isl_stream_push_char(s, c);
		if (c != -1)
			isl_stream_ungetc(s, c);
		isl_stream_push_char(s, '\0');
		if (!strcasecmp(s->buffer, "exists"))
			tok->type = ISL_TOKEN_EXISTS;
		else {
			tok->type = ISL_TOKEN_IDENT;
			tok->u.s = strdup(s->buffer);
		}
		return tok;
	}
	if (c == '>') {
		int c;
		if ((c = isl_stream_getc(s)) == '=') {
			tok = isl_token_new(s->ctx, line, col, old_line != line);
			if (!tok)
				return NULL;
			tok->type = ISL_TOKEN_GE;
			return tok;
		}
		if (c != -1)
			isl_stream_ungetc(s, c);
	}
	if (c == '<') {
		int c;
		if ((c = isl_stream_getc(s)) == '=') {
			tok = isl_token_new(s->ctx, line, col, old_line != line);
			if (!tok)
				return NULL;
			tok->type = ISL_TOKEN_LE;
			return tok;
		}
		if (c != -1)
			isl_stream_ungetc(s, c);
	}
	if (c == '&') {
		tok = isl_token_new(s->ctx, line, col, old_line != line);
		if (!tok)
			return NULL;
		tok->type = ISL_TOKEN_AND;
		if ((c = isl_stream_getc(s)) != '&' && c != -1)
			isl_stream_ungetc(s, c);
		return tok;
	}

	tok = isl_token_new(s->ctx, line, col, old_line != line);
	if (!tok)
		return NULL;
	tok->type = ISL_TOKEN_UNKNOWN;
	return tok;
error:
	isl_token_free(tok);
	return NULL;
}

int isl_stream_eat(struct isl_stream *s, int type)
{
	struct isl_token *tok;

	tok = isl_stream_next_token(s);
	if (!tok)
		return -1;
	if (tok->type == type) {
		isl_token_free(tok);
		return 0;
	}
	isl_stream_error(s, tok, "expecting other token");
	isl_stream_push_token(s, tok);
	return -1;
}

void isl_stream_free(struct isl_stream *s)
{
	if (!s)
		return;
	free(s->buffer);
	if (s->n_token != 0) {
		struct isl_token *tok = isl_stream_next_token(s);
		isl_stream_error(s, tok, "unexpected token");
		isl_token_free(tok);
	}
	isl_ctx_deref(s->ctx);
	free(s);
}

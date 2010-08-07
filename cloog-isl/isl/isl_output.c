#include <isl_set.h>

static void print_constraint_polylib(struct isl_basic_set *bset,
	int ineq, int n,
	FILE *out, int indent, const char *prefix, const char *suffix)
{
	int i;
	unsigned dim = isl_basic_set_n_dim(bset);
	unsigned nparam = isl_basic_set_n_param(bset);
	isl_int *c = ineq ? bset->ineq[n] : bset->eq[n];

	fprintf(out, "%*s%s", indent, "", prefix ? prefix : "");
	fprintf(out, "%d", ineq);
	for (i = 0; i < dim; ++i) {
		fprintf(out, " ");
		isl_int_print(out, c[1+nparam+i], 5);
	}
	for (i = 0; i < bset->n_div; ++i) {
		fprintf(out, " ");
		isl_int_print(out, c[1+nparam+dim+i], 5);
	}
	for (i = 0; i < nparam; ++i) {
		fprintf(out, " ");
		isl_int_print(out, c[1+i], 5);
	}
	fprintf(out, " ");
	isl_int_print(out, c[0], 5);
	fprintf(out, "%s\n", suffix ? suffix : "");
}

static void print_constraints_polylib(struct isl_basic_set *bset,
	FILE *out, int indent, const char *prefix, const char *suffix)
{
	int i;

	for (i = 0; i < bset->n_eq; ++i)
		print_constraint_polylib(bset, 0, i, out,
					indent, prefix, suffix);
	for (i = 0; i < bset->n_ineq; ++i)
		print_constraint_polylib(bset, 1, i, out,
					indent, prefix, suffix);
}

static void isl_basic_set_print_polylib(struct isl_basic_set *bset, FILE *out,
	int indent, const char *prefix, const char *suffix)
{
	unsigned total = isl_basic_set_total_dim(bset);
	fprintf(out, "%*s%s", indent, "", prefix ? prefix : "");
	fprintf(out, "%d %d", bset->n_eq + bset->n_ineq, 1 + total + 1);
	fprintf(out, "%s\n", suffix ? suffix : "");
	print_constraints_polylib(bset, out, indent, prefix, suffix);
}

static void isl_set_print_polylib(struct isl_set *set, FILE *out, int indent)
{
	int i;

	fprintf(out, "%*s", indent, "");
	fprintf(out, "%d\n", set->n);
	for (i = 0; i < set->n; ++i) {
		fprintf(out, "\n");
		isl_basic_set_print_polylib(set->p[i], out, indent, NULL, NULL);
	}
}

void isl_basic_set_print(struct isl_basic_set *bset, FILE *out, int indent,
	const char *prefix, const char *suffix, unsigned output_format)
{
	if (!bset)
		return;
	if (output_format == ISL_FORMAT_POLYLIB)
		isl_basic_set_print_polylib(bset, out, indent, prefix, suffix);
	else if (output_format == ISL_FORMAT_POLYLIB_CONSTRAINTS)
		print_constraints_polylib(bset, out, indent, prefix, suffix);
	else
		isl_assert(bset->ctx, 0, return);
}

void isl_set_print(struct isl_set *set, FILE *out, int indent,
	unsigned output_format)
{
	if (!set)
		return;
	if (output_format == ISL_FORMAT_POLYLIB)
		isl_set_print_polylib(set, out, indent);
	else
		isl_assert(set->ctx, 0, return);
}

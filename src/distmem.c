/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution. 
 *
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "pluto.h"
#include "program.h"
#include "transforms.h"


#define __APPENDFILENAME "packunpack.c"

void get_data_tag(PlutoProg *prog, PlutoAccess *acc, char **data_tag)
{
	int i;
	char *acc_name = acc->name;
	for (i=0; i<prog->num_data; i++) {
		if (!strcmp(prog->data_names[i], acc_name)) {
			sprintf(*data_tag, "%d123", i+1);
			break;
		}
	}
}


/* Get list of <Stmt, acc> lists - each list corresponds to
 * statements with accesses to the same data; *num_data gives the number of
 * data elements; nstmts_per_acc[i] is the number of stmts in list[i]
 * Consider statements passed via 'stmts'
 * Return: num_data - number of variables, i.e., number of lists
 * nstmts_per_acc[]: number of stmt, acc pairs in each list
 * */
struct stmt_access_pair ***get_read_write_access_with_stmts(Stmt **stmts,
		int nstmts, int *num_data,
		int **nstmts_per_acc)
{
	int i, j, k, curr_num;
	int *num_stmts_per_acc;
	struct stmt_access_pair ***acc_stmts;

	curr_num = 0;
	num_stmts_per_acc = NULL;
	acc_stmts = NULL;

	for (i=0; i<nstmts; i++)  {
		Stmt *stmt = stmts[i];
		struct stmt_access_pair *new;
		for (k=0; k<stmt->nreads; k++){

			if(is_access_scalar(stmt->reads[k]))
				continue;

			new = malloc(sizeof(struct stmt_access_pair));
			new->stmt = stmt;
			new->acc = stmt->reads[k];

			for (j=0; j<curr_num; j++)  {
				if (!strcmp(stmt->reads[k]->name, acc_stmts[j][0]->acc->name)) {
					/* Add to end of array */
					acc_stmts[j] = (struct stmt_access_pair **) realloc(acc_stmts[j],
							(num_stmts_per_acc[j]+1)*sizeof(struct stmt_access_pair*));
					acc_stmts[j][num_stmts_per_acc[j]] = new;

					num_stmts_per_acc[j]++;
					break;
				}
			}
			if (j==curr_num)    {
				/* New data variable */
				acc_stmts = (struct stmt_access_pair ***)realloc(acc_stmts,
						(curr_num+1)*sizeof(struct stmt_access_pair **));
				acc_stmts[curr_num] = (struct stmt_access_pair **) malloc(
						sizeof(struct stmt_access_pair*));
				acc_stmts[curr_num][0] = new;

				num_stmts_per_acc = (int *)realloc(num_stmts_per_acc,
						(curr_num+1)*sizeof(int));
				num_stmts_per_acc[curr_num] = 1;
				curr_num++;
			}
		}

		stmt = stmts[i];
		new = malloc(sizeof(struct stmt_access_pair));
		new->stmt = stmt;
		new->acc = stmt->writes[0];

		for (j=0; j<curr_num; j++)  {
			if (!strcmp(stmt->writes[0]->name, acc_stmts[j][0]->acc->name)) {
				/* Add to end of array */
				acc_stmts[j] = (struct stmt_access_pair **) realloc(acc_stmts[j],
						(num_stmts_per_acc[j]+1)*sizeof(struct stmt_access_pair*));
				acc_stmts[j][num_stmts_per_acc[j]] = new;

				num_stmts_per_acc[j]++;
				break;
			}
		}
		if (j==curr_num)    {
			/* New data variable */
			acc_stmts = (struct stmt_access_pair ***)realloc(acc_stmts,
					(curr_num+1)*sizeof(struct stmt_access_pair **));
			acc_stmts[curr_num] = (struct stmt_access_pair **) malloc(
					sizeof(struct stmt_access_pair*));
			acc_stmts[curr_num][0] = new;

			num_stmts_per_acc = (int *)realloc(num_stmts_per_acc,
					(curr_num+1)*sizeof(int));
			num_stmts_per_acc[curr_num] = 1;
			curr_num++;
		}

	}

	*num_data = curr_num;
	*nstmts_per_acc = num_stmts_per_acc;

	return acc_stmts;
}

/* Get list of <Stmt, read acc> lists - each list corresponds to 
 * statements with read accesses to the same data; *num_data gives the number of 
 * data elements; nstmts_per_acc[i] is the number of stmts in list[i] 
 * Consider statements passed via 'stmts'
 * Return: num_data - number of variables, i.e., number of lists
 * nstmts_per_acc[]: number of stmt, acc pairs in each list
 * */
// not being used currently
struct stmt_access_pair ***get_read_access_with_stmts(Stmt **stmts, 
		int nstmts, int *num_data,
		int **nstmts_per_acc)
{
	int i, j, k, curr_num;
	int *num_stmts_per_acc;
	struct stmt_access_pair ***racc_stmts;

	curr_num = 0;
	num_stmts_per_acc = NULL;
	racc_stmts = NULL;

	for (i=0; i<nstmts; i++)  {
		Stmt *stmt = stmts[i];
		for (k=0; k<stmt->nreads; k++){

			if(is_access_scalar(stmt->reads[k]))
				continue;

			struct stmt_access_pair *new = malloc(sizeof(struct stmt_access_pair));
			new->stmt = stmt;
			new->acc = stmt->reads[k];

			for (j=0; j<curr_num; j++)  {
				if (!strcmp(stmt->reads[k]->name, racc_stmts[j][0]->acc->name)) {
					/* Add to end of array */
					racc_stmts[j] = (struct stmt_access_pair **) realloc(racc_stmts[j],
							(num_stmts_per_acc[j]+1)*sizeof(struct stmt_access_pair*));
					racc_stmts[j][num_stmts_per_acc[j]] = new;

					num_stmts_per_acc[j]++;
					break;
				}
			}
			if (j==curr_num)    {
				/* New data variable */
				racc_stmts = (struct stmt_access_pair ***)realloc(racc_stmts,
						(curr_num+1)*sizeof(struct stmt_access_pair **));
				racc_stmts[curr_num] = (struct stmt_access_pair **) malloc(
						sizeof(struct stmt_access_pair*));
				racc_stmts[curr_num][0] = new;

				num_stmts_per_acc = (int *)realloc(num_stmts_per_acc,
						(curr_num+1)*sizeof(int));
				num_stmts_per_acc[curr_num] = 1;
				curr_num++;
			}
		}
	}

	*num_data = curr_num;
	*nstmts_per_acc = num_stmts_per_acc;

	return racc_stmts;
}


/* Get list of <Stmt, write acc> lists - each list corresponds to 
 * statements with write accesses to the same data; *num_data gives the number of 
 * data elements; nstmts_per_acc[i] is the number of stmts in list[i] 
 * Consider statements passed via 'stmts'
 * Return: num_data - number of variables, i.e., number of lists
 * nstmts_per_acc[]: number of stmt, acc pairs in each list
 * */
struct stmt_access_pair ***get_write_access_with_stmts(Stmt **stmts, 
		int nstmts, int *num_data,
		int **nstmts_per_acc)
{
	int i, j, k, curr_num;
	int *num_stmts_per_acc;
	struct stmt_access_pair ***wacc_stmts;

	curr_num = 0;
	num_stmts_per_acc = NULL;
	wacc_stmts = NULL;

	for (i=0; i<nstmts; i++) {
		Stmt *stmt = stmts[i];
		for (k=0; k<stmt->nwrites; k++) {
			struct stmt_access_pair *new = malloc(sizeof(struct stmt_access_pair));
			new->stmt = stmt;
			new->acc = stmt->writes[k];

			for (j=0; j<curr_num; j++)  {
				if (!strcmp(stmt->writes[k]->name, wacc_stmts[j][0]->acc->name)) {
					/* Add to end of array */
					wacc_stmts[j] = (struct stmt_access_pair **) realloc(wacc_stmts[j],
							(num_stmts_per_acc[j]+1)*sizeof(struct stmt_access_pair*));
					wacc_stmts[j][num_stmts_per_acc[j]] = new;

					num_stmts_per_acc[j]++;
					break;
				}
			}
			if (j==curr_num)    {
				/* New data variable */
				wacc_stmts = (struct stmt_access_pair ***)realloc(wacc_stmts,
						(curr_num+1)*sizeof(struct stmt_access_pair **));
				wacc_stmts[curr_num] = (struct stmt_access_pair **) malloc(
						sizeof(struct stmt_access_pair*));
				wacc_stmts[curr_num][0] = new;

				num_stmts_per_acc = (int *)realloc(num_stmts_per_acc,
						(curr_num+1)*sizeof(int));
				num_stmts_per_acc[curr_num] = 1;
				curr_num++;
			}
		}
	}

	*num_data = curr_num;
	*nstmts_per_acc = num_stmts_per_acc;

	return wacc_stmts;
}


static void pluto_mark_statements(PlutoProg *prog)
{
	int i;

	if (options->distmem) {
		pluto_prog_add_param(prog, "my_rank", prog->npar);
		pluto_constraints_add_lb(prog->context, prog->npar-1, 0);
	}

	for (i=0; i<prog->nstmts; i++) {
		Stmt *stmt = prog->stmts[i];
		if (stmt->type == LW_COPY_IN) {
			/* Add my_rank == 0; write-out set is copied-back only by master proc */
			pluto_constraints_set_var(stmt->domain, stmt->domain->ncols-2, 0);
		}
	}
}


/* Get access string */
char *reconstruct_access(PlutoAccess *acc)
{
	int ndims, i;
	ndims = acc->mat->nrows;

	char *access;

	access = (char *) malloc(strlen(acc->name)+ndims*(2+4)+1);

	strcpy(access, acc->name);

	for (i=0; i<ndims; i++) {
		strcat(access, "[");
		char tmp[5];
		sprintf(tmp, "d%d", i+1);
		strcat(access, tmp);
		strcat(access, "]");
	}

	return access;
}


void free_char_array_buffers(char** buffer, int size) {

	int i = 0;
	for(i = 0; i < size; ++i) {
		free(buffer[i]);
	}

	free(buffer);

	return;
}

Stmt *get_new_anchor_stmt(Stmt **stmts, int nstmts)
{
	int i, j;

	int max_dim_index, max_dim;
	max_dim = stmts[0]->dim;
	max_dim_index = 0;
	for (i=1; i<nstmts; i++) {
		if (stmts[i]->dim > max_dim) {
			max_dim = stmts[i]->dim;
			max_dim_index = i;
		}
	}

	Stmt *anchor_stmt = stmts[max_dim_index];
	PlutoConstraints *domain = pluto_constraints_dup(anchor_stmt->domain);
	for (i=0; i<nstmts; i++) {
		PlutoConstraints *new_domain = pluto_constraints_dup(stmts[i]->domain);
		int ncols = new_domain->ncols;
		for (j=0; j<domain->ncols-ncols; j++) {
			pluto_constraints_add_dim(new_domain, stmts[i]->dim, NULL);
		}
		pluto_constraints_unionize(domain, new_domain);
		pluto_constraints_free(new_domain);
	}

	Stmt *new_anchor_stmt = pluto_create_stmt(anchor_stmt->dim, domain, anchor_stmt->trans,
			anchor_stmt->iterators, anchor_stmt->text, anchor_stmt->type);
	new_anchor_stmt->id = anchor_stmt->id;
	pluto_constraints_free(domain);
	return new_anchor_stmt;
}

char *pluto_dist_gen_declarations(char *arr_name, PlutoProg *prog){


	Array* arr = pluto_get_corrs_array(arr_name, prog);
	assert(arr != NULL);

	PlutoConstraints *data_tiles = pluto_dist_get_required_data_tiles(arr->array_bounds,
			0, arr_name, prog);

	/* project out copy level dims
	 */
	//	pluto_constraints_project_out(data_tiles, 0, arr->copy_level);
	print_polylib_visual_sets("set", data_tiles);

	/* Number of points in the projected out domain gives the
	 * total number of data tiles
	 */
	assert(data_tiles->ncols ==  arr->num_tiled_loops + prog->npar + 1);

	char* decl = (char *) malloc(1024*50* sizeof(char));
	decl[0] = 0;

	char *index = (char *)malloc(1024*5* sizeof(char));
	char *size_text = (char *)malloc(1024*5 *sizeof(char));
	size_text[0] = '\0';
	index[0] = '\0';
	int i;

	for (i = 0; i < arr->num_tiled_loops; ++i) {

		char *buf_size =
				get_parametric_bounding_box(data_tiles,  i,
						1, prog->npar, (const char **)prog->params);
		//
		//			get_parametric_bounding_box(data_tiles,arr->first_tile_dim - arr->copy_level + i,
		//				1, prog->npar, (const char **)prog->params);

		sprintf(decl+strlen(decl), "int %s_size_%d = %s;\n" , arr->text, arr->first_tile_dim + i, buf_size);
		free(buf_size);
	}

	for (i = 0; i < arr->num_tiled_loops; ++i) {
		sprintf(index + strlen(index),"[%s_size_%d]",arr->text,arr->first_tile_dim + i);

		sprintf(size_text + strlen(size_text), "%s_size_%d", arr->text, arr->first_tile_dim + i);

		if(i != arr->num_tiled_loops -1)
			sprintf(size_text + strlen(size_text), " * ");
	}

	sprintf(decl+strlen(decl), "%s **", arr->data_type);
	sprintf(decl+strlen(decl), "%s_trans = (%s **)malloc((%s) * sizeof(%s));\n", arr->text,
			arr->data_type, size_text, arr->data_type);

	sprintf(decl+strlen(decl), "BufferManager *buff_mang_%s = ", arr->text);
	sprintf(decl+strlen(decl), "init_buffer_manager(%s);\n", size_text);
	//
	//	sprintf(decl+strlen(decl), "int *");
	//	sprintf(decl+strlen(decl), "%s_ref = (int *)malloc((%s) * sizeof(int));\n", arr->text, size_text);
	//
	//	sprintf(decl+strlen(decl), "BufferList *list_%s = NULL;\n", arr->text);

	pluto_constraints_free(data_tiles);

	return decl;
}

//TODO: Update arr iterators
void pluto_dist_update_arr_domain(struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num){

	int i, src_copy_level;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	//    char *acc_name = anchor_stmt->reads[0]->name;
	char *acc_name = acc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	PlutoConstraints *all_data= NULL;
	for (i=0; i<num_accs; i++)  {
		Stmt *stmt = acc_stmts[i]->stmt;
		PlutoAccess *acc = acc_stmts[i]->acc;

		PlutoConstraints *srcdomain = pluto_get_new_domain(stmt);

		PlutoConstraints *data_in_one =  pluto_compute_region_data(stmt,
				srcdomain, acc, src_copy_level, prog);

		if (all_data == NULL) all_data = pluto_constraints_dup(data_in_one);
		else all_data = pluto_constraints_unionize(all_data, data_in_one);
		pluto_constraints_free(data_in_one);
	}

	IF_DEBUG(printf("all data accessed set\n"););
	print_polylib_visual_sets("all_data", all_data);

	Array *arr = pluto_get_corrs_array(acc_name, prog);

	assert(arr->parmetric_domain != NULL);
	assert(arr->copy_level!= NULL);
	arr->parmetric_domain[loop_num] = pluto_constraints_dup(all_data);
	arr->copy_level[loop_num] = src_copy_level;


	pluto_constraints_project_out(all_data, 0, src_copy_level);

	if(arr->array_bounds == NULL){
		arr->array_bounds = pluto_constraints_dup(all_data);
	}
	else
		pluto_constraints_unionize(arr->array_bounds, all_data);

	pluto_constraints_free(all_data);

	return;
}

char *get_data_tile_ref_count_stmt_text(Array *arr, PlutoProg *prog){

	char *stmt_text = (char *)malloc(1024 * sizeof(char));
	stmt_text[0] = '\0';

	sprintf(stmt_text+ strlen(stmt_text), "decrement_ref_count(buff_mang_%s, ", arr->text);
	pluto_dist_print_tiles_index_string(stmt_text, arr, prog);
	sprintf(stmt_text+ strlen(stmt_text), ");");

	return stmt_text;
}

/*
 * Code generation for the data tile ref count decrement stmts
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * acc_stmts:  all <statements, acc> writing to this data variable
 * This function is called per data variable
 * This function can be called by any flow-based communication scheme
 * caller is responsible for freeing read_in_stmts returned
 */
Stmt **gen_data_tile_ref_count_update_code(struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num)
{
	int j, src_copy_level;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	Stmt *anchor_stmt = acc_stmts[0]->stmt;
	//    char *acc_name = anchor_stmt->reads[0]->name;
	char *acc_name = acc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);
	assert(arr->parmetric_domain[loop_num]!=NULL);

	/* Tile the data set to compute the required data tiles to be allocated
	 */
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(arr->parmetric_domain[loop_num], src_copy_level, acc_name,
					prog);

	print_polylib_visual_sets("datatile", data_tiles);
	// ===== from here on starts the code generation ===============================

	char *stmt_text = get_data_tile_ref_count_stmt_text(arr, prog);

	Stmt *alloc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			stmt_text, DATA_DIST_MANG, loop_num);

	/* Add dimensions that actually scan the data tile space */
	for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
		pluto_stmt_add_dim(alloc_stmt, alloc_stmt->dim,
				alloc_stmt->trans->nrows, arr->iterators[j], H_LOOP, prog);
	}

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	//    pluto_constraints_intersect(alloc_stmt->domain, data_tiles);
	alloc_stmt->domain = data_tiles;

	print_polylib_visual_sets("alloc_stmt", alloc_stmt->domain);
	Stmt **data_alloc_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	data_alloc_stmts[0] = alloc_stmt;

	free(stmt_text);

	return data_alloc_stmts;
}


char* get_tile_pi_check_ref_count_stmt_text(Array *arr, int src_copy_level, int loop_num,
		PlutoProg *prog){

	int i;
	char *stmt_text = (char *)malloc(1024 * sizeof(char));
	stmt_text[0] = '\0';

	if(options->distmem){
		sprintf(stmt_text , "if(pi_%d(",loop_num);

		for (i=0; i<src_copy_level; i++) {
			sprintf(stmt_text+strlen(stmt_text), "t%d,",i+1);
		}
		for (i = 0; i < prog->npar; ++i) {
			sprintf(stmt_text+strlen(stmt_text), "%s,",prog->params[i]);
		}

		sprintf(stmt_text+strlen(stmt_text), "nprocs) == my_rank)\t");
	}

	sprintf(stmt_text+strlen(stmt_text), "data_tile_ref_count_init_%s_%d(", arr->text, loop_num );
	for (i=0; i<src_copy_level; i++) {
		sprintf(stmt_text+strlen(stmt_text), "t%d",i+1);
		if(i!=src_copy_level - 1)
			sprintf(stmt_text+strlen(stmt_text), ", ");
	}

	print_data_dist_parm_call(stmt_text, arr);
	sprintf(stmt_text+strlen(stmt_text), ");");

	return stmt_text;
}

//generate code for write out data retention
Stmt **gen_write_out_tiles_ref_count_code (struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num)
{
	int i, src_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	acc_nrows = acc_stmts[0]->acc->mat->nrows;
	Stmt *anchor_stmt = acc_stmts[0]->stmt;
	//    char *acc_name = anchor_stmt->reads[0]->name;
	char *acc_name = acc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);

	PlutoConstraints *write_out = NULL;
	for (i=0; i<num_accs; i++)  {
		PlutoConstraints *write_out_one;
		//write_out_one = compute_last_writes(wacc_stmts[k], src_copy_level, prog);
		write_out_one = compute_write_out(acc_stmts[i], src_copy_level, prog);
		if (write_out == NULL) write_out = pluto_constraints_dup(write_out_one);
		else write_out = pluto_constraints_unionize(write_out, write_out_one);
		pluto_constraints_free(write_out_one);
	}

	IF_DEBUG(printf("Write out set\n"););
	IF_DEBUG(pluto_constraints_print(stdout, write_out));

	assert(arr->parmetric_domain[loop_num]!=NULL);
	/* Tile the data set to compute the required data tiles to be allocated
	 */
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(arr->parmetric_domain[loop_num], src_copy_level, acc_name, prog);


	assert(data_tiles->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	pluto_constraints_project_out(data_tiles, src_copy_level, acc_nrows);

	char *stmt_text = get_tile_pi_check_ref_count_stmt_text(arr, src_copy_level, loop_num, prog);

	Stmt *alloc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			stmt_text, DATA_DIST_INIT, loop_num);

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	pluto_constraints_intersect(alloc_stmt->domain, data_tiles);

	Stmt **data_alloc_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	data_alloc_stmts[0] = alloc_stmt;

	free(stmt_text);

	return data_alloc_stmts;
}

char* get_data_tile_ref_count_init_stmt_text(Array *arr, int src_copy_level, int loop_num,
		PlutoProg *prog){

	char *stmt_text = (char *)malloc(1024 * sizeof(char));
	stmt_text[0] = '\0';

	sprintf(stmt_text +strlen(stmt_text), "increment_ref_count(buff_mang_%s, ",arr->text);

	pluto_dist_print_tiles_index_string(stmt_text, arr, prog);

	sprintf(stmt_text+ strlen(stmt_text), ");");

	return stmt_text;
}

/*
 * Code generation for the total data tile ref counts
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * acc_stmts:  all <statements, acc> writing to this data variable
 * This function is called per data variable
 * This function can be called by any flow-based communication scheme
 * caller is responsible for freeing read_in_stmts returned
 */
Stmt **gen_data_tile_ref_count_code (struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num)
{
	int  src_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	acc_nrows = acc_stmts[0]->acc->mat->nrows;
	Stmt *anchor_stmt = acc_stmts[0]->stmt;
	//    char *acc_name = anchor_stmt->reads[0]->name;
	char *acc_name = acc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);
	assert(arr->parmetric_domain[loop_num] != NULL);
	assert(arr->copy_level[loop_num] == src_copy_level);

	/* Tile the data set to compute the required data tiles to be allocated
	 */
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(arr->parmetric_domain[loop_num], src_copy_level, acc_name, prog);


	assert(data_tiles->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	pluto_constraints_project_out(data_tiles, src_copy_level, acc_nrows);

	char *stmt_text = get_tile_pi_check_ref_count_stmt_text(arr, src_copy_level, loop_num, prog);

	Stmt *alloc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			stmt_text, DATA_DIST_INIT, loop_num);

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	pluto_constraints_intersect(alloc_stmt->domain, data_tiles);

	Stmt **data_alloc_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	data_alloc_stmts[0] = alloc_stmt;

	free(stmt_text);

	return data_alloc_stmts;
}

/*
 * Optimized communication code generation for the
 * data tile alloc for main kernal
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * acc_stmts:  all <statements, acc> writing to this data variable
 * This function is called per data variable
 * This function can be called by any flow-based communication scheme
 * caller is responsible for freeing read_in_stmts returned
 */
Stmt **gen_data_tile_alloc_code (struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num)
{
	int j, src_copy_level;
	src_copy_level = copy_level[loop_num];


	assert(num_accs >= 1);
	Stmt *anchor_stmt = acc_stmts[0]->stmt;
	//    char *acc_name = anchor_stmt->reads[0]->name;
	char *acc_name = acc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);
	assert(arr->parmetric_domain[loop_num]!=NULL);

	/* Tile the data set to compute the required data tiles to be allocated
	 */
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(arr->parmetric_domain[loop_num], src_copy_level, acc_name, prog);

	// ===== from here on starts the code generation ===============================

	char* alloc_stmt_text = pluto_dist_malloc_stmt_text(acc_name, prog, 1);

	Stmt *alloc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			alloc_stmt_text, DATA_DIST_MANG, loop_num);

	/* Add dimensions that actually scan the data tile space */
	for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
		pluto_stmt_add_dim(alloc_stmt, alloc_stmt->dim,
				alloc_stmt->trans->nrows, arr->iterators[j], H_LOOP, prog);
	}

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	pluto_constraints_intersect(alloc_stmt->domain, data_tiles);

	Stmt **data_alloc_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	data_alloc_stmts[0] = alloc_stmt;

	free(alloc_stmt_text);

	return data_alloc_stmts;
}

Stmt *gen_comm_data_alloc_code(Stmt *anchor_stmt, PlutoConstraints *constraints,int src_copy_level, char *acc_name,
		PlutoProg *prog, int loop_num){

	int j =0;
	Array *arr = pluto_get_corrs_array(acc_name, prog);

	assert(arr!= NULL);

	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(constraints, src_copy_level, acc_name,
					prog);

	char **new_iter = (char **)malloc(arr->num_tiled_loops * sizeof(char*));
	int count = 0;
	for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
		new_iter[count++] = strdup(arr->iterators[j]);
	}
	char* alloc_stmt_text = pluto_dist_malloc_stmt_text(acc_name, prog, 1);

	Stmt *alloc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			alloc_stmt_text, DATA_DIST_MANG , loop_num);

	/* Add dimensions that actually scan the data tile space */
	for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
		pluto_stmt_add_dim(alloc_stmt, alloc_stmt->dim,
				alloc_stmt->trans->nrows, arr->iterators[j], H_LOOP, prog);
	}

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	//	alloc_stmt->domain = pluto_constraints_dup(data_tiles);
	pluto_constraints_intersect(alloc_stmt->domain, data_tiles);

	return alloc_stmt;

}

char *pluto_dist_copy_back_stmt_text(Array *arr, PlutoProg *prog){

	char *text = (char *)malloc(1024 * sizeof(char));
	text[0] = '\0';

	sprintf(text+strlen(text), "%s", arr->text);
	int i = 0;

	for (i = 0; i < arr->dim_orig; ++i) {
		sprintf(text+strlen(text), "[%s]", arr->iterators[i]);
	}

	sprintf(text+strlen(text), ";");

	char * mod_text = pluto_dist_modify_stmt_text(text, 1, prog);

	text[0] = '\0';

	sprintf(text+strlen(text), "%s", arr->text);
	for (i = 0; i < arr->dim_orig; ++i) {
		sprintf(text+strlen(text), "[%s]", arr->iterators[i]);
	}

	sprintf(text + strlen(text), " = %s",mod_text );
	free(mod_text);

	return text;
}

char *pluto_dist_copy_back_func_call_text(int loop_num, int src_copy_level, Array *arr, PlutoProg *prog){
	char *text = (char *)malloc(1024 * sizeof(char));
	text[0] = '\0';

	int i =0;
	if (options->distmem) {
		sprintf(text , "if(pi_%d(",loop_num);

		for (i=0; i<src_copy_level; i++) {
			sprintf(text+strlen(text), "t%d,",i+1);
		}
		for (i = 0; i < prog->npar; ++i) {
			sprintf(text+strlen(text), "%s,",prog->params[i]);
		}

		sprintf(text+strlen(text), "nprocs) == my_rank)\t");
	}

	sprintf(text+strlen(text), "data_tile_copy_back_%s_%d(", arr->text, loop_num);
	for (i=0; i<src_copy_level; i++) {
		sprintf(text+strlen(text), "t%d",i+1);
		if(i != src_copy_level -1)
			sprintf(text+strlen(text), ",");
	}

	if(options->variables_not_global)
		sprintf(text+strlen(text), ", %s", arr->text);

	print_data_dist_parm_call(text, arr);

	sprintf(text+strlen(text), ");");

	return text;
}

Stmt **gen_copy_back_code(struct stmt_access_pair **racc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num)
{
	int i, src_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	acc_nrows = racc_stmts[0]->acc->mat->nrows;
	Stmt *anchor_stmt = racc_stmts[0]->stmt;
	//    char *acc_name = anchor_stmt->reads[0]->name;
	char *acc_name = racc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	PlutoConstraints *all_data= NULL;
	for (i=0; i<num_accs; i++)  {
		Stmt *stmt = racc_stmts[i]->stmt;
		PlutoAccess *acc = racc_stmts[i]->acc;

		PlutoConstraints *srcdomain = pluto_get_new_domain(stmt);

		PlutoConstraints *data_in_one =  pluto_compute_region_data(stmt,
				srcdomain, acc, src_copy_level, prog);

		if (all_data == NULL) all_data = pluto_constraints_dup(data_in_one);
		else all_data = pluto_constraints_unionize(all_data, data_in_one);
		pluto_constraints_free(data_in_one);
	}

	pluto_constraints_project_out(all_data, src_copy_level, acc_nrows);

	Array *arr = pluto_get_corrs_array(acc_name, prog);

	// ===== from here on starts the code generation =====================================================

	char *copy_back_stmt_text = pluto_dist_copy_back_func_call_text(loop_num, src_copy_level, arr, prog);

	Stmt *copy_back_stmt = create_helper_stmt(anchor_stmt, src_copy_level, copy_back_stmt_text,
			DATA_DIST_COPY , loop_num);

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	pluto_constraints_intersect(copy_back_stmt->domain, all_data);

	Stmt **copy_back_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	copy_back_stmts[0] = copy_back_stmt;

	free(copy_back_stmt_text);

	return copy_back_stmts;
}


Stmt **gen_null_init_code(Array *arr, PlutoProg *prog, int loop_num)
{
	int j;

	char *acc_name = arr->text;

	assert(arr->array_bounds != NULL);

	PlutoConstraints *all_tiles =
			pluto_dist_get_required_data_tiles(arr->array_bounds, 0, acc_name, prog);

	// ===== from here on starts the code generation =====================================================

	char *ptr_init_stmt_text = pluto_dist_ptr_init_stmt_text(acc_name, prog);

	Stmt *ptr_init_stmt = create_helper_stmt(prog->stmts[0], 0, ptr_init_stmt_text, DATA_SETUP, loop_num);

	/* Add dimensions that actually scan the data tile space */
	for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
		pluto_stmt_add_dim(ptr_init_stmt, ptr_init_stmt->dim,
				ptr_init_stmt->trans->nrows, arr->iterators[j], H_LOOP, prog);
	}

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	//   pluto_constraints_intersect(ptr_init_stmt->domain, all_tiles);
	ptr_init_stmt->domain = pluto_constraints_dup(all_tiles);

	Stmt **null_init_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	null_init_stmts[0] = ptr_init_stmt;

	free(ptr_init_stmt_text);

	return null_init_stmts;
}

int starts_with(const char *pre, const char *str)
{
	int lenpre = strlen(pre),
			lenstr = strlen(str);
	return lenstr < lenpre ? 0 : strncmp(pre, str, lenpre) == 0;
}


char *pluto_dist_copy_init_text(Array *arr, PlutoProg *prog){
	char *text = (char *)malloc(1024 * 5 * sizeof(char));
	text[0] = '\0';

	sprintf(text+strlen(text), "%s", arr->text);
	int i = 0;

	for (i = 0; i < arr->dim_orig; ++i) {
		sprintf(text+strlen(text), "[d%d]", i+1);
	}

	sprintf(text+strlen(text), " = ");

	char *mod_text = pluto_dist_modify_stmt_text(text, 1, prog);
	free(text);

	sprintf(mod_text+strlen(mod_text), "%s", arr->text);
	for (i = 0; i < arr->dim_orig; ++i) {
		sprintf(mod_text+strlen(mod_text), "[d%d]", i+1);
	}

	sprintf(mod_text + strlen(mod_text), ";");

	return mod_text;

}
char *pluto_dist_read_init_text(char *arr_name, PlutoProg *prog){
	FILE *fp = NULL;

	size_t linesiz=0;
	char* linebuf = (char *)malloc(1024 * 5 * sizeof(char));
	linebuf[0] = '\0';
	ssize_t linelen=0;

	fp = fopen ("array_init", "rt");
	if(fp == NULL) return linebuf;

	while ((linelen=getline(&linebuf, &linesiz, fp)>0)) {

		if(starts_with(arr_name, linebuf)){
			//remove the trailing new line
			linebuf[strlen(linebuf) - 1] = '\0';
			return linebuf;
		}
	}

	fclose(fp);

	char *mod_text = pluto_dist_modify_stmt_text(linebuf, 1, prog);
	free(linebuf);

	return mod_text;

}

char* get_tile_pi_check_read_in_text(Array *arr, int src_copy_level, int loop_num,
		PlutoProg *prog){

	int i;
	char *stmt_text = (char *)malloc(1024 * sizeof(char));
	stmt_text[0] = '\0';

	if(options->distmem){
		sprintf(stmt_text , "if(pi_%d(",loop_num);

		for (i=0; i<src_copy_level; i++) {
			sprintf(stmt_text+strlen(stmt_text), "t%d,",i+1);
		}
		for (i = 0; i < prog->npar; ++i) {
			sprintf(stmt_text+strlen(stmt_text), "%s,",prog->params[i]);
		}

		sprintf(stmt_text+strlen(stmt_text), "nprocs) == my_rank){\t");
	}

	sprintf(stmt_text+strlen(stmt_text), "read_in_data_tile_alloc_%s_%d(", arr->text, loop_num );
	for (i=0; i<src_copy_level; i++) {
		sprintf(stmt_text+strlen(stmt_text), "t%d",i+1);
		if(i!=src_copy_level - 1)
			sprintf(stmt_text+strlen(stmt_text), ", ");
	}

	print_data_dist_parm_call(stmt_text, arr);
	sprintf(stmt_text+strlen(stmt_text), ");");

	sprintf(stmt_text+strlen(stmt_text), "read_in_copy_%s_%d(", arr->text, loop_num);


	for (i=0; i<src_copy_level; i++) {
		sprintf(stmt_text+strlen(stmt_text), "t%d",i+1);
		if(i!=src_copy_level - 1)
			sprintf(stmt_text+strlen(stmt_text), ", ");
	}

	if(options->variables_not_global){
		sprintf(stmt_text+strlen(stmt_text), " __ifndef USE_LOCAL_ARRAYS");
		sprintf(stmt_text+strlen(stmt_text), ", %s", arr->text);
		sprintf(stmt_text+strlen(stmt_text), " __endif");
	}

	print_data_dist_parm_call(stmt_text, arr);
	sprintf(stmt_text+strlen(stmt_text), ");");
	if(options->distmem)
		sprintf(stmt_text+strlen(stmt_text), "}");
	return stmt_text;
}

Stmt **gen_read_in_code (PlutoConstraints *read_in, struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num)
{
	int src_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	acc_nrows = acc_stmts[0]->acc->mat->nrows;
	Stmt *anchor_stmt = acc_stmts[0]->stmt;
	char *acc_name = acc_stmts[0]->acc->name;

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);
	assert(arr->parmetric_domain[loop_num] != NULL);
	assert(arr->copy_level[loop_num] == src_copy_level);

	/* Tile the data set to compute the required data tiles to be allocated
	 */
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(read_in, src_copy_level, acc_name, prog);

	assert(data_tiles->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	pluto_constraints_project_out(data_tiles, src_copy_level, acc_nrows);

	char *stmt_text = get_tile_pi_check_read_in_text(arr, src_copy_level, loop_num, prog);

	Stmt *alloc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			stmt_text, DATA_DIST_INIT, loop_num);

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	pluto_constraints_intersect(alloc_stmt->domain, data_tiles);

	Stmt **read_in_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	read_in_stmts[0] = alloc_stmt;

	free(stmt_text);

	return read_in_stmts;
}

/*
 * Optimized communication code generation for the read-in set
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * racc_stmts:  all <statements, racc> writing to this data variable
 * This function is called per data variable
 * This function can be called by any flow-based communication scheme
 * caller is responsible for freeing read_in_stmts returned
 */
Stmt **gen_read_in_code_old(struct stmt_access_pair **racc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int loop_num)
{
	int i,j, src_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	acc_nrows = racc_stmts[0]->acc->mat->nrows;
	Stmt *anchor_stmt = racc_stmts[0]->stmt;
	//    char *acc_name = anchor_stmt->reads[0]->name;
	char *acc_name = racc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	PlutoConstraints *read_in = NULL;
	for (i=0; i<num_accs; i++)  {
		PlutoConstraints *read_in_one;
		//write_out_one = compute_last_writes(wacc_stmts[k], src_copy_level, prog);
		read_in_one = compute_read_in(racc_stmts[i], src_copy_level, prog);
		if (read_in == NULL) read_in = pluto_constraints_dup(read_in_one);
		else read_in = pluto_constraints_unionize(read_in , read_in_one);
		pluto_constraints_free(read_in_one);
	}

	IF_DEBUG(printf("read in set\n"););
	IF_DEBUG(pluto_constraints_print(stdout, read_in));

	/* Tile the read in set to allocate the required data tiles
	 */
	PlutoConstraints *read_in_tiles =
			pluto_dist_get_required_data_tiles(read_in, src_copy_level, acc_name, prog);

	Array *arr = pluto_get_corrs_array(acc_name, prog);

	// ===== from here on starts the code generation =====================================================


	/* Read the array initialization stmt text from file
	 */
	char* init_stmt_text = (char *)malloc(2048 * sizeof(char));

	init_stmt_text[0] = '\0';

	sprintf(init_stmt_text+strlen(init_stmt_text), " __ifndef USE_LOCAL_ARRAYS");
	sprintf(init_stmt_text+strlen(init_stmt_text), " %s ",  pluto_dist_copy_init_text(arr, prog));
	sprintf(init_stmt_text+strlen(init_stmt_text), " __else");
	char *text= pluto_dist_read_init_text(arr->text, prog);
	char *mod_stmt_text = pluto_dist_modify_stmt_text(text,1, prog);

	if(mod_stmt_text[strlen(mod_stmt_text) -1] == '\n')
		mod_stmt_text[strlen(mod_stmt_text) -1] = 0;

	sprintf(init_stmt_text+strlen(init_stmt_text), " %s ", mod_stmt_text );
	sprintf(init_stmt_text+strlen(init_stmt_text), "__endif");

	/*
    if(options->verify_output)
		init_stmt_text = pluto_dist_copy_init_text(arr, prog);
    else{
		char *text= pluto_dist_read_init_text(arr->text, prog);

		 init_stmt_text = pluto_dist_modify_stmt_text(text,1, prog);
		 free(text);
    }
	 */

	char* alloc_stmt_text = pluto_dist_malloc_stmt_text(acc_name, prog, 1);

	Stmt *init_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			init_stmt_text, DATA_DIST_COPY, loop_num);

	Stmt *alloc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			alloc_stmt_text, DATA_DIST_INIT, loop_num);

	/* Add dimensions that actually scan the data tile space */
	for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
		pluto_stmt_add_dim(alloc_stmt, alloc_stmt->dim,
				alloc_stmt->trans->nrows, arr->iterators[j], H_LOOP, prog);
	}

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	pluto_constraints_intersect(alloc_stmt->domain, read_in_tiles);

	/* Add dimensions that actually scan the data space */
	//    if(options->verify_output){
	//		for (i=0; i < acc_nrows; i++) {
	//			pluto_stmt_add_dim(init_stmt, init_stmt->dim,
	//					init_stmt->trans->nrows, arr->iterators[i], H_LOOP, prog);
	//		}
	//    }
	//    else {
	for (i=0; i < acc_nrows; i++) {
		char iter[5];
		sprintf(iter, "d%d", i+1);
		pluto_stmt_add_dim(init_stmt, init_stmt->dim,
				init_stmt->trans->nrows, iter, H_LOOP, prog);
	}
	//    }

	/* Add constraints on copy loops - completes domain of alloc_stmt */
	pluto_constraints_intersect(init_stmt->domain, read_in);

	//correct the trans function of init stmt
	PlutoMatrix *trans_new = pluto_matrix_dup(init_stmt->trans);
	for (i = 0; i < arr->dim_orig; ++i) {
		for (j = 0; j < arr->dim_orig + prog->npar + 1; ++j) {
			trans_new->val[src_copy_level+i][src_copy_level + j] =
					arr->trans_orig->val[i][j];

		}
	}

	pluto_matrix_free(init_stmt->trans);
	init_stmt->trans = trans_new;

	Stmt **read_in_stmts = (Stmt **) malloc(2*sizeof(Stmt *));
	read_in_stmts[0] = alloc_stmt;
	read_in_stmts[1] = init_stmt;

	free(init_stmt_text);
	free(alloc_stmt_text);

	return read_in_stmts;
}

/* 
 * Optimized communication code generation for the write-out set
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * wacc_stmts:  all <statements, wacc> writing to this data variable
 * This function is called per data variable
 * This function can be called by any flow-based communication scheme
 * caller is responsible for freeing write_out_stmts returned
 */  
Stmt **gen_write_out_code(struct stmt_access_pair **wacc_stmts, int num_accs,
		PlutoProg *prog, Stmt *anchor_stmt, int *copy_level,
		int outer_dist_loop_level, int loop_num, FILE *headerfp)
{
	int i, src_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	char *access = reconstruct_access(wacc_stmts[0]->acc);
	acc_nrows = wacc_stmts[0]->acc->mat->nrows;
	char *acc_name = wacc_stmts[0]->acc->name;

	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	if (src_copy_level >= 1)    {
		assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
	}

	PlutoConstraints *write_out = NULL;
	for (i=0; i<num_accs; i++)  {
		PlutoConstraints *write_out_one;
		//write_out_one = compute_last_writes(wacc_stmts[k], src_copy_level, prog);
		write_out_one = compute_write_out(wacc_stmts[i], src_copy_level, prog);
		if (write_out == NULL) write_out = pluto_constraints_dup(write_out_one);
		else write_out = pluto_constraints_unionize(write_out, write_out_one);
		pluto_constraints_free(write_out_one);
	}

	IF_DEBUG(printf("Write out set\n"););
	IF_DEBUG(pluto_constraints_print(stdout, write_out));


	// ===== from here on starts the code generation =====================================================

	/***************************************************************************************************/

	char *lw_buf_size = get_parametric_bounding_box(write_out, src_copy_level, acc_nrows, prog->npar, (const char **)prog->params);
	char *lw_recv_buf_size = malloc(1280);
	strcpy(lw_recv_buf_size, lw_buf_size);

	if (src_copy_level >= 1) {
		PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
		char *total_extent = malloc(1024);
		strcpy(total_extent, "1");
		/* Assumes load-balanced distribution */
		for (i=outer_dist_loop_level; i<src_copy_level; i++) {
			char *extent;
			/* Just the first one in flow_out is enough (rest all should give the
			 * same since they are all under the same parallel loop and each
			 * iteration of the parallel loop writes to distinct data) */
			get_parametric_extent_const(anchor_stmt_new_dom, i, prog->npar,
					(const char **)prog->params, &extent, NULL);
			/* The + nprocs is needed since the send buffer size should be larger
			 * when some processors have more iterations than others - some processors
			 * have an extra iteration in each dimension in the worst case */
			sprintf(total_extent+strlen(total_extent), "*(%s+nprocs)", extent);
			free(extent);
		}
		sprintf(lw_buf_size+strlen(lw_buf_size), "*floorf((%s)/(float)nprocs)",
				total_extent);
		free(total_extent);
		pluto_constraints_free(anchor_stmt_new_dom);
	}
	IF_DEBUG(printf("Last writer buffer size for %s: %s\n", acc_name, lw_buf_size););

	if (src_copy_level >= 1)    {
		assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
				prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);
	}

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	char args[1024];
	strcpy(args,"");
	sprintf(args+strlen(args), "t%d", 1);
	/* make it src_copy_level+1 since we use an extra dimension to separate
	 * statements */
	for (i=1; i<src_copy_level; i++) {
		sprintf(args+strlen(args), ",t%d", i+1);
	}

	char decl_args[1024];
	strcpy(decl_args,"");
	sprintf(decl_args+strlen(decl_args), "int ts%d", 1);
	for (i=1; i<src_copy_level; i++) {
		sprintf(decl_args+strlen(decl_args), ",int ts%d", i+1);
	}

	char params[1024];
	strcpy(params, "");
	if (prog->npar>=1) {
		sprintf(params+strlen(params), "%s", prog->params[0]);
		for (i=1; i<prog->npar; i++) {
			sprintf(params+strlen(params), ",%s", prog->params[i]);
		}
	}

	char **iters;
	iters = malloc(acc_nrows * sizeof(char *));
	for (i=0; i < acc_nrows; i++) {
		iters[i] = malloc(5);
		sprintf(iters[i], "d%d", i+1);
	}

	char *lwbufname = concat("lw_buf_", acc_name);
	char *lw_count_name = concat("lw_count_", acc_name);
	Array *arr = pluto_get_corrs_array(acc_name, prog);

	sprintf(prog->decls+strlen(prog->decls), "if (fopen(\".test\", \"r\")) {");
	sprintf(prog->decls+strlen(prog->decls), "%s\
            = (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &lw_buf_size_%s);\n", 
            lwbufname, lwbufname, lw_buf_size, acc_name);

	char *lw_holder_text = malloc(512);
	strcpy(lw_holder_text, "");
	sprintf(lw_holder_text+strlen(lw_holder_text), "if (pi_%d(%s,%s,nprocs) == my_rank)",
			loop_num, args, params);
	sprintf(lw_holder_text+strlen(lw_holder_text), "%s = writeout_pack_%s_%d(%s,%s,%s",
			lw_count_name, acc_name, loop_num, args, lwbufname, lw_count_name);
	if (options->variables_not_global) {
		sprintf(lw_holder_text+strlen(lw_holder_text), ",%s", acc_name);
	}
	print_data_dist_parm_call(lw_holder_text, arr);
	sprintf(lw_holder_text+strlen(lw_holder_text),");");

	char *lw_stmt_text = malloc(strlen(lwbufname) + strlen("[")+
			strlen(lw_count_name) + strlen("++] = ") + strlen(access) + 1);
	sprintf(lw_stmt_text, "%s[%s++] = %s", lwbufname, lw_count_name, access);

	fprintf(packfp, "int writeout_pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, lwbufname, lw_count_name);
	fprintf(headerfp, "int writeout_pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, lwbufname, lw_count_name);
	if (options->variables_not_global) {
		fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
		fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
	}
	add_data_dist_parm_decl(headerfp, acc_name, prog);
	fprintf(headerfp,");\n");
	generate_pack_or_unpack(packfp, prog, write_out, lw_stmt_text, IN_FUNCTION, iters, src_copy_level, acc_nrows, acc_name, lw_count_name);

	/* Create statement stub to get hold of the outer loops for the copy
	 * stmt; rest of the loops to be added are the actual copy loops */

	Stmt *lw_copy_stmt = create_helper_stmt(anchor_stmt, src_copy_level, lw_stmt_text, LW_COPY_OUT, loop_num);
	Stmt *lw_holder = create_helper_stmt(anchor_stmt, src_copy_level, lw_holder_text, LW_COPY_OUT, loop_num);

	//IF_DEBUG(pluto_stmt_print(stdout, flow_copy_stmt););
	//IF_DEBUG(pluto_stmt_print(stdout, lw_copy_stmt););

	/* Add dimensions that actually scan the data space */
	for (i=0; i < acc_nrows; i++) {
		char iter[5];
		sprintf(iter, "d%d", i+1);
		pluto_stmt_add_dim(lw_copy_stmt, lw_copy_stmt->dim, lw_copy_stmt->trans->nrows, iter, H_LOOP, prog);
	}

	/* Add constraints on copy loops - completes domain of copy stmt */
	pluto_constraints_intersect(lw_copy_stmt->domain, write_out);

	if(options->data_dist){
		//		Array *arr = pluto_get_corrs_array(acc_name, prog);
		//		int j;
		//		//correct the trans function of init stmt
		//		PlutoMatrix *trans_new = pluto_matrix_dup(lw_copy_stmt->trans);
		//		for (i = 0; i < arr->dim_orig; ++i) {
		//			for (j = 0; j < arr->dim_orig + prog->npar + 1; ++j) {
		//				trans_new->val[src_copy_level+i][src_copy_level + j] =
		//						arr->trans_orig->val[i][j];
		//
		//			}
		//		}
		//
		//		pluto_matrix_free(lw_copy_stmt->trans);
		//		lw_copy_stmt->trans = trans_new;
	}

	char *lwrecvbufname = concat("lw_recv_buf_", acc_name);
	char *lw_displsname = concat("displs_lw_", acc_name);
	char *lw_recv_counts_name = concat("lw_recv_counts_", acc_name);

	char *comm_text_w = malloc(1024);

	strcpy(comm_text_w, "");

	sprintf(comm_text_w+strlen(comm_text_w), "\
        assert((nprocs == 1) || (%s <= %s[1]));\
        MPI_Gather(&%s, 1, MPI_INT,\
        %s, 1, MPI_INT, 0, MPI_COMM_WORLD);", 
        lw_count_name, lw_displsname, lw_count_name, lw_recv_counts_name);
	sprintf(comm_text_w+strlen(comm_text_w), "MPI_Gatherv(%s, %s, MPI_DOUBLE,\
        %s, %s, %s, MPI_DOUBLE, 0, MPI_COMM_WORLD); %s = 0;",
        lwbufname, lw_count_name, lwrecvbufname, lw_recv_counts_name,
        lw_displsname, lw_count_name);
	sprintf(comm_text_w+strlen(comm_text_w), "\
       for (__p=0; __p<nprocs; __p++) curr_%s[__p] = %s[__p];", 
       lw_displsname, lw_displsname);

	/* Create statement stub to get hold of the outer loops for the copy
	 * stmt; rest of the loops to be added are the actual copy loops */
	Stmt *comm_stmt_w = create_helper_stmt(anchor_stmt, outer_dist_loop_level, comm_text_w, LW_COMM_CALL, loop_num);

	sprintf(prog->decls+strlen(prog->decls),
			"%s = (double *) polyrt_max_alloc(%s, nprocs*lw_buf_size_%s, &lw_recv_buf_size_%s);\n",
			lwrecvbufname, lwrecvbufname, acc_name, acc_name);

	sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
        %s[__p] = __p*lw_buf_size_%s/sizeof(double);}\n\n", lw_displsname, acc_name);

	sprintf(prog->decls+strlen(prog->decls), "}");

	char lw_proc_stmt_text[1024];
	sprintf(lw_proc_stmt_text, "proc = pi_%d(", loop_num);
	for (i=0; i<src_copy_level+prog->npar; i++)    {
		if (i>=1)  sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), ", ");
		if (i<=src_copy_level-1) sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), "t%d", i+1);
		else sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), "%s", prog->params[i-src_copy_level]);
	}
	sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), ", nprocs);");
	sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text),
			"if ((my_rank != proc) && (%s[proc] > 0)) {", lw_recv_counts_name);
	sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text),
			"curr_%s[proc] = writeout_unpack_%s_%d(%s,%s,curr_%s[proc]",
			lw_displsname, acc_name, loop_num, args, lwrecvbufname, lw_displsname);
	if (options->variables_not_global) {
		sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), ",%s", acc_name);
	}
	print_data_dist_parm_call(lw_proc_stmt_text, arr);
	sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text),");}");

	char *lw_copyback_text = malloc(strlen(access) + strlen(" = lw_recv_buf_")
			+ strlen(acc_name) + strlen("[") + strlen(lw_displsname) +
			strlen("[proc]+ count++];") + 1 + strlen(" __ifndef USE_LOCAL_ARRAYS   __endif"));
	sprintf(lw_copyback_text, " __ifndef USE_LOCAL_ARRAYS %s = %s[%s++]; __endif",
			access, lwrecvbufname, lw_displsname);

	Stmt *lw_copyback_stmt = create_helper_stmt(anchor_stmt, src_copy_level, lw_copyback_text, LW_COPY_IN, loop_num);
	Stmt *lw_copyback_proc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			lw_proc_stmt_text, LW_COPY_IN, loop_num);

	/* Add dimensions that actually scan the data space */
	for (i=0; i < acc_nrows; i++) {
		char iter[5];
		sprintf(iter, "d%d", i+1);
		pluto_stmt_add_dim(lw_copyback_stmt, lw_copyback_stmt->dim,	lw_copyback_stmt->trans->nrows, iter, H_LOOP, prog);
	}
	sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text),");}");

	fprintf(packfp, "int writeout_unpack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, lwrecvbufname, lw_displsname);
	fprintf(headerfp, "int writeout_unpack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, lwrecvbufname, lw_displsname);
	if (options->variables_not_global) {
		fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
		fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
	}
	add_data_dist_parm_decl(headerfp, acc_name, prog);
	fprintf(headerfp,");\n");
	generate_pack_or_unpack(packfp, prog, write_out, lw_copyback_text, LW_COPY_IN, iters, src_copy_level, acc_nrows,
			acc_name, lw_displsname);
	if(options->data_dist){
		//		Array *arr = pluto_get_corrs_array(acc_name, prog);
		//		int j;
		//		//correct the trans function of init stmt
		//		PlutoMatrix *trans_new = pluto_matrix_dup(lw_copyback_stmt->trans);
		//		for (i = 0; i < arr->dim_orig; ++i) {
		//			for (j = 0; j < arr->dim_orig + prog->npar + 1; ++j) {
		//				trans_new->val[src_copy_level+i][src_copy_level + j] =
		//						arr->trans_orig->val[i][j];
		//
		//			}
		//		}
		//
		//		pluto_matrix_free(lw_copyback_stmt->trans);
		//		lw_copyback_stmt->trans = trans_new;
	}

	//    Stmt *lw_copyback_proc_stmt = create_helper_stmt(anchor_stmt, src_copy_level, lw_proc_stmt_text, LW_COPY_IN, loop_num);
	// IF_DEBUG(pluto_stmt_print(stdout, lw_copyback_proc_stmt););
	/***************************************************************************************************/

	// write out stmts organized as follows (used later during code generation)
	// pack-guard stmt, comm stmt, unpack-guard stmt
#define NUM_WRITE_OUT_STMTS 3
	Stmt **write_out_stmts = (Stmt **) malloc(NUM_WRITE_OUT_STMTS*sizeof(Stmt *));
	write_out_stmts[0] = lw_holder;
	write_out_stmts[1] = comm_stmt_w;
	write_out_stmts[2] = lw_copyback_proc_stmt;

	fclose(packfp);
	for (i=0; i < acc_nrows; i++) {
		free(iters[i]);
	}
	free(iters);

	pluto_constraints_free(write_out);
	free(lw_stmt_text);
	free(lw_holder_text);
	free(lw_copyback_text);
	free(lw_buf_size);
	free(lw_recv_buf_size);
	free(lw_count_name);
	free(lwbufname);
	free(lwrecvbufname);
	free(lw_displsname);
	free(lw_recv_counts_name);
	free(comm_text_w);
	free(access);

	return write_out_stmts;

}

void print_data_dist_parm_call(char* str, Array *arr){

	int s;
	if(options->data_dist && arr != NULL){

		sprintf(str+strlen(str),", %s_trans", arr->text);
		for (s = 1; s < arr->num_tiled_loops; ++s) {
			sprintf(str+strlen(str),", %s_size_%d", arr->text, s + arr->first_tile_dim);
		}

		sprintf(str+strlen(str),",buff_mang_%s", arr->text);
	}

	return;
}


void print_data_dist_parm_call_from_access(char* str, PlutoAccess *access, PlutoProg *prog){

	Array *arr = pluto_get_corrs_array(access->name, prog);

	int s;
	if(options->data_dist && arr != NULL){

		sprintf(str+strlen(str),", %s_trans", arr->text);
		for (s = 1; s < arr->num_tiled_loops; ++s) {
			sprintf(str+strlen(str),", %s_size_%d", arr->text, s + arr->first_tile_dim);
		}

		sprintf(str+strlen(str),",buff_mang_%s", arr->text);
	}

	return;
}

void fprint_data_dist_parm_call(FILE *fp, char *acc_name, PlutoProg *prog){

	int s;
	if(options->data_dist){

		Array *arr = pluto_get_corrs_array(acc_name, prog);
        if(arr == NULL) return;
		fprintf(fp,", %s_trans", arr->text);
		for (s = 1; s < arr->num_tiled_loops; ++s) {
			fprintf(fp,", %s_size_%d", arr->text, s + arr->first_tile_dim);
		}

		fprintf(fp,", buff_mang_%s", arr->text);
	}

	return;
}
void add_data_dist_parm_decl(FILE *fp, char *acc_name, PlutoProg *prog){


	if(options->data_dist) {
		Array *arr = pluto_get_corrs_array(acc_name, prog);
		int s =0;
		fprintf(fp,", %s",arr->data_type);
		for (s = 0; s < 2; ++s) {
			fprintf(fp,"*");
		}
		fprintf(fp," %s_trans",arr->text);

		for (s = 1; s < arr->num_tiled_loops; ++s) {
			fprintf(fp,",int %s_size_%d", arr->text, s + arr->first_tile_dim);
		}

		fprintf(fp,", BufferManager *buff_mang_%s",arr->text);

	}

	return;
}

void generate_pack_or_unpack(FILE *packfp, PlutoProg *prog,
		PlutoConstraints *constraints, char *stmttext, PlutoStmtType stmttype, char **iters,
		int src_copy_level, int acc_nrows, char* acc_name, char *send_count)
{

	add_data_dist_parm_decl(packfp, acc_name, prog);

	fprintf(packfp,"){\n");
	fprintf(packfp, "\nint recv_proc;\n");

	if (options->dynschedule && options->commopt_foifi) {
		fprintf(packfp, "\nint _i;\n");
		fprintf(packfp, "\nint current_index = *p_current_index;\n");
	}

	FILE *packcloogfp;
	PlutoProg *pack;
	PlutoMatrix *packtrans;
	int i, j, total_level = src_copy_level + acc_nrows;
	PlutoConstraints *tdpoly;

	pack = pluto_prog_alloc();

	if(options->data_dist){
		pack->arrays = prog->arrays;
		pack->narrays = prog->narrays;
	}

	for (i=0; i<src_copy_level; i++) {
		char param[6];
		sprintf(param, "ts%d",i+1);
		pluto_prog_add_param(pack, param, pack->npar);
	}

	for (i=0; i<prog->npar; i++) {
		pluto_prog_add_param(pack, prog->params[i], pack->npar);
	}

	for (i=0;i<src_copy_level; i++) {
		pluto_prog_add_hyperplane(pack,0,H_LOOP);
	}

	tdpoly = pluto_constraints_dup(constraints);

	// move acc_nrows loop iterators to the beginning/top
	for (i=0; i<acc_nrows; i++) {
		pluto_constraints_add_dim(tdpoly, 0, NULL);
	}
	for (i=0; i<acc_nrows; i++) {
		pluto_constraints_interchange_cols(tdpoly, i, i+total_level);
	}
	for (i=0; i<acc_nrows; i++) {
		pluto_constraints_remove_dim(tdpoly, total_level);
	}

	packtrans = pluto_matrix_identity(acc_nrows);
	for (i=0; i<src_copy_level+prog->npar+1; i++) {
		pluto_matrix_add_col(packtrans, packtrans->ncols);
	}

	pluto_add_stmt(pack, tdpoly, packtrans, iters, stmttext, stmttype);

	if (options->data_dist && stmttype != LW_COPY_IN) {
		Array *arr = pluto_get_corrs_array(acc_name, prog);

		PlutoConstraints *data_tiles =
				pluto_dist_get_required_data_tiles(constraints, src_copy_level, acc_name,
						prog);

		acc_nrows = arr->num_tiled_loops;

		// move acc_nrows loop iterators to the beginning/top
		for (i=0; i<acc_nrows; i++) {
			pluto_constraints_add_dim(data_tiles, 0, NULL);
		}
		for (i=0; i<acc_nrows; i++) {
			pluto_constraints_interchange_cols(data_tiles, i, i+total_level);
		}
		for (i=0; i<acc_nrows; i++) {
			pluto_constraints_remove_dim(data_tiles, total_level);
		}

		int num_dims = data_tiles->ncols - prog->npar - 1;
		if(num_dims < arr->num_tiled_loops)
			num_dims = arr->num_tiled_loops;

		char **new_iter = (char **)malloc(num_dims * sizeof(char*));
		int count = 0;
		for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
			new_iter[count++] = strdup(arr->iterators[j]);
		}
		for(j=count; j <num_dims; j++){
			new_iter[j] = strdup("temp1");
		}
		char* alloc_stmt_text = pluto_dist_malloc_stmt_text(acc_name, prog, 1);

		pluto_add_stmt(pack, data_tiles, packtrans, new_iter, alloc_stmt_text, stmttype);
		//
		//		Stmt *alloc_stmt = create_helper_stmt(pack->stmts[0], 0,
		//				alloc_stmt_text, DATA_DIST_MANG);
		//
		//	   /* Add dimensions that actually scan the data tile space */
		//		for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
		//		   pluto_stmt_add_dim(alloc_stmt, alloc_stmt->dim,
		//			alloc_stmt->trans->nrows, arr->iterators[j], H_LOOP, prog);
		//	   }
		//
		//		/* Add constraints on copy loops - completes domain of alloc_stmt */
		//		alloc_stmt->domain = pluto_constraints_dup(data_tiles);
		//
		//		pluto_add_given_stmt(pack, alloc_stmt);
		//
		pluto_prog_add_hyperplane(pack, 0, H_SCALAR);

		assert(pack->nstmts == 2);

		//Seperate out data alloc and pack stmts
		Stmt *s = pack->stmts[0];
		s->trans->val[0][s->trans->ncols - 1] = 2;

	}

	packcloogfp = fopen("packunpack.cloog", "w+");
	assert(packcloogfp != NULL);
	pluto_gen_cloog_file(packcloogfp, pack);
	rewind(packcloogfp);
	fflush(packcloogfp);

	generate_declarations(pack, packfp);
	pluto_gen_cloog_code(pack, 1, pack->num_hyperplanes, packcloogfp, packfp);
	if (options->dynschedule && options->commopt_foifi) {
		fprintf(packfp, "\n*p_current_index = current_index;\n");
	}

	fprintf(packfp, "\nreturn %s;\n", send_count);
	for (j = 0; j < pack->nstmts; ++j) {
		fprintf(packfp, "#undef S%d\n", j+1);
	}
	fprintf(packfp, "\n}\n\n");

	pluto_matrix_free(packtrans);
	pluto_constraints_free(tdpoly);
	fclose(packcloogfp);
	pluto_prog_free(pack);
}

Stmt *gen_comm_code_count_recvs(struct stmt_access_pair ***wacc_stmts, int *num_stmts_per_wacc, int num_data,
		PlutoProg *prog, Stmt *anchor_stmt, int *copy_level, int loop_num)
{
	int i, src_copy_level;
	src_copy_level = copy_level[loop_num];
	assert(src_copy_level>=1);

	char args[1024];
	strcpy(args,"");
	sprintf(args+strlen(args), "t%d", 1);
	/* make it src_copy_level+1 since we use an extra dimension to separate
	 * statements */
	for (i=1; i<src_copy_level; i++) {
		sprintf(args+strlen(args), ",t%d", i+1);
	}

	char params[1024];
	strcpy(params, "");
	if (prog->npar>=1) {
		sprintf(params+strlen(params), "%s", prog->params[0]);
		for (i=1; i<prog->npar; i++) {
			sprintf(params+strlen(params), ",%s", prog->params[i]);
		}
	}

	char *count_recvs_stmt_text = malloc(8192);
	sprintf(count_recvs_stmt_text, "proc = pi_%d(%s,%s, nprocs); ", loop_num, args, params);
	sprintf(count_recvs_stmt_text+strlen(count_recvs_stmt_text), "if (my_rank != proc) { ");
	for (i=0; i<num_data; i++) {
		char *acc_name = wacc_stmts[i][0]->acc->name;
		sprintf(count_recvs_stmt_text+strlen(count_recvs_stmt_text),
				"_num_recvs_%s += is_receiver_%s_%d(%s,%s,my_rank,nprocs);",
				acc_name, acc_name, loop_num, args, params);
	}
	sprintf(count_recvs_stmt_text+strlen(count_recvs_stmt_text), "}");

	Stmt *count_recvs_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			count_recvs_stmt_text, COPY_IN, loop_num);

	return count_recvs_stmt;
}

Stmt *gen_comm_code_async_recv(struct stmt_access_pair ***wacc_stmts, int *num_stmts_per_wacc, int num_data,
		PlutoProg *prog, Stmt *anchor_stmt, int *copy_level, int loop_num)
{
	int i, src_copy_level;
	src_copy_level = copy_level[loop_num];
	assert(src_copy_level>=1);

	char args[1024];
	strcpy(args,"");
	sprintf(args+strlen(args), "t%d", 1);
	/* make it src_copy_level+1 since we use an extra dimension to separate
	 * statements */
	for (i=1; i<src_copy_level; i++) {
		sprintf(args+strlen(args), ",t%d", i+1);
	}

	char params[1024];
	strcpy(params, "");
	if (prog->npar>=1) {
		sprintf(params+strlen(params), "%s", prog->params[0]);
		for (i=1; i<prog->npar; i++) {
			sprintf(params+strlen(params), ",%s", prog->params[i]);
		}
	}

	char *async_recv_stmt_text = malloc(8192);
	if (options->timereport)
		sprintf(async_recv_stmt_text, "IF_TIME(t_unpack_start = rtclock());");
	else
		strcpy(async_recv_stmt_text, "");
	for (i=0; i<num_data; i++) {
		char *data_tag;
		data_tag = malloc(10);
		get_data_tag(prog, wacc_stmts[i][0]->acc, &data_tag);
		char *acc_name = wacc_stmts[i][0]->acc->name;
		char *max_elements_name = concat("max_num_elements_", acc_name);
		sprintf(async_recv_stmt_text+strlen(async_recv_stmt_text),
				"assert(\"increase __MAX_NUM_RECVS\" && (_num_recvs_%s <= __MAX_NUM_RECVS));",
				acc_name);
		sprintf(async_recv_stmt_text+strlen(async_recv_stmt_text),
				"for (__p=0; __p<_num_recvs_%s; __p++) { \
			MPI_Irecv(recv_buf_%s[__p], %s, MPI_DOUBLE,\
			MPI_ANY_SOURCE, %s, MPI_COMM_WORLD, &recv_reqs_%s[__p]); }",
			acc_name, acc_name, max_elements_name, data_tag, acc_name);
		free(data_tag);
	}
	if (options->timereport)
		sprintf(async_recv_stmt_text+strlen(async_recv_stmt_text), "IF_TIME(t_unpack += rtclock() - t_unpack_start);");

	Stmt *async_recv_stmt = create_helper_stmt(anchor_stmt, src_copy_level-1,
			async_recv_stmt_text, COMM_CALL, loop_num);

	return async_recv_stmt;
}

Stmt **gen_tasks_code(Stmt **loop_stmts, int nstmts, int *copy_level, PlutoProg *prog, int num_data,
		int loop_num, int nloops, int *pi_mappings, char *tasks_loops_decl, FILE *outfp, FILE *headerfp)
{
	int l, i, j;
	Stmt *anchor_stmt = get_new_anchor_stmt(loop_stmts, nstmts);

	generate_outgoing(loop_stmts, nstmts, copy_level, prog, loop_num,
			pi_mappings, tasks_loops_decl, outfp, headerfp);
	generate_incoming(loop_stmts, nstmts, copy_level, prog, loop_num,
			pi_mappings, outfp, headerfp);

	PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
	int src_copy_level = copy_level[loop_num];
	char *extent;
	char *lbexpr;
	char *dim_sizes = malloc((strlen("max_num_tasks_loop_dim*") + 10)*src_copy_level + 1);
	strcpy(dim_sizes,"");
	int num_dims = 0;
	int dims[src_copy_level];
	for (i=0; i<src_copy_level; i++) {
		if (pluto_is_hyperplane_loop(anchor_stmt, i)) {
			/* Just the first one in flow_out is enough (rest all should give the
			 * same since they are all under the same parallel loop and each
			 * iteration of the parallel loop writes to distinct data) */
			get_parametric_extent_const(anchor_stmt_new_dom, i, prog->npar,
					(const char **)prog->params, &extent, &lbexpr);
			fprintf(headerfp, "extern int max_num_tasks_loop%d_dim%d;\n", loop_num, num_dims);
			fprintf(headerfp, "int max_num_tasks_loop%d_dim%d;\n", loop_num, num_dims);
			fprintf(headerfp, "extern int lb_tasks_loop%d_dim%d;\n", loop_num, num_dims);
			fprintf(headerfp, "int lb_tasks_loop%d_dim%d;\n", loop_num, num_dims);
			sprintf(prog->decls+strlen(prog->decls), "max_num_tasks_loop%d_dim%d = ceil(%s);\n", loop_num, num_dims, extent);
			sprintf(prog->decls+strlen(prog->decls), "lb_tasks_loop%d_dim%d = floor(%s);\n", loop_num, num_dims, lbexpr);
			sprintf(dim_sizes+strlen(dim_sizes), "max_num_tasks_loop%d_dim%d*", loop_num, num_dims);
			free(extent);
			free(lbexpr);
			dims[num_dims++] = i;
		}
	}
	if (options->dynschedule_graph) {
		sprintf(prog->decls+strlen(prog->decls), "tbb::flow::continue_node<tbb::flow::continue_msg> **tasks_loop%d;\n", loop_num);
		sprintf(prog->decls+strlen(prog->decls), "tasks_loop%d = (tbb::flow::continue_node<tbb::flow::continue_msg> **)malloc(%ssizeof(tbb::flow::continue_node<tbb::flow::continue_msg> *));\n", loop_num, dim_sizes);
	}
	else {
		fprintf(headerfp, "extern int *tasks_loop%d;\n", loop_num);
		fprintf(headerfp, "int *tasks_loop%d;\n", loop_num);
		sprintf(prog->decls+strlen(prog->decls), "tasks_loop%d = (int *)malloc(%ssizeof(int));\n", loop_num, dim_sizes);
	}
	free(dim_sizes);
	pluto_constraints_free(anchor_stmt_new_dom);

	char *args = malloc(256);
	strcpy(args,"");
	sprintf(args+strlen(args), "t%d", 1);
	/* make it src_copy_level+1 since we use an extra dimension to separate
	 * statements */
	for (i=1; i<src_copy_level; i++) {
		sprintf(args+strlen(args), ",t%d", i+1);
	}

	char *params = malloc (1024);
	strcpy(params, "");
	if (prog->npar>=1) {
		sprintf(params+strlen(params), "%s", prog->params[0]);
		for (i=1; i<prog->npar; i++) {
			sprintf(params+strlen(params), ",%s", prog->params[i]);
		}
	}

	char *indices = malloc(1024);
	strcpy(indices,"[");
	for (i=0; i<num_dims; i++) {
		sprintf(indices+strlen(indices), "(t%d-lb_tasks_loop%d_dim%d)",
				dims[i]+1, loop_num, i);
		for (j=i+1; j<num_dims; j++) {
			sprintf(indices+strlen(indices), "*max_num_tasks_loop%d_dim%d",
					loop_num, j);
		}
		strcat(indices,"+");
	}
	strcat(indices,"0]");

	char *tasks_loops_args = NULL;
	tasks_loops_args = malloc(256);
	strcpy(tasks_loops_args, "");
	for (l=0; l<nloops; l++) {
		sprintf(tasks_loops_args+strlen(tasks_loops_args), ",tasks_loop%d", l);
	}

	char *init_all_tasks_text = malloc(4096);
	strcpy(init_all_tasks_text, "");

	char *add_outgoing_edges_text = NULL;

	if (options->distmem) {
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"if (pi_%d(%s,%s,nprocs) == my_rank) { ",
				loop_num, args, params);

		// owned tasks
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text), "_num_tasks_to_execute++;");
		// initialize firing count of owned tasks
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"tasks_loop%d%s = get_num_remote_src_tasks_%d(%s,%s,my_rank,nprocs) + \
                get_num_local_src_tasks_%d(%s,%s,my_rank,nprocs);",
                loop_num, indices, loop_num, args, params, loop_num, args, params);
		// enqueue if firing count is 0
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"if (tasks_loop%d%s == 0) { \
                remote_dep_tasks = get_num_remote_dep_tasks_%d(%s,%s,my_rank,nprocs); \
                local_dep_tasks = get_num_local_dep_tasks_%d(%s,%s,my_rank,nprocs); \
                if (__is_multi_partitioned[%d] || __is_block_cyclic[%d]) affinity = -1; \
                else affinity = pi_threads_%d(%s,%s,my_rank,nprocs,num_threads); \
                __Task task(%d,%s,remote_dep_tasks,local_dep_tasks,affinity); \
                pqueue.push(task); ",
                loop_num, indices, loop_num, args, params, loop_num, args, params,
                loop_num, loop_num, loop_num, args, params, loop_num, args);
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"IF_DYNSCHEDULER_MORE_DEBUG_PRINT(\
                fprintf(__debug_print_fp, \"added node %%d loop %d task ",
                loop_num);
		for (i=0; i<src_copy_level; i++) {
			sprintf(init_all_tasks_text+strlen(init_all_tasks_text),  "%%d ");
		}
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),  "\
                remote_dep_tasks %%d local_dep_tasks %%d ready tasks %%lu\\n\", \
                my_rank, %s, remote_dep_tasks, local_dep_tasks, pqueue.size())); \
                IF_DEBUG_FLUSH(fflush(__debug_print_fp)); }",
                args);

		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"} else if (is_receiver_%d(%s,%s,my_rank,nprocs)) { ",
				loop_num, args, params);

		// remote tasks
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text), "_num_tasks_to_unpack++;");
		// initialize firing count of remote tasks (for unpack)
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"tasks_loop%d%s = %d + \
                get_num_local_src_tasks_%d(%s,%s,my_rank,nprocs);",
                loop_num, indices, num_data, loop_num, args, params);
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"IF_DYNSCHEDULER_MORE_DEBUG_PRINT(\
                fprintf(__debug_print_fp, \"node %%d loop %d task ",
                loop_num);
		for (i=0; i<src_copy_level; i++) {
			sprintf(init_all_tasks_text+strlen(init_all_tasks_text),  "%%d ");
		}
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),  "\
                initial count %%d\\n\", \
                my_rank, %s, tasks_loop%d%s)); \
                IF_DEBUG_FLUSH(fflush(__debug_print_fp));",
                args, loop_num, indices);

		sprintf(init_all_tasks_text+strlen(init_all_tasks_text), "}");
	}
	else if (options->dynschedule_graph) {
		int num_data;
		PlutoAccess **accs;
		accs = pluto_get_all_distinct_vars(prog->stmts, prog->nstmts, &num_data);

		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"tasks_loop%d%s = \
                new tbb::flow::continue_node<tbb::flow::continue_msg>(*__graph,__Task(%d,%s", 
                loop_num, indices, loop_num, args);
		if (options->variables_not_global) {
			for (i=0; i<num_data; i++) {
				sprintf(init_all_tasks_text+strlen(init_all_tasks_text), ",%s", accs[i]->name);
			}
		}
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text), "));");
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"if (get_num_local_src_tasks_%d(%s,%s) == 0) \
                { tbb::flow::make_edge(*__start, *tasks_loop%d%s); }",
                loop_num, args, params, loop_num, indices);

		add_outgoing_edges_text = malloc(1024);
		sprintf(add_outgoing_edges_text,
				"add_outgoing_edges_%d(%s,%s%s,tasks_loop%d%s)",
				loop_num, args, params, tasks_loops_args, loop_num, indices);

		free(accs);
	}
	else {
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text), "_num_tasks_to_execute++;");
		// initialize firing count of owned tasks
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"tasks_loop%d%s = \
                get_num_local_src_tasks_%d(%s,%s);",
                loop_num, indices, loop_num, args, params);
		// enqueue if firing count is 0
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"if (tasks_loop%d%s == 0) { \
                local_dep_tasks = get_num_local_dep_tasks_%d(%s,%s); \
                affinity = pi_%d(%s,%s,num_threads); \
                __Task task(%d,%s,local_dep_tasks,affinity); \
                pqueue.push(task); ",
                loop_num, indices, loop_num, args, params,
                loop_num, args, params, loop_num, args);
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),
				"IF_DYNSCHEDULER_MORE_DEBUG_PRINT(\
                fprintf(__debug_print_fp, \"added node %%d loop %d task ",
                loop_num);
		for (i=0; i<src_copy_level; i++) {
			sprintf(init_all_tasks_text+strlen(init_all_tasks_text),  "%%d ");
		}
		sprintf(init_all_tasks_text+strlen(init_all_tasks_text),  "\
                local_dep_tasks %%d ready tasks %%lu\\n\", \
                my_rank, %s, local_dep_tasks, pqueue.size())); \
                IF_DEBUG_FLUSH(fflush(__debug_print_fp)); }",
                args);
	}

	Stmt *init_all_tasks_stmt = create_helper_stmt(anchor_stmt, src_copy_level, init_all_tasks_text, ALL_TASKS, loop_num);
	free(init_all_tasks_text);
	for (i=0; i<src_copy_level; i++) { // required for generating iterator declarations
		init_all_tasks_stmt->hyp_types[i] = loop_stmts[0]->hyp_types[i];
	}

	Stmt *add_outgoing_edges_stmt = NULL;
	if (options->dynschedule_graph) {
		add_outgoing_edges_stmt = create_helper_stmt(anchor_stmt, src_copy_level, add_outgoing_edges_text, ALL_TASKS, loop_num);
		free(add_outgoing_edges_text);
		for (i=0; i<src_copy_level; i++) { // required for generating iterator declarations
			add_outgoing_edges_stmt->hyp_types[i] = loop_stmts[0]->hyp_types[i];
		}
	}

	Stmt **tasks_stmts = NULL;
	if (options->dynschedule_graph) {
#define NUM_GEN_GRAPH_STMTS 2
		tasks_stmts = (Stmt **) malloc(NUM_GEN_GRAPH_STMTS*sizeof(Stmt*));
		tasks_stmts[0] = init_all_tasks_stmt;
		tasks_stmts[1] = add_outgoing_edges_stmt;
	}
	else {
#define NUM_GEN_TASKS_STMTS 1
		tasks_stmts = (Stmt **) malloc(NUM_GEN_TASKS_STMTS*sizeof(Stmt*));
		tasks_stmts[0] = init_all_tasks_stmt;
	}

	free(args);
	free(params);
	free(indices);
	free(tasks_loops_args);

	return tasks_stmts;
}

void gen_compute_task_cloog_code(PlutoProg *prog, int loop_num, Stmt **loop_stmts, int nstmts,
		int src_copy_level, FILE *outfp, FILE *headerfp)
{
	int i;

	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_all_distinct_vars(loop_stmts, nstmts, &num_data);

	PlutoProg *compute_task = NULL;
	compute_task = pluto_prog_alloc();

	for (i=0; i<prog->npar; i++) {
		pluto_prog_add_param(compute_task, prog->params[i], compute_task->npar);
	}
	if(options->data_dist){
		compute_task->arrays = prog->arrays;
		compute_task->narrays = prog->narrays;
	}

	for (i=0; i<nstmts; i++) {
		Stmt *stmt = loop_stmts[i];
		pluto_add_stmt(compute_task, stmt->domain, stmt->trans, stmt->iterators, stmt->text, ORIG_IN_FUNCTION);
		compute_task->stmts[compute_task->nstmts-1]->first_tile_dim = stmt->first_tile_dim;
		compute_task->stmts[compute_task->nstmts-1]->last_tile_dim = stmt->last_tile_dim;
	}

	IF_DEBUG(pluto_prog_print(stdout, compute_task));
	if (compute_task->nstmts >= 1) {
		assert(compute_task->stmts[0]->trans->nrows == compute_task->num_hyperplanes);
	}

	compute_task->num_parameterized_loops = src_copy_level;

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	fprintf(packfp, "void compute_task_%d(", loop_num);
	fprintf(headerfp, "void compute_task_%d(", loop_num);
	fprintf(packfp, "int t1");
	fprintf(headerfp, "int t1");
	for (i=1; i<src_copy_level; i++) {
		fprintf(packfp, ",int t%d", i+1);
		fprintf(headerfp, ",int t%d", i+1);
	}
	if (options->variables_not_global && !options->data_dist) {
		for (i=0; i<num_data; i++) {
			char *acc_name = accs[i]->name;
			fprintf(packfp, ",__DECLARATION_OF_%s", acc_name);
			fprintf(headerfp, ",__DECLARATION_OF_%s", acc_name);
		}
	}
	if(options->data_dist){
        accs = pluto_get_accs(loop_stmts, nstmts, &num_data);
		for (i=0; i<num_data; i++) {
			char *acc_name = accs[i]->name;
			add_data_dist_parm_decl(packfp, acc_name, prog);
			add_data_dist_parm_decl(headerfp, acc_name, prog);
		}
	}
	fprintf(packfp, ")\n{\n");
	fprintf(headerfp, ");\n");

	FILE *cloogfp = fopen("compute_task.cloog", "w+");
	pluto_gen_cloog_file(cloogfp, compute_task);
	rewind(cloogfp);

	pluto_mark_vec_stmts(cloogfp, compute_task);
	//    for(i=0;i<compute_task->nstmts;i++){
	//    	compute_task->stmts[i]->inner_loop_vec = 1;
	//    }

	generate_declarations(compute_task, packfp);
	pluto_gen_cloog_code(compute_task, -1, -1, cloogfp, packfp);
	fclose(cloogfp);

	for (i=0; i<compute_task->nstmts; i++) {
		fprintf(packfp, "#undef S%d\n", i+1);
		fprintf(packfp, "#undef decl_S%d\n", i+1);
		fprintf(packfp, "#undef orig_S%d\n", i+1);
	}
	fprintf(packfp, "}\n\n");

	fclose(packfp);

	pluto_prog_free(compute_task);

	free(accs);
}


void data_tile_prog_add_parm(PlutoProg *data_tile_prog, PlutoProg* prog, int src_copy_level){

	int i;

	for (i=0; i<src_copy_level; i++) {
		char param[6];

		sprintf(param, "ts%d",i+1);
		pluto_prog_add_param(data_tile_prog, param, data_tile_prog->npar);
	}

	for (i=0; i<prog->npar; i++) {
		pluto_prog_add_param(data_tile_prog, prog->params[i], data_tile_prog->npar);
	}

}

void data_tile_domain_remove_copy_level_dims(PlutoConstraints *data_tiles, int src_copy_level, PlutoProg *prog){

	int total_level = data_tiles->ncols - prog->npar - 1;
	int j;

	for (j=0; j<src_copy_level; j++) {
		pluto_constraints_add_dim(data_tiles,total_level, NULL);
	}
	for (j=0; j<src_copy_level; j++) {
		pluto_constraints_interchange_cols(data_tiles, j, j+total_level);
	}
	for (j=0; j<src_copy_level; j++) {
		pluto_constraints_remove_dim(data_tiles, 0);
	}

}

PlutoMatrix *get_data_tile_trans_func(int src_copy_level, int acc_nrows, PlutoProg *prog){

	int j;
	PlutoMatrix *trans = pluto_matrix_identity(acc_nrows);
	for (j=0; j<src_copy_level+prog->npar+1; j++) {
		pluto_matrix_add_col(trans, trans->ncols);
	}

	return trans;
}

char **get_data_tile_iterators(Array *arr, int acc_nrows){

	char** iters = malloc(acc_nrows * sizeof(char *));
	int j, k;

	assert(acc_nrows == arr->last_tile_dim - arr->first_tile_dim + 1);

	for(j=arr->last_tile_dim, k=0; j>= arr->first_tile_dim; --j, k++){
		iters[k] = arr->iterators[j];
	}

	return iters;

}

char **get_data_scan_iterators(Array *arr, int acc_nrows){

	char** iters = malloc(acc_nrows * sizeof(char *));
	int j;

	for(j=0;j<acc_nrows;j++){
		iters[j] = malloc(50);
		iters[j][0] = 0;
		sprintf(iters[j], "d%d", j+1);
	}

	return iters;

}

void gen_func_cloog_code(char *func_name, FILE *headerfp, int loop_num, int src_copy_level, int num_data, PlutoAccess **accs,
		PlutoProg *data_tile_prog, PlutoProg *prog, int writeout_func){

	int i;
	FILE *packfp = fopen("packunpack.c", "a");
	assert(packfp != NULL);

	fprintf(packfp, "void %s_%d(", func_name, loop_num);
	fprintf(headerfp, "void %s_%d(", func_name, loop_num);
	if(src_copy_level != 0){
		fprintf(packfp, "int ts1");
		fprintf(headerfp, "int ts1");
		for (i=1; i<src_copy_level; i++) {
			fprintf(packfp, ",int ts%d", i+1);
			fprintf(headerfp, ",int ts%d", i+1);
		}
	}
	else {
		fprintf(packfp, "int my_rank, int nprocs");
		fprintf(headerfp, "int my_rank, int nprocs");
	}

	if (options->variables_not_global) {
		fprintf(packfp," __ifndef USE_LOCAL_ARRAYS ");
		fprintf(headerfp," __ifndef USE_LOCAL_ARRAYS ");

		for (i=0; i<num_data; i++) {
			char *acc_name = accs[i]->name;

			if(is_access_scalar(accs[i])) continue;

			fprintf(packfp, ",__DECLARATION_OF_%s", acc_name);
			fprintf(headerfp, ",__DECLARATION_OF_%s", acc_name);
		}

		fprintf(packfp, " __endif ");
		fprintf(headerfp, " __endif ");
	}

	for (i=0; i<num_data; i++) {
		char *acc_name = accs[i]->name;

		if(is_access_scalar(accs[i])) continue;

		add_data_dist_parm_decl(packfp, acc_name, prog);
		add_data_dist_parm_decl(headerfp, acc_name, prog);
	}


	fprintf(packfp, ")\n{\n");
	fprintf(headerfp, ");\n");

	if(writeout_func)
		fprintf(packfp, "##ifndef USE_LOCAL_ARRAYS\n");

	char cloog_func_name[256] = "";

	sprintf(cloog_func_name, "%s.cloog",func_name);
	FILE *cloogfp = fopen(cloog_func_name, "w+");
	pluto_gen_cloog_file(cloogfp, data_tile_prog);
	rewind(cloogfp);
	// pluto_gen_cloog_code(compute_task, cloogfp, stdout);
	// rewind(cloogfp);
	generate_declarations(data_tile_prog, packfp);
	pluto_gen_cloog_code(data_tile_prog, -1, -1, cloogfp, packfp);
	fclose(cloogfp);

	if(writeout_func)
		fprintf(packfp, "##endif\n");

	for (i=0; i<data_tile_prog->nstmts; i++) {
		fprintf(packfp, "#undef S%d\n", i+1);
	}
	fprintf(packfp, "}\n\n");

	fclose(packfp);

}

PlutoConstraints *get_data_tile_alloc_constraints(Array *arr, int loop_num, int src_copy_level, char *acc_name,
		int acc_nrows, PlutoProg *prog){

	assert(arr->parmetric_domain[loop_num]!=NULL);
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(arr->parmetric_domain[loop_num], src_copy_level, acc_name, prog);

	assert(data_tiles->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	data_tile_domain_remove_copy_level_dims(data_tiles, src_copy_level, prog);
	return data_tiles;
}

void gen_data_tile_alloc_cloog_code(PlutoProg *prog, int loop_num, Stmt **loop_stmts, int nstmts,
		int src_copy_level, FILE *outfp, FILE *headerfp)
{
	int i;

	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_accs(loop_stmts, nstmts, &num_data);

	PlutoProg *data_tile_prog = NULL;
	data_tile_prog = pluto_prog_alloc();

	data_tile_prog_add_parm(data_tile_prog, prog, src_copy_level);

	char *acc_name;
	int acc_nrows = 0;
	for(i=0;i<num_data;i++){
		acc_name = accs[i]->name;
		acc_nrows = accs[i]->mat->nrows;

		if(is_access_scalar(accs[i])) continue;

		Array *arr = pluto_get_corrs_array(acc_name, prog);
		assert(arr!=NULL);

		PlutoConstraints *data_tiles = get_data_tile_alloc_constraints(arr,loop_num, src_copy_level, acc_name,
				acc_nrows, prog);

		char* alloc_stmt_text = pluto_dist_malloc_stmt_text(acc_name, prog, 1);

		PlutoMatrix *trans = get_data_tile_trans_func(src_copy_level, acc_nrows, prog);

		char **iters = get_data_tile_iterators(arr, acc_nrows);

		pluto_add_stmt(data_tile_prog, data_tiles, trans, iters, alloc_stmt_text
				, DATA_DIST_MANG);

		pluto_matrix_free(trans);
		free(iters);
	}

	pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

	IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

	int temp = options->variables_not_global;
	options->variables_not_global = 0;
	gen_func_cloog_code("data_tile_alloc", headerfp, loop_num, src_copy_level, num_data, accs, data_tile_prog, prog, 0);
	options->variables_not_global = temp;

	pluto_prog_free(data_tile_prog);

	free(accs);
}

PlutoConstraints *get_data_tile_ref_count_update_constraints(Array *arr, int loop_num, int src_copy_level, char *acc_name,
		int acc_nrows, PlutoProg *prog){

	assert(arr->parmetric_domain[loop_num]!=NULL);
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(arr->parmetric_domain[loop_num], src_copy_level, acc_name, prog);

	assert(data_tiles->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	data_tile_domain_remove_copy_level_dims(data_tiles, src_copy_level, prog);
	return data_tiles;
}

PlutoConstraints *get_data_tile_ref_count_init_constraints(Array *arr, int loop_num, int src_copy_level, char *acc_name,
		int acc_nrows, PlutoProg *prog){

    assert(arr->parmetric_domain[loop_num]!=NULL);
	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(arr->parmetric_domain[loop_num], src_copy_level, acc_name, prog);

	assert(data_tiles->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	data_tile_domain_remove_copy_level_dims(data_tiles, src_copy_level, prog);
	return data_tiles;
}

void gen_data_tile_ref_count_update_cloog_code(PlutoProg *prog, int loop_num, Stmt **loop_stmts, int nstmts,
		int src_copy_level, FILE *outfp, FILE *headerfp)
{
	int i;

	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_accs(loop_stmts, nstmts, &num_data);

	PlutoProg *data_tile_prog = NULL;
	data_tile_prog = pluto_prog_alloc();

	data_tile_prog_add_parm(data_tile_prog, prog, src_copy_level);

	char *acc_name;
	int acc_nrows = 0;
	for(i=0;i<num_data;i++){
		acc_name = accs[i]->name;
		acc_nrows = accs[i]->mat->nrows;

		Array *arr = pluto_get_corrs_array(acc_name, prog);
		assert(arr!=NULL);

		PlutoConstraints *data_tiles = get_data_tile_ref_count_update_constraints(arr,loop_num,
				src_copy_level, acc_name, acc_nrows, prog);

		char* ref_count_update_stmt_text = get_data_tile_ref_count_stmt_text(arr, prog);

		PlutoMatrix *trans = get_data_tile_trans_func(src_copy_level, acc_nrows, prog);

		char **iters = get_data_tile_iterators(arr, acc_nrows);

		pluto_add_stmt(data_tile_prog, data_tiles, trans, iters, ref_count_update_stmt_text
				, DATA_DIST_MANG);

		pluto_matrix_free(trans);
		free(iters);
	}

	pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

	IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

	int temp = options->variables_not_global;
	options->variables_not_global = 0;
	gen_func_cloog_code("data_tile_ref_count_update", headerfp, loop_num, src_copy_level, num_data, accs, data_tile_prog, prog, 0);
	options->variables_not_global = temp;

	pluto_prog_free(data_tile_prog);

	free(accs);
}

PlutoConstraints *get_read_in_constraints(struct stmt_access_pair **racc_stmts, int num_accs,
		int src_copy_level, int loop_num,  PlutoProg *prog){

	//char *acc_name = racc_stmts[0]->acc->name;
	//	Array *arr = pluto_get_corrs_array(acc_name, prog);


	PlutoConstraints *read_in = NULL;
	//    read_in = arr->parmetric_domain[loop_num];

	int i;
	for (i=0; i<num_accs; i++)  {
		PlutoConstraints *read_in_one;
		//write_out_one = compute_last_writes(wacc_stmts[k], src_copy_level, prog);
		read_in_one = compute_read_in(racc_stmts[i], src_copy_level, prog);
		if (read_in == NULL) read_in = pluto_constraints_dup(read_in_one);
		else read_in = pluto_constraints_unionize(read_in , read_in_one);
		pluto_constraints_free(read_in_one);
	}

	return read_in;
}

PlutoConstraints *get_read_in_data_tile_alloc_constraints(PlutoConstraints *read_in,
		Array *arr, int loop_num, int src_copy_level, char *acc_name,
		int acc_nrows, PlutoProg *prog){

	PlutoConstraints *data_tiles =
			pluto_dist_get_required_data_tiles(read_in, src_copy_level, acc_name, prog);

	assert(data_tiles->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	data_tile_domain_remove_copy_level_dims(data_tiles, src_copy_level, prog);
	return data_tiles;
}

void gen_read_in_data_tile_alloc_cloog_code(struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoConstraints *read_in, PlutoProg *prog, int loop_num, int src_copy_level, FILE *headerfp)
{


	char *acc_name;
	int acc_nrows = 0;

	PlutoProg *data_tile_prog = NULL;
	data_tile_prog = pluto_prog_alloc();

	data_tile_prog->arrays = prog->arrays;
	data_tile_prog->narrays = prog->narrays;

	data_tile_prog_add_parm(data_tile_prog, prog, src_copy_level);

	PlutoAccess *acc = acc_stmts[0]->acc;
	acc_name = acc->name;
	acc_nrows = acc->mat->nrows;

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);

	PlutoConstraints *data_tiles = get_read_in_data_tile_alloc_constraints(read_in,
			arr,loop_num, src_copy_level, acc_name, acc_nrows, prog);

	char* alloc_stmt_text = pluto_dist_malloc_stmt_text(acc_name, prog, 1);

	PlutoMatrix *trans = get_data_tile_trans_func(src_copy_level, acc_nrows, prog);

	char **iters = get_data_tile_iterators(arr, acc_nrows);

	pluto_add_stmt(data_tile_prog, data_tiles, trans, iters, alloc_stmt_text
			, DATA_DIST_INIT);

	pluto_matrix_free(trans);
	free(iters);

	pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

	IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

	char func_name[512] = "read_in_data_tile_alloc_";
	sprintf(func_name+strlen(func_name), "%s", acc_name);

	int temp = options->variables_not_global;
	options->variables_not_global = 0;
	gen_func_cloog_code(func_name, headerfp, loop_num, src_copy_level, 1 , &acc, data_tile_prog, prog, 0);
	options->variables_not_global = temp;

	pluto_prog_free(data_tile_prog);

}

PlutoConstraints *get_read_in_copy_constraints(PlutoConstraints *read_in,
		Array *arr, int loop_num, int src_copy_level, char *acc_name,
		int acc_nrows, PlutoProg *prog){

	PlutoConstraints *data_tiles = pluto_constraints_dup(read_in);

	data_tile_domain_remove_copy_level_dims(data_tiles, src_copy_level, prog);
	return data_tiles;
}

char* get_read_in_copy_stmt_text(Array *arr, PlutoProg *prog){

	char* init_stmt_text = (char *)malloc(2048 * sizeof(char));

	init_stmt_text[0] = '\0';

	sprintf(init_stmt_text+strlen(init_stmt_text), " __ifndef USE_LOCAL_ARRAYS");
	sprintf(init_stmt_text+strlen(init_stmt_text), " %s ",  pluto_dist_copy_init_text(arr, prog));
	sprintf(init_stmt_text+strlen(init_stmt_text), " __else");
	char *text= pluto_dist_read_init_text(arr->text, prog);
	char *mod_stmt_text = pluto_dist_modify_stmt_text(text,1, prog);

	if(mod_stmt_text[strlen(mod_stmt_text) -1] == '\n')
		mod_stmt_text[strlen(mod_stmt_text) -1] = 0;

	sprintf(init_stmt_text+strlen(init_stmt_text), " %s ", mod_stmt_text );
	sprintf(init_stmt_text+strlen(init_stmt_text), "__endif");

	return init_stmt_text;
}

void gen_read_in_copy_cloog_code(struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoConstraints *read_in, PlutoProg *prog, int loop_num,  int src_copy_level, FILE *headerfp)
{


	char *acc_name;
	int acc_nrows = 0;

	PlutoProg *data_tile_prog = NULL;
	data_tile_prog = pluto_prog_alloc();

	data_tile_prog->arrays = prog->arrays;
	data_tile_prog->narrays = prog->narrays;

	data_tile_prog_add_parm(data_tile_prog, prog, src_copy_level);

	PlutoAccess *acc = acc_stmts[0]->acc;
	acc_name = acc->name;
	acc_nrows = acc->mat->nrows;

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);

	PlutoConstraints *data_tiles = get_read_in_copy_constraints(read_in,
			arr,loop_num, src_copy_level, acc_name, acc_nrows, prog);

	char* copy_stmt_text = get_read_in_copy_stmt_text(arr, prog);

	PlutoMatrix *trans = get_data_tile_trans_func(src_copy_level, acc_nrows, prog);

	char **iters = get_data_scan_iterators(arr, acc_nrows);

	pluto_add_stmt(data_tile_prog, data_tiles, trans, iters, copy_stmt_text
			, DATA_DIST_COPY);

	pluto_matrix_free(trans);
	free(iters);

	pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

	IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

	char func_name[512] = "read_in_copy_";
	sprintf(func_name+strlen(func_name), "%s", acc_name);

	gen_func_cloog_code(func_name, headerfp, loop_num, src_copy_level, 1 , &acc, data_tile_prog, prog, 0);

	pluto_prog_free(data_tile_prog);

}
void gen_data_tile_ref_count_init_cloog_code(PlutoProg *prog, int loop_num, Stmt **loop_stmts, int nstmts,
		int src_copy_level, FILE *headerfp)
{
	int i;

	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_accs(loop_stmts, nstmts, &num_data);

	char *acc_name;
	int acc_nrows = 0;
	for(i=0;i<num_data;i++){


		if(is_access_scalar(accs[i])) continue;

		PlutoProg *data_tile_prog = NULL;
		data_tile_prog = pluto_prog_alloc();

		data_tile_prog->arrays = prog->arrays;
		data_tile_prog->narrays = prog->narrays;

		data_tile_prog_add_parm(data_tile_prog, prog, src_copy_level);

		acc_name = accs[i]->name;
		acc_nrows = accs[i]->mat->nrows;

		Array *arr = pluto_get_corrs_array(acc_name, prog);
		assert(arr!=NULL);

		PlutoConstraints *data_tiles = get_data_tile_ref_count_init_constraints(arr,loop_num,
				src_copy_level, acc_name, acc_nrows, prog);

		char* ref_count_update_stmt_text = get_data_tile_ref_count_init_stmt_text(arr, src_copy_level, loop_num, prog);

		PlutoMatrix *trans = get_data_tile_trans_func(src_copy_level, acc_nrows, prog);

		char **iters = get_data_tile_iterators(arr, acc_nrows);

		pluto_add_stmt(data_tile_prog, data_tiles, trans, iters, ref_count_update_stmt_text
				, DATA_DIST_INIT);

		pluto_matrix_free(trans);
		free(iters);

		pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

		IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

		char func_name[512] = "data_tile_ref_count_init_";
		sprintf(func_name+strlen(func_name), "%s", acc_name);

		int temp = options->variables_not_global;
		options->variables_not_global = 0;
		gen_func_cloog_code(func_name, headerfp, loop_num, src_copy_level, 1 , &accs[i], data_tile_prog, prog, 0);
		options->variables_not_global = temp;

		pluto_prog_free(data_tile_prog);

	}

	free(accs);
}

PlutoConstraints *get_data_tile_copy_back_constraints(struct stmt_access_pair **wacc_stmts, int num_accs,
		int src_copy_level, PlutoProg *prog){

	PlutoConstraints *all_data= NULL;
	int i;

	for (i=0; i<num_accs; i++)  {
		Stmt *stmt = wacc_stmts[i]->stmt;
		PlutoAccess *acc = wacc_stmts[i]->acc;

		PlutoConstraints *srcdomain = pluto_get_new_domain(stmt);

		PlutoConstraints *data_in_one =  pluto_compute_region_data(stmt,
				srcdomain, acc, src_copy_level, prog);

		if (all_data == NULL) all_data = pluto_constraints_dup(data_in_one);
		else all_data = pluto_constraints_unionize(all_data, data_in_one);
		pluto_constraints_free(data_in_one);
	}

	assert(all_data->ncols == src_copy_level + wacc_stmts[0]->acc->mat->nrows + prog->npar + 1);

	data_tile_domain_remove_copy_level_dims(all_data, src_copy_level, prog);
	return all_data;
}

void gen_data_tile_copy_back_cloog_code(PlutoProg *prog, int loop_num,struct stmt_access_pair **wacc_stmts, int num_accs,
		int src_copy_level, FILE *headerfp)
{

	int num_data = 1;
	char *acc_name = wacc_stmts[0]->acc->name;

	int acc_nrows = wacc_stmts[0]->acc->mat->nrows;

	PlutoProg *data_tile_prog = NULL;
	data_tile_prog = pluto_prog_alloc();

	data_tile_prog->arrays = prog->arrays;
	data_tile_prog->narrays = prog->narrays;

	data_tile_prog_add_parm(data_tile_prog, prog, src_copy_level);

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(arr!=NULL);

	PlutoConstraints *data_tiles = get_data_tile_copy_back_constraints(wacc_stmts, num_accs,
			src_copy_level, prog);

	char *copy_back_stmt_text = pluto_dist_copy_back_stmt_text(arr, prog);

	PlutoMatrix *trans = get_data_tile_trans_func(src_copy_level, acc_nrows, prog);


	pluto_add_stmt(data_tile_prog, data_tiles, trans, arr->iterators, copy_back_stmt_text
			, DATA_DIST_COPY);

	pluto_matrix_free(trans);

	pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

	IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

	char func_name[512] = "data_tile_copy_back_";
	sprintf(func_name+strlen(func_name), "%s", acc_name);

	gen_func_cloog_code(func_name, headerfp, loop_num, src_copy_level, num_data, &(wacc_stmts[0]->acc) , data_tile_prog, prog, 1);

	pluto_prog_free(data_tile_prog);


}
void gen_init_tasks_cloog_code(PlutoProg *prog, Stmt ***tasks_stmts, int nloops,
		char *tasks_loops_decl, FILE *headerfp)
{
	int i, j;

	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_all_distinct_vars(prog->stmts, prog->nstmts, &num_data);

	PlutoProg *init_tasks = NULL;
	init_tasks = pluto_prog_alloc();

	for (i=0; i<prog->npar; i++) {
		pluto_prog_add_param(init_tasks, prog->params[i], init_tasks->npar);
	}

	int num_stmts_per_loop = (options->dynschedule_graph) ? NUM_GEN_GRAPH_STMTS : NUM_GEN_TASKS_STMTS;
	for (j=0; j<num_stmts_per_loop; j++) {
		for (i=0; i<nloops; i++) {
			pluto_add_given_stmt(init_tasks, tasks_stmts[i][j]);
		}
	}

	pluto_pad_stmt_transformations(init_tasks);
	IF_DEBUG(pluto_prog_print(stdout, init_tasks));
	if (init_tasks->nstmts >= 1) {
		assert(init_tasks->stmts[0]->trans->nrows == init_tasks->num_hyperplanes);
	}

	pluto_separate_stmts(init_tasks, init_tasks->stmts, init_tasks->nstmts, 0, 0);

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	if (options->dynschedule) {
		fprintf(packfp, "void init_tasks(");
		fprintf(headerfp, "void init_tasks(");
	}
	else {
		fprintf(packfp, "void create_dag(");
		fprintf(headerfp, "void create_dag(");
	}
	if (options->distmem) {
		fprintf(packfp, "int my_rank,int nprocs,int num_threads%s,long int *p_num_tasks_to_execute,long int *p_num_tasks_to_unpack", tasks_loops_decl);
		fprintf(headerfp, "int my_rank,int nprocs,int num_threads%s,long int *p_num_tasks_to_execute,long int *p_num_tasks_to_unpack", tasks_loops_decl);
	}
	else if (options->dynschedule) {
		fprintf(packfp, "int num_threads,long int *p_num_tasks_to_execute%s", tasks_loops_decl);
		fprintf(headerfp, "int num_threads,long int *p_num_tasks_to_execute%s", tasks_loops_decl);
	}
	else { // options->dynschedule_graph
		fprintf(packfp, "tbb::flow::graph *__graph, tbb::flow::broadcast_node<tbb::flow::continue_msg> *__start%s", tasks_loops_decl);
		fprintf(headerfp, "tbb::flow::graph *__graph, tbb::flow::broadcast_node<tbb::flow::continue_msg> *__start%s", tasks_loops_decl);
		if (options->variables_not_global) {
			for (i=0; i<num_data; i++) {
				fprintf(packfp, ", __DECLARATION_OF_%s", accs[i]->name);
				fprintf(headerfp, ", __DECLARATION_OF_%s", accs[i]->name);
			}
		}
	}
	fprintf(packfp, ")\n{\n");
	fprintf(headerfp, ");\n");

	if (options->distmem) {
		fprintf(packfp, "int remote_dep_tasks, local_dep_tasks, affinity;\n");
	}
	else if (options->dynschedule) {
		fprintf(packfp, "int local_dep_tasks, affinity;\n");
	}
	if (options->dynschedule) {
		fprintf(packfp, "long int _num_tasks_to_execute = 0;\n");
	}
	if (options->distmem) {
		fprintf(packfp, "long int _num_tasks_to_unpack = 0;\n");
	}

	FILE *cloogfp = fopen("init_tasks.cloog", "w+");
	pluto_gen_cloog_file(cloogfp, init_tasks);
	rewind(cloogfp);
	generate_declarations(init_tasks, packfp);
	pluto_gen_cloog_code(init_tasks, -1, -1, cloogfp, packfp);
	fclose(cloogfp);

	for (i=0; i<init_tasks->nstmts; i++) {
		fprintf(packfp, "#undef S%d\n", i+1);
	}

	if (options->dynschedule) {
		fprintf(packfp, "*p_num_tasks_to_execute = _num_tasks_to_execute;\n");
	}
	if (options->distmem) {
		fprintf(packfp, "*p_num_tasks_to_unpack = _num_tasks_to_unpack;\n");
	}
	fprintf(packfp, "}\n\n");

	fclose(packfp);

	pluto_prog_free(init_tasks);

	free(accs);
}

void gen_pack_send_text_code(PlutoProg *prog, Stmt ***copy_comm_stmts, struct stmt_access_pair ***wacc_stmts,
		int loop_num, int num_data, int num_comm_stmts, int copy_level, FILE *headerfp, Ploop *loops)
{
	int i, j;

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	fprintf(packfp, "void pack_send_%d(", loop_num);
	fprintf(headerfp, "void pack_send_%d(", loop_num);
	fprintf(packfp, "int my_rank,int nprocs");
	fprintf(headerfp, "int my_rank,int nprocs");
	for (i=0; i<copy_level; i++) {
		fprintf(packfp, ",int t%d", i+1);
		fprintf(headerfp, ",int t%d", i+1);
	}

	if(options->data_dist){
		int num_rw_data;
		PlutoAccess **accs;
		accs = pluto_get_accs(loops->stmts, loops->nstmts, &num_rw_data);

		for (i=0; i<num_rw_data; i++) {
			char *acc_name = accs[i]->name;
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ",__DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ",__DECLARATION_OF_%s", acc_name);
			}
			if(options->data_dist){
				add_data_dist_parm_decl(packfp, acc_name, prog);
				add_data_dist_parm_decl(headerfp, acc_name, prog);
			}

			/*
			fprintf(packfp, ",double **send_buf_%s", acc_name);
			fprintf(headerfp, ",double **send_buf_%s", acc_name);
			fprintf(packfp, ",MPI_Request *send_reqs_%s", acc_name);
			fprintf(headerfp, ",MPI_Request *send_reqs_%s", acc_name);
			fprintf(packfp, ",MPI_Status *send_stats_%s", acc_name);
			fprintf(headerfp, ",MPI_Status *send_stats_%s", acc_name);
			 */
			fprintf(packfp, ",std::vector<double *>& send_buf_%s", acc_name);
			fprintf(headerfp, ",std::vector<double *>& send_buf_%s", acc_name);
			fprintf(packfp, ",std::vector<MPI_Request>& send_reqs_%s", acc_name);
			fprintf(headerfp, ",std::vector<MPI_Request>& send_reqs_%s", acc_name);
		}

		free(accs);
	}
	else {

		for (i=0; i<num_data; i++) {
			char *acc_name = wacc_stmts[i][0]->acc->name;
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ",__DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ",__DECLARATION_OF_%s", acc_name);
			}
			fprintf(packfp, ",std::vector<double *>& send_buf_%s", acc_name);
			fprintf(headerfp, ",std::vector<double *>& send_buf_%s", acc_name);
			fprintf(packfp, ",std::vector<MPI_Request>& send_reqs_%s", acc_name);
			fprintf(headerfp, ",std::vector<MPI_Request>& send_reqs_%s", acc_name);
		}
	}
	if (options->timereport) {
		fprintf(packfp, ",double *p_total_count");
		fprintf(headerfp, ",double *p_total_count");
	}
	fprintf(packfp, ")\n{\n");
	fprintf(headerfp, ");\n");

	fprintf(packfp, "int __p, receiver_list[nprocs];\n");
	if (options->commopt_fop || options->commopt_foifi) {
		for (i=0; i<num_data; i++) {
			char *acc_name = wacc_stmts[i][0]->acc->name;
			fprintf(packfp, "int send_counts_%s[nprocs];\n", acc_name);
			fprintf(packfp, "for (__p=0;__p<nprocs;__p++) send_counts_%s[__p] = 0;\n", acc_name);
			fprintf(packfp, "int index_%s[nprocs];\n", acc_name);
			fprintf(packfp, "int current_index_%s = 0;\n", acc_name);
		}
		if (options->fop_unicast_runtime) {
			fprintf(packfp, "int distinct_recv;\n");
		}
	}
	else {
		for (i=0; i<num_data; i++) {
			char *acc_name = wacc_stmts[i][0]->acc->name;
			fprintf(packfp, "int send_count_%s = 0;\n", acc_name);
			fprintf(packfp, "int index_%s;\n", acc_name);
		}
		fprintf(packfp, "int _i;\n");
	}
	fprintf(packfp, "int flag;\n");
	if (options->timereport)
		fprintf(packfp, "double __total_count = 0;\n");

	for (i=0; i<num_data; i++) {
		// assuming last statement is unpack statement
		for (j=0; j<num_comm_stmts-1; j++) {
			fprintf(packfp, "%s;\n", copy_comm_stmts[i][j]->text);
		}
	}

	if (options->timereport)
		fprintf(packfp, "*p_total_count += __total_count;\n");

	fprintf(packfp, "}\n\n");

	fclose(packfp);
}

void gen_unpack_text_code(PlutoProg *prog, Stmt ***copy_comm_stmts, struct stmt_access_pair ***wacc_stmts,
		int loop_num, int num_data, int num_comm_stmts, int copy_level, FILE *headerfp)
{
	int i;

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	for (i=0; i<num_data; i++) {
		char *acc_name = wacc_stmts[i][0]->acc->name;
		fprintf(packfp, "void unpack_buf_%s_%d(", acc_name, loop_num);
		fprintf(headerfp, "void unpack_buf_%s_%d(", acc_name, loop_num);
		fprintf(packfp, "int my_rank,int nprocs");
		fprintf(headerfp, "int my_rank,int nprocs");
		if (options->variables_not_global && !options->data_dist) {
			fprintf(packfp, ",__DECLARATION_OF_%s", acc_name);
			fprintf(headerfp, ",__DECLARATION_OF_%s", acc_name);
		}
		if(options->data_dist){
			add_data_dist_parm_decl(packfp, acc_name, prog);
			add_data_dist_parm_decl(headerfp, acc_name, prog);
		}
		fprintf(packfp, ",double *recv_buf_%s)\n{\n", acc_name);
		fprintf(headerfp, ",double *recv_buf_%s);\n", acc_name);

		if (options->commopt_fop || options->commopt_foifi) {
			fprintf(packfp, "int curr_displs_%s;\n", acc_name);
			if (options->fop_unicast_runtime) {
				fprintf(packfp, "int distinct_recv;\n");
			}
		}
		else {
			fprintf(packfp, "int displs_%s;\n", acc_name);
		}

		fprintf(packfp, "%s;\n", copy_comm_stmts[i][num_comm_stmts-1]->text);

		fprintf(packfp, "}\n\n");
	}

	fclose(packfp);
}

void gen_write_out_cloog_code(PlutoProg *prog, PlutoProg *write_out_prog, FILE *headerfp)
{
	int i;

	int num_data;
	PlutoAccess **waccs;
	waccs = pluto_get_all_distinct_write_vars(prog, &num_data);

	pluto_mark_statements(write_out_prog);
	pluto_detect_scalar_dimensions(write_out_prog);
	pluto_pad_stmt_transformations(write_out_prog);
	IF_DEBUG(pluto_prog_print(stdout, write_out_prog));
	if (write_out_prog->nstmts >= 1) {
		assert(write_out_prog->stmts[0]->trans->nrows == write_out_prog->num_hyperplanes);
	}

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	fprintf(packfp, "void write_out(");
	fprintf(headerfp, "void write_out(");
	fprintf(packfp, "int my_rank,int nprocs");
	fprintf(headerfp, "int my_rank,int nprocs");
	for (i=0; i<num_data; i++) {
		char *name = waccs[i]->name;
		fprintf(packfp, ",double *lw_buf_%s", name);
		fprintf(headerfp, ",double *lw_buf_%s", name);
		fprintf(packfp, ",double *lw_recv_buf_%s", name);
		fprintf(headerfp, ",double *lw_recv_buf_%s", name);
		fprintf(packfp, ",int *displs_lw_%s", name);
		fprintf(headerfp, ",int *displs_lw_%s", name);
	}
	if (options->variables_not_global) {
		for (i=0; i<num_data; i++) {
			fprintf(packfp, ", __DECLARATION_OF_%s", waccs[i]->name);
			fprintf(headerfp, ", __DECLARATION_OF_%s", waccs[i]->name);
		}
	}
	for (i=0; i<num_data; i++) {
		add_data_dist_parm_decl(packfp, waccs[i]->name, prog);
		add_data_dist_parm_decl(headerfp, waccs[i]->name, prog);
	}
	fprintf(packfp, ")\n{\n");
	fprintf(headerfp, ");\n");

	if(options->data_dist)
		fprintf(packfp, "##ifndef USE_LOCAL_ARRAYS\n");
	//    fprintf(packfp, "double t_writeout_start = 0.0;");
	fprintf(packfp, "int __p, proc;\n");
	fprintf(packfp, "int _lb_dist, _ub_dist");
	for (i=0; i<prog->num_hyperplanes; i++) {
		fprintf(packfp, ", lbd_t%d, ubd_t%d", i+1, i+1);
	}
	fprintf(packfp, ";\n");
	for (i=0; i<num_data; i++) {
		char *name = waccs[i]->name;
		fprintf(packfp, "int lw_recv_counts_%s[nprocs];\n", name);
		fprintf(packfp, "int curr_displs_lw_%s[nprocs];\n", name);
		fprintf(packfp, "int lw_count_%s = 0;\n", name);
	}

	FILE *cloogfp = fopen("write_out.cloog", "w+");
	pluto_gen_cloog_file(cloogfp, write_out_prog);
	rewind(cloogfp);
	generate_declarations(write_out_prog, packfp);
	pluto_gen_cloog_code(write_out_prog, -1, -1, cloogfp, packfp);
	fclose(cloogfp);

	for (i=0; i<write_out_prog->nstmts; i++) {
		fprintf(packfp, "#undef S%d\n", i+1);
	}

	if(options->data_dist)
		fprintf(packfp, "##endif\n");
	fprintf(packfp, "}\n\n");

	fclose(packfp);

	free(waccs);

	pluto_prog_free(write_out_prog);
}

void gen_dynschedule_graph_main_text_code(PlutoProg *prog, Ploop **loops, int nloops,
		int copy_level[nloops], FILE *outfp)
{
	int l, i;

	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_all_distinct_vars(prog->stmts, prog->nstmts, &num_data);

	char *tasks_loops_args = malloc(256);
	strcpy(tasks_loops_args, "");
	for (l=0; l<nloops; l++) {
		sprintf(tasks_loops_args+strlen(tasks_loops_args), ",tasks_loop%d", l);
	}

	fprintf(outfp, "\n\n\n\t/* Start of dynamic scheduling code */\n");
	if (options->timereport)
		fprintf(outfp, "\tIF_TIME(t_create_dag = rtclock());\n");
	fprintf(outfp, "\tcreate_dag(&__graph,&__start%s", tasks_loops_args);
	if (options->variables_not_global) {
		for (i=0; i<num_data; i++) {
			fprintf(outfp, ", %s", accs[i]->name);
		}
	}
	fprintf(outfp, ");\n");
	if (options->timereport)
		fprintf(outfp, "\tIF_TIME(t_create_dag = rtclock() - t_create_dag);\n");
	fprintf(outfp, "\t__start.try_put(tbb::flow::continue_msg());\n");
	fprintf(outfp, "\t__graph.wait_for_all();\n");
	fprintf(outfp, "\t/* End of dynamic scheduling code */\n\n\n");

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	fprintf(packfp, "void __Task::operator()(tbb::flow::continue_msg) const\n");
	fprintf(packfp, "{\n");
	fprintf(packfp, "\tswitch(loop){\n");
	for (l=0; l<nloops; l++) {
		fprintf(packfp, "\t\tcase %d:\n", l);

		int num_rw_data;
		PlutoAccess **laccs;
		laccs = pluto_get_all_distinct_vars(loops[l]->stmts, loops[l]->nstmts, &num_rw_data);

		fprintf(packfp, "\t\t\tcompute_task_%d(", l);
		fprintf(packfp, "id1");
		for (i=1; i<copy_level[l]; i++) {
			fprintf(packfp, ",id%d", i+1);
		}
		if (options->variables_not_global) {
			for (i=0; i<num_rw_data; i++) {
				char *acc_name = laccs[i]->name;
				fprintf(packfp, ",%s", acc_name);
			}
		}
		fprintf(packfp, ");\n");

		free(laccs);
	}
	fprintf(packfp, "\t\tdefault:\n");
	fprintf(packfp, "\t\t\t;\n");
	fprintf(packfp, "\t}\n");
	fprintf(packfp, "}\n");

	fclose(packfp);

	free(tasks_loops_args);
	free(accs);
}

void gen_dynschedule_main_text_code(PlutoProg *prog, Ploop **loops, int nloops,
		int copy_level[nloops], FILE *outfp)
{
	fprintf(outfp, "\n\n\n/*Start of dynamic scheduling code*/\n\n\n");

	fprintf(outfp, "##if defined(__DYNSCHEDULER_DEBUG_PRINT) || defined(__DYNSCHEDULER_MORE_DEBUG_PRINT) || defined(__DEBUG_FLUSH)\n");
	fprintf(outfp, "\tchar *debug_print_name = (char *)malloc(128*sizeof(char));\n");
	fprintf(outfp, "\tstrcpy(debug_print_name, \"debug_print\");\n");
	fprintf(outfp, "\tsprintf(debug_print_name+strlen(debug_print_name), \"_node%%d\", my_rank);\n");
	fprintf(outfp, "\t__debug_print_fp = fopen(debug_print_name, \"w\");\n");
	fprintf(outfp, "\tfree(debug_print_name);\n");
	fprintf(outfp, "##endif\n");

	int i, j, k, l, ii;

	int num_data;
	PlutoAccess **waccs;
	waccs = pluto_get_all_distinct_write_vars(prog, &num_data);
	char *acc_names[num_data];
	char *max_elements_name[num_data];
	char *data_tag[num_data];
	for (i=0; i<num_data; i++) {
		acc_names[i] = waccs[i]->name;
		max_elements_name[i] = concat("max_num_elements_", acc_names[i]);
		if (options->distmem) {
			data_tag[i] = malloc(10);
			get_data_tag(prog, waccs[i], &data_tag[i]);
		} else {
			data_tag[i] = NULL;
		}
	}

	int num_write_data[nloops];
	int *num_stmts_per_wacc[nloops]; // indexed by data variable
	struct stmt_access_pair ***wacc_stmts[nloops]; // indexed by data variable
	for (l=0; l<nloops; l++) {
		wacc_stmts[l] = get_write_access_with_stmts(loops[l]->stmts,
				loops[l]->nstmts, &num_write_data[l], &num_stmts_per_wacc[l]);
	}

	int max_copy_level = copy_level[0];
	for (l=1; l<nloops; l++) {
		if (copy_level[l] > max_copy_level) max_copy_level = copy_level[l];
	}

	char *params = malloc(256);
	strcpy(params, "");
	if (prog->npar>=1) {
		sprintf(params+strlen(params), "%s", prog->params[0]);
		for (i=1; i<prog->npar; i++) {
			sprintf(params+strlen(params), ",%s", prog->params[i]);
		}
	}
	// my_rank is already a param

	char *tasks_loops_args = NULL;
	tasks_loops_args = malloc(256);
	strcpy(tasks_loops_args, "");
	for (l=0; l<nloops; l++) {
		sprintf(tasks_loops_args+strlen(tasks_loops_args), ",tasks_loop%d", l);
	}

	if (options->timereport)
		fprintf(outfp, "\tIF_TIME(t_tasks_create = rtclock());\n");
	if (options->distmem) {
		fprintf(outfp, "\tinit_tasks(my_rank,nprocs,_num_threads%s,&_num_tasks_to_execute,&_num_tasks_to_unpack);\n", tasks_loops_args);
	}
	else {
		fprintf(outfp, "\tinit_tasks(_num_threads,&_num_tasks_to_execute%s);\n", tasks_loops_args);
	}
	if (options->timereport)
		fprintf(outfp, "\tIF_TIME(t_tasks_create = rtclock() - t_tasks_create);\n");
	fprintf(outfp, "##ifdef __DYNSCHEDULER_MORE_DEBUG_PRINT\n");
	for (i=0; i<num_data; i++) {
		for (l=0; l<nloops; l++) {
			fprintf(outfp, "\tstd::vector<int> prev_firing_count_%s_%d;\n", waccs[i]->name, l);
		}
	}
	fprintf(outfp, 	"\tfor (__p=0; __p<__MAX_NUM_RECVS; __p++) {\n");
	for (i=0; i<num_data; i++) {
		for (l=0; l<nloops; l++) {
			fprintf(outfp, "\t\tprev_firing_count_%s_%d.push_back(-1);\n", waccs[i]->name, l);
		}
	}
	fprintf(outfp, "\t}\n");
	fprintf(outfp, "##endif\n");

	if (options->distmem) {
		if (options->timereport)
			fprintf(outfp, "\tIF_TIME(t_unpack_start = rtclock());\n");
		fprintf(outfp, "\tfor (__p=0; __p<__MAX_NUM_RECVS; __p++) {\n");
		fprintf(outfp, "\t\tint __error;\n");
		for (i=0; i<num_data; i++) {
			fprintf(outfp,
					"\t\t__error = MPI_Irecv(recv_buf_%s[__p], %s, MPI_DOUBLE,\
                MPI_ANY_SOURCE, %s, MPI_COMM_WORLD, &recv_reqs_%s[__p]);\n",
                acc_names[i], max_elements_name[i], data_tag[i], acc_names[i]);
			fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
			fprintf(outfp, "\t\tif (__error) {\n");
			fprintf(outfp, "\t\t\tchar *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len;\n");
			fprintf(outfp, "\t\t\tMPI_Error_string(__error, msg, &len);\n");
			fprintf(outfp, "\t\t\tfprintf(__debug_print_fp, \"node %%d %s recv buf %%d error %%d %%s\\n\", my_rank, __p, __error, msg);\n", acc_names[i]);
			fprintf(outfp, "\t\t\tfree(msg);\n");
			fprintf(outfp, "\t\t}\n");
			fprintf(outfp, "##endif\n");
		}
		fprintf(outfp, "\t}\n");
		if (options->timereport)
			fprintf(outfp, "\tIF_TIME(t_unpack += rtclock() - t_unpack_start);\n");
		fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
		fprintf(outfp, "\tfprintf(__debug_print_fp, \"node %%d num tasks to execute %%ld num tasks to unpack %%ld max recv buffers %%d");
		for (i=0; i<num_data; i++) {
			fprintf(outfp, " %s max elements %%d", waccs[i]->name);
		}
		fprintf(outfp, "\\n\", my_rank, _num_tasks_to_execute, _num_tasks_to_unpack, __MAX_NUM_RECVS");
		for (i=0; i<num_data; i++) {
			fprintf(outfp, ", max_num_elements_%s", waccs[i]->name);
		}
		fprintf(outfp, ");\n");
		fprintf(outfp, "\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
		fprintf(outfp, "##endif\n");
	}

	fprintf(outfp, "\n##pragma omp parallel private(__tid)\n{\n");
	fprintf(outfp, "\t__tid = omp_get_thread_num();\n");
	fprintf(outfp, "\tint __flag;\n");
	fprintf(outfp, "\tint __loop;\n");
	fprintf(outfp, "\tint task_id1");
	for (i=1; i<max_copy_level; i++) {
		fprintf(outfp, ",task_id%d", i+1);
	}
	fprintf(outfp, ";\n");
	fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
	fprintf(outfp, "\tfprintf(__debug_print_fp, \"node %%d thread %%d\\n\", my_rank, __tid);\n");
	fprintf(outfp, "\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
	fprintf(outfp, "##endif\n");
	if (options->timereport)
		fprintf(outfp, "\tIF_TIME(t_wait_comp_start[__tid] = rtclock());\n");

	if (options->distmem) {
		fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
		fprintf(outfp, "\tint __busy_wait_count = 0;\n");
		fprintf(outfp, "##endif\n");
		fprintf(outfp, "\n##pragma omp single nowait\n{ /* receive-unpack */\n");
		if (options->timereport) {
			fprintf(outfp, "\t__tid_receiver = __tid;\n");
			fprintf(outfp, "\tIF_TIME(t_wait_unpack_start = rtclock());\n");
		}
		fprintf(outfp, "\twhile(_num_tasks_to_unpack > 0) {\n");
		fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
		fprintf(outfp, "##pragma omp atomic\n");
		fprintf(outfp, "\t\t__busy_wait_count++;\n");
		fprintf(outfp, "##endif\n");
		fprintf(outfp, "\t\tfor (__p=0; (__p<recv_buf_%s.size())", acc_names[0]);
		for (i=1; i<num_data; i++) {
			fprintf(outfp, 	" || (__p<recv_buf_%s.size())", acc_names[i]);
		}
		fprintf(outfp, 	"; __p++) {\n");
		for (i=0; i<num_data; i++) {
			fprintf(outfp, "\t\tif (__p<recv_buf_%s.size()) {\n", acc_names[i]);
			fprintf(outfp, "\t\t\tif (recv_reqs_%s[__p] != MPI_REQUEST_NULL) { \n", acc_names[i]);
			fprintf(outfp, "##ifndef __DYNSCHEDULER_DEBUG_PRINT\n");
			fprintf(outfp, "\t\t\t\tMPI_Test(&recv_reqs_%s[__p],&recv_completed_%s[__p],MPI_STATUS_IGNORE);\n",
					acc_names[i], acc_names[i]);
			fprintf(outfp, "##else\n");
			fprintf(outfp, "\t\t\t\tMPI_Status recv_stats_%s;\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\trecv_stats_%s.MPI_ERROR = MPI_Test(&recv_reqs_%s[__p],&recv_completed_%s[__p],&recv_stats_%s);\n",
					acc_names[i], acc_names[i], acc_names[i], acc_names[i]);
			fprintf(outfp, "\t\t\t\tif (recv_completed_%s[__p]) { \n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\t__loop = (int)recv_buf_%s[__p][0];\n", acc_names[i]);
			for (j=0; j<max_copy_level; j++) {
				fprintf(outfp, "\t\t\t\t\ttask_id%d = (int)recv_buf_%s[__p][%d];\n", j+1, acc_names[i], j+1);
			}
			fprintf(outfp, "\t\t\t\t\tint size; \n");
			fprintf(outfp, "\t\t\t\t\tMPI_Get_count(&recv_stats_%s, MPI_DOUBLE, &size); \n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\tif (recv_stats_%s.MPI_ERROR != MPI_SUCCESS) {\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\t\tchar *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len;\n");
			fprintf(outfp, "\t\t\t\t\t\tMPI_Error_string(recv_stats_%s.MPI_ERROR, msg, &len);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\t\tfprintf(__debug_print_fp, \"received %s node %%d from node %%d buf %%d size %%d loop %%d task", acc_names[i]);
			for (j=0; j<max_copy_level; j++) {
				fprintf(outfp, " %%d");
			}
			fprintf(outfp, " error %%d %%s\\n\", my_rank, recv_stats_%s.MPI_SOURCE, __p, size, __loop", acc_names[i]);
			for (j=0; j<max_copy_level; j++) {
				fprintf(outfp, ", task_id%d", j+1);
			}
			fprintf(outfp, ", recv_stats_%s.MPI_ERROR, msg);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\t\tfree(msg);\n");
			fprintf(outfp, "\t\t\t\t\t} else {\n");
			fprintf(outfp, "\t\t\t\t\t\tfprintf(__debug_print_fp, \"received %s node %%d from node %%d buf %%d size %%d loop %%d task", acc_names[i]);
			for (j=0; j<max_copy_level; j++) {
				fprintf(outfp, " %%d");
			}
			fprintf(outfp, "\\n\", my_rank, recv_stats_%s.MPI_SOURCE, __p, size, __loop", acc_names[i]);
			for (j=0; j<max_copy_level; j++) {
				fprintf(outfp, ", task_id%d", j+1);
			}
			fprintf(outfp, ");\n");
			fprintf(outfp, "\t\t\t\t\t}\n");
			fprintf(outfp, "\t\t\t\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
			fprintf(outfp, "\t\t\t\t} \n");
			fprintf(outfp, "##endif\n");
			fprintf(outfp, "\t\t\t}\n");
			fprintf(outfp, "\t\t\tif (recv_completed_%s[__p]) { \n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t__loop = (int)recv_buf_%s[__p][0];\n", acc_names[i]);
			for (j=0; j<max_copy_level; j++) {
				fprintf(outfp, "\t\t\t\ttask_id%d = (int)recv_buf_%s[__p][%d];\n", j+1, acc_names[i], j+1);
			}
			fprintf(outfp, "\t\t\t\tswitch(__loop){\n");
			for (l=0; l<nloops; l++) {
				for (j=0; j<num_write_data[l];j++) {
					if (!strcmp(wacc_stmts[l][j][0]->acc->name, acc_names[i])) {
						Stmt *anchor_stmt = loops[l]->stmts[0];
						int dims[copy_level[l]];
						int num_dims = 0;
						for (k=0; k<copy_level[l]; k++) {
							if (anchor_stmt->hyp_types[k] != H_SCALAR) {
								dims[num_dims++] = k;
							}
						}

						char *indices = malloc(1024);
						strcpy(indices,"[");
						for (ii=0; ii<num_dims; ii++) {
							sprintf(indices+strlen(indices), "(task_id%d-lb_tasks_loop%d_dim%d)",
									dims[ii]+1, l, ii);
							for (k=ii+1; k<num_dims; k++) {
								sprintf(indices+strlen(indices), "*max_num_tasks_loop%d_dim%d",
										l, k);
							}
							strcat(indices,"+");
						}
						strcat(indices,"0]");

						char *args = malloc(512);
						strcpy(args,"");
						sprintf(args+strlen(args), "task_id%d", 1);
						/* make it src_copy_level+1 since we use an extra dimension to separate
						 * statements */
						for (ii=1; ii<copy_level[l]; ii++) {
							sprintf(args+strlen(args), ",task_id%d", ii+1);
						}

						fprintf(outfp, "\t\t\t\t\tcase %d:\n", l);

						fprintf(outfp, "##ifdef __DYNSCHEDULER_MORE_DEBUG_PRINT\n");
						fprintf(outfp, "\t\t\t\t\t\tif (prev_firing_count_%s_%d[__p] != tasks_loop%d%s) {\n", acc_names[i], l, l, indices);
						fprintf(outfp, "\t\t\t\t\t\t\tprev_firing_count_%s_%d[__p] = tasks_loop%d%s;\n", acc_names[i], l, l, indices);
						fprintf(outfp, "\t\t\t\t\t\t\tfprintf(__debug_print_fp, \"waiting to unpack %s node %%d buf %%d loop %%d task", acc_names[i]);
						for (j=0; j<max_copy_level; j++) {
							fprintf(outfp, " %%d");
						}
						fprintf(outfp, " count %%d\\n\", my_rank, __p, __loop");
						for (j=0; j<max_copy_level; j++) {
							fprintf(outfp, ", task_id%d", j+1);
						}
						fprintf(outfp, ", tasks_loop%d%s);\n", l, indices);
						fprintf(outfp, "\t\t\t\t\t\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
						fprintf(outfp, "\t\t\t\t\t\t} \n");
						fprintf(outfp, "##endif\n");

						fprintf(outfp, "\t\t\t\t\t\tif (tasks_loop%d%s <= %d) {\n", l, indices, num_write_data[l]);

						if (options->timereport) {
							fprintf(outfp, "\t\t\t\t\t\t\tIF_TIME(t_wait_unpack += rtclock() - t_wait_unpack_start);\n");
							fprintf(outfp, "\t\t\t\t\t\t\tIF_TIME(t_unpack_start = rtclock());\n");
						}
						fprintf(outfp, "\t\t\t\t\t\t\tif (tasks_loop%d%s == %d) {\n", l, indices, num_write_data[l]);
						for (ii=0; ii<num_write_data[l]; ii++) {
							char *acc_name = wacc_stmts[l][ii][0]->acc->name;
							fprintf(outfp, "\t\t\t\t\t\t\t\ttasks_loop%d%s -=\
                                    !is_receiver_%s_%d(%s,%s,nprocs);\n",
                                    l, indices, acc_name, l, args, params);
						}
						fprintf(outfp, "\t\t\t\t\t\t\t}\n");

						fprintf(outfp, "\t\t\t\t\t\t\tunpack_buf_%s_%d(", acc_names[i], l);
						fprintf(outfp, "my_rank,nprocs");
						if (options->variables_not_global && !options->data_dist) {
							fprintf(outfp, ",%s", acc_names[i]);
						}
						if(options->data_dist){
							fprint_data_dist_parm_call(outfp, acc_names[i], prog);
						}
						fprintf(outfp, ",recv_buf_%s[__p]);\n", acc_names[i]);

						fprintf(outfp, "\t\t\t\t\t\t\ttasks_loop%d%s--;\n", l, indices);

						fprintf(outfp, "\t\t\t\t\t\t\trecv_completed_%s[__p] = 0;\n", acc_names[i]);

						fprintf(outfp,
								"\t\t\t\t\t\t\tint __error = MPI_Irecv(recv_buf_%s[__p], %s, MPI_DOUBLE,\
                            MPI_ANY_SOURCE, %s, MPI_COMM_WORLD, &recv_reqs_%s[__p]);\n",
                            acc_names[i], max_elements_name[i], data_tag[i], acc_names[i]);
						fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
						fprintf(outfp, "\t\t\t\t\t\t\tif (__error) {\n");
						fprintf(outfp, "\t\t\t\t\t\t\t\tchar *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len;\n");
						fprintf(outfp, "\t\t\t\t\t\t\t\tMPI_Error_string(__error, msg, &len);\n");
						fprintf(outfp, "\t\t\t\t\t\t\t\tfprintf(__debug_print_fp, \"node %%d %s recv buf %%d error %%d %%s\\n\", my_rank, __p, __error, msg);\n", acc_names[i]);
						fprintf(outfp, "\t\t\t\t\t\t\t\tfree(msg);\n");
						fprintf(outfp, "\t\t\t\t\t\t\t}\n");
						fprintf(outfp, "##pragma omp atomic\n");
						fprintf(outfp, "\t\t\t\t\t\t\t__busy_wait_count *= 0;\n");
						fprintf(outfp, "\t\t\t\t\t\t\tfprintf(__debug_print_fp, \"unpacked %s node %%d buf %%d loop %%d task", acc_names[i]);
						for (j=0; j<max_copy_level; j++) {
							fprintf(outfp, " %%d");
						}
						fprintf(outfp, " count %%d\\n\", my_rank, __p, __loop");
						for (j=0; j<max_copy_level; j++) {
							fprintf(outfp, ", task_id%d", j+1);
						}
						fprintf(outfp, ", tasks_loop%d%s);\n", l, indices);
						fprintf(outfp, "\t\t\t\t\t\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
						fprintf(outfp, "##ifdef __DYNSCHEDULER_MORE_DEBUG_PRINT\n");
						fprintf(outfp, "\t\t\t\t\t\t\tprev_firing_count_%s_%d[__p] = -1;\n", acc_names[i], l);
						fprintf(outfp, "##endif\n");
						fprintf(outfp, "##endif\n");

						fprintf(outfp, "\t\t\t\t\t\t\tif (tasks_loop%d%s == 0) {\n", l, indices);
						fprintf(outfp, "\t\t\t\t\t\t\t\tremote_update_dep_tasks_%d(", l);
						for (k=0; k<copy_level[l]+prog->npar; k++)    {
							if (k!=0) fprintf(outfp, ",");
							if (k<=copy_level[l]-1) fprintf(outfp, "task_id%d", k+1);
							else fprintf(outfp, "%s", prog->params[k-copy_level[l]]);
						}
						// my_rank is already added as a parameter
						fprintf(outfp, ",nprocs,_num_threads%s);\n", tasks_loops_args);
						fprintf(outfp, "\t\t\t\t\t\t\t\t_num_tasks_to_unpack--;\n");
						fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
						fprintf(outfp, "\t\t\t\t\t\t\t\tfprintf(__debug_print_fp, \"remote updated %s node %%d loop %%d task", acc_names[i]);
						for (j=0; j<max_copy_level; j++) {
							fprintf(outfp, " %%d");
						}
						fprintf(outfp, " num tasks left to unpack %%ld\\n\", my_rank, __loop");
						for (j=0; j<max_copy_level; j++) {
							fprintf(outfp, ", task_id%d", j+1);
						}
						fprintf(outfp, ", _num_tasks_to_unpack);\n");
						fprintf(outfp, "\t\t\t\t\t\t\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
						fprintf(outfp, "##endif\n");
						fprintf(outfp, "\t\t\t\t\t\t\t}\n");
						if (options->timereport) {
							fprintf(outfp, "\t\t\t\t\t\t\tIF_TIME(t_unpack += rtclock() - t_unpack_start);\n");
							fprintf(outfp, "\t\t\t\t\t\t\tIF_TIME(t_wait_unpack_start = rtclock());\n");
						}

						fprintf(outfp, "\t\t\t\t\t\t}\n");

						fprintf(outfp, "\t\t\t\t\t\tbreak;\n");

						free(indices);
						free(args);
						break;
					}
				}
			}
			fprintf(outfp, "\t\t\t\t\tdefault:\n");
			fprintf(outfp, "\t\t\t\t\t\tfprintf(stderr, \"ERROR: %s data cannot be sent by a task of loop %%d\\n\", __loop);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\t\tassert(0);\n");
			fprintf(outfp, "\t\t\t\t}\n");
			fprintf(outfp, "\t\t\t}\n");
			fprintf(outfp, "\t\t}\n");
		}
		fprintf(outfp, "\t\t}\n");

		for (i=0; i<num_data; i++) {
			fprintf(outfp, "\t\tfor (__p=0; __p<recv_buf_%s.size(); __p++) {\n", acc_names[i]);
			fprintf(outfp, "\t\t\tif (recv_reqs_%s[__p] != MPI_REQUEST_NULL) break;\n", acc_names[i]);
			fprintf(outfp, "\t\t}\n");
			fprintf(outfp, "\t\tif (__p == recv_buf_%s.size()) {\n", acc_names[i]);
			fprintf(outfp, "\t\t\tfor (__p=0; __p<__MAX_NUM_RECVS; __p++) {\n");
			fprintf(outfp, "\t\t\t\tdouble *buf = NULL;\n");
			fprintf(outfp, "\t\t\t\tsize_t buf_size = 0;\n");
			fprintf(outfp, "\t\t\t\tbuf = (double *) polyrt_max_alloc(buf, sizeof(double)*%s, &buf_size);\n", max_elements_name[i]);
			fprintf(outfp, "\t\t\t\trecv_buf_%s.push_back(buf);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\trecv_reqs_%s.push_back(MPI_REQUEST_NULL);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\trecv_completed_%s.push_back(0);\n", acc_names[i]);
			fprintf(outfp, "##ifdef __DYNSCHEDULER_MORE_DEBUG_PRINT\n");
			for (l=0; l<nloops; l++) {
				fprintf(outfp, "\t\t\t\tprev_firing_count_%s_%d.push_back(-1);\n", acc_names[i], l);
			}
			fprintf(outfp, "##endif\n");
			fprintf(outfp,
					"\t\t\t\tMPI_Irecv(buf, %s, MPI_DOUBLE,\
                MPI_ANY_SOURCE, %s, MPI_COMM_WORLD, &recv_reqs_%s[recv_reqs_%s.size()-1]);\n",
                max_elements_name[i], data_tag[i], acc_names[i], acc_names[i]);
			fprintf(outfp, "\t\t\t}\n");
			fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
			fprintf(outfp, "\t\t\tfprintf(__debug_print_fp, \"node %%d recv thread %%d %s max recv buffers %%lu\\n\", my_rank, __tid, recv_buf_%s.size());\n", acc_names[i], acc_names[i]);
			fprintf(outfp, "\t\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
			fprintf(outfp, "##endif\n");
			fprintf(outfp, "\t\t}\n");
		}
		free(params);

		fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
		fprintf(outfp, "\t\tif (__busy_wait_count > 1000000) {\n");
		fprintf(outfp, "\t\t\tfprintf(stderr, \"\\nNO PROGRESS in Node %%d\\n\", my_rank);\n");
		fprintf(outfp, "\t\t\tfflush(stderr);\n");
		fprintf(outfp, "\t\t\texit(1);\n");
		fprintf(outfp, "\t\t}\n");
		fprintf(outfp, "##endif\n");
	}

	int section;
	// section == 0 -> receiver thread
	// section == 1 -> other threads
	for (section=0; section<2; section++) {
		if (section == 0) {
			if (!options->distmem) continue;
			fprintf(outfp, "\n##ifndef __DYNSCHEDULER_DEDICATED_RECEIVER\n");
		}
		fprintf(outfp, "\n{ /* compute-pack-send */\n");
		fprintf(outfp, "\twhile (_num_tasks_to_execute > 0) {\n");
		fprintf(outfp, "\t\t__Task task;\n");
		fprintf(outfp, "\t\t__flag = pqueue.try_pop(task);\n");
		fprintf(outfp, "\t\tif(__flag==0) { ");
		if (options->timereport)
			fprintf(outfp, "IF_TIME(waiting_pops[__tid]++); ");
		fprintf(outfp, "if (pqueue.size() == 0) {\n");
		fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
		fprintf(outfp, "\t\t\tint __i;\n");
		for (i=0; i<num_data; i++) {
			fprintf(outfp, "\t\t\tfor (__i=0;__i<send_reqs_%s[__tid].size();__i++) {\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\tif (send_reqs_%s[__tid][__i] != MPI_REQUEST_NULL) {\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\tint flag, __error = MPI_Test(&send_reqs_%s[__tid][__i], &flag, MPI_STATUS_IGNORE);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\tif (__error) {\n");
			fprintf(outfp, "\t\t\t\t\t\tchar *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len;\n");
			fprintf(outfp, "\t\t\t\t\t\tMPI_Error_string(__error, msg, &len);\n");
			fprintf(outfp, "\t\t\t\t\t\tfprintf(__debug_print_fp, \"node %%d %s send test buf %%d error %%d %%s\\n\", my_rank, __i, __error, msg);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\t\tfree(msg);\n");
			fprintf(outfp, "\t\t\t\t\t}\n");
			fprintf(outfp, "\t\t\t\t\tif (flag) {\n");
			fprintf(outfp, "\t\t\t\t\t\tfprintf(__debug_print_fp, \"node %%d %s send completed buf %%d\\n\", my_rank, __i);\n", acc_names[i]);
			fprintf(outfp, "\t\t\t\t\t}\n");
			fprintf(outfp, "\t\t\t\t}\n");
			fprintf(outfp, "\t\t\t}\n");
		}
		fprintf(outfp, "\t\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
		fprintf(outfp, "##endif\n");
		if (section == 0) {
			fprintf(outfp, "\t\t\tbreak;\n");
			fprintf(outfp, "\t\t} else {\n");
			fprintf(outfp, "\t\t\tcontinue;\n");
			fprintf(outfp, "\t\t} }\n");
		} else {
			fprintf(outfp, "\t\t}\n");
			fprintf(outfp, "\t\tcontinue;\n");
			fprintf(outfp, "\t\t}\n");
		}
		fprintf(outfp, "\t\t__loop = task.loop;\n");
		for (i=1; i<=max_copy_level; i++) {
			fprintf(outfp, "\t\ttask_id%d = task.id%d;\n", i, i);
		}
		fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
		fprintf(outfp, "\t\tlong int captured;\n");
		fprintf(outfp, "##pragma omp atomic capture\n");
		fprintf(outfp, "\t\tcaptured = --_num_tasks_to_execute;\n");
		fprintf(outfp, "\t\tfprintf(__debug_print_fp, \"computing node %%d thread %%d loop %%d task");
		for (j=0; j<max_copy_level; j++) {
			fprintf(outfp, " %%d");
		}
		fprintf(outfp, " ready tasks %%lu\\n\", my_rank, __tid, __loop");
		for (j=0; j<max_copy_level; j++) {
			fprintf(outfp, ", task_id%d", j+1);
		}
		fprintf(outfp, ", pqueue.size());\n");
		fprintf(outfp, "\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
		fprintf(outfp, "##else\n");
		fprintf(outfp, "##pragma omp atomic\n");
		fprintf(outfp, "\t\t_num_tasks_to_execute--;\n");
		fprintf(outfp, "##endif\n");
		fprintf(outfp, "\t\tswitch(__loop){\n");
		for (l=0; l<nloops; l++) {
			Stmt *anchor_stmt = loops[l]->stmts[0];
			int dims[copy_level[l]];
			int num_dims = 0;
			for (k=0; k<copy_level[l]; k++) {
				if (anchor_stmt->hyp_types[k] != H_SCALAR) {
					dims[num_dims++] = k;
				}
			}

			char *indices = malloc(1024);
			strcpy(indices,"[");
			for (ii=0; ii<num_dims; ii++) {
				sprintf(indices+strlen(indices), "(task_id%d-lb_tasks_loop%d_dim%d)",
						dims[ii]+1, l, ii);
				for (k=ii+1; k<num_dims; k++) {
					sprintf(indices+strlen(indices), "*max_num_tasks_loop%d_dim%d",
							l, k);
				}
				strcat(indices,"+");
			}
			strcat(indices,"0]");

			fprintf(outfp, "\t\t\tcase %d:\n", l);

			int num_rw_data;
			PlutoAccess **accs;
			accs = pluto_get_all_distinct_vars(loops[l]->stmts, loops[l]->nstmts, &num_rw_data);

			if (options->timereport) {
				fprintf(outfp, "\t\t\t\tIF_TIME(t_wait_comp[__tid] += rtclock() - t_wait_comp_start[__tid]);\n");

				if(options->data_dist){
					fprintf(outfp, "\t\t\t\tIF_TIME(t_data_mang_start[__tid] = rtclock());\n");
					fprintf(outfp, "\t\t\t\tdata_tile_alloc_%d(", l);
					fprintf(outfp, "task_id1");
					for (i=1; i<copy_level[l]; i++) {
						fprintf(outfp, ",task_id%d", i+1);
					}
					if (options->variables_not_global && !options->data_dist) {
						for (i=0; i<num_rw_data; i++) {
							char *acc_name = accs[i]->name;
							fprintf(outfp, ",%s", acc_name);
						}
					}
					for (i=0; i<num_rw_data; i++) {
						char *acc_name = accs[i]->name;
						fprint_data_dist_parm_call(outfp, acc_name, prog);
					}

					fprintf(outfp, ");\n");
					fprintf(outfp, "\t\t\t\tIF_TIME(t_data_mang[__tid] += rtclock() - t_data_mang_start[__tid]);\n");

				}



				fprintf(outfp, "\t\t\t\tIF_TIME(t_comp_start[__tid] = rtclock());\n");
			}
			fprintf(outfp, "\t\t\t\tcompute_task_%d(", l);
			fprintf(outfp, "task_id1");
			for (i=1; i<copy_level[l]; i++) {
				fprintf(outfp, ",task_id%d", i+1);
			}
			if (options->variables_not_global && !options->data_dist) {
				for (i=0; i<num_rw_data; i++) {
					char *acc_name = accs[i]->name;
					fprintf(outfp, ",%s", acc_name);
				}
			}
			if(options->data_dist){
				for (i=0; i<num_rw_data; i++) {
					char *acc_name = accs[i]->name;
					fprint_data_dist_parm_call(outfp, acc_name, prog);
				}
			}
			fprintf(outfp, ");\n");
			if (options->timereport)
				fprintf(outfp, "\t\t\t\tIF_TIME(t_comp[__tid] += rtclock() - t_comp_start[__tid]);\n");


			if(options->data_dist){
				fprintf(outfp, "\t\t\t\tIF_TIME(t_data_mang_start[__tid] = rtclock());\n");
				fprintf(outfp, "\t\t\t\tdata_tile_ref_count_update_%d(", l);
				fprintf(outfp, "task_id1");
				for (i=1; i<copy_level[l]; i++) {
					fprintf(outfp, ",task_id%d", i+1);
				}
				if (options->variables_not_global && !options->data_dist) {
					for (i=0; i<num_rw_data; i++) {
						char *acc_name = accs[i]->name;
						fprintf(outfp, ",%s", acc_name);
					}
				}
				for (i=0; i<num_rw_data; i++) {
					char *acc_name = accs[i]->name;
					fprint_data_dist_parm_call(outfp, acc_name, prog);
				}

				fprintf(outfp, ");\n");
				fprintf(outfp, "\t\t\t\tIF_TIME(t_data_mang[__tid] += rtclock() - t_data_mang_start[__tid]);\n");

			}

			if (options->distmem) {
				if (options->timereport)
					fprintf(outfp, "\t\t\t\tIF_TIME(t_pack_start[__tid] = rtclock());\n");
				fprintf(outfp, "\t\t\t\tpack_send_%d(", l);
				fprintf(outfp, "my_rank,nprocs");
				for (i=0; i<copy_level[l]; i++) {
					fprintf(outfp, ",task_id%d", i+1);
				}

				if(options->data_dist){
					for (i=0; i<num_rw_data; i++) {
						char *acc_name = accs[i]->name;
						if (options->variables_not_global && !options->data_dist) fprintf(outfp, ",%s", acc_name);
						fprint_data_dist_parm_call(outfp, acc_name, prog);
						fprintf(outfp, ",send_buf_%s[__tid]", acc_name);
						fprintf(outfp, ",send_reqs_%s[__tid]", acc_name);
//						fprintf(outfp, ",send_stats_%s[__tid]", acc_name);
					}
				}
				else{
					for (i=0; i<num_write_data[l]; i++) {
						char *acc_name = wacc_stmts[l][i][0]->acc->name;
						if (options->variables_not_global) fprintf(outfp, ",%s", acc_name);
						fprintf(outfp, ",send_buf_%s[__tid]", acc_name);
						fprintf(outfp, ",send_reqs_%s[__tid]", acc_name);
					}
				}
				if (options->timereport)
					fprintf(outfp, ",&__total_count[__tid]");
				fprintf(outfp, ");\n");
				if (options->timereport)
					fprintf(outfp, "\t\t\t\tIF_TIME(t_pack[__tid] += rtclock() - t_pack_start[__tid]);\n");
			}

			if (options->timereport)
				fprintf(outfp, "\t\t\t\tIF_TIME(t_tasks_manage_start[__tid] = rtclock());\n");
			fprintf(outfp, "\t\t\t\tlocal_update_dep_tasks_%d(", l);
			for (i=0; i<copy_level[l]+prog->npar; i++)    {
				if (i!=0) fprintf(outfp, ",");
				if (i<=copy_level[l]-1) fprintf(outfp, "task_id%d", i+1);
				else fprintf(outfp, "%s", prog->params[i-copy_level[l]]);
			}
			if (options->distmem) {
				// my_rank is already added as a parameter
				fprintf(outfp, ",nprocs");
			}
			fprintf(outfp, ",_num_threads");
			fprintf(outfp, "%s);\n", tasks_loops_args);
			if (options->timereport) {
				fprintf(outfp, "\t\t\t\tIF_TIME(t_tasks_manage[__tid] += rtclock() - t_tasks_manage_start[__tid]);\n");

				fprintf(outfp, "\t\t\t\tIF_TIME(t_wait_comp_start[__tid] = rtclock());\n");
			}
			fprintf(outfp, "\t\t\t\tbreak;\n");

			free(accs);
		}
		fprintf(outfp, "\t\t\tdefault:\n");
		fprintf(outfp, "\t\t\t\tassert(0);\n");
		fprintf(outfp, "\t\t}\n");
		fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
		fprintf(outfp, "\t\tfprintf(__debug_print_fp, \"computed node %%d thread %%d loop %%d task");
		for (j=0; j<max_copy_level; j++) {
			fprintf(outfp, " %%d");
		}
		fprintf(outfp, " tasks left %%ld\\n\", my_rank, __tid, __loop");
		for (j=0; j<max_copy_level; j++) {
			fprintf(outfp, ", task_id%d", j+1);
		}
		fprintf(outfp, ", captured);\n");
		fprintf(outfp, "\t\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
		fprintf(outfp, "##pragma omp atomic\n");
		fprintf(outfp, "\t\t__busy_wait_count *= 0;\n");
		fprintf(outfp, "##endif\n");
		fprintf(outfp, "\t}\n");
		fprintf(outfp, "}\n");
		if (section == 0) {
			fprintf(outfp, "\n##endif\n\n");
			fprintf(outfp, "\t}\n");
			fprintf(outfp, "}\n");
		}
	}

	fprintf(outfp, "\n##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
	fprintf(outfp, "\tfprintf(__debug_print_fp, \"node %%d thread %%d done\\n\", my_rank, __tid);\n");
	fprintf(outfp, "\tIF_DEBUG_FLUSH(fflush(__debug_print_fp));\n");
	fprintf(outfp, "##endif\n");

	fprintf(outfp, "\n}\n\n");

	if (options->distmem) {
		if (options->timereport)
			fprintf(outfp, "\tIF_TIME(t_unpack_start = rtclock());\n");
		for (i=0; i<num_data; i++) {
			fprintf(outfp, 	"\tfor (__p=0; __p<recv_buf_%s.size(); __p++) {\n", acc_names[i]);
			fprintf(outfp,
					"\t\tif (recv_reqs_%s[__p] != MPI_REQUEST_NULL) \
                MPI_Request_free(&recv_reqs_%s[__p]);\n",
                acc_names[i], acc_names[i]);
			fprintf(outfp, "\t}\n");
		}
		if (options->timereport)
			fprintf(outfp, "\tIF_TIME(t_unpack += rtclock() - t_unpack_start);\n");
	}

	fprintf(outfp, "##ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
	fprintf(outfp, "\tfclose(__debug_print_fp);\n");
	fprintf(outfp, "##endif\n");

	fprintf(outfp, "\n\n\n/*End of dynamic scheduling code*/\n\n\n");

	free(tasks_loops_args);

	for (l=1; l<nloops; l++) {
		for (i=0; i<num_write_data[l]; i++) {
			free(wacc_stmts[l][i]);
		}
		free(wacc_stmts[l]);
		free(num_stmts_per_wacc[l]);
	}
	for (i=0; i<num_data; i++) {
		if (data_tag[i] != NULL) free(data_tag[i]);
	}
	free(waccs);
}

void gen_dynschedule_graph_header_text_code(PlutoProg *prog, int *copy_level, int nloops, FILE *headerfp)
{
	int i, l;

	int max_copy_level = copy_level[0];
	for (l=1; l<nloops; l++) {
		if (copy_level[l] > max_copy_level) max_copy_level = copy_level[l];
	}

	int distinct_copy_level[nloops];
	int num_distinct_copy_level = 0;
	for (l=0; l<nloops; l++) {
		for (i=0; i<num_distinct_copy_level; i++) {
			if (copy_level[l] == distinct_copy_level[i]) break;
		}
		if (i == num_distinct_copy_level) {
			distinct_copy_level[num_distinct_copy_level++] = copy_level[l];
		}
	}

	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_all_distinct_vars(prog->stmts, prog->nstmts, &num_data);

	fprintf(headerfp, "#include \"tbb/flow_graph.h\"\n\n");

	fprintf(headerfp, "class __Task {");
	fprintf(headerfp, "public:");
	fprintf(headerfp, "int loop;");
	for (i=1; i<=max_copy_level; i++) {
		fprintf(headerfp, "int id%d;", i);
	}
	if (options->variables_not_global) {
		for (i=0; i<num_data; i++) {
			fprintf(headerfp, "__DECLARATION_OF_POINTER_TO_%s(%s);", accs[i]->name, accs[i]->name);
		}
	}
	fprintf(headerfp, " __Task() {};");
	for (l=0; l<num_distinct_copy_level; l++) {
		fprintf(headerfp, " __Task(const int l");
		for (i=1; i<=distinct_copy_level[l]; i++) {
			fprintf(headerfp, ", const int _id%d", i);
		}
		if (options->variables_not_global) {
			for (i=0; i<num_data; i++) {
				fprintf(headerfp, ", __DECLARATION_OF_POINTER_TO_%s(_%s)", accs[i]->name, accs[i]->name);
			}
		}
		fprintf(headerfp, ") :");
		fprintf(headerfp, "loop(l)");
		for (i=1; i<=distinct_copy_level[l]; i++) {
			fprintf(headerfp, ",id%d(_id%d)", i, i);
		}
		for (; i<=max_copy_level; i++) {
			fprintf(headerfp, ",id%d(0)", i);
		}
		if (options->variables_not_global) {
			for (i=0; i<num_data; i++) {
				fprintf(headerfp, ",%s(_%s)", accs[i]->name, accs[i]->name);
			}
		}
		fprintf(headerfp, " {};");
	}
	fprintf(headerfp, " __Task(const __Task &task) :");
	fprintf(headerfp, "loop(task.loop)");
	for (i=1; i<=max_copy_level; i++) {
		fprintf(headerfp, ",id%d(task.id%d)", i, i);
	}
	if (options->variables_not_global) {
		for (i=0; i<num_data; i++) {
			fprintf(headerfp, ",%s(task.%s)", accs[i]->name, accs[i]->name);
		}
	}
	fprintf(headerfp, " {};");
	fprintf(headerfp, " ~__Task() {};");
	fprintf(headerfp, " void operator()(tbb::flow::continue_msg) const;");
	fprintf(headerfp, " };\n\n");

	free(accs);
}

void gen_dynschedule_header_text_code(int *copy_level, int nloops, FILE *headerfp)
{
	int i, l;

	int max_copy_level = copy_level[0];
	for (l=1; l<nloops; l++) {
		if (copy_level[l] > max_copy_level) max_copy_level = copy_level[l];
	}

	int distinct_copy_level[nloops];
	int num_distinct_copy_level = 0;
	for (l=0; l<nloops; l++) {
		for (i=0; i<num_distinct_copy_level; i++) {
			if (copy_level[l] == distinct_copy_level[i]) break;
		}
		if (i == num_distinct_copy_level) {
			distinct_copy_level[num_distinct_copy_level++] = copy_level[l];
		}
	}

	if (options->distmem) {
		fprintf(headerfp, "#include <mpi.h>\n");
	}
	fprintf(headerfp, "#include \"tbb/concurrent_priority_queue.h\"\n\n");
	fprintf(headerfp, "#include <omp.h>\n\n");

	fprintf(headerfp, "#ifdef __DEBUG_FLUSH\n");
	fprintf(headerfp, "#define IF_DEBUG_FLUSH(foo) foo\n");
	fprintf(headerfp, "#else\n");
	fprintf(headerfp, "#define IF_DEBUG_FLUSH(foo) \n");
	fprintf(headerfp, "#endif\n\n");

	fprintf(headerfp, "#ifdef __DYNSCHEDULER_MORE_DEBUG_PRINT\n");
	fprintf(headerfp, "#define IF_DYNSCHEDULER_MORE_DEBUG_PRINT(foo) foo\n");
	fprintf(headerfp, "#else\n");
	fprintf(headerfp, "#define IF_DYNSCHEDULER_MORE_DEBUG_PRINT(foo) \n");
	fprintf(headerfp, "#endif\n\n");

	fprintf(headerfp, "#ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
	fprintf(headerfp, "#define IF_DYNSCHEDULER_DEBUG_PRINT(foo) foo\n");
	fprintf(headerfp, "#else\n");
	fprintf(headerfp, "#define IF_DYNSCHEDULER_DEBUG_PRINT(foo) \n");
	fprintf(headerfp, "#endif\n\n");

	fprintf(headerfp, "extern FILE *__debug_print_fp;\n\n");
	fprintf(headerfp, "FILE *__debug_print_fp;\n");

	fprintf(headerfp, "class __Task {\n");
	fprintf(headerfp, "\tpublic:\n");
	fprintf(headerfp, "\t\tint loop;\n");
	for (i=1; i<=max_copy_level; i++) {
		fprintf(headerfp, "\t\tint id%d;\n", i);
	}
	if (options->distmem) {
		fprintf(headerfp, "\t\tint remote_dep_tasks;\n");
	}
	fprintf(headerfp, "\t\tint local_dep_tasks;\n");
	fprintf(headerfp, "\t\tint thread_affinity;\n");
	fprintf(headerfp, "\t\t__Task() {}\n");
	for (l=0; l<num_distinct_copy_level; l++) {
		fprintf(headerfp, "\t\t__Task(const int l");
		for (i=1; i<=distinct_copy_level[l]; i++) {
			fprintf(headerfp, ", const int _id%d", i);
		}
		if (options->distmem) {
			fprintf(headerfp, ", const int remote");
		}
		fprintf(headerfp, ", const int local, const int affinity) :\n");
		fprintf(headerfp, "\t\t\tloop(l),\n");
		for (i=1; i<=distinct_copy_level[l]; i++) {
			fprintf(headerfp, "\t\t\tid%d(_id%d),\n", i, i);
		}
		for (; i<=max_copy_level; i++) {
			fprintf(headerfp, "\t\t\tid%d(0),\n", i);
		}
		if (options->distmem) {
			fprintf(headerfp, "\t\t\tremote_dep_tasks(remote),\n");
		}
		fprintf(headerfp, "\t\t\tlocal_dep_tasks(local),\n");
		fprintf(headerfp, "\t\t\tthread_affinity(affinity) {}\n");
	}
	fprintf(headerfp, "\t\t__Task(const __Task &task) :\n");
	fprintf(headerfp, "\t\t\tloop(task.loop),\n");
	for (i=1; i<=max_copy_level; i++) {
		fprintf(headerfp, "\t\t\tid%d(task.id%d),\n", i, i);
	}
	if (options->distmem) {
		fprintf(headerfp, "\t\t\tremote_dep_tasks(task.remote_dep_tasks),\n");
	}
	fprintf(headerfp, "\t\t\tlocal_dep_tasks(task.local_dep_tasks),\n");
	fprintf(headerfp, "\t\t\tthread_affinity(task.thread_affinity) {}\n");
	fprintf(headerfp, "\t\t~__Task() {}\n");
	fprintf(headerfp, "};\n\n");

	fprintf(headerfp, "class __Compare {\n");
	fprintf(headerfp, "\tpublic:\n");
	fprintf(headerfp, "\t\t// returns true only if task2 has higher priority than task1\n");
	fprintf(headerfp, "\t\tbool operator()(const __Task &task1, const __Task &task2) const {\n");
	fprintf(headerfp, "#ifndef __DYNSCHEDULER_NO_PRIORITY\n");
	if (options->distmem) {
		fprintf(headerfp, "\t\t\tif (task2.remote_dep_tasks>task1.remote_dep_tasks) return true;\n");
		fprintf(headerfp, "\t\t\tif (task1.remote_dep_tasks>task2.remote_dep_tasks) return false;\n");
	}
	fprintf(headerfp, "\t\t\tif (task2.local_dep_tasks>task1.local_dep_tasks) return true;\n");
	fprintf(headerfp, "\t\t\tif (task1.local_dep_tasks>task2.local_dep_tasks) return false;\n");
	fprintf(headerfp, "\t\t\tint thread_id = omp_get_thread_num();\n");
	fprintf(headerfp, "\t\t\tif ((thread_id == task2.thread_affinity) && (thread_id != task1.thread_affinity)) return true;\n");
	fprintf(headerfp, "\t\t\tif ((thread_id == task1.thread_affinity) && (thread_id != task2.thread_affinity)) return false;\n");
	fprintf(headerfp, "\t\t\t// lexicographic compare \n");
	fprintf(headerfp, "\t\t\tif ( (task2.id1 > task1.id1) ");
	for (i=2; i<=max_copy_level; i++) {
		fprintf(headerfp, "\n\t\t\t\t|| ( (task2.id1 == task1.id1) ");
		for (l=2; l<i; l++) {
			fprintf(headerfp, "&& (task2.id%d == task1.id%d) ", l, l);
		}
		fprintf(headerfp, "&& (task2.id%d > task1.id%d) ", i, i);
		fprintf(headerfp, ") ");
	}
	fprintf(headerfp, ")\n\t\t\t\t return true;\n");
	fprintf(headerfp, "#endif\n");
	fprintf(headerfp, "\t\t\treturn false;\n");
	fprintf(headerfp, "\t\t}\n");
	fprintf(headerfp, "};\n\n");

	fprintf(headerfp, "extern tbb::concurrent_priority_queue<__Task, __Compare> pqueue;\n\n");
	fprintf(headerfp, "tbb::concurrent_priority_queue<__Task, __Compare> pqueue;\n");
}

/*
 * Optimized communication code generation wih dependence spliting
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * wacc_stmts:  all <statements, wacc> writing to this data variable
 * This function is called per data variable
 */
Stmt **gen_comm_code_opt_fop(int data_id, struct stmt_access_pair **wacc_stmts, int num_accs,
		int nloops, int num_data, PlutoProg *prog, Stmt *anchor_stmt, int *copy_level, int outer_dist_loop_level,
		int loop_num, int *pi_mappings, int* num_comm_stmts, FILE *sigmafp, FILE *headerfp)
{
	int i, k,s, src_copy_level, max_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];
	assert(src_copy_level>=1);
	max_copy_level = copy_level[0];
	for (i=1; i<nloops; i++) {
		if (max_copy_level < copy_level[i]) {
			max_copy_level = copy_level[i];
		}
	}

	assert(num_accs >= 1);
	char *access = reconstruct_access(wacc_stmts[0]->acc);
	acc_nrows = wacc_stmts[0]->acc->mat->nrows;
	char *acc_name = wacc_stmts[0]->acc->name;
	char *data_tag;
	data_tag = malloc(10);
	get_data_tag(prog, wacc_stmts[0]->acc, &data_tag);

	Array *arr = pluto_get_corrs_array(acc_name, prog);
	IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

	/* To be inside a loop: can't foresee other use */
	assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));

	/* Sender-side copying */
	PlutoConstraintsList *atomic_flowouts = pluto_constraints_list_alloc((PlutoConstraints *) NULL);
	PlutoConstraints *flow_out = NULL; // used only for get parametric extent/bounding box
	for (k=0; k<num_accs; k++)  {
		PlutoConstraints *flow_out_one;
		flow_out_one = compute_flow_out(wacc_stmts[k], src_copy_level, copy_level, prog, pi_mappings);
		if (flow_out == NULL) flow_out = pluto_constraints_dup(flow_out_one);
		else{
			pluto_constraints_unionize(flow_out, flow_out_one);
		}
		pluto_constraints_free(flow_out_one);

		compute_flow_out_partitions(wacc_stmts[k], src_copy_level, copy_level, prog, atomic_flowouts, pi_mappings);
	}

	PlutoConstraintsList *curr;
	if (outer_dist_loop_level == (src_copy_level-1)) { // implies !multi_level_distribution
		PlutoConstraintsList *prev = NULL;
		curr = atomic_flowouts;
		PlutoConstraintsList *const_cst = NULL;
		while(curr != NULL) {
			PlutoConstraints *cst = curr->constraints;
			int const_bounds = 1;
			while (cst != NULL) {
				for (i=src_copy_level; i<src_copy_level+acc_nrows; i++) {
					if (get_const_bound_difference(cst, i) == -1) {
						const_bounds = 0;
						break;
					}
				}
				if (const_bounds == 0) break;
				cst = cst->next;
			}
			if (const_bounds == 1) {
				if (const_cst == NULL) {
					const_cst = curr;
				}
				else {
					pluto_constraints_unionize(const_cst->constraints, curr->constraints);
					PlutoDepList *dep_list = curr->deps;
					while (dep_list != NULL) {
						pluto_deps_list_append(const_cst->deps, dep_list->dep);
						dep_list = dep_list->next;
					}

					prev->next = curr->next;
					curr->next = NULL;
					pluto_constraints_list_free(curr);
					curr = prev->next;
					continue;
				}
			}
			prev = curr;
			curr = curr->next;
		}
	}

	split_deps_acc_flowout(atomic_flowouts, src_copy_level , acc_nrows, prog);

	generate_sigma_dep_split(wacc_stmts, num_accs, copy_level, prog, atomic_flowouts, loop_num,
			pi_mappings, sigmafp, headerfp);

	// if opencl is not enabled, generate the MPI code for distributed cluster
	//
		// ===== from here on starts the code generation =====================================================

		/***************************************************************************************************/

		char *send_buf_size =
				get_parametric_bounding_box(flow_out, src_copy_level, acc_nrows,
						prog->npar, (const char **)prog->params);

		if (!(options->dynschedule)) {
			PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
			char *total_extent = malloc(1024);
			strcpy(total_extent, "1");
			/* Assumes load-balanced distribution in each dimension */
			for (i=outer_dist_loop_level; i<src_copy_level; i++) {
				char *extent;
				/* Just the first one in flow_out is enough (rest all should give the
				 * same since they are all under the same parallel loop and each
				 * iteration of the parallel loop writes to distinct data) */
				get_parametric_extent_const(anchor_stmt_new_dom, i, prog->npar,
						(const char **)prog->params, &extent, NULL);
				/* The + nprocs is needed since the send buffer size should be larger
				 * when some processors have more iterations than others - some processors
				 * have an extra iteration in each dimension in the worst case */
				sprintf(total_extent+strlen(total_extent), "*(%s+nprocs)", extent);
				free(extent);
			}
			sprintf(send_buf_size+strlen(send_buf_size),
					"*floorf((%s)/(float)nprocs)", total_extent);
			free(total_extent);
			pluto_constraints_free(anchor_stmt_new_dom);
		}
		IF_DEBUG(printf("Send buffer size for %s: %s\n", acc_name, send_buf_size););

		char *sendbufname = concat("send_buf_", acc_name);
		char *send_counts_name = concat("send_counts_", acc_name);
		char *max_elements_name = concat("max_num_elements_", acc_name);

		if (!options->dynschedule) {
            sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
                %s[__p] = 0;\n}\n", send_counts_name);
		}

		if (options->dynschedule) { // one more for the loop number
			sprintf(prog->decls+strlen(prog->decls), "\
            %s = max(%s, %s + %d);\n", max_elements_name, max_elements_name, send_buf_size, max_copy_level + 1);
			sprintf(prog->decls+strlen(prog->decls), "\
                for (__tid=0; __tid<_num_threads; __tid++) { \
                for (__p=0; __p<__MAX_NUM_SENDS; __p++) {\n\
                %s[__tid][__p] = (double *) polyrt_max_alloc(%s[__tid][__p], \
                sizeof(double)*%s, &send_buf_size_%s[__tid][__p]);\n}\n}\n",
                sendbufname, sendbufname, max_elements_name, acc_name);
		}
		else {
			sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
				%s[__p]\
				= (double *) polyrt_max_alloc(%s[__p], sizeof(double)*(%s), &send_buf_size_%s[__p]);\n}",
				sendbufname, sendbufname, send_buf_size, acc_name);
		}

		char *flow_copy_text = malloc(strlen(sendbufname) + strlen("[")+
				strlen(send_counts_name) + strlen("++] = ") + strlen(access) + 1);
		sprintf(flow_copy_text, "%s[%s++] = %s", sendbufname, send_counts_name, access);

		char *displsname = concat("displs_", acc_name);
		char *currdisplsname = concat("curr_displs_", acc_name);
		char *recvbufname = concat("recv_buf_", acc_name);
		char *recv_counts_name = concat("recv_counts_", acc_name);

		if (options->dynschedule) {
			sprintf(prog->decls+strlen(prog->decls),
					"for (__p=0; __p<__MAX_NUM_RECVS; __p++) { \
                %s[__p] = (double *) polyrt_max_alloc(%s[__p], sizeof(double)*%s, &recv_buf_size_%s[__p]); }\n",
                recvbufname, recvbufname, max_elements_name, acc_name);
		}
		else {
			sprintf(prog->decls+strlen(prog->decls),
					"%s = (double *) polyrt_max_alloc(%s, nprocs*send_buf_size_%s[0], &recv_buf_size_%s);\n",
					recvbufname, recvbufname, acc_name, acc_name);
		}
		if (!options->dynschedule) {
			sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
                %s[__p] = __p*send_buf_size_%s[__p]/sizeof(double);}\n\n", displsname, acc_name);
		}

		char args[1024];
		strcpy(args,"");
		sprintf(args+strlen(args), "t%d", 1);
		/* make it src_copy_level+1 since we use an extra dimension to separate
		 * statements */
		for (i=1; i<src_copy_level; i++) {
			sprintf(args+strlen(args), ",t%d", i+1);
		}

		char passed_args[1024];
		strcpy(passed_args,"");
		sprintf(passed_args+strlen(passed_args), "ts%d", 1);
		for (i=1; i<src_copy_level; i++) {
			sprintf(passed_args+strlen(passed_args), ",ts%d", i+1);
		}

		char decl_args[1024];
		strcpy(decl_args,"");
		sprintf(decl_args+strlen(decl_args), "int ts%d", 1);
		for (i=1; i<src_copy_level; i++) {
			sprintf(decl_args+strlen(decl_args), ",int ts%d", i+1);
		}

		char params[1024];
		strcpy(params, "");
		if (prog->npar>=1) {
			sprintf(params+strlen(params), "%s", prog->params[0]);
			for (i=1; i<prog->npar; i++) {
				sprintf(params+strlen(params), ",%s", prog->params[i]);
			}
		}

		char **iters;
		iters = malloc(acc_nrows * sizeof(char *));
		for (i=0; i < acc_nrows; i++) {
			iters[i] = malloc(5);
			sprintf(iters[i], "d%d", i+1);
		}

		char *comm_text = malloc(2048);

		if (options->dynschedule) {
			strcpy(comm_text, "");
			sprintf(comm_text+strlen(comm_text),  "\
				for (__p=0; __p<nprocs; __p++) {\
				if(%s[__p] > 0) {",
				send_counts_name);
			sprintf(comm_text+strlen(comm_text), "%s[index_%s[__p]][0] = %d;", sendbufname, acc_name, loop_num);
			for (i=0; i<src_copy_level; i++) {
				sprintf(comm_text+strlen(comm_text), "%s[index_%s[__p]][%d] = t%d;", sendbufname, acc_name, i+1, i+1);
			}
			sprintf(comm_text+strlen(comm_text),  "\
		        IF_DYNSCHEDULER_DEBUG_PRINT(\
		        fprintf(__debug_print_fp, \"sending node %%d buf %%d loop %d task ",
		        loop_num);
			for (i=0; i<src_copy_level; i++) {
				sprintf(comm_text+strlen(comm_text),  "%%d ");
			}
			sprintf(comm_text+strlen(comm_text),  "\
		        %s size %%d to node %%d\\n\", \
		        my_rank, index_%s[__p], %s, %s[__p], __p)); \
                IF_DEBUG_FLUSH(fflush(__debug_print_fp));",
                acc_name, acc_name, args, send_counts_name);
			if (options->timereport)
				sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s[__p]);", send_counts_name);
			sprintf(comm_text+strlen(comm_text),  "\
				int __error = MPI_Isend(%s[index_%s[__p]], %s[__p], MPI_DOUBLE,\
					__p, %s, MPI_COMM_WORLD, &send_reqs_%s[index_%s[__p]]);\
				%s[__p] = 0;\
                IF_DYNSCHEDULER_DEBUG_PRINT(\
                if (__error) {\
                char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                MPI_Error_string(__error, msg, &len); \
                fprintf(__debug_print_fp, \"node %%d %s send buf %%d error %%d %%s\\n\", my_rank, index_%s[__p], __error, msg); } )\
                }}",
                sendbufname, acc_name, send_counts_name, data_tag, acc_name, acc_name, send_counts_name, acc_name, acc_name);
		}
		else {
			if (options->timereport)
				sprintf(comm_text,  "IF_TIME(t_comm_start = rtclock());");
			else
				strcpy(comm_text, "");
			sprintf(comm_text+strlen(comm_text), "\
				MPI_Alltoall(%s, 1, MPI_INT,\
					%s, 1, MPI_INT, MPI_COMM_WORLD);",
					send_counts_name, recv_counts_name);
			sprintf(comm_text+strlen(comm_text),  "\
				req_count=0;\
				for (__p=0; __p<nprocs; __p++) {\
				if (%s[__p] >= 1) {",
				send_counts_name);
			if (options->timereport)
				sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s[__p]);", send_counts_name);
			sprintf(comm_text+strlen(comm_text),  "\
				MPI_Isend(%s[__p], %s[__p], MPI_DOUBLE,\
					__p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}",
					sendbufname, send_counts_name);
			sprintf(comm_text+strlen(comm_text),  "for (__p=0; __p<nprocs; __p++) {\
				if(%s[__p] >= 1) {\
				MPI_Irecv(%s+%s[__p], %s[__p], MPI_DOUBLE,\
					__p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}\
				MPI_Waitall(req_count, reqs, stats);\
				for (__p=0; __p<nprocs; __p++) {\
				%s[__p] = 0;}",
				recv_counts_name, recvbufname, displsname, recv_counts_name,
				send_counts_name);
			sprintf(comm_text+strlen(comm_text), "for (__p=0; __p<nprocs; __p++) {\
				%s[__p] = %s[__p]; }", currdisplsname, displsname);
			if (options->timereport)
				sprintf(comm_text+strlen(comm_text),  "IF_TIME(t_comm += rtclock() - t_comm_start);");
		}

		/* Receiver-side copy */
		/* With send/recv -based more exact communication scheme, recv side copy
		 * does not include all the sender-side copies; processes that do not send
		 * data are to be skipped */
		assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
				prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);

		char *flow_copyback_text = malloc(strlen(access) + strlen(" = ") + strlen(recvbufname)
				+ strlen("[") + strlen(currdisplsname) + strlen("++]") + 1);
		sprintf(flow_copyback_text, "%s = %s[%s++]",
				access, recvbufname, currdisplsname);

		/* Create statement stub to get hold of the outer loops for the copy
		 * stmt; rest of the loops to be added are the actual copy loops */
		Stmt *comm_stmt = create_helper_stmt(anchor_stmt, outer_dist_loop_level, comm_text, COMM_CALL, loop_num);

		FILE *packfp = fopen(__APPENDFILENAME, "a");
		assert(packfp != NULL);

		char *pack_stmt_text = malloc(40960);
		assert(pack_stmt_text != NULL);
		strcpy(pack_stmt_text, "");
		char *unpack_stmt_text = malloc(40960);
		assert(unpack_stmt_text != NULL);
		if (options->dynschedule) {
			strcpy(unpack_stmt_text, "");
			// one more for the loop number
			sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
					"%s = %d; ", currdisplsname, max_copy_level + 1);
			if (options->fop_unicast_runtime) {
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "int __p = pi_%d(", loop_num);
				for (i=0; i<src_copy_level; i++) {
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
							"(int)%s[%d],",
							recvbufname, i+1);
				}
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "%s, nprocs); ", params);
			}
		}
		else {
			sprintf(unpack_stmt_text, "proc = pi_%d(%s,%s, nprocs); ", loop_num, args, params);
			sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "if ((my_rank != proc) && (%s[proc] > 0)) { ", recv_counts_name);
		}
		k = 0;
		curr = atomic_flowouts;
		while(curr != NULL) {
			if ((curr->constraints == NULL) || pluto_constraints_is_empty(curr->constraints)) {
				curr = curr->next;
				continue;
			}

			int j, l;
			int broadcast = 0;

			char acc_name_k[512];
			sprintf(acc_name_k, "%s_%d", acc_name, k+1);

			if (options->fop_unicast_runtime) {
				PlutoConstraints *foifi_sets[nloops];
				PlutoConstraints *receivers[nloops];
				PlutoDepList *dep_list = curr->deps;
				for (l=0; l<nloops; l++) {
					foifi_sets[l] = NULL;
					receivers[l] = NULL;
				}
				while (dep_list != NULL) {
					Dep *dep = dep_list->dep;
					Stmt *dest = prog->stmts[dep->dest];
					int dependent_loop = pi_mappings[dest->id];
					if (dependent_loop == -1) {
						// destination statement will be executed by all processors
						broadcast = 1;
						break;
					}

					int dest_copy_level = copy_level[dependent_loop];

					PlutoConstraints *fi = compute_flow_in_of_dep(dep, dest_copy_level, prog, 1);
					PlutoConstraints *rt = get_receiver_tiles_of_dep(dep, src_copy_level, dest_copy_level, prog, 1);

					if (foifi_sets[dependent_loop] == NULL)
						foifi_sets[dependent_loop] = pluto_constraints_dup(fi);
					else
						pluto_constraints_unionize(foifi_sets[dependent_loop], fi);

					if (receivers[dependent_loop] == NULL)
						receivers[dependent_loop] = pluto_constraints_dup(rt);
					else
						pluto_constraints_unionize(receivers[dependent_loop], rt);

					pluto_constraints_free(fi);
					pluto_constraints_free(rt);
					dep_list = dep_list->next;
				}
				if (!broadcast) {
					char sigma_check_text[1024];
					sprintf(sigma_check_text, "distinct_recv = sigma_check_%s_%d(%s,%s", acc_name_k, loop_num, args, params);

					sprintf(pack_stmt_text+strlen(pack_stmt_text), "%s, my_rank, nprocs); if (distinct_recv == 1) {", sigma_check_text);
					if (options->dynschedule) {
						sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "distinct_recv = sigma_check_%s_%d(", acc_name_k, loop_num);
						for (i=0; i<src_copy_level; i++) {
							sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
									"(int)%s[%d],",
									recvbufname, i+1);
						}
						sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "%s,", params);
					}
					else {
						sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "%s,", sigma_check_text);
					}
					if (options->dynschedule) {
						sprintf(unpack_stmt_text+strlen(unpack_stmt_text), " __p, nprocs); if (distinct_recv == 1) {");
					}
					else {
						sprintf(unpack_stmt_text+strlen(unpack_stmt_text), " proc, nprocs); if (distinct_recv == 1) {");
					}

					for (l=0; l<nloops; l++) {
						PlutoConstraints *foifi = foifi_sets[l];
						if (foifi == NULL) continue;

						int dest_copy_level = copy_level[l];
						assert(dest_copy_level>=1);
						int total_copy_level = src_copy_level + dest_copy_level;
						assert(foifi->ncols == dest_copy_level + acc_nrows + prog->npar + 1);

						PlutoConstraints *fo = pluto_constraints_dup(curr->constraints);
						// add target tile iterators for flow-out set
						for (j=0; j<dest_copy_level; j++) {
							pluto_constraints_add_dim(fo,src_copy_level, NULL);
						}

						// add source tile iterators for flow-in set
						for (j=0; j<src_copy_level; j++) {
							pluto_constraints_add_dim(foifi,0, NULL);
						}

						pluto_constraints_intersect(foifi, fo);

						if (pluto_constraints_is_empty(foifi)) {
							pluto_constraints_free(fo);
							pluto_constraints_free(foifi);
							pluto_constraints_free(receivers[l]);
							continue;
						}

						PlutoConstraints *receiver_tiles = receivers[l];
						assert(receiver_tiles->ncols == total_copy_level + prog->npar + 1);
						assert(!pluto_constraints_is_empty(receiver_tiles));

						char dest_args[1024];
						strcpy(dest_args,"");
						sprintf(dest_args+strlen(dest_args), "t%d", src_copy_level+1);
						for (i=src_copy_level+1; i<total_copy_level; i++) {
							sprintf(dest_args+strlen(dest_args), ",t%d", i+1);
						}

						char decl_dest_args[1024];
						strcpy(decl_dest_args,"");
						sprintf(decl_dest_args+strlen(decl_dest_args), "int ts%d", src_copy_level+1);
						for (i=src_copy_level+1; i<total_copy_level; i++) {
							sprintf(decl_dest_args+strlen(decl_dest_args), ",int ts%d", i+1);
						}

						char **dest_iters;
						dest_iters = malloc(dest_copy_level * sizeof(char *));
						for (i=0; i<dest_copy_level; i++) {
							dest_iters[i] = malloc(5);
							sprintf(dest_iters[i], "t%d", i+src_copy_level+1);
						}

						char pi_text[1024];
						sprintf(pi_text, "recv_proc = pi_%d(%s,%s, nprocs)", l, dest_args, params);

						if (options->dynschedule) {
							sprintf(pack_stmt_text+strlen(pack_stmt_text),
									"pack_foifi_%s_%d_%d(%s,%s,%s, my_rank, nprocs, send_reqs_%s, index_%s, &current_index_%s",
									acc_name_k, loop_num, l, args, sendbufname, send_counts_name, acc_name, acc_name, acc_name);
						}
						else {
							sprintf(pack_stmt_text+strlen(pack_stmt_text), "pack_foifi_%s_%d_%d(%s,%s,%s, my_rank, nprocs",
									acc_name_k, loop_num, l, args, sendbufname, send_counts_name);
						}
						if (options->variables_not_global && !options->data_dist) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);
						print_data_dist_parm_call(pack_stmt_text, arr);

						sprintf(pack_stmt_text+strlen(pack_stmt_text), ");");

						char flow_cg_text[2048];
						sprintf(flow_cg_text, "%s; if (recv_proc != my_rank) { ",
								pi_text);
						if (options->dynschedule) {
							sprintf(flow_cg_text+strlen(flow_cg_text), "\
                                if (%s[recv_proc] == 0) { \
                                for (;current_index<send_reqs.size();current_index++) { \
                                if (send_reqs[current_index] == MPI_REQUEST_NULL) break; \
                                int __error = MPI_Test(&send_reqs[current_index], &flag, MPI_STATUS_IGNORE); \
                                IF_DYNSCHEDULER_DEBUG_PRINT(\
                                if (flag) {\
                                fprintf(__debug_print_fp, \"node %%d %s send completed buf %%d\\n\", my_rank, current_index); } \
                                if (__error) {\
                                char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                                MPI_Error_string(__error, msg, &len); \
                                fprintf(__debug_print_fp, \"node %%d %s send test buf %%d error %%d %%s\\n\", my_rank, current_index, __error, msg); } )\
                                if (flag) break; } \
                                if (current_index==send_reqs.size()) { \
                                double *buf = NULL; \
                                size_t buf_size = 0; \
                                buf = (double *) polyrt_max_alloc(buf, sizeof(double)*max_num_elements_%s, &buf_size); \
                                %s.push_back(buf); \
                                IF_DYNSCHEDULER_DEBUG_PRINT(\
                                fprintf(__debug_print_fp, \"node %%d %s max send buffers %%d\\n\", my_rank, current_index+1)); \
                                send_reqs.push_back(MPI_REQUEST_NULL); }\
                                index[recv_proc] = current_index++;",
                                send_counts_name, acc_name, acc_name, acc_name, sendbufname, acc_name);
							sprintf(flow_cg_text+strlen(flow_cg_text), "\
								%s[recv_proc] = %d; }",
								send_counts_name, max_copy_level + 1); // one more for the loop number
						}
						if (options->dynschedule) {
							sprintf(flow_cg_text+strlen(flow_cg_text), "\
                                %s[recv_proc] = pack_foifi_recv_%s_%d_%d(%s,%s,%s[index[recv_proc]],%s[recv_proc]",
                                send_counts_name, acc_name_k, loop_num, l, passed_args, dest_args, sendbufname, send_counts_name);
						}
						else {
							sprintf(flow_cg_text+strlen(flow_cg_text), "\
                                %s[recv_proc] = pack_foifi_recv_%s_%d_%d(%s,%s,%s[recv_proc],%s[recv_proc]",
                                send_counts_name, acc_name_k, loop_num, l, passed_args, dest_args, sendbufname, send_counts_name);
						}
						if (options->variables_not_global && !options->data_dist) sprintf(flow_cg_text+strlen(flow_cg_text), ",%s", acc_name);

						print_data_dist_parm_call(flow_cg_text, arr);

						sprintf(flow_cg_text+strlen(flow_cg_text), "); }");

						if (options->dynschedule) {
							fprintf(packfp, "int pack_foifi_%s_%d_%d(%s,std::vector<double *>&%s,int *%s, \
								int my_rank, int nprocs, std::vector<MPI_Request>& send_reqs",
								acc_name_k, loop_num, l, decl_args, sendbufname, send_counts_name);
							fprintf(headerfp, "int pack_foifi_%s_%d_%d(%s,std::vector<double *>&%s,int *%s, \
								int my_rank, int nprocs, std::vector<MPI_Request>& send_reqs",
								acc_name_k, loop_num, l, decl_args, sendbufname, send_counts_name);
							if (options->dynschedule) {
								fprintf(packfp,",int *index, int *p_current_index");
								fprintf(headerfp,",int *index, int *p_current_index");
							}
						}
						else {
							fprintf(packfp, "int pack_foifi_%s_%d_%d(%s,double **%s,int *%s, \
								int my_rank, int nprocs",
								acc_name_k, loop_num, l, decl_args, sendbufname, send_counts_name);
							fprintf(headerfp, "int pack_foifi_%s_%d_%d(%s,double **%s,int *%s, \
								int my_rank, int nprocs",
								acc_name_k, loop_num, l, decl_args, sendbufname, send_counts_name);
						}
						if (options->variables_not_global && !options->data_dist) {
							fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
							fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
						}
						add_data_dist_parm_decl(headerfp, acc_name, prog);
						fprintf(headerfp,");\n");

						generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_cg_text, IN_FUNCTION, dest_iters, src_copy_level,
								dest_copy_level, acc_name, send_counts_name);

						fprintf(packfp, "int pack_foifi_recv_%s_%d_%d(%s,%s,double *%s,int %s",
								acc_name_k, loop_num, l, decl_args, decl_dest_args, sendbufname, send_counts_name);
						fprintf(headerfp, "int pack_foifi_recv_%s_%d_%d(%s,%s,double *%s,int %s",
								acc_name_k, loop_num, l, decl_args, decl_dest_args, sendbufname, send_counts_name);
						if (options->variables_not_global && !options->data_dist) {
							fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
							fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
						}
						add_data_dist_parm_decl(headerfp, acc_name, prog);
						fprintf(headerfp,");\n");
						generate_pack_or_unpack(packfp, prog, foifi, flow_copy_text, IN_FUNCTION, iters, total_copy_level,
								acc_nrows, acc_name, send_counts_name);

						if (options->dynschedule) {
							sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
									"%s = unpack_foifi_%s_%d_%d(",
									currdisplsname, acc_name_k, loop_num, l);
							for (i=0; i<src_copy_level; i++) {
								sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
										"(int)%s[%d],",
										recvbufname, i+1);
							}
							sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
									"%s,%s, my_rank, nprocs",
									recvbufname, currdisplsname);
						}
						else {
							sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
									"%s[proc] = unpack_foifi_%s_%d_%d(%s,%s,%s[proc], my_rank, nprocs",
									currdisplsname, acc_name_k, loop_num, l, args, recvbufname, currdisplsname);
						}
						if (options->variables_not_global && !options->data_dist) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);

						print_data_dist_parm_call(unpack_stmt_text, arr);

						sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ");");

						char flow_copyback_guard_text[1024];
						sprintf(flow_copyback_guard_text, "%s; if (recv_proc == my_rank) \
                            { %s = unpack_foifi_recv_%s_%d_%d(%s,%s,%s,%s", 
                            pi_text, currdisplsname, acc_name_k, loop_num, l, passed_args, dest_args, recvbufname, currdisplsname);
						if (options->variables_not_global && !options->data_dist) sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), ",%s", acc_name);
						if(options->data_dist){

							sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text),", %s_trans", arr->text);
							for (s = 1; s < arr->num_tiled_loops; ++s) {
								sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text),", %s_size_%d", arr->text, s + arr->first_tile_dim);
							}

							sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text),", %s_ref", arr->text);
							sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text),", list_%s", arr->text);
						}

						sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), "); }");

						fprintf(packfp, "int unpack_foifi_%s_%d_%d(%s,double *%s,int %s, int my_rank, int nprocs",
								acc_name_k, loop_num, l, decl_args, recvbufname, currdisplsname);
						fprintf(headerfp, "int unpack_foifi_%s_%d_%d(%s,double *%s,int %s, int my_rank, int nprocs",
								acc_name_k, loop_num, l, decl_args, recvbufname, currdisplsname);
						if (options->variables_not_global && !options->data_dist) {
							fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
							fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
						}
						add_data_dist_parm_decl(headerfp, acc_name, prog);
						fprintf(headerfp,");\n");
						generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_copyback_guard_text, IN_FUNCTION, dest_iters,
								src_copy_level, dest_copy_level, acc_name, currdisplsname);

						fprintf(packfp, "int unpack_foifi_recv_%s_%d_%d(%s,%s,double *%s,int %s",
								acc_name_k, loop_num, l, decl_args, decl_dest_args, recvbufname, currdisplsname);
						fprintf(headerfp, "int unpack_foifi_recv_%s_%d_%d(%s,%s,double *%s,int %s",
								acc_name_k, loop_num, l, decl_args, decl_dest_args, recvbufname, currdisplsname);
						if (options->variables_not_global && !options->data_dist) {
							fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
							fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
						}
						add_data_dist_parm_decl(headerfp, acc_name, prog);
						fprintf(headerfp,");\n");
						generate_pack_or_unpack(packfp, prog, foifi, flow_copyback_text, IN_FUNCTION, iters, total_copy_level, acc_nrows
								, acc_name, currdisplsname);

						for (i=0; i<dest_copy_level; i++) {
							free(dest_iters[i]);
						}
						free(dest_iters);

						pluto_constraints_free(fo);
						pluto_constraints_free(foifi);
						pluto_constraints_free(receiver_tiles);
					}

					sprintf(pack_stmt_text+strlen(pack_stmt_text), "} else {");
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "} else {");
				}
			}

			char sigma_text[1024];
			sprintf(sigma_text, "sigma_%s_%d(%s,%s", acc_name_k, loop_num, args, params);

			// !!!roshan should add the below keyword to parallelize per-receiver copy loop
			// this currently leads to slow down; should be investigated
			/*__omp_par_for_guided \*/
			if (options->dynschedule) {
				sprintf(pack_stmt_text+strlen(pack_stmt_text),
						"for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; } \
				%s, my_rank, nprocs, receiver_list); \
				for (__p=0; __p<nprocs; __p++) { \
				if (receiver_list[__p] != 0) { ",
				sigma_text);
			}
			else {
				sprintf(pack_stmt_text+strlen(pack_stmt_text),
						"for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; } \
				%s, my_rank, nprocs, receiver_list); \
				for (__p=0; __p<nprocs; __p++) { \
				if (receiver_list[__p] != 0) { ",
				sigma_text);
			}
			if (options->dynschedule) {
				sprintf(pack_stmt_text+strlen(pack_stmt_text), "\
                    if (%s[__p] == 0) { \
                    for (;current_index_%s<send_reqs_%s.size();current_index_%s++) { \
                    if (send_reqs_%s[current_index_%s] == MPI_REQUEST_NULL) break; \
                    int __error = MPI_Test(&send_reqs_%s[current_index_%s], &flag, MPI_STATUS_IGNORE); \
                    IF_DYNSCHEDULER_DEBUG_PRINT(\
                    if (flag) {\
                    fprintf(__debug_print_fp, \"node %%d %s send completed buf %%d\\n\", my_rank, current_index_%s); } \
                    if (__error) {\
                    char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                    MPI_Error_string(__error, msg, &len); \
                    fprintf(__debug_print_fp, \"node %%d %s send test buf %%d error %%d %%s\\n\", my_rank, current_index_%s, __error, msg); } )\
                    IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                    if (flag) break; } \
                    if (current_index_%s==send_reqs_%s.size()) { \
                    double *buf = NULL; \
                    size_t buf_size = 0; \
                    buf = (double *) polyrt_max_alloc(buf, sizeof(double)*max_num_elements_%s, &buf_size); \
                    %s.push_back(buf); \
                    IF_DYNSCHEDULER_DEBUG_PRINT(\
                    fprintf(__debug_print_fp, \"node %%d %s max send buffers %%d\\n\", my_rank, current_index_%s+1)); \
                    IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                    send_reqs_%s.push_back(MPI_REQUEST_NULL); }\
                    index_%s[__p] = current_index_%s++;",
                    send_counts_name, acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, 
                    acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, 
                    sendbufname, acc_name, acc_name, acc_name, acc_name, acc_name);
				sprintf(pack_stmt_text+strlen(pack_stmt_text), "\
					%s[__p] = %d; }",
					send_counts_name, max_copy_level + 1); // one more for the loop number
				sprintf(pack_stmt_text+strlen(pack_stmt_text),
						"%s[__p] = pack_%s_%d(%s, %s[index_%s[__p]], %s[__p]",
						send_counts_name, acc_name_k, loop_num, args, sendbufname, acc_name, send_counts_name);
			}
			else {
				sprintf(pack_stmt_text+strlen(pack_stmt_text),
						"%s[__p] = pack_%s_%d(%s, %s[__p], %s[__p]",
						send_counts_name, acc_name_k, loop_num, args, sendbufname, send_counts_name);
			}
			if (options->variables_not_global && !options->data_dist) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);

			print_data_dist_parm_call(pack_stmt_text, arr);

			sprintf(pack_stmt_text+strlen(pack_stmt_text), "); } }");

			if (options->fop_unicast_runtime && !broadcast) sprintf(pack_stmt_text+strlen(pack_stmt_text), "}");

			if(curr!= NULL && curr->deps == NULL){
				curr->constraints = pluto_constraints_empty(src_copy_level + acc_nrows + prog->npar +1);
			}
			assert(curr->constraints->ncols == src_copy_level + acc_nrows + prog->npar + 1);

			fprintf(packfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name_k, loop_num, decl_args, sendbufname, send_counts_name);
			fprintf(headerfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name_k, loop_num, decl_args, sendbufname, send_counts_name);
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
			}
			add_data_dist_parm_decl(headerfp, acc_name, prog);
			fprintf(headerfp,");\n");
			generate_pack_or_unpack(packfp, prog, curr->constraints, flow_copy_text, IN_FUNCTION, iters, src_copy_level, acc_nrows,
					acc_name, send_counts_name);

			if (options->dynschedule) {
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
						"if (is_receiver_%s_%d(",
						acc_name_k, loop_num);
				for (i=0; i<src_copy_level; i++) {
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
							"(int)%s[%d],",
							recvbufname, i+1);
				}
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
						"%s,my_rank,nprocs) != 0) { \
				%s = unpack_%s_%d(",
				params, currdisplsname, acc_name_k, loop_num);
				for (i=0; i<src_copy_level; i++) {
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
							"(int)%s[%d],",
							recvbufname, i+1);
				}
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
						"%s,%s",
						recvbufname, currdisplsname);
			}
			else {
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
						"if (is_receiver_%s_%d(%s,%s,my_rank,nprocs) != 0) { \
				%s[proc] = unpack_%s_%d(%s,%s,%s[proc]",
				acc_name_k, loop_num, args, params, currdisplsname, acc_name_k, loop_num, args, recvbufname, currdisplsname);
			}
			if (options->variables_not_global && !options->data_dist) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);

			print_data_dist_parm_call(unpack_stmt_text, arr);

			sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "); }");

			if (options->fop_unicast_runtime && !broadcast) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "}");

			fprintf(packfp, "int unpack_%s_%d(%s,double *%s,int %s", acc_name_k, loop_num, decl_args, recvbufname, currdisplsname);
			fprintf(headerfp, "int unpack_%s_%d(%s,double *%s,int %s", acc_name_k, loop_num, decl_args, recvbufname, currdisplsname);
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
			}
			add_data_dist_parm_decl(headerfp, acc_name, prog);
			fprintf(headerfp,");\n");
			generate_pack_or_unpack(packfp, prog, curr->constraints, flow_copyback_text, IN_FUNCTION, iters, src_copy_level,
					acc_nrows, acc_name,currdisplsname);

			k++;
			curr = curr->next;
		}
		if (!options->dynschedule) {
            sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "}");
		}

		/* Create statement stub to get hold of the outer loops for the copy
		 * stmt; rest of the loops to be added are the actual copy loops */
		Stmt *pack_stmt = create_helper_stmt(anchor_stmt, src_copy_level, 
                pack_stmt_text, COPY_OUT, loop_num);
		Stmt *unpack_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
                unpack_stmt_text, COPY_IN, loop_num);

		// comm stmts organized as follows
		// pack stmt, comm stmt, unpack stmt
		*num_comm_stmts = 3;
		Stmt **copy_comm_stmts = (Stmt **) malloc((*num_comm_stmts)*sizeof(Stmt *));
		copy_comm_stmts[0] = pack_stmt;
		copy_comm_stmts[1] = comm_stmt;
		copy_comm_stmts[2] = unpack_stmt;

		for (i=0; i<acc_nrows; i++) {
			free(iters[i]);
		}
		free(iters);

		free(displsname);
		free(currdisplsname);
		free(sendbufname);
		free(send_counts_name);
		free(recv_counts_name);
		free(pack_stmt_text);
		free(unpack_stmt_text);
		free(flow_copy_text);
		free(flow_copyback_text);
		free(comm_text);
		free(send_buf_size);
		free(data_tag);
		free(access);
		fclose(packfp);

		pluto_constraints_list_free(atomic_flowouts);

		return copy_comm_stmts;
}

void free_stmt_array_buffers(Stmt** buffer, int size) {

	int i = 0;
	for(i = 0; i < size; ++i) {
		free(buffer[i]);
	}

	free(buffer);

	return;
}

/*
 * Optimized communication code generation using FOIFI scheme
 * copy_level: number of outer loops to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop is normally the parallel loop)
 * wacc_stmts: all <statements, write access> of this data variable in this loop
 * num_accs: number of write accesses of this data variable in this loop
 * flow_in: flow_in set of this data variable in each loop - indexed by loop
 * This function is called per data variable
 */
Stmt **gen_comm_code_opt_foifi(int data_id, struct stmt_access_pair **wacc_stmts, int num_accs,
		int nloops, int num_data, PlutoProg *prog, Stmt *anchor_stmt, int *copy_level, int outer_dist_loop_level,
		int loop_num, int *pi_mappings, int *num_comm_stmts, FILE *headerfp)
{
	int i, j, k, l, src_copy_level, max_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];
	assert(src_copy_level>=1);
	max_copy_level = copy_level[0];
	for (i=1; i<nloops; i++) {
		if (max_copy_level < copy_level[i]) {
			max_copy_level = copy_level[i];
		}
	}


	assert(num_accs >= 1);
	char *access = reconstruct_access(wacc_stmts[0]->acc);
	acc_nrows = wacc_stmts[0]->acc->mat->nrows;
	char *acc_name = wacc_stmts[0]->acc->name;
	char *data_tag;
	Array *arr = pluto_get_corrs_array(acc_name, prog);
	data_tag = malloc(10);
	get_data_tag(prog, wacc_stmts[0]->acc, &data_tag);

	/* To be inside a loop: can't foresee other use */
	assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
			prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);

	PlutoConstraints *foifi_sets[nloops];
	for (l=0; l<nloops; l++) {
		foifi_sets[l] = NULL;
	}

	int broadcast = 0;
	/* Sender-side copying */
	PlutoConstraints *flow_out = NULL; // used only for get parametric extent/bounding box
	for (k=0; k<num_accs; k++)  {
		PlutoConstraints *flow_out_one;
		flow_out_one = compute_flow_out(wacc_stmts[k], src_copy_level, copy_level, prog, pi_mappings);
		if (flow_out == NULL) flow_out = pluto_constraints_dup(flow_out_one);
		else{
			flow_out = pluto_constraints_unionize(flow_out, flow_out_one);
		}
		pluto_constraints_free(flow_out_one);

		if (broadcast) continue;

		for (i=0; i<prog->ndeps; i++)   {
			Dep *dep = prog->deps[i];

			if (options->data_dist) {
				/* Only RAW and RAR deps matter when data is distributed*/
				if (dep->type != OSL_DEPENDENCE_RAW && dep->type != OSL_DEPENDENCE_RAR) continue;
			}
			else {
				/* Only RAW deps matter */
				if (dep->type != OSL_DEPENDENCE_RAW) continue;
			}

			// if (dep->dirvec[copy_level+1] == DEP_ZERO) continue;

			/* If the dependence doesn't originate from this access */
			if (dep->src_acc != wacc_stmts[k]->acc) continue;

			assert(dep->dest_acc != NULL);

			Stmt *dest = prog->stmts[dep->dest];
			int dependent_loop = pi_mappings[dest->id];
			if (dependent_loop == -1) {
				// destination statement will be executed by all processors
				broadcast = 1;
				break;
			}

			int dest_copy_level = copy_level[dependent_loop];

			PlutoConstraints *fo = compute_flow_out_of_dep(dep, src_copy_level, copy_level, prog, 0, NULL, pi_mappings);
			// add target tile iterators for flow-out set
			for (j=0; j<dest_copy_level; j++) {
				pluto_constraints_add_dim(fo,src_copy_level, NULL);
			}

			PlutoConstraints *fi = compute_flow_in_of_dep(dep, dest_copy_level, prog, 0);
			// add source tile iterators for flow-in set
			for (j=0; j<src_copy_level; j++) {
				pluto_constraints_add_dim(fi,0, NULL);
			}

			PlutoConstraints *foifi = pluto_constraints_intersection(fo, fi);

			if (foifi_sets[dependent_loop] == NULL)
				foifi_sets[dependent_loop] = pluto_constraints_dup(foifi);
			else
				pluto_constraints_unionize(foifi_sets[dependent_loop], foifi);

			pluto_constraints_free(fo);
			pluto_constraints_free(fi);
			pluto_constraints_free(foifi);
		}
	}

	// === from here on we begin the code generation common for all dependent loops===================================
	//

	// !!!roshan may need to increase send buffer size to an arbitrarily large number
	// due to the possible amount of duplication in the FOIFI scheme
	char *send_buf_size =
			get_parametric_bounding_box(flow_out, src_copy_level, acc_nrows,
					prog->npar, (const char **)prog->params);
	sprintf(send_buf_size+strlen(send_buf_size), "*__FOIFI_MAX_DUP_FACTOR");

	if (!(options->dynschedule)) {
		PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
		char *total_extent = malloc(1024);
		strcpy(total_extent, "1");
		/* Assumes load-balanced distribution */
		for (i=outer_dist_loop_level; i<src_copy_level; i++) {
			char *extent;
			/* Just the first one in flow_out is enough (rest all should give the
			 * same since they are all under the same parallel loop and each
			 * iteration of the parallel loop writes to distinct data) */
			get_parametric_extent_const(anchor_stmt_new_dom, i, prog->npar,
					(const char **)prog->params, &extent, NULL);
			/* The + nprocs is needed since the send buffer size should be larger
			 * when some processors have more iterations than others - some processors
			 * have an extra iteration in each dimension in the worst case */
			sprintf(total_extent+strlen(total_extent), "*(%s+nprocs)", extent);
			free(extent);
		}
		sprintf(send_buf_size+strlen(send_buf_size),
				"*floorf((%s)/(float)nprocs)", total_extent);
		free(total_extent);
		pluto_constraints_free(anchor_stmt_new_dom);
	}
	IF_DEBUG(printf("Send buffer size for %s: %s\n", acc_name, send_buf_size););

	char *sendbufname = concat("send_buf_", acc_name);
	char *send_counts_name = concat("send_counts_", acc_name);
	char *max_elements_name = concat("max_num_elements_", acc_name);

	if (!options->dynschedule) {
        sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
                %s[__p] = 0;\n}\n", send_counts_name);
	}

	if (options->dynschedule) { // one more for the loop number
		sprintf(prog->decls+strlen(prog->decls), "\
            %s = max(%s, %s + %d);\n", max_elements_name, max_elements_name, send_buf_size, max_copy_level+1);
		sprintf(prog->decls+strlen(prog->decls), "\
                for (__tid=0; __tid<_num_threads; __tid++) { \
                for (__p=0; __p<__MAX_NUM_SENDS; __p++) {\n\
                %s[__tid][__p] = (double *) polyrt_max_alloc(%s[__tid][__p], \
                sizeof(double)*%s, &send_buf_size_%s[__tid][__p]);\n}\n}\n",
                sendbufname, sendbufname, max_elements_name, acc_name);
	}
	else {
		sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
				%s[__p]\
				= (double *) polyrt_max_alloc(%s[__p], sizeof(double)*(%s), &send_buf_size_%s[__p]);\n}",
				sendbufname, sendbufname, send_buf_size, acc_name);
	}

	char *flow_copy_text = malloc(strlen(sendbufname) + strlen("[")+
			strlen(send_counts_name) + strlen("++] = ") + strlen(access) + 1);
	sprintf(flow_copy_text, "%s[%s++] = %s", sendbufname, send_counts_name, access);

	char *displsname = concat("displs_", acc_name);
	char *currdisplsname = concat("curr_displs_", acc_name);
	char *recvbufname = concat("recv_buf_", acc_name);
	char *recv_counts_name = concat("recv_counts_", acc_name);

	char args[1024];
	strcpy(args,"");
	sprintf(args+strlen(args), "t%d", 1);
	/* make it src_copy_level+1 since we use an extra dimension to separate
	 * statements */
	for (i=1; i<src_copy_level; i++) {
		sprintf(args+strlen(args), ",t%d", i+1);
	}

	char passed_args[1024];
	strcpy(passed_args,"");
	sprintf(passed_args+strlen(passed_args), "ts%d", 1);
	for (i=1; i<src_copy_level; i++) {
		sprintf(passed_args+strlen(passed_args), ",ts%d", i+1);
	}

	char decl_args[1024];
	strcpy(decl_args,"");
	sprintf(decl_args+strlen(decl_args), "int ts%d", 1);
	for (i=1; i<src_copy_level; i++) {
		sprintf(decl_args+strlen(decl_args), ",int ts%d", i+1);
	}

	char params[1024];
	strcpy(params, "");
	if (prog->npar>=1) {
		sprintf(params+strlen(params), "%s", prog->params[0]);
		for (i=1; i<prog->npar; i++) {
			sprintf(params+strlen(params), ",%s", prog->params[i]);
		}
	}

	char **iters;
	iters = malloc(acc_nrows * sizeof(char *));
	for (i=0; i < acc_nrows; i++) {
		iters[i] = malloc(5);
		sprintf(iters[i], "d%d", i+1);
	}

	char *comm_text = malloc(2048);

	if (options->dynschedule) {
		strcpy(comm_text, "");
		sprintf(comm_text+strlen(comm_text),  "\
				for (__p=0; __p<nprocs; __p++) {\
				if(%s[__p] > 0) {",
				send_counts_name);
		sprintf(comm_text+strlen(comm_text), "%s[index_%s[__p]][0] = %d;", sendbufname, acc_name, loop_num);
		for (i=0; i<src_copy_level; i++) {
			sprintf(comm_text+strlen(comm_text), "%s[index_%s[__p]][%d] = t%d;", sendbufname, acc_name, i+1, i+1);
		}
		sprintf(comm_text+strlen(comm_text),  "\
		        IF_DYNSCHEDULER_DEBUG_PRINT(\
		        fprintf(__debug_print_fp, \"sending node %%d buf %%d loop %d task ",
		        loop_num);
		for (i=0; i<src_copy_level; i++) {
			sprintf(comm_text+strlen(comm_text),  "%%d ");
		}
		sprintf(comm_text+strlen(comm_text),  "\
		        %s size %%d to node %%d\\n\", \
		        my_rank, index_%s[__p], %s, %s[__p], __p)); \
                IF_DEBUG_FLUSH(fflush(__debug_print_fp));",
                acc_name, acc_name, args, send_counts_name);
		if (options->timereport)
			sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s[__p]);", send_counts_name);
		sprintf(comm_text+strlen(comm_text),  "\
				MPI_Isend(%s[index_%s[__p]], %s[__p], MPI_DOUBLE,\
					__p, %s, MPI_COMM_WORLD, &send_reqs_%s[index_%s[__p]]);\
				%s[__p] = 0;\
                IF_DYNSCHEDULER_DEBUG_PRINT(\
                if (__error) {\
                char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                MPI_Error_string(__error, msg, &len); \
                fprintf(__debug_print_fp, \"node %%d %s send buf %%d error %%d %%s\\n\", my_rank, index_%s[__p], __error, msg); } )\
                }}",
                sendbufname, acc_name, send_counts_name, data_tag, acc_name, acc_name, send_counts_name, acc_name, acc_name);
	}
	else {
		/* Message passing code (MPI calls) */
		if (options->timereport)
			sprintf(comm_text,  "IF_TIME(t_comm_start = rtclock());");
		else
			strcpy(comm_text, "");
		sprintf(comm_text+strlen(comm_text), "\
				MPI_Alltoall(%s, 1, MPI_INT,\
					%s, 1, MPI_INT, MPI_COMM_WORLD);",
					send_counts_name, recv_counts_name);
		sprintf(comm_text+strlen(comm_text),  "\
				req_count=0;\
				for (__p=0; __p<nprocs; __p++) {\
				if (%s[__p] >= 1) {\
				assert(\"increase __FOIFI_MAX_DUP_FACTOR\" && (%s[__p]*8 <= send_buf_size_%s[__p]));",
				send_counts_name, send_counts_name, acc_name);
		if (options->timereport)
			sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s[__p]);", send_counts_name);
		sprintf(comm_text+strlen(comm_text),  "\
				MPI_Isend(%s[__p], %s[__p], MPI_DOUBLE,\
					__p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}",
					sendbufname, send_counts_name);
		sprintf(comm_text+strlen(comm_text),  "for (__p=0; __p<nprocs; __p++) {\
				if(%s[__p] >= 1) {\
				MPI_Irecv(%s+%s[__p], %s[__p], MPI_DOUBLE,\
					__p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}\
				MPI_Waitall(req_count, reqs, stats);\
				for (__p=0; __p<nprocs; __p++) {\
				%s[__p] = 0;}",
				recv_counts_name, recvbufname, displsname, recv_counts_name,
				send_counts_name);
		sprintf(comm_text+strlen(comm_text), "for (__p=0; __p<nprocs; __p++) {\
				%s[__p] = %s[__p]; }", currdisplsname, displsname);
		if (options->timereport)
			sprintf(comm_text+strlen(comm_text),  "IF_TIME(t_comm += rtclock() - t_comm_start);");
	}

	/* Create statement stub to get hold of the outer loops for the copy
	 * stmt; rest of the loops to be added are the actual copy loops */
	Stmt *comm_stmt = create_helper_stmt(anchor_stmt, outer_dist_loop_level, comm_text, COMM_CALL, loop_num);

	/* Receiver-side copy */
	/* With send/recv -based more exact communication scheme, recv side copy
	 * does not include all the sender-side copies; processes that do not send
	 * data are to be skipped */
	if (options->dynschedule) {
		sprintf(prog->decls+strlen(prog->decls),
				"for (__p=0; __p<__MAX_NUM_RECVS; __p++) { \
                %s[__p] = (double *) polyrt_max_alloc(%s[__p], sizeof(double)*%s, &recv_buf_size_%s[__p]); }\n",
                recvbufname, recvbufname, max_elements_name, acc_name);
	}
	else {
		sprintf(prog->decls+strlen(prog->decls),
				"%s = (double *) polyrt_max_alloc(%s, nprocs*send_buf_size_%s[0], &recv_buf_size_%s);\n",
				recvbufname, recvbufname, acc_name, acc_name);
	}
	if (!options->dynschedule) {
		sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
                %s[__p] = __p*send_buf_size_%s[__p]/sizeof(double);}\n\n", displsname, acc_name);
	}

	char *flow_copyback_text = malloc(strlen(access) + strlen(" = ") + strlen(recvbufname)
			+ strlen("[") + strlen(currdisplsname) + strlen("++]") + 1);
	sprintf(flow_copyback_text, "%s = %s[%s++]",
			access, recvbufname, currdisplsname);

	FILE *packfp = fopen(__APPENDFILENAME, "a");
	assert(packfp != NULL);

	char *pack_stmt_text = malloc(20480);
	strcpy(pack_stmt_text,"");
	char *unpack_stmt_text = malloc(20480);
	if (options->dynschedule) {
		strcpy(unpack_stmt_text, "");
		// one more for the loop number
		sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
				"%s = %d; ", currdisplsname, max_copy_level + 1);
	}
	else {
		sprintf(unpack_stmt_text, "proc = pi_%d(%s,%s, nprocs); ", loop_num, args, params);
		sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "if ((my_rank != proc) && (%s[proc] > 0)) { ", recv_counts_name);
	}
	if (broadcast) {
		if ((flow_out != NULL) && !pluto_constraints_is_empty(flow_out)) {
			sprintf(pack_stmt_text+strlen(pack_stmt_text),
					"for (__p=0; __p<nprocs; __p++) \
                    { if (__p != my_rank) { ");
			if (options->dynschedule) {
				sprintf(pack_stmt_text+strlen(pack_stmt_text), "\
                        if (%s[__p] == 0) { \
                        for (;current_index_%s<send_reqs_%s.size();current_index_%s++) { \
                        if (send_reqs_%s[current_index_%s] == MPI_REQUEST_NULL) break; \
                        int __error = MPI_Test(&send_reqs_%s[current_index_%s], &flag, MPI_STATUS_IGNORE); \
                        IF_DYNSCHEDULER_DEBUG_PRINT(\
                        if (flag) {\
                        fprintf(__debug_print_fp, \"node %%d %s send completed buf %%d\\n\", my_rank, current_index_%s); } \
                        if (__error) {\
                        char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                        MPI_Error_string(__error, msg, &len); \
                        fprintf(__debug_print_fp, \"node %%d %s send test buf %%d error %%d %%s\\n\", my_rank, current_index_%s, __error, msg); } )\
                        IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                        if (flag) break; } \
                        if (current_index_%s==send_reqs_%s.size()) { \
                        double *buf = NULL; \
                        size_t buf_size = 0; \
                        buf = (double *) polyrt_max_alloc(buf, sizeof(double)*max_num_elements_%s, &buf_size); \
                        %s.push_back(buf); \
                        IF_DYNSCHEDULER_DEBUG_PRINT(\
                        fprintf(__debug_print_fp, \"node %%d %s max send buffers %%d\\n\", my_rank, current_index_%s+1)); \
                        IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                        send_reqs_%s.push_back(MPI_REQUEST_NULL); }\
                        index_%s[__p] = current_index_%s++;",
                        send_counts_name, acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, 
                        acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, acc_name, 
                        sendbufname, acc_name, acc_name, acc_name, acc_name, acc_name);
				sprintf(pack_stmt_text+strlen(pack_stmt_text), "\
						%s[__p] = %d; }",
						send_counts_name, max_copy_level + 1); // one more for the loop number
				sprintf(pack_stmt_text+strlen(pack_stmt_text),
						"%s[__p] = pack_%s_%d(%s,%s[index_%s[__p]],%s[__p]",
						send_counts_name, acc_name, loop_num, args, sendbufname, acc_name, send_counts_name);
			}
			else {
				sprintf(pack_stmt_text+strlen(pack_stmt_text),
						"%s[__p] = pack_%s_%d(%s,%s[__p],%s[__p]",
						send_counts_name, acc_name, loop_num, args, sendbufname, send_counts_name);
			}
			if (options->variables_not_global && !options->data_dist) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);
			print_data_dist_parm_call(pack_stmt_text, arr);

			sprintf(pack_stmt_text+strlen(pack_stmt_text), "); } }");

			fprintf(packfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, sendbufname, send_counts_name);
			fprintf(headerfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, sendbufname, send_counts_name);
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
			}
			add_data_dist_parm_decl(headerfp, acc_name, prog);
			fprintf(headerfp,");\n");
			generate_pack_or_unpack(packfp, prog, flow_out, flow_copy_text, IN_FUNCTION, iters, src_copy_level, acc_nrows,
					acc_name, send_counts_name);

			if (options->dynschedule) {
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
						"%s = unpack_%s_%d(",
						currdisplsname, acc_name, loop_num);
				for (i=0; i<src_copy_level; i++) {
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
							"(int)%s[%d],",
							recvbufname, i+1);
				}
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
						"%s,%s",
						recvbufname, currdisplsname);
			}
			else {
				sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
						"%s[proc] = unpack_%s_%d(%s,%s,%s[proc]",
						currdisplsname, acc_name, loop_num, args, recvbufname, currdisplsname);
			}
			if (options->variables_not_global && !options->data_dist) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);

			print_data_dist_parm_call(unpack_stmt_text, arr);

			sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ");");

			fprintf(packfp, "int unpack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, recvbufname, currdisplsname);
			fprintf(headerfp, "int unpack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, recvbufname, currdisplsname);
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
			}
			add_data_dist_parm_decl(headerfp, acc_name, prog);
			fprintf(headerfp,");\n");
			generate_pack_or_unpack(packfp, prog, flow_out, flow_copyback_text, IN_FUNCTION, iters, src_copy_level, acc_nrows,
					acc_name, "count");
		}
	}
	else {
		for (l=0; l<nloops; l++) {
			PlutoConstraints *foifi = foifi_sets[l];
			if (foifi == NULL) continue;
			if (pluto_constraints_is_empty(foifi)) {
				pluto_constraints_free(foifi);
				continue;
			}

			int dest_copy_level = copy_level[l];
			assert(dest_copy_level>=1);
			int total_copy_level = src_copy_level + dest_copy_level;
			assert(foifi->ncols == total_copy_level + acc_nrows + prog->npar + 1);

			PlutoConstraints *receiver_tiles = get_receiver_tiles(wacc_stmts, num_accs, src_copy_level, dest_copy_level, prog, l, pi_mappings);
			assert(receiver_tiles->ncols == total_copy_level + prog->npar + 1);
			assert(!pluto_constraints_is_empty(receiver_tiles));

			IF_DEBUG(printf("Data flow out intersection flow in set for %s and dependent loop %d\n", acc_name, l););
			IF_DEBUG(pluto_constraints_print(stdout, foifi));

			// === from here on we begin the code generation for each dependent loop==========================================
			//
				char dest_args[1024];
				strcpy(dest_args,"");
				sprintf(dest_args+strlen(dest_args), "t%d", src_copy_level+1);
				for (i=src_copy_level+1; i<total_copy_level; i++) {
					sprintf(dest_args+strlen(dest_args), ",t%d", i+1);
				}

				char decl_dest_args[1024];
				strcpy(decl_dest_args,"");
				sprintf(decl_dest_args+strlen(decl_dest_args), "int ts%d", src_copy_level+1);
				for (i=src_copy_level+1; i<total_copy_level; i++) {
					sprintf(decl_dest_args+strlen(decl_dest_args), ",int ts%d", i+1);
				}

				char **dest_iters;
				dest_iters = malloc(dest_copy_level * sizeof(char *));
				for (i=0; i<dest_copy_level; i++) {
					dest_iters[i] = malloc(5);
					sprintf(dest_iters[i], "t%d", i+src_copy_level+1);
				}

				char pi_text[1024];
				sprintf(pi_text, "recv_proc = pi_%d(%s,%s, nprocs)", l, dest_args, params);

				if (options->dynschedule) {
					sprintf(pack_stmt_text+strlen(pack_stmt_text),
							"pack_%s_%d_%d(%s,%s,%s, my_rank, nprocs, send_reqs_%s, index_%s, &current_index_%s",
							acc_name, loop_num, l, args, sendbufname, send_counts_name, acc_name, acc_name, acc_name);
				}
				else {
					sprintf(pack_stmt_text+strlen(pack_stmt_text),
							"pack_%s_%d_%d(%s,%s,%s, my_rank, nprocs",
							acc_name, loop_num, l, args, sendbufname, send_counts_name);
				}
				if (options->variables_not_global && !options->data_dist) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);
				print_data_dist_parm_call(pack_stmt_text, arr);

				sprintf(pack_stmt_text+strlen(pack_stmt_text), ");");

				char flow_cg_text[2048];
				sprintf(flow_cg_text, "%s; if (recv_proc != my_rank) { ",
						pi_text);
				if (options->dynschedule) {
					sprintf(flow_cg_text+strlen(flow_cg_text), "\
                        if (%s[recv_proc] == 0) { \
                        for (;current_index<send_reqs.size();current_index++) { \
                        if (send_reqs[current_index] == MPI_REQUEST_NULL) break; \
                        int __error = MPI_Test(&send_reqs[current_index], &flag, MPI_STATUS_IGNORE); \
                        IF_DYNSCHEDULER_DEBUG_PRINT(\
                        if (flag) {\
                        fprintf(__debug_print_fp, \"node %%d %s send completed buf %%d\\n\", my_rank, current_index); } \
                        if (__error) {\
                        char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                        MPI_Error_string(__error, msg, &len); \
                        fprintf(__debug_print_fp, \"node %%d %s send test buf %%d error %%d %%s\\n\", my_rank, current_index, __error, msg); } )\
                        IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                        if (flag) break; } \
                        if (current_index==send_reqs.size()) { \
                        double *buf = NULL; \
                        size_t buf_size = 0; \
                        buf = (double *) polyrt_max_alloc(buf, sizeof(double)*max_num_elements_%s, &buf_size); \
                        %s.push_back(buf); \
                        IF_DYNSCHEDULER_DEBUG_PRINT(\
                        fprintf(__debug_print_fp, \"node %%d %s max send buffers %%d\\n\", my_rank, current_index+1)); \
                        IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                        send_reqs.push_back(MPI_REQUEST_NULL); }\
                        index[recv_proc] = current_index++;",
                        send_counts_name, acc_name, acc_name, acc_name, sendbufname, acc_name);
					sprintf(flow_cg_text+strlen(flow_cg_text), "\
                        %s[recv_proc] = %d; }",
                        send_counts_name, max_copy_level + 1); // one more for the loop number
				}
				if (options->dynschedule) {
					sprintf(flow_cg_text+strlen(flow_cg_text), "\
                        %s[recv_proc] = pack_recv_%s_%d_%d(%s,%s,%s[index[recv_proc]],%s[recv_proc]",
                        send_counts_name, acc_name, loop_num, l, passed_args, dest_args, sendbufname, send_counts_name);
				}
				else {
					sprintf(flow_cg_text+strlen(flow_cg_text), "\
                        %s[recv_proc] = pack_recv_%s_%d_%d(%s,%s,%s[recv_proc],%s[recv_proc]",
                        send_counts_name, acc_name, loop_num, l, passed_args, dest_args, sendbufname, send_counts_name);
				}
				if (options->variables_not_global && !options->data_dist) sprintf(flow_cg_text+strlen(flow_cg_text), ",%s", acc_name);

				print_data_dist_parm_call(flow_cg_text, arr);

				sprintf(flow_cg_text+strlen(flow_cg_text), "); }");

				if (options->dynschedule) {
					fprintf(packfp, "int pack_%s_%d_%d(%s,std::vector<double *>&%s,int *%s, \
						int my_rank, int nprocs, std::vector<MPI_Request>& send_reqs",
						acc_name, loop_num, l, decl_args, sendbufname, send_counts_name);
					fprintf(headerfp, "int pack_%s_%d_%d(%s,std::vector<double *>&%s,int *%s, \
						int my_rank, int nprocs, std::vector<MPI_Request>& send_reqs",
						acc_name, loop_num, l, decl_args, sendbufname, send_counts_name);
					if (options->dynschedule) {
						fprintf(packfp,",int *index,int *p_current_index");
						fprintf(headerfp,",int *index,int *p_current_index");
					}
				}
				else {
					fprintf(packfp, "int pack_%s_%d_%d(%s,double **%s,int *%s, \
						int my_rank, int nprocs",
						acc_name, loop_num, l, decl_args, sendbufname, send_counts_name);
					fprintf(headerfp, "int pack_%s_%d_%d(%s,double **%s,int *%s, \
						int my_rank, int nprocs",
						acc_name, loop_num, l, decl_args, sendbufname, send_counts_name);
				}
				if (options->variables_not_global && !options->data_dist) {
					fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
					fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
				}
				add_data_dist_parm_decl(headerfp, acc_name, prog);
				fprintf(headerfp,");\n");

				generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_cg_text, IN_FUNCTION, dest_iters, src_copy_level, dest_copy_level,
						acc_name, send_counts_name);

				fprintf(packfp, "int pack_recv_%s_%d_%d(%s,%s,double *%s,int %s",
						acc_name, loop_num, l, decl_args, decl_dest_args, sendbufname, send_counts_name);
				fprintf(headerfp, "int pack_recv_%s_%d_%d(%s,%s,double *%s,int %s",
						acc_name, loop_num, l, decl_args, decl_dest_args, sendbufname, send_counts_name);
				if (options->variables_not_global && !options->data_dist) {
					fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
					fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
				}
				add_data_dist_parm_decl(headerfp, acc_name, prog);
				fprintf(headerfp,");\n");
				generate_pack_or_unpack(packfp, prog, foifi, flow_copy_text, IN_FUNCTION, iters, total_copy_level, acc_nrows,
						acc_name, send_counts_name);

				if (options->dynschedule) {
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
							"%s = unpack_%s_%d_%d(",
							currdisplsname, acc_name, loop_num, l);
					for (i=0; i<src_copy_level; i++) {
						sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
								"(int)%s[%d],",
								recvbufname, i+1);
					}
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
							"%s,%s, my_rank, nprocs",
							recvbufname, currdisplsname);
				}
				else {
					sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "\
						%s[proc] = unpack_%s_%d_%d(%s,%s,%s[proc], my_rank, nprocs",
						currdisplsname, acc_name, loop_num, l, args, recvbufname, currdisplsname);
				}
				if (options->variables_not_global && !options->data_dist) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);

				print_data_dist_parm_call(unpack_stmt_text, arr);

				sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ");");

				char flow_copyback_guard_text[2048];
				sprintf(flow_copyback_guard_text, "%s; if (recv_proc == my_rank) \
                    { %s = unpack_recv_%s_%d_%d(%s,%s,%s,%s", 
                    pi_text, currdisplsname, acc_name, loop_num, l, passed_args, dest_args, recvbufname, currdisplsname);
				if (options->variables_not_global && !options->data_dist) sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), ",%s", acc_name);

				print_data_dist_parm_call(flow_copyback_guard_text, arr);

				sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), "); }");

				fprintf(packfp, "int unpack_%s_%d_%d(%s,double *%s,int %s, int my_rank, int nprocs", acc_name, loop_num, l, decl_args, recvbufname, currdisplsname);
				fprintf(headerfp, "int unpack_%s_%d_%d(%s,double *%s,int %s, int my_rank, int nprocs", acc_name, loop_num, l, decl_args, recvbufname, currdisplsname);
				if (options->variables_not_global && !options->data_dist) {
					fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
					fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
				}
				add_data_dist_parm_decl(headerfp, acc_name, prog);
				fprintf(headerfp,");\n");
				fprintf(packfp, "\nint recv_proc;\n");
				generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_copyback_guard_text, IN_FUNCTION, dest_iters, src_copy_level, dest_copy_level,
						acc_name, "count");

				fprintf(packfp, "int unpack_recv_%s_%d_%d(%s,%s,double *%s,int %s", acc_name, loop_num, l, decl_args, decl_dest_args, recvbufname, currdisplsname);
				fprintf(headerfp, "int unpack_recv_%s_%d_%d(%s,%s,double *%s,int %s", acc_name, loop_num, l, decl_args, decl_dest_args, recvbufname, currdisplsname);
				if (options->variables_not_global && !options->data_dist) {
					fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
					fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
				}
				add_data_dist_parm_decl(headerfp, acc_name, prog);
				fprintf(headerfp,");\n");
				generate_pack_or_unpack(packfp, prog, foifi, flow_copyback_text, IN_FUNCTION, iters, total_copy_level, acc_nrows,
						acc_name, "count");

				for (i=0; i<dest_copy_level; i++) {
					free(dest_iters[i]);
				}
				free(dest_iters);

				pluto_constraints_free(foifi);
				pluto_constraints_free(receiver_tiles);

		}
	}
	if (!options->dynschedule) {
        sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "}");
	}

	/* Create statement stub to get hold of the outer loops for the copy
	 * stmt; rest of the loops to be added are the actual copy loops */
	Stmt *pack_stmt = create_helper_stmt(anchor_stmt, src_copy_level, 
            pack_stmt_text, COPY_OUT, loop_num);
	Stmt *unpack_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
            unpack_stmt_text, COPY_IN, loop_num);

	// comm stmts organized as follows
	// pack stmt, comm stmt, unpack stmt
	*num_comm_stmts = 3;
	Stmt **copy_comm_stmts = (Stmt **) malloc((*num_comm_stmts)*sizeof(Stmt *));
	copy_comm_stmts[0] = pack_stmt;
	copy_comm_stmts[1] = comm_stmt;
	copy_comm_stmts[2] = unpack_stmt;

	for (i=0; i<acc_nrows; i++) {
		free(iters[i]);
	}
	free(iters);

	free(pack_stmt_text);
	free(unpack_stmt_text);
	free(flow_copy_text);
	free(flow_copyback_text);
	free(comm_text);

	free(displsname);
	free(sendbufname);
	free(send_counts_name);
	free(recv_counts_name);
	free(send_buf_size);
	free(data_tag);
	fclose(packfp);

	pluto_constraints_free(flow_out);
	free(access);
	return copy_comm_stmts;
}


/*
 * Optimized communication code generation
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop is normally the parallel loop)
 * wacc_stmts:  all <statements, wacc> writing to this data variable
 * This function is called per data variable
 */
Stmt **gen_comm_code_opt(int data_id, struct stmt_access_pair **wacc_stmts, int num_accs,
		int nloops, int num_data, PlutoProg *prog, Stmt *anchor_stmt, int *copy_level, int outer_dist_loop_level,
		int loop_num, int *pi_mappings, int *num_comm_stmts, FILE *headerfp)
{
	int i, k, src_copy_level, max_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];
	assert(src_copy_level>=1);
	max_copy_level = copy_level[0];
	for (i=1; i<nloops; i++) {
		if (max_copy_level < copy_level[i]) {
			max_copy_level = copy_level[i];
		}
	}

	assert(num_accs >= 1);
	char *access = reconstruct_access(wacc_stmts[0]->acc);
	acc_nrows = wacc_stmts[0]->acc->mat->nrows;
	char *data_tag;
	data_tag = malloc(10);
	get_data_tag(prog, wacc_stmts[0]->acc, &data_tag);

	char *acc_name = wacc_stmts[0]->acc->name;
	Array *arr = pluto_get_corrs_array(acc_name, prog);
	assert(num_accs >= 1);


	// printf("%d accesses: %s\n", num_accs, acc_name);

	/* To be inside a loop: can't foresee other use */
	assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
			prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);

	/* Sender-side copying */
	PlutoConstraints *flow_out = NULL;
	for (k=0; k<num_accs; k++)  {
		PlutoConstraints *flow_out_one;
		flow_out_one = compute_flow_out(wacc_stmts[k], src_copy_level, copy_level, prog, pi_mappings);
		if (flow_out == NULL) flow_out = pluto_constraints_dup(flow_out_one);
		else{
			flow_out = pluto_constraints_unionize(flow_out, flow_out_one);
		}
		pluto_constraints_free(flow_out_one);
	}
	IF_DEBUG(print_polylib_visual_sets("data", flow_out));

	// generate_tau(wacc_stmts, num_accs, copy_level, prog);

	IF_DEBUG(printf("Data flow out set for %s\n", acc_name););
	IF_DEBUG(pluto_constraints_print(stdout, flow_out));
	assert(flow_out->ncols == src_copy_level + acc_nrows + prog->npar + 1);

	// ==== from here onwards we start the code generation ============================================
	//

		char *send_buf_size =
				get_parametric_bounding_box(flow_out, src_copy_level, acc_nrows,
						prog->npar, (const char **)prog->params);

		if (!(options->dynschedule)) {
			PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
			char *total_extent = malloc(1024);
			strcpy(total_extent, "1");
			/* Assumes load-balanced distribution */
			for (i=outer_dist_loop_level; i<src_copy_level; i++) {
				char *extent;
				/* Just the first one in flow_out is enough (rest all should give the
				 * same since they are all under the same parallel loop and each
				 * iteration of the parallel loop writes to distinct data) */
				get_parametric_extent_const(anchor_stmt_new_dom, i, prog->npar,
						(const char **)prog->params, &extent, NULL);
				/* The + nprocs is needed since the send buffer size should be larger
				 * when some processors have more iterations than others - some processors
				 * have an extra iteration in each dimension in the worst case */
				sprintf(total_extent+strlen(total_extent), "*(%s+nprocs)", extent);
				free(extent);
			}
			sprintf(send_buf_size+strlen(send_buf_size),
					"*floorf((%s)/(float)nprocs)", total_extent);
			free(total_extent);
			pluto_constraints_free(anchor_stmt_new_dom);
		}
		IF_DEBUG(printf("Send buffer size for %s: %s\n", acc_name, send_buf_size););

		FILE *packfp = fopen(__APPENDFILENAME, "a");
		assert(packfp != NULL);

		char args[1024];
		strcpy(args,"");
		sprintf(args+strlen(args), "t%d", 1);
		/* make it src_copy_level+1 since we use an extra dimension to separate
		 * statements */
		for (i=1; i<src_copy_level; i++) {
			sprintf(args+strlen(args), ",t%d", i+1);
		}

		char decl_args[1024];
		strcpy(decl_args,"");
		sprintf(decl_args+strlen(decl_args), "int ts%d", 1);
		for (i=1; i<src_copy_level; i++) {
			sprintf(decl_args+strlen(decl_args), ",int ts%d", i+1);
		}

		char params[1024];
		strcpy(params, "");
		if (prog->npar>=1) {
			sprintf(params+strlen(params), "%s", prog->params[0]);
			for (i=1; i<prog->npar; i++) {
				sprintf(params+strlen(params), ",%s", prog->params[i]);
			}
		}

		char **iters;
		iters = malloc(acc_nrows * sizeof(char *));
		for (i=0; i < acc_nrows; i++) {
			iters[i] = malloc(5);
			sprintf(iters[i], "d%d", i+1);
		}

		/* Sender-side copy */
		char *sendbufname = concat("send_buf_", acc_name);
		char *send_counts_name = concat("send_counts_", acc_name);
		char *send_count_name = concat("send_count_", acc_name);
		char *max_elements_name = concat("max_num_elements_", acc_name);

		if (options->dynschedule) {
			// one more for the loop number
			sprintf(prog->decls+strlen(prog->decls), "\
            %s = max(%s, %s + %d);\n", max_elements_name, max_elements_name, send_buf_size, max_copy_level+1);
			sprintf(prog->decls+strlen(prog->decls), "\
                for (__tid=0; __tid<_num_threads; __tid++) { \
                for (_i=0; _i<__MAX_NUM_SENDS; _i++) {\n\
                %s[__tid][_i] = (double *) polyrt_max_alloc(%s[__tid][_i], sizeof(double)*%s, &send_buf_size_%s[__tid][_i]);\n}\n}\n",
                sendbufname, sendbufname, max_elements_name, acc_name);
		}
		else {
			sprintf(prog->decls+strlen(prog->decls), "%s\
				= (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &send_buf_size_%s);\n",
				sendbufname, sendbufname, send_buf_size, acc_name);
		}

		char *flow_copy_text = malloc(strlen(sendbufname) + strlen("[")+
				strlen(send_count_name) + strlen("++] = ") + strlen(access) + 1);
		sprintf(flow_copy_text, "%s[%s++] = %s", sendbufname, send_count_name, access);

		char *sigma_text = malloc(1024);
		sprintf(sigma_text, "sigma_%s_%d(%s,%s", acc_name, loop_num, args, params);

		char *sigma_sender_text = malloc(1024);
		strcpy(sigma_sender_text, sigma_text);
        sprintf(sigma_sender_text+strlen(sigma_sender_text), ", my_rank, nprocs, receiver_list)");

		char *flow_cg_text = malloc(2048);
		strcpy(flow_cg_text, "");

		if ((flow_out != NULL) && !pluto_constraints_is_empty(flow_out)) {
            sprintf(flow_cg_text+strlen(flow_cg_text), "\
                for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; } \
                %s; \
                for (__p=0; __p<nprocs; __p++) { if (receiver_list[__p] != 0) { ",
                sigma_sender_text);
			if (options->dynschedule) {
				sprintf(flow_cg_text+strlen(flow_cg_text), "\
                for (_i=0; _i<%s.size(); _i++) {\
                if (send_reqs_%s[_i] == MPI_REQUEST_NULL) break; \
                int __error = MPI_Testall(nprocs, &send_reqs_%s[_i*nprocs], &flag, MPI_STATUSES_IGNORE); \
                IF_DYNSCHEDULER_DEBUG_PRINT(\
                if (flag) {\
                fprintf(__debug_print_fp, \"node %%d %s send completed buf %%d\\n\", my_rank, _i); } \
                if (__error) {\
                char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                MPI_Error_string(__error, msg, &len); \
                fprintf(__debug_print_fp, \"node %%d %s send test buf %%d error %%d %%s\\n\", my_rank, _i, __error, msg); } )\
                IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                if (flag) break; } \
                if (_i == %s.size()) { \
                double *buf = NULL; \
                size_t buf_size = 0; \
                buf = (double *) polyrt_max_alloc(buf, sizeof(double)*max_num_elements_%s, &buf_size); \
                %s.push_back(buf); \
                IF_DYNSCHEDULER_DEBUG_PRINT(\
                fprintf(__debug_print_fp, \"node %%d %s max send buffers %%d\\n\", my_rank, _i+1)); \
                IF_DEBUG_FLUSH(fflush(__debug_print_fp));\
                for (__p=0; __p<nprocs; __p++) send_reqs_%s.push_back(MPI_REQUEST_NULL); }\
                index_%s = _i;",
                sendbufname, acc_name, acc_name, acc_name, acc_name, 
                sendbufname, acc_name, sendbufname, acc_name, acc_name, acc_name);
				sprintf(flow_cg_text+strlen(flow_cg_text), "\
					%s = %d;",
					send_count_name, max_copy_level + 1); // one more for the loop number
				sprintf(flow_cg_text+strlen(flow_cg_text), "\
					%s = pack_%s_%d(%s,%s[index_%s],%s",
					send_count_name, acc_name, loop_num, args, sendbufname, acc_name, send_count_name);
			}
			else {
				sprintf(flow_cg_text+strlen(flow_cg_text), "\
					%s = pack_%s_%d(%s,%s,%s",
					send_count_name, acc_name, loop_num, args, sendbufname, send_count_name);

			}
			if (options->variables_not_global && !options->data_dist) sprintf(flow_cg_text+strlen(flow_cg_text), ",%s", acc_name);

			print_data_dist_parm_call(flow_cg_text, arr);

			sprintf(flow_cg_text+strlen(flow_cg_text), "); break;} }");

			fprintf(packfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, sendbufname, send_count_name);
			fprintf(headerfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, sendbufname, send_count_name);
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
			}
			add_data_dist_parm_decl(headerfp, acc_name, prog);
			fprintf(headerfp,");\n");
			generate_pack_or_unpack(packfp, prog, flow_out, flow_copy_text, IN_FUNCTION, iters, src_copy_level, acc_nrows,
					acc_name, send_count_name);
		}

		/* Create statement stub to get hold of the outer loops for the copy
		 * stmt; rest of the loops to be added are the actual copy loops */
		Stmt *flow_copy_guard = create_helper_stmt(anchor_stmt, src_copy_level, flow_cg_text, COPY_OUT, loop_num);

		/* sigma, clear statements */
		char *clear_text = "for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; }";

		Stmt *send_recv_list_clear = create_helper_stmt(anchor_stmt, outer_dist_loop_level, clear_text, COMM_CALL, loop_num);
		Stmt *sigma_stmt = create_helper_stmt(anchor_stmt, src_copy_level, sigma_sender_text, SIGMA, loop_num);

		char *displsname = concat("displs_", acc_name);
		char *recvbufname = concat("recv_buf_", acc_name);
		char *recv_counts_name = concat("recv_counts_", acc_name);

		char *comm_text = malloc(2048);

		/* Message passing code (MPI calls) */
		if (options->dynschedule) {
			strcpy(comm_text, "");
			sprintf(comm_text+strlen(comm_text),  "\
				for (__p=0; __p<nprocs; __p++) {\
				if(receiver_list[__p] != 0) {");
			sprintf(comm_text+strlen(comm_text), "%s[index_%s][0] = %d;", sendbufname, acc_name, loop_num);
			for (i=0; i<src_copy_level; i++) {
				sprintf(comm_text+strlen(comm_text), "%s[index_%s][%d] = t%d;", sendbufname, acc_name, i+1, i+1);
			}
			sprintf(comm_text+strlen(comm_text),  "\
		        IF_DYNSCHEDULER_DEBUG_PRINT(\
		        fprintf(__debug_print_fp, \"sending node %%d buf %%d loop %d task ",
		        loop_num);
			for (i=0; i<src_copy_level; i++) {
				sprintf(comm_text+strlen(comm_text),  "%%d ");
			}
			sprintf(comm_text+strlen(comm_text),  "\
		        %s size %%d to node %%d\\n\", \
		        my_rank, index_%s[__p], %s, %s, __p)); \
                IF_DEBUG_FLUSH(fflush(__debug_print_fp));",
                acc_name, acc_name, args, send_count_name);
			if (options->timereport)
				sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s);", send_count_name);
			sprintf(comm_text+strlen(comm_text),  "\
				MPI_Isend(%s[index_%s], %s, MPI_DOUBLE,\
					__p, %s, MPI_COMM_WORLD, &send_reqs_%s[index_%s*nprocs+__p]);\
                IF_DYNSCHEDULER_DEBUG_PRINT(\
                if (__error) {\
                char *msg=(char*)malloc((MPI_MAX_ERROR_STRING+1)*sizeof(char)); int len; \
                MPI_Error_string(__error, msg, &len); \
                fprintf(__debug_print_fp, \"node %%d %s send buf %%d error %%d %%s\\n\", my_rank, index_%s, __error, msg); } )\
                }}",
                sendbufname, acc_name, send_count_name, data_tag, acc_name, acc_name, acc_name, acc_name);
		}
		else {
			if (options->timereport)
				sprintf(comm_text,  "IF_TIME(t_comm_start = rtclock());");
			else
				strcpy(comm_text, "");
			sprintf(comm_text+strlen(comm_text), "\
				for (__p=0; __p<nprocs; __p++) {\
				%s[__p] = receiver_list[__p]? %s: 0;\
				}\
				MPI_Alltoall(%s, 1, MPI_INT,\
					%s, 1, MPI_INT, MPI_COMM_WORLD);",
					send_counts_name, send_count_name, send_counts_name, recv_counts_name);
			sprintf(comm_text+strlen(comm_text),  "\
				req_count=0;\
				for (__p=0; __p<nprocs; __p++) {\
				if(%s[__p] >= 1) {",
				send_counts_name);
			if (options->timereport)
				sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s);", send_count_name);
			sprintf(comm_text+strlen(comm_text),  "\
				MPI_Isend(%s, %s, MPI_DOUBLE,\
					__p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}",
					sendbufname, send_count_name);
			sprintf(comm_text+strlen(comm_text),  "for (__p=0; __p<nprocs; __p++) {\
				if(%s[__p] >= 1) {\
				MPI_Irecv(%s+%s[__p], %s[__p], MPI_DOUBLE,\
					__p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}\
				MPI_Waitall(req_count, reqs, stats);\
				%s = 0; \
                for (__p=0; __p<nprocs; __p++) curr_%s[__p] = 0;", 
                recv_counts_name, recvbufname, displsname, recv_counts_name,
                send_count_name, displsname);
			if (options->timereport)
				sprintf(comm_text+strlen(comm_text),  "IF_TIME(t_comm += rtclock() - t_comm_start);");
		}

		/*sprintf(text, "MPI_Allgather(send_buf, send_count, MPI_DOUBLE,\
      recv_buf, max_recv_count, MPI_DOUBLE, MPI_COMM_WORLD)");*/

		/* Create statement stub to get hold of the outer loops for the copy
		 * stmt; rest of the loops to be added are the actual copy loops */
		Stmt *comm_stmt = create_helper_stmt(anchor_stmt, outer_dist_loop_level, comm_text, COMM_CALL, loop_num);
		//IF_DEBUG(pluto_stmt_print(stdout, comm_stmt););

		/* Receiver-side copy */
		/* With send/recv -based more exact communication scheme, recv side copy
		 * does not include all the sender-side copies; processes that do not send
		 * data are to be skipped */
		if (options->dynschedule) {
			sprintf(prog->decls+strlen(prog->decls),
					"for (__p=0; __p<__MAX_NUM_RECVS; __p++) { \
                %s[__p] = (double *) polyrt_max_alloc(%s[__p], sizeof(double)*%s, &recv_buf_size_%s[__p]); }\n",
                recvbufname, recvbufname, max_elements_name, acc_name);
		}
		else {
			sprintf(prog->decls+strlen(prog->decls),
					"%s = (double *) polyrt_max_alloc(%s, nprocs*send_buf_size_%s, &recv_buf_size_%s);\n",
					recvbufname, recvbufname, acc_name, acc_name);
			sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
				%s[__p] = __p*send_buf_size_%s/sizeof(double);}\n\n", displsname, acc_name);
		}

		char *flow_copyback_text = malloc(strlen(access) + strlen(" = ") + strlen(recvbufname)
				+ strlen("[") + strlen(displsname) + strlen("[proc]+ count++];") + 1);
		sprintf(flow_copyback_text, "%s = %s[%s + count++];",
				access, recvbufname, displsname);

		char proc_stmt_text[1024];
		strcpy(proc_stmt_text, "");

		if ((flow_out != NULL) && !pluto_constraints_is_empty(flow_out)) {
			if (options->dynschedule) {
				// one more for the loop number
				sprintf(proc_stmt_text+strlen(proc_stmt_text),
						"%s = %d; \
					%s = unpack_%s_%d(",
					displsname, max_copy_level + 1, displsname, acc_name, loop_num);
				for (i=0; i<src_copy_level; i++) {
					sprintf(proc_stmt_text+strlen(proc_stmt_text),
							"(int)%s[%d],",
							recvbufname, i+1);
				}
				sprintf(proc_stmt_text+strlen(proc_stmt_text),
						"%s,0,%s",
						recvbufname, displsname);
				if (options->variables_not_global && !options->data_dist) sprintf(proc_stmt_text+strlen(proc_stmt_text), ",%s", acc_name);
				if(options->data_dist)
					print_data_dist_parm_call(proc_stmt_text, arr);
				sprintf(proc_stmt_text+strlen(proc_stmt_text), ");");
			}
			else {
				sprintf(proc_stmt_text+strlen(proc_stmt_text), "proc = pi_%d(%s,%s, nprocs); ", loop_num, args, params);
				sprintf(proc_stmt_text+strlen(proc_stmt_text),
						"if ((proc != my_rank) && (%s[proc] > 0)) { \
					for (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; } \
					%s, proc, nprocs, receiver_list); ", recv_counts_name, sigma_text);
				sprintf(proc_stmt_text+strlen(proc_stmt_text),
						"for (__p=0; __p<nprocs; __p++) { if (receiver_list[__p] != 0) \
					{ curr_%s[proc] = unpack_%s_%d(%s,%s,%s[proc],curr_%s[proc]",
					displsname, acc_name, loop_num, args, recvbufname, displsname, displsname);
				if (options->variables_not_global && !options->data_dist) sprintf(proc_stmt_text+strlen(proc_stmt_text), ",%s", acc_name);
				print_data_dist_parm_call(proc_stmt_text, arr);
				sprintf(proc_stmt_text+strlen(proc_stmt_text), "); break; } } }");
			}

			fprintf(packfp, "int unpack_%s_%d(%s,double *%s,int %s,int count", acc_name, loop_num, decl_args, recvbufname, displsname);
			fprintf(headerfp, "int unpack_%s_%d(%s,double *%s,int %s,int count", acc_name, loop_num, decl_args, recvbufname, displsname);
			if (options->variables_not_global && !options->data_dist) {
				fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
				fprintf(headerfp, ", __DECLARATION_OF_%s", acc_name);
			}
			add_data_dist_parm_decl(headerfp, acc_name, prog);
			fprintf(headerfp,");\n");
			generate_pack_or_unpack(packfp, prog, flow_out, flow_copyback_text, IN_FUNCTION, iters, src_copy_level, acc_nrows
					, acc_name, "count");
		}

		/* Create statement stub to get hold of the outer loops for the copy
		 * stmt; rest of the loops to be added are the actual copy loops */
		Stmt *copyback_proc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
					proc_stmt_text, COPY_IN, loop_num);

		// comm stmts organized as follows
		if (options->dynschedule) {
			// pack stmt, comm stmt, unpack stmt
			*num_comm_stmts = 3;
		}
		else {
			// pack stmt, clear stmt, sigma stmt, comm stmt, unpack stmt
			*num_comm_stmts = 5;
		}
		Stmt **copy_comm_stmts = (Stmt **) malloc((*num_comm_stmts)*sizeof(Stmt *));
		copy_comm_stmts[0] = flow_copy_guard;
		if (options->dynschedule) {
			copy_comm_stmts[1] = comm_stmt;
			copy_comm_stmts[2] = copyback_proc_stmt;
		}
		else {
			copy_comm_stmts[1] = send_recv_list_clear;
			copy_comm_stmts[2] = sigma_stmt;
			copy_comm_stmts[3] = comm_stmt;
			copy_comm_stmts[4] = copyback_proc_stmt;
		}

		for (i=0; i<acc_nrows; i++) {
			free(iters[i]);
		}
		free(iters);

		pluto_constraints_free(flow_out);
		free(sigma_text);
		free(flow_copy_text);
		free(flow_cg_text);
		free(flow_copyback_text);
		free(comm_text);

		free(displsname);
		free(sendbufname);
		free(send_count_name);
		free(send_counts_name);
		free(recv_counts_name);
		free(send_buf_size);
		free(data_tag);
		free(access);
		fclose(packfp);

		return copy_comm_stmts;

}

// returns for each domain parallel loop or band of loops:
// outer_dist_loop_level : the level/dimension of the outermost distributed/parallel loop
// copy_level : the level/dimension of a statement within the innermost distributed/parallel loop
void init_copy_level(PlutoProg *prog, Ploop **loops, int nloops, int *copy_level, int *outer_dist_loop_level)
{
	int l, inner_dist_loop_level;
	for (l=0; l<nloops; l++) {
		if (options->multi_level_distribution) {
			Stmt *stmt = loops[l]->stmts[0];
			assert(stmt->last_tile_dim >= loops[l]->depth);
#ifdef distributing_all_inter_tile_loops
			for (inner_dist_loop_level=stmt->last_tile_dim; inner_dist_loop_level>=loops[l]->depth; inner_dist_loop_level--) {
				if (pluto_is_hyperplane_loop(stmt, inner_dist_loop_level)) break;
			}
#else
			for (inner_dist_loop_level=loops[l]->depth; inner_dist_loop_level<stmt->last_tile_dim; inner_dist_loop_level++) {
				if (pluto_is_hyperplane_loop(stmt, inner_dist_loop_level+1)) {
					inner_dist_loop_level++;
					break; // for now, distributing only 2 loops at the most
				}
			}
#endif
			if (options->distmem) {
				int i, num_inner_loops = 0;
				Ploop **inner_loops = pluto_get_loops_under(loops[l]->stmts, loops[l]->nstmts, loops[l]->depth,
						prog, &num_inner_loops);
				int is_distribution_set = 0;
				for (i=0; i<num_inner_loops; i++) {
					if (inner_loops[i]->depth <= inner_dist_loop_level) {
						if (!pluto_loop_is_parallel(prog, inner_loops[i])) {
							if (options->dynschedule) {
								sprintf(prog->decls+strlen(prog->decls), "__is_multi_partitioned[%d] = 1;\n", l);
								is_distribution_set = 1;
							}
							else {
								inner_dist_loop_level = loops[l]->depth;
							}
							break;
						}
					}
				}
				if (options->dynschedule && !is_distribution_set) {
					// !!!roshan can use some loop analysis to see if block cyclic is beneficial
					sprintf(prog->decls+strlen(prog->decls), "##ifdef __USE_BLOCK_CYCLIC\n");
					sprintf(prog->decls+strlen(prog->decls), "__is_block_cyclic[%d] = 1;\n", l);
					sprintf(prog->decls+strlen(prog->decls), "##endif\n");
				}
				pluto_loops_free(inner_loops,num_inner_loops);
			}
		}
		else {
			inner_dist_loop_level = loops[l]->depth;
		}
		copy_level[l] = inner_dist_loop_level + 1;
		if (outer_dist_loop_level != NULL) outer_dist_loop_level[l] = loops[l]->depth;
		if (options->blockcyclic) {
			/* FIXME: copy_level[l] in the presence of blockcyclic
			 * will depend on the order in which these parallel loops were tiled
			 * for block cylic partitioning (will work if they were in the
			 * increasing order of depths: i.e., loops in ploops array were in
			 * increasing order of their depths
			 */
			copy_level[l] += l;
			if (outer_dist_loop_level != NULL) outer_dist_loop_level[l] += l;
		}
	}
}

void pluto_dynschedule_common_codegen(PlutoProg *prog, FILE *sigmafp, FILE *outfp, FILE *headerfp)
{
	int nloops=0;
	Ploop **loops = pluto_get_dom_parallel_loops(prog, &nloops);
	/* Loops should be in the increasing order of their depths */
	qsort(loops, nloops, sizeof(Ploop *), pluto_loop_compar);

	int *copy_level = (int *) malloc(nloops*sizeof(int));
	int *outer_dist_copy_level = (int *) malloc(nloops*sizeof(int));

	init_copy_level(prog, loops, nloops, copy_level, outer_dist_copy_level);

	int l;
	for (l=0; l<nloops; l++) {
		Ploop *loop = loops[l];
		if(options->data_dist) {
			//due to my_rank as a parameter, data tile domains doesnt assume my_rank parameter
			if(options->distmem)
				prog->npar--;

			gen_data_tile_alloc_cloog_code(prog, l, loop->stmts, loop->nstmts, copy_level[l], sigmafp, headerfp);
			gen_data_tile_ref_count_update_cloog_code(prog, l, loop->stmts, loop->nstmts, copy_level[l], sigmafp, headerfp);

			if(options->distmem)
				prog->npar++;
		}
		gen_compute_task_cloog_code(prog, l, loop->stmts, loop->nstmts, copy_level[l], sigmafp, headerfp);
	}

	if(options->data_dist && !options->distmem){

		pluto_dist_generate_copy_back_code(loops, nloops,copy_level,
				headerfp, prog);
	}
	if (options->dynschedule_graph) {
		gen_dynschedule_graph_main_text_code(prog, loops, nloops, copy_level, outfp);
	}
	else {
		gen_dynschedule_main_text_code(prog, loops, nloops, copy_level, outfp);
	}

	free(copy_level);

	pluto_loops_free(loops,nloops);
}

int pluto_dynschedule_graph_codegen(PlutoProg *prog, FILE *sigmafp, FILE *outfp, FILE *headerfp)
{
	fprintf(outfp, "#include <assert.h>\n\n");
	fprintf(outfp, "#include <omp.h>\n\n");
	fprintf(outfp, "#include \"tbb/task_scheduler_init.h\"\n\n");

	if (options->timereport) {
		fprintf(outfp, "\tdouble t_local, t_create_dag;\n");
	}

	fprintf(outfp, "\tint _num_threads;\n##pragma omp parallel\n{_num_threads = omp_get_num_threads();}\n");
	fprintf(outfp, "tbb::task_scheduler_init init(_num_threads);\n");
	fprintf(outfp, "tbb::flow::graph __graph;\n");
	fprintf(outfp, "tbb::flow::broadcast_node<tbb::flow::continue_msg> __start(__graph);\n");
	fprintf(outfp, "%s", prog->decls);

	if (options->timereport) {
		fprintf(outfp, "\tIF_TIME(t_local = rtclock());\n");
	}

	pluto_dynschedule_common_codegen(prog, sigmafp, outfp, headerfp);

	if (options->timereport) {
		fprintf(outfp, "##ifdef TIME\n");
		fprintf(outfp, "\tt_local = rtclock() - t_local;\n");
		fprintf(outfp, "\tchar *buffer = (char *)malloc(4096 * sizeof(char));\n");
		fprintf(outfp, "\tstrcpy(buffer, \"\");\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"DAG creation time: %%0.6lf s\\n\", t_create_dag);\n\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Total time: %%0.6lf s\\n\", t_local);\n\n");
		fprintf(outfp, "\tfprintf(stdout, \"%%s\\n\", buffer);\n");
		fprintf(outfp, "\tfree(buffer);\n");
		fprintf(outfp, "##endif\n");
	}

	return 0;
}

int pluto_dynschedule_codegen(PlutoProg *prog, FILE *sigmafp, FILE *outfp, FILE *headerfp)
{
	fprintf(outfp, "#include \"polyrt.h\"\n");

	fprintf(outfp, "#include <assert.h>\n");
	fprintf(outfp, "#include <omp.h>\n\n");

	fprintf(outfp, "\tint _num_threads;\n##pragma omp parallel\n{_num_threads = omp_get_num_threads();}\n");
	fprintf(outfp, "\tint __tid;\n");
	fprintf(outfp, "\tlong int _num_tasks_to_execute;\n");

	if (options->timereport) {
		fprintf(outfp, "\tint waiting_pops_start[_num_threads], waiting_pops[_num_threads], \
                waiting_pops_threads;\n");
		fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) waiting_pops[__tid] = 0;\n");
		fprintf(outfp, "\tdouble t_local;\n");
		fprintf(outfp, "\tdouble t_tasks_create;\n");
		fprintf(outfp, "\tdouble t_tasks_manage_start[_num_threads], t_tasks_manage[_num_threads], \
                t_tasks_manage_threads;\n");
		fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_tasks_manage[__tid] = 0;\n");
		fprintf(outfp, "\tdouble t_comp_start[_num_threads], t_comp[_num_threads], t_comp_threads_min, t_comp_threads;\n");
		fprintf(outfp, "\tdouble t_comp_mean, t_comp_sd;\n");
		fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_comp[__tid] = 0;\n");
		fprintf(outfp, "\tdouble t_wait_comp_start[_num_threads], t_wait_comp[_num_threads], t_wait_comp_threads;\n");
		fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_wait_comp[__tid] = 0;\n");

		if(options->data_dist){
			fprintf(outfp, "\tdouble t_data_mang_start[_num_threads], t_data_mang[_num_threads], t_data_mang_threads;\n");
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_data_mang[__tid] = 0;\n");
		}

		if(options->data_dist){
			fprintf(outfp, "\tdouble t_writeout_start, t_writeout = 0.0;\n");
			fprintf(outfp, "\tdouble t_data_init_start, t_data_init = 0.0;\n");

		}
	}

	fprintf(outfp, "%s", prog->decls);

	fprintf(outfp, "\tpolyrt_init_grid_size();\n");

    int i;
	if(options->compute_pi){
		fprintf(outfp, "compute_pi(");
		fprintf(outfp,"%d", options->num_inital_partitions);
		for(i=0;i<prog->npar;i++){
			//			if(i!=0)
			fprintf(outfp,",");
			fprintf(outfp," %s ", prog->params[i]);
		}
		fprintf(outfp,"); \n return ;\n");
	}

	if(options->data_dist){
		fprintf(outfp, "\tIF_TIME(t_data_init_start = rtclock());\n");
		fprintf(outfp, "data_tile_init_0(0, 0");

		int num = 0;
		PlutoAccess **access = pluto_get_all_accs(prog, &num);

		if(options->variables_not_global){
			fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
			for(i=0;i<num;i++){
				fprintf(outfp,", %s ", access[i]->name);
			}
			fprintf(outfp," __endif ");
		}
		for(i=0;i<num;i++){
			fprint_data_dist_parm_call(outfp, access[i]->name, prog);
		}
		fprintf(outfp, ");\n");
		fprintf(outfp, "\tIF_TIME(t_data_init += rtclock() - t_data_init_start)\n");
		free(access);
	}

    if (options->timereport) {
        fprintf(outfp, "\tIF_TIME(t_local = rtclock());\n");
    }

	pluto_dynschedule_common_codegen(prog, sigmafp, outfp, headerfp);

    if (options->timereport) {
		fprintf(outfp, "\tIF_TIME(t_local = rtclock() - t_local);\n");
    }

	if (options->data_dist) {
		if (options->timereport) {
			fprintf(outfp, "\tIF_TIME(t_writeout = rtclock());\n");
		}

		int num_data;
		PlutoAccess **waccs;
		waccs = pluto_get_all_waccs(prog, &num_data);
		fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
		fprintf(outfp, "\tif (fopen(\".test\", \"r\")) {\n ");
		if(options->data_dist){
			fprintf(outfp, "copy_back_writeout_0(0, 0");

			if(options->variables_not_global){
				fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
				for(i=0;i<num_data;i++){
					fprintf(outfp,", %s ", waccs[i]->name);
				}
				fprintf(outfp," __endif ");
			}

			for (i=0; i<num_data; i++) {
				fprint_data_dist_parm_call(outfp, waccs[i]->name,prog);
			}
			fprintf(outfp, ");\n");
		}

		fprintf(outfp, "\n}");
		fprintf(outfp," __endif ");
		free(waccs);
		if (options->timereport) {
			fprintf(outfp, "\tIF_TIME(t_writeout = rtclock() - t_writeout);\n");
		}
	}

	if(options->data_dist){
		PlutoAccess **waccs;
		int num = 0;
		if(options->data_dist)
			waccs = pluto_get_all_accs(prog, &num);
		else
			waccs = pluto_get_all_waccs(prog, &num);

		for (i = 0; i < num ; ++i) {
			fprintf(outfp, "free_buffer_mang(buff_mang_%s);\n",waccs[i]->name);
		}

		free(waccs);
	}

	if (options->timereport) {
		fprintf(outfp, "##ifdef TIME\n");
		fprintf(outfp, "\twaiting_pops_threads = waiting_pops[0];\n\
                for (__tid=1; __tid<_num_threads; __tid++) if (waiting_pops_threads < waiting_pops[__tid]) \
                waiting_pops_threads = waiting_pops[__tid];\n");
		fprintf(outfp, "\tt_tasks_manage_threads = t_tasks_manage[0];\n\
                for (__tid=1; __tid<_num_threads; __tid++) if (t_tasks_manage_threads < t_tasks_manage[__tid]) \
                t_tasks_manage_threads = t_tasks_manage[__tid];\n");
		fprintf(outfp, "\tt_comp_threads_min = t_comp[0];\n\
                for (__tid=1; __tid<_num_threads; __tid++) if (t_comp_threads_min > t_comp[__tid]) \
                t_comp_threads_min = t_comp[__tid];\n");
		fprintf(outfp, "\tt_comp_threads = t_comp[0];\n\
                for (__tid=1; __tid<_num_threads; __tid++) if (t_comp_threads < t_comp[__tid]) \
                t_comp_threads = t_comp[__tid];\n");
		fprintf(outfp, "\tt_comp_mean = 0.0;\n\
                for (__tid=0; __tid<_num_threads; __tid++) t_comp_mean += t_comp[__tid];\n\
                t_comp_mean = t_comp_mean/_num_threads;\n");
		fprintf(outfp, "\tt_comp_sd = 0.0;\n\
                for (__tid=0; __tid<_num_threads; __tid++) t_comp_sd += (t_comp[__tid]-t_comp_mean)*(t_comp[__tid]-t_comp_mean);\n\
                t_comp_sd = sqrt(t_comp_sd/_num_threads);\n");
		fprintf(outfp, "\tt_wait_comp_threads = t_wait_comp[0];\n\
                for (__tid=1; __tid<_num_threads; __tid++) if (t_wait_comp_threads < t_wait_comp[__tid]) \
                t_wait_comp_threads = t_wait_comp[__tid];\n");
        if (options->data_dist) {
            fprintf(outfp, "\tt_data_mang_threads = t_data_mang[0];\n\
                    for (__tid=1; __tid<_num_threads; __tid++) if (t_data_mang_threads < t_data_mang[__tid]) \
                    t_data_mang_threads = t_data_mang[__tid];\n");
        }
		fprintf(outfp, "\tchar *buffer = (char *)malloc(4096 * sizeof(char));\n");
		fprintf(outfp, "\tstrcpy(buffer, \"\");\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Minimum computation time across all threads: %%0.6lf s\\n\", t_comp_threads_min);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Maximum computation time across all threads: %%0.6lf s\\n\", t_comp_threads);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Mean computation time across all threads: %%0.6lf s\\n\", t_comp_mean);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Standard deviation in computation time across all threads: %%0.6lf s\\n\", t_comp_sd);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Load imbalance factor (standard-deviation/mean): %%0.2lf %%%%\\n\", 100*t_comp_sd/t_comp_mean);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Maximum 'waiting for computation' time across all threads: %%0.6lf s\\n\", t_wait_comp_threads);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Maximum number of waiting/unsuccessful POPs across all threads: %%d\\n\", waiting_pops_threads);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Maximum tasks management time across all threads: %%0.6lf s\\n\", t_tasks_manage_threads);\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Tasks creation time: %%0.6lf s\\n\", t_tasks_create);\n");
        if (options->data_dist) {
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Data init time: %%0.6lf s\\n\", t_data_init);\n");
            fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Data mang time: %%0.6lf s\\n\", t_data_mang_threads);\n");
            fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Copy back time: %%0.6lf s\\n\", t_writeout);\n");
            fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Sum of split times: %%0.6lf s\\n\", t_comp_threads+t_wait_comp_threads+t_tasks_manage_threads+t_tasks_create+t_data_init+t_data_mang_threads+t_writeout);\n");
        } else {
            fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Sum of split times: %%0.6lf s\\n\", t_comp_threads+t_wait_comp_threads+t_tasks_manage_threads+t_tasks_create);\n");
        }
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Total time: %%0.6lf s\\n\", t_local);\n");
		fprintf(outfp, "\tfprintf(stdout, \"%%s\\n\", buffer);\n");
		fprintf(outfp, "\tfree(buffer);\n");
		fprintf(outfp, "##endif\n");
	}

	return 0;
}

int pluto_distmem_codegen(PlutoProg *prog, FILE *cloogfp, FILE *sigmafp, FILE *outfp, FILE *headerfp)
{
	int i;

	fprintf(outfp, "#include <string.h>\n");
	fprintf(outfp, "#include \"polyrt.h\"\n");

	if(options->data_dist)
		fprintf(outfp, "#include \"buffer_manager.h\"\n");

	fprintf(outfp, "#include \"pi_defs.h\"\n");
	fprintf(outfp, "#include <assert.h>\n");

	if (options->commopt_foifi) {
		fprintf(outfp, "##ifndef __FOIFI_MAX_DUP_FACTOR\n");
		fprintf(outfp, "##define __FOIFI_MAX_DUP_FACTOR 1\n");
		fprintf(outfp, "##endif\n\n");
	}
	fprintf(outfp, "#include <omp.h>\n\n");
	fprintf(outfp, "#include <mpi.h>\n\n");
	fprintf(outfp, "#define MPI \n\n");
	fprintf(outfp, "#include <limits.h>\n\n");

	pluto_mark_vec_stmts(cloogfp, prog);
	generate_declarations(prog, outfp);

	fprintf(outfp, "int _num_threads;\n##pragma omp parallel\n{_num_threads = omp_get_num_threads();}\n");

	if(options->data_dist){
		fprintf(outfp, "\tdouble t_data_init_start, t_data_init= 0.0, t_globaldata_init= 0.0;\n");
		if(!options->dynschedule)
			fprintf(outfp, "\tdouble t_data_mang_start, t_data_mang= 0.0, t_globaldata_mang= 0.0;\n");
	}

	fprintf(outfp, "\n##ifndef GLOBAL_MY_RANK\n\tint my_rank;\n##endif\n");
	fprintf(outfp, "\n##ifndef __MAX_NUM_RECVS\n\t##define __MAX_NUM_RECVS (nprocs*_num_threads*16)\n##endif\n");
	fprintf(outfp, "\n##ifndef __MAX_NUM_SENDS\n\t##define __MAX_NUM_SENDS (nprocs*_num_threads*8)\n##endif\n");
	fprintf(outfp, "\tint nprocs, my_start, my_end, _i,\
            __p, proc, recv_proc, distinct_recv,\
    		__tid;\n");
	fprintf(outfp, "\tlong int _num_tasks_to_execute, _num_tasks_to_unpack;\n");
	fprintf(outfp, "\tint count;\n");
	fprintf(outfp, "int req_count;\n");
	fprintf(outfp, "\tint _lb_dist, _ub_dist");
	for (i=0; i<prog->num_hyperplanes; i++) {
		fprintf(outfp, ", lbd_t%d, ubd_t%d", i+1, i+1);
	}
	fprintf(outfp, ";\n");
	if (options->dynschedule) {
		fprintf(outfp, "\tint thread_support_provided;\n");
		fprintf(outfp, "\tMPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &thread_support_provided);\n");
		fprintf(outfp, "\tassert(thread_support_provided == MPI_THREAD_MULTIPLE);\n");
	}
	else {
		fprintf(outfp, "\tMPI_Init(NULL, NULL);\n");
	}
	fprintf(outfp, "\tMPI_Comm_rank(MPI_COMM_WORLD, &my_rank);\n");
	fprintf(outfp, "\tMPI_Comm_size(MPI_COMM_WORLD, &nprocs);\n");

	if (!options->dynschedule) {
        fprintf(outfp, "\tint receiver_list[nprocs]; \n");
        fprintf(outfp, "\tfor (__p=0; __p<nprocs; __p++) { receiver_list[__p] = 0; } \n");
        fprintf(outfp, "MPI_Request %s[2*nprocs];\n", "reqs");
        fprintf(outfp, "MPI_Status %s[2*nprocs];\n\n", "stats");
	}

	if (options->timereport) {
		if (options->dynschedule) {
			fprintf(outfp, "\tint __tid_receiver;\n");
			fprintf(outfp, "\tint waiting_pops[_num_threads], \
                    waiting_pops_threads, waiting_pops_global;\n");
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) waiting_pops[__tid] = 0;\n");
			fprintf(outfp, "\tdouble t_tasks_create, t_globaltasks_create;\n");
			fprintf(outfp, "\tdouble t_tasks_manage_start[_num_threads], t_tasks_manage[_num_threads], \
                    t_tasks_manage_threads, t_globaltasks_manage;\n");
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_tasks_manage[__tid] = 0;\n");
			fprintf(outfp, "\tdouble t_comp_start[_num_threads], t_comp[_num_threads], \
                    t_comp_threads_min, t_globalcomp_min, t_comp_threads, t_globalcomp;\n");
			fprintf(outfp, "\tdouble t_comp_mean, t_comp_sd, t_globalcomp_mean, t_globalcomp_sd;\n");
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_comp[__tid] = 0;\n");
			fprintf(outfp, "\tdouble t_wait_comp_start[_num_threads], t_wait_comp[_num_threads], t_wait_comp_threads, t_globalwait_comp;\n");
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_wait_comp[__tid] = 0;\n");
			if(options->data_dist){
				fprintf(outfp, "\tdouble t_data_mang_start[_num_threads], t_data_mang[_num_threads], t_data_mang_threads, t_globaldata_mang;\n");
				fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_data_mang[__tid] = 0;\n");
			}
			fprintf(outfp, "\tdouble t_pack_start[_num_threads], t_pack[_num_threads], t_pack_threads, t_globalpack;\n");
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) t_pack[__tid] = 0;\n");
		}
		else {
			fprintf(outfp, "\tdouble t_comm_start, t_comm = 0.0, t_globalcomm = 0.0;\n");
			fprintf(outfp, "\tdouble t_comp_mean, t_comp_sd;\n");
			fprintf(outfp, "\tdouble t_comp_start, t_comp = 0.0, t_globalcomp = 0.0;\n");
			fprintf(outfp, "\tdouble t_pack_start, t_pack = 0.0, t_globalpack = 0.0;\n");
		}
		fprintf(outfp, "\tdouble t_writeout = 0.0;\n");
		fprintf(outfp, "\tdouble t_unpack_start, t_unpack = 0.0, t_globalunpack = 0.0;\n");
		if (options->dynschedule) {
			fprintf(outfp, "\tdouble t_wait_unpack_start, t_wait_unpack = 0.0, t_globalwait_unpack = 0.0;\n");
		}
		fprintf(outfp, "\tdouble t_local = 0.0, t_global = 0.0;\n");
		if (options->dynschedule) {
			fprintf(outfp, "\tdouble __total_count[_num_threads];\n");
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) __total_count[__tid] = 0;\n");
			fprintf(outfp, "\tdouble __total_count_threads = 0;\n");
		}
		else {
			fprintf(outfp, "\tdouble __total_count = 0;\n");
		}
		fprintf(outfp, "\tdouble __total_count_all = 0;\n");
	}

	fprintf(outfp, "%s", prog->decls);

	fprintf(outfp, "\tpolyrt_init_grid_size();\n");

	if(options->compute_pi){
		fprintf(outfp, "compute_pi(");

		fprintf(outfp,"%d", options->num_inital_partitions);
		for(i=0;i<prog->npar;i++){
			//			if(i!=0)
			fprintf(outfp,",");
			fprintf(outfp," %s ", prog->params[i]);
		}
		fprintf(outfp,");\n");
	}

	if(options->data_dist){
		fprintf(outfp, "\tIF_TIME(t_data_init_start = rtclock());\n");
		fprintf(outfp, "data_tile_init_0(my_rank, nprocs");

		int num = 0;
		PlutoAccess **access = pluto_get_all_accs(prog, &num);

		if(options->variables_not_global){
			fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
			for(i=0;i<num;i++){
				fprintf(outfp,", %s ", access[i]->name);
			}
			fprintf(outfp," __endif ");
		}
		for(i=0;i<num;i++){
			fprint_data_dist_parm_call(outfp, access[i]->name, prog);
		}
		fprintf(outfp, ");\n");
		fprintf(outfp, "\tIF_TIME(t_data_init += rtclock() - t_data_init_start)\n");
		free(access);
	}

	if (options->timereport) {
		fprintf(outfp, "\tIF_TIME(t_local = rtclock());\n");
	}

	if (options->dynschedule) {
		pluto_dynschedule_common_codegen(prog, sigmafp, outfp, headerfp);
	}
	else {
		pluto_gen_cloog_code(prog, -1, -1, cloogfp, outfp);
	}

    if (options->timereport) {
		fprintf(outfp, "\tIF_TIME(t_local = rtclock() - t_local);\n");
    }

	if (options->verify_output) {
		if (options->timereport) {
			fprintf(outfp, "\tIF_TIME(t_writeout = rtclock());\n");
		}

		int num_data;
		PlutoAccess **waccs;
		waccs = pluto_get_all_distinct_write_vars(prog, &num_data);
		fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
		fprintf(outfp, "\tif (fopen(\".test\", \"r\")) {\n ");
		if(options->data_dist){
			fprintf(outfp, "copy_back_writeout_0(my_rank, nprocs");

			if(options->variables_not_global){
				fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
				for(i=0;i<num_data;i++){
					fprintf(outfp,", %s ", waccs[i]->name);
				}
				fprintf(outfp," __endif ");
			}

			for (i=0; i<num_data; i++) {
				fprint_data_dist_parm_call(outfp, waccs[i]->name,prog);
			}
			fprintf(outfp, ");\n");
		}
		fprintf(outfp, "write_out(");
		fprintf(outfp, "my_rank,nprocs");
		for (i=0; i<num_data; i++) {
			fprintf(outfp, ",lw_buf_%s", waccs[i]->name);
			fprintf(outfp, ",lw_recv_buf_%s", waccs[i]->name);
			fprintf(outfp, ",displs_lw_%s", waccs[i]->name);
		}
		if (options->variables_not_global) {
			for (i=0; i<num_data; i++) {
				fprintf(outfp, ",%s", waccs[i]->name);
			}
		}
		if(options->data_dist){
			for (i=0; i<num_data; i++)
				fprint_data_dist_parm_call(outfp, waccs[i]->name, prog);

		}

		fprintf(outfp, ");\n}");
		fprintf(outfp," __endif ");
		free(waccs);
		if (options->timereport) {
			fprintf(outfp, "\tIF_TIME(t_writeout = rtclock() - t_writeout);\n");
		}
	}

	if(options->data_dist){
		PlutoAccess **waccs;
		int num = 0;
		if(options->data_dist)
			waccs = pluto_get_all_accs(prog, &num);
		else
			waccs = pluto_get_all_waccs(prog, &num);

		for (i = 0; i < num ; ++i) {
			fprintf(outfp, "free_buffer_mang(buff_mang_%s);\n",waccs[i]->name);
		}

		free(waccs);
	}

	if (options->timereport) {
		fprintf(outfp, "##ifdef TIME\n");

		if (options->dynschedule) {
			fprintf(outfp, "\tfor (__tid=0; __tid<_num_threads; __tid++) __total_count_threads += __total_count[__tid];\n");
			fprintf(outfp, "\twaiting_pops_threads = waiting_pops[0];\n\
                    for (__tid=1; __tid<_num_threads; __tid++) if (waiting_pops_threads < waiting_pops[__tid]) \
                    waiting_pops_threads = waiting_pops[__tid];\n");
			fprintf(outfp, "\tt_tasks_manage_threads = t_tasks_manage[0];\n\
                    for (__tid=1; __tid<_num_threads; __tid++) if (t_tasks_manage_threads < t_tasks_manage[__tid]) \
                    t_tasks_manage_threads = t_tasks_manage[__tid];\n");
			fprintf(outfp, "\tif (__tid_receiver != 0) t_comp_threads_min = t_comp[0];\n\
                    else t_comp_threads_min = t_comp[1];\n\
                    for (__tid=1; __tid<_num_threads; __tid++) if ((__tid != __tid_receiver) && (t_comp_threads_min > t_comp[__tid])) \
                    t_comp_threads_min = t_comp[__tid];\n");
			fprintf(outfp, "\tt_comp_threads = t_comp[0];\n\
                    for (__tid=1; __tid<_num_threads; __tid++) if (t_comp_threads < t_comp[__tid]) \
                    t_comp_threads = t_comp[__tid];\n");
			fprintf(outfp, "\tt_comp_mean = 0.0;\n\
                    for (__tid=0; __tid<_num_threads; __tid++) if (__tid != __tid_receiver) t_comp_mean += t_comp[__tid];\n\
                    t_globalcomp_mean = t_comp_mean;\n\
                    t_comp_mean = t_comp_mean/(_num_threads-1);\n");
			fprintf(outfp, "\tt_comp_sd = 0.0;\n\
                    for (__tid=0; __tid<_num_threads; __tid++) if (__tid != __tid_receiver) t_comp_sd += (t_comp[__tid]-t_comp_mean)*(t_comp[__tid]-t_comp_mean);\n\
                    t_comp_sd = sqrt(t_comp_sd/(_num_threads-1));\n");
			fprintf(outfp, "\tt_wait_comp_threads = t_wait_comp[0];\n\
                    for (__tid=1; __tid<_num_threads; __tid++) if (t_wait_comp_threads < t_wait_comp[__tid]) \
                    t_wait_comp_threads = t_wait_comp[__tid];\n");
			if(options->data_dist){
				fprintf(outfp, "\tt_data_mang_threads = t_data_mang[0];\n\
						for (__tid=1; __tid<_num_threads; __tid++) if (t_data_mang_threads < t_data_mang[__tid]) \
						t_data_mang_threads = t_data_mang[__tid];\n");
			}
			fprintf(outfp, "\tt_pack_threads = t_pack[0];\n\
                    for (__tid=1; __tid<_num_threads; __tid++) if (t_pack_threads < t_pack[__tid]) \
                    t_pack_threads = t_pack[__tid];\n");
		}
		fprintf(outfp, "\tchar *buffer = (char *)malloc(4096 * sizeof(char));\n");
		fprintf(outfp, "\tstrcpy(buffer, \"\");\n");
        fprintf(outfp, "\tif (fopen(\".detailed_time_report\", \"r\")) {\n");
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Node %%d\\n\", my_rank);\n");
		if(options->verify_output)
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Write-out time: %%0.6lf s\\n\", t_writeout);\n\n");
		if (options->dynschedule) {
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Values communicated: %%0.0lf\\n\", __total_count_threads);\n");
		}
		else {
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Values communicated: %%0.0lf\\n\", __total_count);\n");
		}
		if(options->data_dist){
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Data init time: %%0.6lf s\\n\", t_data_init);\n\n");
			if(options->dynschedule) {
				fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Data mang time: %%0.6lf s\\n\", t_data_mang_threads);\n\n");
            } else {
				fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Data mang time: %%0.6lf s\\n\", t_data_mang);\n\n");
            }
		}

		if (options->dynschedule) {
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Minimum computation time across all threads: %%0.6lf s\\n\", t_comp_threads_min);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Maximum computation time across all threads: %%0.6lf s\\n\", t_comp_threads);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Mean computation time across all threads: %%0.6lf s\\n\", t_comp_mean);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Standard deviation in computation time across all threads: %%0.6lf s\\n\", t_comp_sd);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Load imbalance factor (standard-deviation/mean): %%0.2lf %%%%\\n\", 100*t_comp_sd/t_comp_mean);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"'Waiting for computation' time: %%0.6lf s\\n\", t_wait_comp_threads);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Tasks creation time: %%0.6lf s\\n\", t_tasks_create);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Tasks management time: %%0.6lf s\\n\", t_tasks_manage_threads);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Packing time: %%0.6lf s\\n\", t_pack_threads);\n");
		}
		else {
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Computation time: %%0.6lf s\\n\", t_comp);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Communication time: %%0.6lf s\\n\", t_comm);\n");
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Packing time: %%0.6lf s\\n\", t_pack);\n");
		}
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Unpacking time: %%0.6lf s\\n\", t_unpack);\n");
		if (options->dynschedule) {
			fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"'Waiting for unpacking' time: %%0.6lf s\\n\", t_wait_unpack);\n");
			// unpacking time isn't considered since ideally, it should be overlapped with the rest
            if(options->data_dist){
                fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Sum of split times: %%0.6lf s\\n\", "
                        "t_comp_threads + t_wait_comp_threads + t_tasks_create + t_tasks_manage_threads + t_pack_threads + t_data_mang_threads);\n\n");
            }
            else {
                fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Sum of split times: %%0.6lf s\\n\", "
                        "t_comp_threads + t_wait_comp_threads + t_tasks_create + t_tasks_manage_threads + t_pack_threads);\n\n");
            }
		}
		else {
			if(options->data_dist){
				fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Sum of split times: %%0.6lf s\\n\", t_comp + t_comm + t_pack + t_unpack + t_data_mang);\n\n");
			}
			else
				fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Sum of split times: %%0.6lf s\\n\", t_comp + t_comm + t_pack + t_unpack);\n\n");
		}
        if (options->data_dist) {
            fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Total time - write-out - data_init time: %%0.6lf s\\n\", t_local);\n\n");
        } else {
            fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Total time minus write-out time: %%0.6lf s\\n\", t_local);\n");
        }
		fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"-------\");\n");
		fprintf(outfp, "\tfprintf(stdout, \"%%s\\n\", buffer);\n");
        fprintf(outfp, "\t}\n");
		if (options->dynschedule) {
			fprintf(outfp, "\tMPI_Reduce(&__total_count_threads, &__total_count_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);\n");
		}
		else {
			fprintf(outfp, "\tMPI_Reduce(&__total_count, &__total_count_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);\n");
		}
		if (options->dynschedule) {
			fprintf(outfp, "\tMPI_Allreduce(MPI_IN_PLACE, &t_globalcomp_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tt_globalcomp_mean = t_globalcomp_mean/((_num_threads-1)*nprocs);\n");
			fprintf(outfp, "\tt_globalcomp_sd = 0.0;\n\
                    for (__tid=0; __tid<_num_threads; __tid++) if (__tid != __tid_receiver) t_globalcomp_sd += (t_comp[__tid]-t_globalcomp_mean)*(t_comp[__tid]-t_globalcomp_mean);\n\
                    MPI_Allreduce(MPI_IN_PLACE, &t_globalcomp_sd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);\n\
                    t_globalcomp_sd = sqrt(t_globalcomp_sd/((_num_threads-1)*nprocs));\n");
            if (options->data_dist) {
                fprintf(outfp, "\tMPI_Reduce(&t_data_init, &t_globaldata_init, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
                fprintf(outfp, "\tMPI_Reduce(&t_data_mang_threads, &t_globaldata_mang, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
            }
			fprintf(outfp, "\tMPI_Reduce(&waiting_pops_threads, &waiting_pops_global, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_comp_threads_min, &t_globalcomp_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_comp_threads, &t_globalcomp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_wait_comp_threads, &t_globalwait_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_tasks_create, &t_globaltasks_create, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_tasks_manage_threads, &t_globaltasks_manage, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_pack_threads, &t_globalpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
		}
		else {
			fprintf(outfp, "\tMPI_Allreduce(&t_comp, &t_comp_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tt_comp_mean = t_comp_mean/nprocs;\n");
			fprintf(outfp, "\tt_comp_sd = (t_comp-t_comp_mean)*(t_comp-t_comp_mean);\n");
			fprintf(outfp, "\tMPI_Allreduce(MPI_IN_PLACE, &t_comp_sd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tt_comp_sd = sqrt(t_comp_sd/nprocs);\n");
            if (options->data_dist) {
                fprintf(outfp, "\tMPI_Reduce(&t_data_init, &t_globaldata_init, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
                fprintf(outfp, "\tMPI_Reduce(&t_data_mang, &t_globaldata_mang, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
            }
			fprintf(outfp, "\tMPI_Reduce(&t_comp, &t_globalcomp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_comm, &t_globalcomm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
			fprintf(outfp, "\tMPI_Reduce(&t_pack, &t_globalpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
		}
		fprintf(outfp, "\tMPI_Reduce(&t_unpack, &t_globalunpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
		if (options->dynschedule) {
			fprintf(outfp, "\tMPI_Reduce(&t_wait_unpack, &t_globalwait_unpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
		}
		fprintf(outfp, "\tMPI_Reduce(&t_local, &t_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
		fprintf(outfp, "\tif (my_rank==0) {\n");
		fprintf(outfp, "\t\tstrcpy(buffer, \"SUMMARY\\n\");\n");
		if(options->verify_output)
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Write-out time spent in master node: %%0.6lf s\\n\", t_writeout);\n");
        if(options->data_dist){
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum data init time spent across all nodes: %%0.6lf s\\n\", t_globaldata_init);\n\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum data mang time spent across all nodes: %%0.6lf s\\n\", t_globaldata_mang);\n\n");
        }
		fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Total communication volume across all nodes: %%0.6lf GB\\n\", __total_count_all*8/((double)1024*1024*1024));\n");
		fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum computation time spent across all nodes: %%0.6lf s\\n\", t_globalcomp);\n");
		if (options->dynschedule) {
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Minimum computation time spent across all nodes: %%0.6lf s\\n\", t_globalcomp_min);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Mean computation time across all nodes: %%0.6lf s\\n\", t_globalcomp_mean);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Standard deviation in computation time across all nodes: %%0.6lf s\\n\", t_globalcomp_sd);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Load imbalance factor (standard-deviation/mean): %%0.2lf %%%%\\n\", 100*t_globalcomp_sd/t_globalcomp_mean);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum 'waiting for computation' time spent across all nodes: %%0.6lf s\\n\", t_globalwait_comp);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum number of waiting/unsuccessful POPs across all nodes: %%d\\n\", waiting_pops_global);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum tasks creation time spent across all nodes: %%0.6lf s\\n\", t_globaltasks_create);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum tasks management time spent across all nodes: %%0.6lf s\\n\", t_globaltasks_manage);\n");
			;}
		else {
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Mean computation time across all nodes (not per-thread): %%0.6lf s\\n\", t_comp_mean);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Standard deviation in computation time across all nodes (not per-thread): %%0.6lf s\\n\", t_comp_sd);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Load imbalance factor (standard-deviation/mean): %%0.2lf %%%%\\n\", 100*t_comp_sd/t_comp_mean);\n");
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum communication time spent across all nodes: %%0.6lf s\\n\", t_globalcomm);\n");
		}
		fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum packing time spent across all nodes: %%0.6lf s\\n\", t_globalpack);\n");
		fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum unpacking time spent across all nodes: %%0.6lf s\\n\", t_globalunpack);\n");
		if (options->dynschedule) {
			fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum 'waiting for unpacking' time spent across all nodes: %%0.6lf s\\n\", t_globalwait_unpack);\n");
		}
		fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum total time spent across all nodes: %%0.6lf s\\n\", t_global);\n");
		fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"-------\");\n");
		fprintf(outfp, "\t\tfprintf(stdout, \"%%s\\n\", buffer);\n");
		fprintf(outfp, "\t}\n");
		fprintf(outfp, "\tfree(buffer);\n");
		fprintf(outfp, "##endif\n");
	}
	fprintf(outfp, "\tMPI_Barrier(MPI_COMM_WORLD);\n");
	fprintf(outfp, "\tMPI_Finalize();\n");

	return 0;
}

void pluto_add_distmem_decls(PlutoProg *prog, FILE *headerfp)
{
	int i, num=0;

	PlutoAccess **waccs;
	if(options->data_dist)
		waccs = pluto_get_all_accs(prog, &num);
	else
		waccs = pluto_get_all_distinct_write_vars(prog, &num);

	prog->orig_accesses = waccs;
	prog->num_accesses = num;

	for (i=0; i<num; i++) {

		assert(waccs != NULL);
		assert(waccs[i] != NULL);
		char *name = waccs[i]->name;

		if (options->dynschedule) {
			sprintf(prog->decls+strlen(prog->decls),
					"std::vector<MPI_Request> send_reqs_%s[_num_threads];\n", name);
			sprintf(prog->decls+strlen(prog->decls),
					"std::vector<MPI_Request> recv_reqs_%s;\n", name);
			sprintf(prog->decls+strlen(prog->decls),
					"std::vector<int> recv_completed_%s;\n", name);
			if (options->commopt_foifi || options->commopt_fop) {
				sprintf(prog->decls+strlen(prog->decls), "\
                        for (__tid=0; __tid<_num_threads; __tid++) {\n\
                        for (__p=0; __p<__MAX_NUM_SENDS; __p++) {\n\
                        send_reqs_%s[__tid].push_back(MPI_REQUEST_NULL);\n}\n}\n", 
                        name);
			}
			else {
				sprintf(prog->decls+strlen(prog->decls), "\
                        for (__tid=0; __tid<_num_threads; __tid++) {\n\
                        for (__p=0; __p<__MAX_NUM_SENDS*nprocs; __p++) {\n\
                        send_reqs_%s[__tid].push_back(MPI_REQUEST_NULL);\n}\n}\n", 
                        name);
			}
			sprintf(prog->decls+strlen(prog->decls), "\
                    for (__p=0; __p<__MAX_NUM_RECVS; __p++) {\n\
                    recv_reqs_%s.push_back(MPI_REQUEST_NULL);\n\
                    recv_completed_%s.push_back(0);\n\
                    }\n", name, name);
		}

		/* lw count name */
		sprintf(prog->decls+strlen(prog->decls),
				"int lw_count_%s = 0;\n", name );
		if (options->dynschedule) {
			sprintf(prog->decls+strlen(prog->decls),
					"int send_counts_%s[__tid][nprocs];\n", name);
			sprintf(prog->decls+strlen(prog->decls),
					"max_num_elements_%s = 0;\n", name);
			fprintf(headerfp, "extern int max_num_elements_%s;\n", name);
			fprintf(headerfp, "int max_num_elements_%s;\n", name);
		}
		else {
			sprintf(prog->decls+strlen(prog->decls),
					"int send_counts_%s[nprocs];\n", name);
		}
		sprintf(prog->decls+strlen(prog->decls),
				"int recv_counts_%s[nprocs];\n", name);
		sprintf(prog->decls+strlen(prog->decls),
				"int lw_recv_counts_%s[nprocs];\n", name);
		sprintf(prog->decls+strlen(prog->decls),
				"int displs_%s[nprocs];\n", name);
		sprintf(prog->decls+strlen(prog->decls),
				"int displs_lw_%s[nprocs];\n", name);
		if (options->commopt_foifi || options->commopt_fop) {
			if (options->dynschedule) {
				sprintf(prog->decls+strlen(prog->decls),
						"std::vector<double *> send_buf_%s[_num_threads];\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"std::vector<size_t> send_buf_size_%s[_num_threads];\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"for (__tid=0; __tid<_num_threads; __tid++) {\n\
						for (__p=0; __p<__MAX_NUM_SENDS; __p++) {\n\
						send_buf_%s[__tid].push_back(NULL);\n\
						send_buf_size_%s[__tid].push_back(0);\n}\n}", name, name);
			}
			else {
				sprintf(prog->decls+strlen(prog->decls),
						"double *send_buf_%s[nprocs];\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"int curr_displs_%s[nprocs];\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"size_t send_buf_size_%s[nprocs];\n", name);
				sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
						send_buf_%s[__p] = NULL;\n\
						send_buf_size_%s[__p] = 0;\n}", name, name);
			}
		}
		else {
			/* send buf count name */
			if (options->dynschedule) {
				sprintf(prog->decls+strlen(prog->decls),
						"std::vector<double *> send_buf_%s[_num_threads];\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"std::vector<size_t> send_buf_size_%s[_num_threads];\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"for (__tid=0; __tid<_num_threads; __tid++) {\n\
                        for (_i=0; _i<__MAX_NUM_SENDS; _i++) {\n\
						send_buf_%s[__tid].push_back(NULL);\n\
						send_buf_size_%s[__tid].push_back(0);\n}\n}",
						name, name);
			}
			else {
				sprintf(prog->decls+strlen(prog->decls),
						"int send_count_%s = 0;\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"int curr_displs_%s[nprocs];\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"double *send_buf_%s = NULL;\n", name);
				sprintf(prog->decls+strlen(prog->decls),
						"size_t send_buf_size_%s = 0;\n", name);
			}
		}
		/* To store previous size */
		if (options->dynschedule) {
			sprintf(prog->decls+strlen(prog->decls),
					"std::vector<double *> recv_buf_%s;\n", name);
			sprintf(prog->decls+strlen(prog->decls),
					"size_t recv_buf_size_%s[__MAX_NUM_RECVS];\n", name);
			sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<__MAX_NUM_RECVS; __p++) {\n\
                    recv_buf_%s.push_back(NULL);\n\
                    recv_buf_size_%s[__p] = 0;\n}", name, name);
		}
		else {
			sprintf(prog->decls+strlen(prog->decls),
					"double *recv_buf_%s = NULL;\n", name);
			sprintf(prog->decls+strlen(prog->decls),
					"size_t recv_buf_size_%s = 0;\n", name);
		}
		sprintf(prog->decls+strlen(prog->decls),
				"double *lw_buf_%s = NULL;\n", name);
		sprintf(prog->decls+strlen(prog->decls),
				"size_t lw_buf_size_%s = 0;\n", name);
		sprintf(prog->decls+strlen(prog->decls),
				"double *lw_recv_buf_%s = NULL;\n", name);
		sprintf(prog->decls+strlen(prog->decls),
				"size_t lw_recv_buf_size_%s = 0;\n", name);
	}

	free(waccs);
}

void init_packfile()
{
	char appendfilename[15];
	strcpy(appendfilename,__APPENDFILENAME);
	FILE *appendfp = fopen(".appendfilename", "w");
	assert(appendfp != NULL);
	fprintf(appendfp, "%s\n", appendfilename);
	fclose(appendfp);
	FILE *packfp = fopen(appendfilename, "w");
	fclose(packfp);
}

void pluto_dist_initalize_array_domains(PlutoProg *prog, int nloops){

	int i, j;
	for (i = 0; i < prog->narrays; ++i) {
		Array *arr = prog->arrays[i];

		arr->parmetric_domain = (PlutoConstraints **)malloc(nloops * sizeof(PlutoConstraints *));
        for (j=0; j<nloops; j++) {
            arr->parmetric_domain[j] = NULL;
        }
		arr->copy_level = (int *)malloc(nloops * sizeof(int));

		arr->array_bounds = NULL;
	}

	return;
}

void init_pi_mappings(PlutoProg *prog, Ploop **loops, int nloops, int *pi_mappings)
{
	int l, i;
	for (i=0; i<prog->nstmts; i++) {
		pi_mappings[i] = -1;
	}

	for (l=0; l<nloops; l++) {
		for (i=0; i<loops[l]->nstmts; i++) {
			pi_mappings[loops[l]->stmts[i]->id] = l;
		}
	}
	// if pi_mappings of a statement is -1, then that statement is not distributed
}

void pluto_dynschedule_common_parallelize(PlutoProg *prog, FILE *sigmafp, FILE *headerfp,
		Ploop **loops, int nloops, int *pi_mappings, int *copy_level)
{
	int l, i;
	int *num_data = (int *) malloc(nloops*sizeof(int));
	Stmt ***tasks_stmts = (Stmt ***) malloc(nloops*sizeof(Stmt **));

	char type[64];
	if (options->dynschedule_graph) {
		strcpy(type, "tbb::flow::continue_node<tbb::flow::continue_msg> **");
	}
	else {
		strcpy(type, "int *");
	}
	char *tasks_loops_decl = malloc(512);
	strcpy(tasks_loops_decl, "");
	for (l=0; l<nloops; l++) {
		sprintf(tasks_loops_decl+strlen(tasks_loops_decl), ",%stasks_loop%d", type, l);
	}

	if (options->dynschedule_graph) {
		gen_dynschedule_graph_header_text_code(prog, copy_level, nloops, headerfp);
	}
	else {
		gen_dynschedule_header_text_code(copy_level, nloops, headerfp);
	}

	for (l=0; l<nloops; l++) {
		Ploop *loop = loops[l];

		int *num_stmts_per_wacc; // indexed by data variable
		struct stmt_access_pair ***wacc_stmts; // indexed by data variable
		wacc_stmts = get_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_data[l], &num_stmts_per_wacc);

		tasks_stmts[l] = gen_tasks_code(loop->stmts, loop->nstmts, copy_level, prog, num_data[l],
				l, nloops, pi_mappings, tasks_loops_decl, sigmafp, headerfp);

		for (i=0; i<num_data[l]; i++) {
			free(wacc_stmts[i]);
		}
		free(wacc_stmts);
		free(num_stmts_per_wacc);
	}

	gen_init_tasks_cloog_code(prog, tasks_stmts, nloops, tasks_loops_decl, headerfp);

	free(tasks_loops_decl);
	for(l=0; l<nloops; l++) {
		free(tasks_stmts[l]);
	}
	free(tasks_stmts);
	free(num_data);
}


/* Dynamic scheduling parallelization - for shared-memory only */
int pluto_dynschedule_parallelize(PlutoProg *prog, FILE *sigmafp, FILE *headerfp, FILE *pifp)
{
	// print_hyperplane_properties(prog);

	FILE *pidefs = NULL;
	pidefs = fopen("pi_defs.h", "w");
	fclose(pidefs);

	fprintf(sigmafp, "\n#ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
	fprintf(sigmafp, "#include <stdio.h>\n");
	fprintf(sigmafp, "#endif\n\n");

	fprintf(sigmafp, "\
#include <math.h>\n\
#include <assert.h>\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");

	int nloops=0;
	Ploop **loops = pluto_get_dom_parallel_loops(prog, &nloops);
	/* Loops should be in the increasing order of their depths */
	qsort(loops, nloops, sizeof(Ploop *), pluto_loop_compar);

	IF_DEBUG(printf("distmem: parallelizing loops\n"););
	IF_DEBUG(pluto_loops_print(loops,nloops););

	int *pi_mappings = malloc(prog->nstmts*sizeof(int));
	init_pi_mappings(prog, loops, nloops, pi_mappings);

	int *copy_level = (int *) malloc(nloops*sizeof(int));
	int *outer_dist_loop_level = (int *) malloc(nloops*sizeof(int));

	init_packfile();

	init_copy_level(prog, loops, nloops, copy_level, outer_dist_loop_level);

	if(options->data_dist)
		pluto_dist_gen_intial_code(nloops, loops, copy_level, prog, headerfp);

	pluto_dynschedule_common_parallelize(prog, sigmafp, headerfp, loops, nloops, pi_mappings, copy_level);

	if (options->dynschedule) {
		fprintf(pifp, "\
#include <math.h>\n\
#include <assert.h>\n\
#include \"pi_defs.h\"\n\
#include \"polyrt.h\"\n\
#include \"buffer_manager.h\"\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");

		fprintf(pifp, "#ifndef __BLOCK_CYCLIC_BLOCK_SIZE\n");
		fprintf(pifp, "#define __BLOCK_CYCLIC_BLOCK_SIZE 1\n");
		fprintf(pifp, "#endif\n\n");

		if(options->compute_pi){
			pluto_dist_gen_compute_pi(nloops, loops, copy_level,
					prog, pi_mappings, pifp, headerfp);
		}

		int l;
		for (l=0; l<nloops; l++) {
			generate_pi(pifp, headerfp, outer_dist_loop_level[l], copy_level[l]-1,
					prog, loops[l], l, NULL, 0);
		}
	}

	free(outer_dist_loop_level);
	free(copy_level);
	free(pi_mappings);

	pluto_mark_statements(prog);

	pluto_loops_free(loops,nloops);

	// !!!roshan not sure if this is required
	// since no dimensions have been added to the program
	pluto_compute_dep_satisfaction_complex(prog);

	return (nloops==0);
}

char *get_add_vertex_weight_stmt_text(int phase_level, int **is_scalar_dim, int loop_num)
{

	char *text = malloc(1024);
	text[0] = '\0';

	sprintf(text+strlen(text), "cur_phase = prev_phase ");
	int i, j;
	for(i=0;i<phase_level;i++){

		if(is_scalar_dim[loop_num][i]) continue;

		sprintf(text+strlen(text), " + ");

		sprintf(text+strlen(text), "t%d ", i+1);
		for(j=i+1;j<phase_level;j++)
			sprintf(text+strlen(text), " * MAX_TILES_PER_DIM");
	}

	sprintf(text+strlen(text), ";");
	sprintf(text+strlen(text), "for(i=0;i<total_num_phases;i++) {vertex_wgt[vertex_count++] = (i==cur_phase)?1:0;}");

	return text;
}

Stmt **gen_vertex_wgt_stmt_code (struct stmt_access_pair **acc_stmts, int num_accs,
		PlutoProg *prog, int *copy_level, int *phase_levels, int **is_scalar_dim, int loop_num)
{
	int src_copy_level;
	src_copy_level = copy_level[loop_num];
	int phase_level = phase_levels[loop_num];

	Stmt *anchor_stmt = acc_stmts[0]->stmt;

	char *stmt_text = get_add_vertex_weight_stmt_text(phase_level, is_scalar_dim, loop_num);

	char *gaurd_stmt_text = malloc(1024);
	sprintf(gaurd_stmt_text, "prev_phase = cur_phase + 1;");

	Stmt *vertex_wights_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			stmt_text, IN_FUNCTION, loop_num);

	Stmt *gaurd_stmt = create_helper_stmt(anchor_stmt, 0,
			gaurd_stmt_text, IN_FUNCTION, loop_num);

	Stmt **vertex_wights_stmts = (Stmt **) malloc(2*sizeof(Stmt *));
	vertex_wights_stmts[0] = vertex_wights_stmt;
	vertex_wights_stmts[1] = gaurd_stmt;

	free(stmt_text);
	free(gaurd_stmt_text);

	return vertex_wights_stmts;
}

char *get_add_edge_weight_stmt_text()
{

	char *text = malloc(1024);
	text[0] = '\0';

	sprintf(text, "for(i=0;i<total_num_phases;i++) {vetex_wgt[count++] = (i==cur_phase)?1:0}");

	return text;
}

Stmt **gen_edge_wgt_stmt_code (struct stmt_access_pair ***acc_stmts, int num_accs,
		Stmt *anchor_stmt, PlutoProg *prog, int *copy_level, int non_scalar_copy_level,
		int **is_scalar_dim, int max_copy_level, int loop_num)
{
	int i,j, src_copy_level;
	src_copy_level = copy_level[loop_num];

	char *text = malloc(1024 * 5);
	text[0] = '\0';

	sprintf(text+strlen(text), "clear_comm_cost(comm_cost); ");

	for (j = 0; j < num_accs; ++j) {
		char *acc_name = acc_stmts[j][0]->acc->name;
		sprintf(text+strlen(text), "edge_wgt_%s_%d(", acc_name, loop_num);

		for(i=0;i<src_copy_level;i++){

			//			if(is_scalar_dim[loop_num][i]) continue;

			if(i==0)
				sprintf(text+strlen(text), " t%d", i+1);
			else
				sprintf(text+strlen(text), ", t%d", i+1);
		}

		for(i=0;i<prog->npar;i++){
			sprintf(text+strlen(text), ", %s",prog->params[i] );
		}
		sprintf(text+strlen(text), ", comm_cost); ");
	}

	sprintf(text+strlen(text), "cur_vertex = get_corrs_node(%d", loop_num);

	int k;
	for(k =0 ; k < src_copy_level;k++){
		if(is_scalar_dim[loop_num][k]) continue;
		sprintf(text+strlen(text), ", t%d ", k+1);
	}
	for(k =0 ; k < max_copy_level - non_scalar_copy_level;k++)
		sprintf(text+strlen(text), ", 0");

	sprintf(text+strlen(text), "); ");

	sprintf(text+strlen(text), "update_edges(cur_vertex, edges, &num_edges, edge_wgts, vadj, &num_vadj, comm_cost); ");
	Stmt *edge_wights_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
			text, IN_FUNCTION, loop_num);

	Stmt **edge_wights_stmts = (Stmt **) malloc(1*sizeof(Stmt *));
	edge_wights_stmts[0] = edge_wights_stmt;

	free(text);

	return edge_wights_stmts;
}

void gen_comm_vol_code( struct stmt_access_pair **wacc_stmts, int num_accs,
		int nloops,  PlutoProg *prog, Stmt *anchor_stmt, int *copy_level,
		int **is_scalar_dim, int loop_num, int *pi_mappings, FILE *headerfp, FILE *outfp)
{
	int i, j, k, l, src_copy_level, acc_nrows;
	src_copy_level = copy_level[loop_num];

	assert(num_accs >= 1);
	acc_nrows = wacc_stmts[0]->acc->mat->nrows;
	char *acc_name = wacc_stmts[0]->acc->name;

	PlutoConstraints *foifi_sets[nloops];
	for (l=0; l<nloops; l++) {
		foifi_sets[l] = NULL;
	}

	int broadcast = 0;
	PlutoConstraints *flow_out = NULL;
	for (k=0; k<num_accs; k++)  {
		PlutoConstraints *flow_out_one;
		flow_out_one = compute_flow_out(wacc_stmts[k], src_copy_level, copy_level, prog, pi_mappings);
		if (flow_out == NULL) flow_out = pluto_constraints_dup(flow_out_one);
		else{
			flow_out = pluto_constraints_unionize(flow_out, flow_out_one);
		}
		pluto_constraints_free(flow_out_one);

		if (broadcast) continue;

		for (i=0; i<prog->ndeps; i++)   {
			Dep *dep = prog->deps[i];

			if (options->data_dist) {
				/* Only RAW and RAR deps matter when data is distributed*/
				if (dep->type != OSL_DEPENDENCE_RAW && dep->type != OSL_DEPENDENCE_RAR) continue;
			}
			else {
				/* Only RAW deps matter */
				if (dep->type != OSL_DEPENDENCE_RAW) continue;
			}

			// if (dep->dirvec[copy_level+1] == DEP_ZERO) continue;

			/* If the dependence doesn't originate from this access */
			//              if (dep->src_acc != wacc_stmts[k]->acc && dep->dest_acc != wacc_stmts[k]->acc) continue;
			if (dep->src_acc != wacc_stmts[k]->acc) continue;

			assert(dep->dest_acc != NULL);

			Stmt *dest = prog->stmts[dep->dest];
			int dependent_loop = pi_mappings[dest->id];
			if (dependent_loop == -1) {
				// destination statement will be executed by all processors
				broadcast = 1;
				break;
			}

			int dest_copy_level = copy_level[dependent_loop];

			PlutoConstraints *fo = compute_flow_out_of_dep(dep, src_copy_level, copy_level, prog, 0, NULL, pi_mappings);
			print_polylib_visual_sets("fo", fo);
			// add target tile iterators for flow-out set
			for (j=0; j<dest_copy_level; j++) {
				pluto_constraints_add_dim(fo,src_copy_level, NULL);
			}

			PlutoConstraints *fi = compute_flow_in_of_dep(dep, dest_copy_level, prog, 0);
			print_polylib_visual_sets("fi", fi);
			// add source tile iterators for flow-in set
			for (j=0; j<src_copy_level; j++) {
				pluto_constraints_add_dim(fi,0, NULL);
			}

			PlutoConstraints *foifi = pluto_constraints_intersection(fo, fi);

			print_polylib_visual_sets(dep->src_acc->name, foifi);

			if (foifi_sets[dependent_loop] == NULL)
				foifi_sets[dependent_loop] = pluto_constraints_dup(foifi);
			else
				pluto_constraints_unionize(foifi_sets[dependent_loop], foifi);

			pluto_constraints_free(fo);
			pluto_constraints_free(fi);
			pluto_constraints_free(foifi);
		}
	}

	char *stmt_text = malloc(1024);
	sprintf(stmt_text, "count++;");

	char *iters[acc_nrows];
	for(i=0;i<acc_nrows;i++){
		iters[i] = malloc(10);
		sprintf(iters[i], "d%d", i+1);
	}

	for(l=0;l<nloops;l++){

		if(foifi_sets[l] == NULL) continue;

		int total_copy_level = copy_level[l] + src_copy_level + acc_nrows;


		for (j=0; j<acc_nrows; j++) {
			pluto_constraints_add_dim(foifi_sets[l],0, NULL);
		}
		for (j=0; j<acc_nrows; j++) {
			pluto_constraints_interchange_cols(foifi_sets[l], j, j+total_copy_level);
		}
		for (j=0; j<acc_nrows; j++) {
			pluto_constraints_remove_dim(foifi_sets[l], total_copy_level);
		}

		PlutoProg *sigma = pluto_prog_alloc();

		for (i=0; i<src_copy_level; i++) {
			//			if(is_scalar_dim[loop_num][i]) continue;

			char param[6];
			sprintf(param, "ts%d",i+1);
			pluto_prog_add_param(sigma, param, sigma->npar);
		}

		for (i=0; i<copy_level[l]; i++) {
			//			if(is_scalar_dim[l][i]) continue;

			char param[6];
			sprintf(param, "td%d",i+1);
			pluto_prog_add_param(sigma, param, sigma->npar);
		}

		for (i=0; i<prog->npar; i++) {
			pluto_prog_add_param(sigma, prog->params[i], sigma->npar);
		}

		total_copy_level = copy_level[l] + src_copy_level;

		PlutoMatrix *trans = get_data_tile_trans_func(total_copy_level, acc_nrows, prog);

		pluto_add_stmt(sigma, foifi_sets[l], trans, iters, stmt_text, IN_FUNCTION);

		pluto_matrix_free(trans);

		pluto_pad_stmt_transformations(sigma);

		pluto_separate_stmts(sigma, sigma->stmts, sigma->nstmts,0,0);

		FILE *cloogfp = NULL;

		fprintf(outfp, "int comm_vol_%s_%d_%d(", acc_name, loop_num, l);
		fprintf(headerfp, "int comm_vol_%s_%d_%d(", acc_name, loop_num, l);

		int first = 0;
		for(i=0;i<src_copy_level;i++){

			//			if(is_scalar_dim[loop_num][i]) continue;

			if(first==0){
				fprintf(outfp,"int ts%d", i+1);
				fprintf(headerfp,"int ts%d", i+1);
				first = 1;
			}
			else{
				fprintf(outfp,",int ts%d", i+1);
				fprintf(headerfp,",int ts%d", i+1);
			}
		}

		for(i=0;i<copy_level[l];i++){

			//			if(is_scalar_dim[l][i]) continue;

			fprintf(outfp,",int td%d", i+1);
			fprintf(headerfp,",int td%d", i+1);
		}

		for(i=0;i<prog->npar;i++){
			fprintf(outfp,",int %s ", prog->params[i] );
			fprintf(headerfp,",int %s ", prog->params[i] );
		}

		fprintf(outfp, ")\n{\n");
		fprintf(headerfp, ");\n");

		fprintf(outfp, "int count = 0;\n");

		cloogfp = fopen("sigma.cloog", "w+");
		pluto_gen_cloog_file(cloogfp, sigma);
		rewind(cloogfp);
		generate_declarations(sigma, outfp);
		pluto_gen_cloog_code(sigma, -1, -1, cloogfp, outfp);
		fclose(cloogfp);

		for (i=0; i<sigma->nstmts; i++) {
			fprintf(outfp, "#undef S%d\n", i+1);
		}


		fprintf(outfp, "printf(\" comm vol of %s_%d_%d  from ",  acc_name, loop_num, l);
		for(i=0;i<src_copy_level;i++){
			fprintf(outfp,"%%d ");
		}

		fprintf(outfp," to ");
		for(i=0;i<copy_level[l];i++){
			fprintf(outfp,"%%d ");
		}
		fprintf(outfp,"  = %%d\\n\"  ");
		for(i=0;i<src_copy_level;i++){
			fprintf(outfp,", ts%d ", i+1);
		}
		for(i=0;i<copy_level[l];i++){
			fprintf(outfp,", td%d ", i+1);
		}

		fprintf(outfp,", count);\n ");

		fprintf(outfp, "return count;\n");
		fprintf(outfp, "}\n\n");
		pluto_prog_free(sigma);

	}
	free(stmt_text);

}


void gen_edge_wgt_func_code(struct stmt_access_pair **wacc_stmts,
		int naccs, int *copy_level, int *scalar_dim_copy_level ,int **is_scalar_dim, int max_copy_level, PlutoProg *prog, int loop_num,
		int nloops,  int *pi_mappings, FILE *outfp, FILE *headerfp)
{


	int i,j, src_copy_level;
	src_copy_level = copy_level[loop_num];
	char *acc_name = wacc_stmts[0]->acc->name;

	PlutoProg *sigma = pluto_prog_alloc();

	for (i=0; i<src_copy_level; i++) {

		//		if(is_scalar_dim[loop_num][i]) continue;

		char param[6];
		sprintf(param, "ts%d",i+1);
		pluto_prog_add_param(sigma, param, sigma->npar);
	}

	for (i=0; i<prog->npar; i++) {
		pluto_prog_add_param(sigma, prog->params[i], sigma->npar);
	}

	int num_dest_loops;
	PlutoConstraints *tdpoly_list[prog->ndeps];
	int dep_loop_num_list[prog->ndeps];
	PlutoMatrix *trans[prog->ndeps];
	char **iters[prog->ndeps];
	char *indices[prog->ndeps];

	int broadcast = 0;

	generate_sigma_common(wacc_stmts, naccs, copy_level, prog, loop_num, pi_mappings, 1,
			&broadcast,
			&num_dest_loops,
			tdpoly_list,
			dep_loop_num_list,
			trans,
			iters,
			indices);

	for (i=0; i<num_dest_loops; i++) {
		PlutoConstraints *tdpoly = tdpoly_list[i];
		int dep_loop_num = dep_loop_num_list[i];
		int dest_copy_level = copy_level[dep_loop_num];

		char *sigma_text = malloc(1024 *4);
		sprintf(sigma_text, "comm_cost");
		sprintf(sigma_text+strlen(sigma_text),"[%d]", dep_loop_num);
		for(j=0;j<dest_copy_level;j++){

			if(is_scalar_dim[dep_loop_num][j]) continue;

			sprintf(sigma_text+strlen(sigma_text),"[%s]", iters[i][j]);
		}
		for(j=0;j<max_copy_level - scalar_dim_copy_level[dep_loop_num];j++){
			sprintf(sigma_text+strlen(sigma_text),"[0]");
		}

		sprintf(sigma_text+strlen(sigma_text), "+= comm_vol_%s_%d_%d(", acc_name, loop_num, dep_loop_num);
		for(j=0;j<src_copy_level;j++){
			//		   if(is_scalar_dim[loop_num][j]) continue;

			if(j==0)
				sprintf(sigma_text+strlen(sigma_text)," ts%d",j+1 );
			else
				sprintf(sigma_text+strlen(sigma_text),", ts%d",j+1 );
		}

		for(j=0;j<dest_copy_level;j++){
			//		   if(is_scalar_dim[dep_loop_num][j]) continue;

			sprintf(sigma_text+strlen(sigma_text),", %s", iters[i][j]);
		}

		for(j=0;j<prog->npar;j++){
			sprintf(sigma_text+strlen(sigma_text),", %s", prog->params[j]);
		}

		sprintf(sigma_text+strlen(sigma_text),");");

		pluto_add_stmt(sigma,tdpoly,trans[i],iters[i],sigma_text, IN_FUNCTION);

		for (j=0; j<dest_copy_level; j++) {
			free(iters[i][j]);
		}

		free(sigma_text);
		free(indices[i]);
		free(iters[i]);
		pluto_matrix_free(trans[i]);
		pluto_constraints_free(tdpoly);
	}

	pluto_pad_stmt_transformations(sigma);
	//   pluto_prog_print(sigma);
	if (sigma->nstmts >= 1) {
		assert(sigma->stmts[0]->trans->nrows == sigma->num_hyperplanes);
	}

	pluto_separate_stmts(sigma, sigma->stmts, sigma->nstmts, 0, 0);

	FILE *cloogfp = NULL;

	fprintf(outfp, "void edge_wgt_%s_%d(", acc_name, loop_num);
	fprintf(headerfp, "void edge_wgt_%s_%d(", acc_name, loop_num);
	for (i=0; i<src_copy_level+prog->npar; i++)    {
		if (i!=0) {
			fprintf(outfp, ", ");
			fprintf(headerfp, ", ");
		}
		if (i<=src_copy_level-1) {
			fprintf(outfp, "int ts%d", i+1);
			fprintf(headerfp, "int ts%d", i+1);
		}
		else {
			fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
			fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
		}
	}
	fprintf(outfp, ", long comm_cost[%d]", nloops);
	fprintf(headerfp, ", long comm_cost[%d]", nloops);

	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[MAX_TILES_PER_DIM]");
		fprintf(headerfp, "[MAX_TILES_PER_DIM]");
	}

	fprintf(outfp, ")\n{\n");
	fprintf(headerfp, ");\n");

	cloogfp = fopen("sigma.cloog", "w+");
	pluto_gen_cloog_file(cloogfp, sigma);
	rewind(cloogfp);
	generate_declarations(sigma, outfp);
	pluto_gen_cloog_code(sigma, -1, -1, cloogfp, outfp);
	fclose(cloogfp);

	for (i=0; i<sigma->nstmts; i++) {
		fprintf(outfp, "#undef S%d\n", i+1);
	}


	fprintf(outfp, "}\n\n");
	pluto_prog_free(sigma);
}

void fprint_clear_comm_cost(FILE *outfp, int max_copy_level, int nloops){

	fprintf(outfp, "void clear_comm_cost(");
	fprintf(outfp, "long comm_cost[%d]", nloops);
	int i = 0;

	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[MAX_TILES_PER_DIM]");
	}
	fprintf(outfp, "){\n");
	fprintf(outfp, "int l ");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, ", t%d", i+1);
	}
	fprintf(outfp, ";\n");
	fprintf(outfp, "for(l = 0; l < %d; l++)\n", nloops);
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "for(t%d = 0; t%d < MAX_TILES_PER_DIM; t%d++)\n", i+1, i+1, i+1);
	}
	fprintf(outfp, "comm_cost[l]");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[t%d]", i+1);
	}
	fprintf(outfp, " = 0;\n return; }");

}

void fprint_get_node_func(FILE *outfp, int max_copy_level, int num_loops, int *copy_level){
	fprintf(outfp, "int get_corrs_node(int l");
	int i, j, k;
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, ", int t%d", i+1);
	}
	fprintf(outfp, "){\n");
	fprintf(outfp, "int i, node = 0, offset = 0, vertices[%d]; \n", num_loops);
	for(i=0;i<num_loops;i++){
		fprintf(outfp, " vertices[%d] = ",i);
		for(j=0;j<copy_level[i];j++){
			if(j !=0 )
				fprintf(outfp, " * ");
			fprintf(outfp, " MAX_TILES_PER_DIM ");
		}
		fprintf(outfp, ";\n");
	}

	fprintf(outfp, "for(i=0;i<l;i++)\n"
			"node += vertices[i];\n");

	fprintf(outfp, "switch(l) {\n ");


	for(k=0;k<num_loops;k++){
		fprintf(outfp, "case %d: \n\t offset = ",k);
		for(i=0;i<copy_level[k];i++){

			if(i!=0)
				fprintf(outfp, " + ");

			fprintf(outfp, " t%d * ( ", i+1);

			if( i == copy_level[k] - 1)
				fprintf(outfp, " 1 ");

			for(j = 0; j < copy_level[k] - i - 1; j++){
				if(j==0)
					fprintf(outfp, "MAX_TILES_PER_DIM ");
				else
					fprintf(outfp, "* MAX_TILES_PER_DIM ");
			}

			fprintf(outfp, ") ");
		}
		fprintf(outfp, "; break;\n ");
	}
	fprintf(outfp, "default: offset = 0; break;\n}\nnode += offset;\n");


	fprintf(outfp, "return node;\n }\n");
}
void fprint_partition_func(FILE *outfp){

	fprintf(outfp, "void partition_print(long num_nodes, long *partition, long obj_val, long ret){"
			"\n int i = 0; \n"
			"remove(\"partition\");\n"
			"FILE *fp = fopen(\"partition\", \"w\");\n"
			"fprintf(fp,\"obj_val = %%d  ret = %%d \\n \", obj_val, ret);\n"
			"for(i=0;i<num_nodes;i++)\n"
			"fprintf(fp, \" #%%d proc = %%d\\n \", i, partition[i]);\n"
			"fclose(fp);\n"
			"return;\n"
			"}\n ");
}

void fprint_debug_graph_func(FILE *outfp){

	fprintf(outfp, "void debug_graph_print(");
	fprintf(outfp, " long *edges, long num_edges, long  *edge_wgts, long *vertex_wgt, long num_vertex_wgt, long total_num_phases, "
			"long *vadj, long num_vadj){ \n int i, j;\n");
	fprintf(outfp, "remove(\"debug_graph\");\n");
	fprintf(outfp, "FILE *fp = fopen(\"debug_graph\", \"w\");\n");
	fprintf(outfp, "fprintf(fp,\"#vertices = %%d\\n vertex weights\\n\", num_vertex_wgt);\n" );
	fprintf(outfp, " \n for(i = 0; i< num_vertex_wgt; i++) { if((i%%total_num_phases) == 0) fprintf(fp, \"\\n\"); \n fprintf(fp, \"%%d \\t\", vertex_wgt[i]);\n }\n " );
	fprintf(outfp, "fprintf(fp,\"\\n#edges = %%d\\n Edge weights\\n\", num_edges);\n" );
	fprintf(outfp, " \n for(i = 0; i< num_vadj-1; i++) { "
			" \n fprintf(fp, \"(%%d)  \", i);"
			"for(j=vadj[i]; j<vadj[i+1]; j++){ "
			"fprintf(fp, \"%%d \\t\", edge_wgts[j]);"
			"}"
			" \n fprintf(fp, \"\\n\");"
			"}\n" );
	fprintf(outfp, "fprintf(fp,\"\\nEdges \\n\");\n" );
	fprintf(outfp, " \n for(i = 0; i< num_vadj-1; i++) {"
			" \n fprintf(fp, \"(%%d)  \", i);"
			" for(j=vadj[i]; j<vadj[i+1]; j++){ "
			"fprintf(fp, \"%%d \\t\", edges[j]);"
			"} "
			"\n fprintf(fp, \"\\n\");"
			"}\n" );

	fprintf(outfp, "fclose(fp);\n return;\n}\n");

	return;
}

void fprint_update_edges_func(FILE *outfp,  int max_copy_level, int nloops){

	fprintf(outfp, "void update_edges(");
	int i = 0;

	fprintf(outfp, "long cur_vertex, long *edges, long  *num_edges, long *edge_wgts, long  *vadj, long *num_vadj, long comm_cost[%d]", nloops);

	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[MAX_TILES_PER_DIM]");
	}

	fprintf(outfp, "){\n");
	fprintf(outfp, "int l, i, j ");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, ", t%d", i+1);
	}
	fprintf(outfp, ";\n");

	fprintf(outfp, " \n for(i = 0; i< (*num_vadj)-1; i++){"
			"  for(j=vadj[i]; j<vadj[i+1]; j++){ "
			"if(cur_vertex == edges[j]){"
			"  edges[(*num_edges)] = i; "
			"  edge_wgts[(*num_edges)++] = edge_wgts[j];"
			"}"
			"}"
			"}\n" );

	fprintf(outfp, "for(l = 0; l < %d; l++)\n", nloops);
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "for(t%d = 0; t%d < MAX_TILES_PER_DIM; t%d++)\n", i+1, i+1, i+1);
	}
	fprintf(outfp, "if(comm_cost[l]");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[t%d]", i+1);
	}
	fprintf(outfp, " != 0) { \n edge_wgts[(*num_edges)] = comm_cost[l]");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[t%d]", i+1);
	}

	fprintf(outfp, ";  edges[(*num_edges)++] = get_corrs_node(l");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, ", t%d", i+1);
	}
	fprintf(outfp, "); \n");

	fprintf(outfp, "\nprintf(\"comm_cost[%%d]");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[%%d]");
	}
	fprintf(outfp, " = %%d \\n \", l ");
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, ", t%d", i+1);
	}
	fprintf(outfp, ", edges[(*num_edges)-1] );\n");

	fprintf(outfp, "}\n ");

	fprintf(outfp, "vadj[(*num_vadj)++] = *num_edges;\n ");

	fprintf(outfp, "return ;");
	fprintf(outfp, "}\n ");
}

char *get_total_num_phases(int nloops, Ploop **loops, int *copy_level, int *outer_dist_loop_level,  int **is_scalar_dim, PlutoProg *prog){

	char *total_phases = malloc(1024);
	total_phases[0] = 0;

	sprintf(total_phases+strlen(total_phases), " 0 ");

	int l;

	for(l=0; l<nloops; l++) {

		int j;
		int first = 0;
		for(j=0;j<outer_dist_loop_level[l];j++){

			if(is_scalar_dim[l][j]) continue;

			if(first==0){
				sprintf(total_phases+strlen(total_phases), "+ ");
				first = 1;
			}
			else
				sprintf(total_phases+strlen(total_phases), " * ");

			sprintf(total_phases+strlen(total_phases), " MAX_TILES_PER_DIM ");
		}

		if(!first)
			sprintf(total_phases+strlen(total_phases), " + 1 ");
	}

	return total_phases;
}

int* get_num_non_scalar_dims(int nloops, Ploop **loops, int *copy_level, int **is_scalar_dim){

	int *non_scalar_copy_level = (int *)malloc(nloops * sizeof(int));

	int i, j;
	for(i=0;i<nloops;i++){
		non_scalar_copy_level[i]= 0;
		is_scalar_dim[i] = malloc(copy_level[i] * sizeof(int));

	}

	for(i=0;i<nloops;i++){
		Stmt *stmt = loops[i]->stmts[0];
		for(j=0;j<copy_level[i];j++){

			if(stmt->hyp_types[j] != H_SCALAR){
				non_scalar_copy_level[i]++;
				is_scalar_dim[i][j] = 0;
			}
			else
				is_scalar_dim[i][j] = 1;
		}
	}

	return non_scalar_copy_level;
}

int pluto_dist_gen_compute_pi(int nloops, Ploop **loops, int *copy_level,
		PlutoProg *prog, int *pi_mappings, FILE *outfp, FILE *headerfp){

	int i, l;

	if(!options->compute_pi)
		return 0;

	int *num_write_data = (int *) malloc(nloops*sizeof(int));
	int *num_read_write_data = (int *) malloc(nloops*sizeof(int));

	int *num_stmts_per_wacc;
	int *num_stmts_per_rwacc;
	struct stmt_access_pair ***wacc_stmts;
	struct stmt_access_pair ***rwacc_stmts;

	Stmt ****vertex_wgt_stmts, ****edge_wgt_stmts;
	vertex_wgt_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	edge_wgt_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));

	int **is_scalar_dim = malloc(nloops * sizeof(int));
	int *non_scalar_copy_level = get_num_non_scalar_dims(nloops, loops, copy_level, is_scalar_dim);

	int max_copy_level = non_scalar_copy_level[0];
	for (i=1; i<nloops; i++) {
		if (max_copy_level < non_scalar_copy_level[i]) {
			max_copy_level = non_scalar_copy_level[i];
		}
	}

	fprintf(outfp, "#include <stdio.h> \n\n");
	fprintf(outfp, "#include <metis.h> \n\n");
	fprintf(outfp, "#define MAX_TILES_PER_DIM %d \n", options->num_tiles_per_dim);
	fprintf(outfp, "#define MAX_NUM_PROCS 32 \n\n");
	fprintf(headerfp, "#define MAX_TILES_PER_DIM %d \n\n", options->num_tiles_per_dim);

	for(l=0;l<nloops;l++){
		fprintf(outfp, "int pi_mappings_%d", l);
		for(i=0;i<non_scalar_copy_level[l];i++){
			fprintf(outfp,"[MAX_NUM_PROCS]");
		}
		fprintf(outfp, ";\n");
	}

	for(l=0; l<nloops; l++) {

		Ploop *loop = loops[l];

		Stmt *anchor_stmt = get_new_anchor_stmt(loop->stmts, loop->nstmts);

		wacc_stmts = get_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_write_data[l], &num_stmts_per_wacc);

		rwacc_stmts = get_read_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_read_write_data[l], &num_stmts_per_rwacc);

		vertex_wgt_stmts[l] = (Stmt ***) malloc(num_write_data[l]*sizeof(Stmt **));
		edge_wgt_stmts[l] = (Stmt ***) malloc(num_read_write_data[l]*sizeof(Stmt **));

		vertex_wgt_stmts[l][0] =
				gen_vertex_wgt_stmt_code(wacc_stmts[0],num_stmts_per_wacc[0],
						prog,copy_level,non_scalar_copy_level, is_scalar_dim, l);

		for(i=0;i<num_read_write_data[l];i++){

			gen_comm_vol_code(rwacc_stmts[i], num_stmts_per_rwacc[i], nloops, prog,
					anchor_stmt, copy_level,
					is_scalar_dim, l,pi_mappings, headerfp, outfp);

			gen_edge_wgt_func_code(rwacc_stmts[i], num_stmts_per_rwacc[i], copy_level, non_scalar_copy_level, is_scalar_dim, max_copy_level,
					prog, l,nloops, pi_mappings, outfp, headerfp);

		}

		edge_wgt_stmts[l][0] =
				gen_edge_wgt_stmt_code(rwacc_stmts,num_read_write_data[l],
						anchor_stmt, prog,copy_level, non_scalar_copy_level[l],
						is_scalar_dim, max_copy_level, l);

		for (i=0; i<num_write_data[l]; i++) {
			free(wacc_stmts[i]);
		}
		free(wacc_stmts);

		for (i=0; i<num_read_write_data[l]; i++) {
			free(rwacc_stmts[i]);
		}
		free(rwacc_stmts);
	}

	PlutoProg *compute_pi = pluto_prog_alloc();
	data_tile_prog_add_parm(compute_pi, prog, 0);

	for (l = 0; l < nloops; ++l) {
		for (i = 0; i < 2; ++i) {
			pluto_add_given_stmt(compute_pi, vertex_wgt_stmts[l][0][i]);
		}

		pluto_add_given_stmt(compute_pi, edge_wgt_stmts[l][0][0]);
	}
	pluto_pad_stmt_transformations(compute_pi);
	pluto_separate_stmts(compute_pi, compute_pi->stmts, compute_pi->nstmts, 0, 0);

	FILE *cloogfp = NULL;

	fprint_get_node_func(outfp, max_copy_level,nloops,non_scalar_copy_level);
	fprint_clear_comm_cost(outfp, max_copy_level, nloops);
	fprint_update_edges_func(outfp, max_copy_level, nloops);
	fprint_debug_graph_func(outfp);
	fprint_partition_func(outfp);

	fprintf(outfp, "void compute_pi(int nprocs");
	fprintf(headerfp, "void compute_pi(int nprocs");
	for (i=0; i<prog->npar; i++)    {

		//		if (i!=0) {
		fprintf(outfp, ", ");
		fprintf(headerfp, ", ");
		//		}

		fprintf(outfp, "int %s", prog->params[i]);
		fprintf(headerfp, "int %s", prog->params[i]);
	}

	fprintf(outfp, ")\n{\n");
	fprintf(headerfp, ");\n");

	int *tile_sizes = malloc(prog->npar * sizeof(int));
	int num_tile_dims = prog->npar;
	pluto_dist_read_tile_sizes(tile_sizes, num_tile_dims);

	for(i=0;i<prog->npar;i++)
		fprintf(outfp, "%s = MAX_TILES_PER_DIM * %d;\n", prog->params[i], tile_sizes[i]);

	cloogfp = fopen("sigma.cloog", "w+");
	pluto_gen_cloog_file(cloogfp, compute_pi);
	rewind(cloogfp);
	generate_declarations(compute_pi, outfp);
	fprintf(outfp, "long comm_cost[%d]", nloops);

	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "[MAX_TILES_PER_DIM]");
	}
	fprintf(outfp, ";\n");
	fprintf(outfp, "long num_edges = 0, num_vadj = 0, cur_vertex = 0;\n");
	fprintf(outfp, "long total_vertices = %d ", nloops);
	for(i=0;i<max_copy_level;i++){
		fprintf(outfp, "* MAX_TILES_PER_DIM ");
	}
	fprintf(outfp, ";\n");
	fprintf(outfp, "long total_edges = total_vertices * total_vertices;\n");
	fprintf(outfp, "long vadj[total_vertices+1];\n");
	fprintf(outfp, "long edges[total_edges];\n");
	fprintf(outfp, "long edge_wgts[total_edges];\n");
	fprintf(outfp, "long i, total_num_phases = %s, vertex_count = 0, cur_phase = 0, prev_phase = 0;\n",
			get_total_num_phases(nloops, loops, copy_level,non_scalar_copy_level, is_scalar_dim, prog));

	fprintf(outfp, "long vertex_wgt[total_vertices * total_num_phases];\n");
	fprintf(outfp, "long partitions[total_vertices], objval = 0, ret = 0, parts = nprocs;\n");

	fprintf(outfp, "vadj[0] = 0; num_vadj = 1;\n");

	pluto_gen_cloog_code(compute_pi, -1, -1, cloogfp, outfp);
	fclose(cloogfp);

	fprintf(outfp, "total_vertices = num_vadj - 1;\n");
	fprintf(outfp, "debug_graph_print(edges, num_edges, edge_wgts, vertex_wgt, vertex_count, "
			"total_num_phases, vadj, num_vadj);\n");

	fprintf(outfp, "ret = METIS_PartGraphKway(&total_vertices, &total_num_phases, vadj, edges,"
			"vertex_wgt, NULL, edge_wgts, &parts, NULL, NULL, NULL, &objval, partitions);\n ");

	fprintf(outfp, "partition_print(total_vertices, partitions, objval, ret);\n");

	fprintf(outfp,"int count = 0");
	for(i=0;i<max_copy_level;i++)
		fprintf(outfp,", p%d", i);
	fprintf(outfp, ";\n\n");

	for(l=0;l<nloops;l++){
		for(i=0;i<non_scalar_copy_level[l];i++){
			fprintf(outfp,"for(p%d=0;p%d<MAX_TILES_PER_DIM;p%d++)\n", i, i, i);
		}
		fprintf(outfp,"pi_mappings_%d", l);
		for(i=0;i<non_scalar_copy_level[l];i++){
			fprintf(outfp, "[p%d]", i);
		}
		fprintf(outfp," = partitions[count++];\n\n");
	}

	//	fprintf(outfp,"polyrt_generalize_partition(pi_mappings, MAX_TILES_PER_DIM, nprocs);\n\n");

	for (i=0; i<compute_pi->nstmts; i++) {
		fprintf(outfp, "#undef S%d\n", i+1);
	}

	fprintf(outfp, "}\n\n");


	pluto_prog_free(compute_pi);

	return 0;
}

int pluto_dist_gen_intial_code(int nloops, Ploop **loops, int *copy_level,
		PlutoProg *prog, FILE *headerfp){

	int i, k, l;


	if(!options->data_dist)
		return 0;

	Stmt ***null_init_stmts, ****read_in_stmts, ****ref_count_initialize_stmts,
	****ref_count_writeout_stmts;

	int *num_write_data, *num_read_data, *num_read_write_data;

	pluto_dist_initalize_array_domains(prog, nloops);

	read_in_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	ref_count_initialize_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	ref_count_writeout_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));

	num_write_data = (int *) malloc(nloops*sizeof(int));
	num_read_data = (int *) malloc(nloops*sizeof(int));
	num_read_write_data = (int *) malloc(nloops*sizeof(int));

	int num_read_in_stmts = 1;
	int num_null_init_stmts = 1;
	int num_ref_count_init_stmts = 1;
	int num_ref_count_write_out_stmts = 1;
	int *num_stmts_per_wacc; // indexed by data variable
	int *num_stmts_per_racc; // indexed by data variable
	int *num_stmts_per_acc; // indexed by data variable

	struct stmt_access_pair ***wacc_stmts; // indexed by data variable
	struct stmt_access_pair ***racc_stmts; // indexed by data variable
	struct stmt_access_pair ***acc_stmts; // indexed by data variable

	int num_inital_stmts = 0;

	for(l=0; l<nloops; l++) {

		Ploop *loop = loops[l];


		acc_stmts = get_read_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_read_write_data[l], &num_stmts_per_acc);
		racc_stmts = get_read_access_with_stmts(loop->stmts,
				loop->nstmts, &num_read_data[l], &num_stmts_per_racc);
		wacc_stmts = get_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_write_data[l], &num_stmts_per_wacc);

		read_in_stmts[l] = (Stmt ***) malloc(num_read_data[l]*sizeof(Stmt **));
		ref_count_initialize_stmts[l] = (Stmt ***) malloc(num_read_write_data[l]*sizeof(Stmt **));
		ref_count_writeout_stmts[l] = (Stmt ***) malloc(num_write_data[l]*sizeof(Stmt **));

		for (i = 0; i < num_read_write_data[l]; ++i) {
			pluto_dist_update_arr_domain( acc_stmts[i],
					num_stmts_per_acc[i],prog,copy_level,l);
		}

		gen_data_tile_ref_count_init_cloog_code(prog, l, loop->stmts, loop->nstmts,
				copy_level[l], headerfp);

		for (i = 0; i < num_read_data[l]; ++i) {
			PlutoConstraints *read_in = get_read_in_constraints(racc_stmts[i],
					num_stmts_per_racc[i],copy_level[l],  l,  prog);

			gen_read_in_data_tile_alloc_cloog_code(racc_stmts[i], num_stmts_per_racc[i], read_in,
					prog, l,copy_level[l], headerfp);

			gen_read_in_copy_cloog_code(racc_stmts[i], num_stmts_per_racc[i], read_in,
					prog, l,copy_level[l], headerfp);

			read_in_stmts[l][i] =
					gen_read_in_code(read_in, racc_stmts[i],num_stmts_per_racc[i],prog,copy_level,l);
		}

		for (i = 0; i < num_read_write_data[l]; ++i) {
			ref_count_initialize_stmts[l][i] =
					gen_data_tile_ref_count_code(acc_stmts[i],num_stmts_per_acc[i],
							prog,copy_level,l);
		}

		for (i = 0; i < num_write_data[l]; ++i) {
			ref_count_writeout_stmts[l][i] =
					gen_write_out_tiles_ref_count_code(wacc_stmts[i],num_stmts_per_wacc[i],
							prog,copy_level,l);
		}

		for (i=0; i<num_write_data[l]; i++) {
			free(wacc_stmts[i]);
		}
		free(wacc_stmts);

		for (i=0; i<num_read_data[l]; i++) {
			free(racc_stmts[i]);
		}
		free(racc_stmts);

		for (i=0; i<num_read_write_data[l]; i++) {
			free(acc_stmts[i]);
		}
		free(acc_stmts);

	}

	null_init_stmts = (Stmt ***) malloc(prog->narrays*sizeof(Stmt **));
	for (i = 0; i < prog->narrays; ++i) {

		Array *arr = prog->arrays[i];
		sprintf(prog->decls+strlen(prog->decls), "%s",
				pluto_dist_gen_declarations(arr->text, prog));

		null_init_stmts[i] =
				gen_null_init_code(arr,prog, 0);
	}

	PlutoProg *data_tile_prog = NULL;
	data_tile_prog = pluto_prog_alloc();

	data_tile_prog->arrays = prog->arrays;
	data_tile_prog->narrays = prog->narrays;

	data_tile_prog_add_parm(data_tile_prog, prog, 0);

	for (i = 0; i < prog->narrays; ++i) {
		for (k = 0; k < num_null_init_stmts; ++k) {
			pluto_add_given_stmt(data_tile_prog, null_init_stmts[i][k]);
			num_inital_stmts++;
		}
	}

	for (l = 0; l < nloops; ++l) {
		for (i = 0; i < num_read_write_data[l]; ++i) {
			for (k = 0; k < num_ref_count_init_stmts; ++k) {
				pluto_add_given_stmt(data_tile_prog, ref_count_initialize_stmts[l][i][k]);
				num_inital_stmts++;
			}
		}

		for (i = 0; i < num_write_data[l]; ++i) {
			for (k = 0; k < num_ref_count_write_out_stmts; ++k) {
				pluto_add_given_stmt(data_tile_prog, ref_count_writeout_stmts[l][i][k]);
				num_inital_stmts++;
			}
		}

		for (i = 0; i < num_read_data[l]; ++i) {
			for (k = 0; k < num_read_in_stmts; ++k) {
				pluto_add_given_stmt(data_tile_prog, read_in_stmts[l][i][k]);
				num_inital_stmts++;
			}
		}
	}

	pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

	IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

	char func_name[512] = "data_tile_init";
	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_accs(prog->stmts, prog->nstmts, &num_data);
	gen_func_cloog_code(func_name, headerfp, 0, 0, num_data,
			accs, data_tile_prog, prog, 0);

	pluto_prog_free(data_tile_prog);

	free(accs);

	return 0;
}

void gen_copy_back_function(Stmt ****copy_back_stmts, int *num_write_data,
		int num_copy_back_stmts, int nloops, PlutoProg *prog, FILE *headerfp){

	if(!options->data_dist ) return;

	if(!options->verify_output) return;

	int i, k, l;
	PlutoProg *data_tile_prog = NULL;
	data_tile_prog = pluto_prog_alloc();

	data_tile_prog->arrays = prog->arrays;
	data_tile_prog->narrays = prog->narrays;

	data_tile_prog_add_parm(data_tile_prog, prog, 0);

	for(l=0;l<nloops;l++){
		for (i=0; i<num_write_data[l]; i++) {
			for (k=0; k<num_copy_back_stmts; k++) {
				pluto_add_given_stmt(data_tile_prog, copy_back_stmts[l][i][k]);
			}
		}
	}

	pluto_separate_stmts(data_tile_prog,data_tile_prog->stmts, data_tile_prog->nstmts, 0, 0);

	IF_DEBUG(pluto_prog_print(stdout, data_tile_prog));

	char func_name[512] = "copy_back_writeout";
	int num_data;
	PlutoAccess **accs;
	accs = pluto_get_waccs(prog->stmts, prog->nstmts, &num_data);
	gen_func_cloog_code(func_name, headerfp, 0, 0, num_data,
			accs, data_tile_prog, prog, 1);

	pluto_prog_free(data_tile_prog);

	free(accs);

	return ;
}

int pluto_sharedmem_data_dist_codegen(PlutoProg *prog,Ploop **loops, int nloops,
		int copy_level[nloops],  FILE *outfp)
{
	fprintf(outfp, "#include \"buffer_manager.h\"\n");
	fprintf(outfp, "%s", prog->decls);

	fprintf(outfp, "\tdouble t_data_init_start, t_data_init= 0.0 ;\n");
	fprintf(outfp, "\tdouble t_comp_start, t_comp= 0.0;\n");
	fprintf(outfp, "\tdouble t_copy_back_start, t_copy_back =  0.0;\n");

	int i;

	int num_data;
	PlutoAccess **waccs;

	fprintf(outfp, "\tIF_TIME(t_data_init_start = rtclock())\n");
	fprintf(outfp, "data_tile_init_0(0, 0");

	int num = 0;
	PlutoAccess **access = pluto_get_all_accs(prog, &num);
	waccs = pluto_get_all_waccs(prog, &num_data);

	if(options->variables_not_global){
		fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
		for(i=0;i<num;i++){
			fprintf(outfp,", %s ", access[i]->name);
		}
		fprintf(outfp, " __endif ");
	}
	for(i=0;i<num;i++){
		fprint_data_dist_parm_call(outfp, access[i]->name, prog);
	}
	fprintf(outfp, ");\n");
	fprintf(outfp, "\tIF_TIME(t_data_init += rtclock() - t_data_init_start)\n");
	fprintf(outfp, "\tIF_TIME(t_comp_start = rtclock())\n");
	free(access);

	int l, k;

//	PlutoProg *execute_task_prog = NULL;
//	execute_task_prog = pluto_prog_alloc();
//
//	execute_task_prog->arrays = prog->arrays;
//	execute_task_prog->narrays = prog->narrays;
//
//	data_tile_prog_add_parm(execute_task_prog, prog, 0);
//	execute_task_prog->hProps = prog->hProps;
//	execute_task_prog->deps = prog->deps;
//	execute_task_prog->ndeps = prog->ndeps;

	int count = 0;
	int *sep_stmts_size = (int *)malloc(nloops * sizeof(int));
	Stmt ***sep_stmts;
	sep_stmts = (Stmt ***) malloc(nloops*sizeof(Stmt **));

	for (l=0; l<nloops; l++) {
		count = 0;
		Stmt *anchor_stmt = loops[l]->stmts[0];
		// int dims[copy_level[l]];
		// int num_dims = 0;
		for (k=0; k<copy_level[l]; k++) {
			if (anchor_stmt->hyp_types[k] != H_SCALAR) {
				// dims[num_dims++] = k;
			}
		}

		sep_stmts_size[l] = 2 + loops[l]->nstmts;
		sep_stmts[l] = (Stmt **)malloc(sep_stmts_size[l]*sizeof(Stmt *));

		char *stmt_text = (char *)malloc(1024 * 10 * sizeof(char));
		stmt_text[0] = '\0';

		int num_rw_data;
		PlutoAccess **accs;
		accs = pluto_get_all_distinct_vars(loops[l]->stmts, loops[l]->nstmts, &num_rw_data);

		sprintf(stmt_text+strlen(stmt_text), "\t\t\t\tdata_tile_alloc_%d(", l);
		sprintf(stmt_text+strlen(stmt_text), "t1");
		for (i=1; i<copy_level[l]; i++) {
			sprintf(stmt_text+strlen(stmt_text), ",t%d", i+1);
		}
		if (options->variables_not_global && !options->data_dist) {
			for (i=0; i<num_rw_data; i++) {
				char *acc_name = accs[i]->name;
				sprintf(stmt_text+strlen(stmt_text), ",%s", acc_name);
			}
		}
		for (i=0; i<num_rw_data; i++) {
			print_data_dist_parm_call_from_access(stmt_text, accs[i], prog);
		}

		sprintf(stmt_text+strlen(stmt_text), ");");

		Stmt *alloc_data_tile = create_helper_stmt(anchor_stmt, copy_level[l],
				stmt_text,  ORIG, l);

		pluto_add_given_stmt(prog, alloc_data_tile);
		sep_stmts[l][count++]=alloc_data_tile;

		for(i=0;i<loops[l]->nstmts; i++)
			sep_stmts[l][count++] = loops[l]->stmts[i];
		/*
		sprintf(stmt_text+strlen(stmt_text), "\t\t\t\tcompute_task_%d(", l);
		sprintf(stmt_text+strlen(stmt_text), "t1");
		for (i=1; i<copy_level[l]; i++) {
			sprintf(stmt_text+strlen(stmt_text), ",t%d", i+1);
		}
		if (options->variables_not_global && !options->data_dist) {
			for (i=0; i<num_rw_data; i++) {
				char *acc_name = accs[i]->name;
				sprintf(stmt_text+strlen(stmt_text), ",%s", acc_name);
			}
		}
		for (i=0; i<num_rw_data; i++) {
			print_data_dist_parm_call_from_access(stmt_text, accs[i], prog);
		}
		sprintf(stmt_text+strlen(stmt_text), ");");

		*/

		stmt_text[0]='\0';
		sprintf(stmt_text+strlen(stmt_text), "\t\t\t\tdata_tile_ref_count_update_%d(", l);
		sprintf(stmt_text+strlen(stmt_text), "t1");
		for (i=1; i<copy_level[l]; i++) {
			sprintf(stmt_text+strlen(stmt_text), ",t%d", i+1);
		}
		if (options->variables_not_global && !options->data_dist) {
			for (i=0; i<num_rw_data; i++) {
				char *acc_name = accs[i]->name;
				sprintf(stmt_text+strlen(stmt_text), ",%s", acc_name);
			}
		}

		for (i=0; i<num_rw_data; i++) {
			print_data_dist_parm_call_from_access(stmt_text, accs[i], prog);
		}

		sprintf(stmt_text+strlen(stmt_text), ");");

		Stmt *update_ref_count = create_helper_stmt(anchor_stmt, copy_level[l],
				stmt_text,  ORIG, l);

		pluto_add_given_stmt(prog, update_ref_count);
		sep_stmts[l][count++]=update_ref_count;

		pluto_separate_stmts(prog, sep_stmts[l],sep_stmts_size[l], copy_level[l], l);
	}


	char cloog_func_name[256] = "execute_tasks.cloog";
//
	FILE *cloogfp = fopen(cloog_func_name, "w+");
	pluto_gen_cloog_file(cloogfp, prog);
	rewind(cloogfp);
	pluto_mark_vec_stmts(cloogfp, prog);
//	// pluto_gen_cloog_code(compute_task, cloogfp, stdout);
//	// rewind(cloogfp);
	generate_declarations(prog, outfp);
	pluto_gen_cloog_code(prog, -1, -1, cloogfp, outfp);
	fclose(cloogfp);


	fprintf(outfp, "\tIF_TIME(t_comp = rtclock() - t_comp_start);\n");

	fprintf(outfp, "\tIF_TIME(t_copy_back_start = rtclock());\n");

	fprintf(outfp, "\tif (fopen(\".test\", \"r\")) {\n ");
	if(options->data_dist){
		fprintf(outfp, "copy_back_writeout_0(0,0");
		if(options->variables_not_global){
			fprintf(outfp," __ifndef USE_LOCAL_ARRAYS ");
			for(i=0;i<num_data;i++){
				fprintf(outfp,", %s ", waccs[i]->name);
			}
			fprintf(outfp," __endif ");
		}

		for (i=0; i<num_data; i++) {
			fprint_data_dist_parm_call(outfp, waccs[i]->name,prog);
		}
		fprintf(outfp, ");\n}\n");
	}

	free(waccs);

	fprintf(outfp, "\tIF_TIME(t_copy_back = rtclock() - t_copy_back_start);\n");

	waccs = pluto_get_all_accs(prog, &num_data);

	for (i = 0; i < num_data ; ++i) {
		fprintf(outfp, "\tfree_buffer_mang(buff_mang_%s);\n",waccs[i]->name);
	}

	fprintf(outfp, "\n\tIF_TIME(fprintf(stdout, \"Computation time: %%0.6lf s\\n\", t_comp));\n\n");
	fprintf(outfp, "\tIF_TIME(fprintf(stdout, \"Data init time: %%0.6lf s\\n\", t_data_init));\n\n");
	fprintf(outfp, "\tIF_TIME(fprintf(stdout, \"copy back time: %%0.6lf s\\n\", t_copy_back));\n\n");

	free(waccs);

    return 0;
}

int pluto_shared_memory_data_dist(PlutoProg *prog, FILE *headerfp, FILE *outfp){

	int nloops=0;
	int l;

	// Stmt ****ref_count_update_stmts;


	Ploop **loops = pluto_get_dom_parallel_loops(prog, &nloops);
	/* Loops should be in the increasing order of their depths */
	qsort(loops, nloops, sizeof(Ploop *), pluto_loop_compar);

	// copy_back_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	// data_alloc_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	// ref_count_update_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));

	int *copy_level;

	copy_level = (int *) malloc(nloops*sizeof(int));
	// num_write_data = (int *) malloc(nloops*sizeof(int));
	// num_read_data = (int *) malloc(nloops*sizeof(int));
	// num_read_write_data = (int *) malloc(nloops*sizeof(int));

	init_packfile();

	init_copy_level(prog, loops, nloops, copy_level, NULL);

	pluto_detect_hyperplane_types(prog);

	pluto_dist_gen_intial_code(nloops, loops, copy_level, prog, headerfp);

	for (l=0; l<nloops; l++) {
		Ploop *loop = loops[l];

		gen_data_tile_alloc_cloog_code(prog, l, loop->stmts, loop->nstmts, copy_level[l], NULL, headerfp);
		gen_data_tile_ref_count_update_cloog_code(prog, l, loop->stmts, loop->nstmts, copy_level[l], NULL, headerfp);
//		gen_compute_task_cloog_code(prog, l, loop->stmts, loop->nstmts, copy_level[l], NULL, headerfp);
	}

	pluto_dist_generate_copy_back_code(loops, nloops,copy_level,
			headerfp, prog);

	pluto_sharedmem_data_dist_codegen(prog, loops, nloops, copy_level, outfp);

//	for (l=0; l<nloops; l++) {
//		Ploop *loop = loops[l];
//
//		for (i=0; i<loop->nstmts; i++) {
//			loop->stmts[i]->ploop_id = l;
//		}
//
//		int *num_stmts_per_wacc; // indexed by data variable
//		int *num_stmts_per_racc; // indexed by data variable
//		int *num_stmts_per_acc; // indexed by data variable
//
//		struct stmt_access_pair ***wacc_stmts; // indexed by data variable
//		struct stmt_access_pair ***racc_stmts; // indexed by data variable
//		struct stmt_access_pair ***acc_stmts; // indexed by data variable
//
//		acc_stmts = get_read_write_access_with_stmts(loop->stmts,
//				loop->nstmts, &num_read_write_data[l], &num_stmts_per_acc);
//		racc_stmts = get_read_access_with_stmts(loop->stmts,
//				loop->nstmts, &num_read_data[l], &num_stmts_per_racc);
//		wacc_stmts = get_write_access_with_stmts(loop->stmts,
//				loop->nstmts, &num_write_data[l], &num_stmts_per_wacc);
//
//		copy_back_stmts[l] = (Stmt ***) malloc(num_write_data[l]*sizeof(Stmt **));
//
//		data_alloc_stmts[l] = (Stmt ***) malloc(num_read_write_data[l]*sizeof(Stmt **));
//		ref_count_update_stmts[l] = (Stmt ***) malloc(num_read_write_data[l]*sizeof(Stmt **));
//
//		for (i = 0; i < num_read_write_data[l]; ++i) {
//			data_alloc_stmts[l][i] =
//					gen_data_tile_alloc_code(acc_stmts[i],num_stmts_per_acc[i],prog,copy_level,l);
//		}
//
//		for (i = 0; i < num_read_write_data[l]; ++i) {
//			ref_count_update_stmts[l][i] =
//					gen_data_tile_ref_count_update_code(acc_stmts[i],num_stmts_per_acc[i],prog,copy_level,l);
//		}
//
//		for (i=0; i<num_write_data[l]; i++)    {
//			copy_back_stmts[l][i] =
//					gen_copy_back_code(wacc_stmts[i],num_stmts_per_wacc[i],prog,copy_level,l);
//
//			gen_data_tile_copy_back_cloog_code(prog, l, wacc_stmts[i],
//					num_stmts_per_wacc[i], copy_level[l], headerfp);
//		}
//
//
//		for (i=0; i<num_write_data[l]; i++) {
//			free(wacc_stmts[i]);
//		}
//		free(wacc_stmts);
//		free(num_stmts_per_wacc);
//
//		for (i=0; i<num_read_data[l]; i++) {
//			free(racc_stmts[i]);
//		}
//		free(racc_stmts);
//		free(num_stmts_per_racc);
//
//		for (i=0; i<num_read_write_data[l]; i++) {
//			free(acc_stmts[i]);
//		}
//		free(acc_stmts);
//		free(num_stmts_per_acc);
//
//	}
//
//	int num_copy_back_stmts = 1;
//	int num_ref_count_update_stmts = 1;
//	int num_alloc_stmts = 1;
//	int *sep_stmts_size = (int *)malloc(nloops * sizeof(int));
//
//	gen_copy_back_function(copy_back_stmts, num_write_data, num_copy_back_stmts,
//			nloops, prog,  headerfp);
//
//	for (l=0; l<nloops; l++) {
//		Ploop *loop = loops[l];
//
//		sep_stmts_size[l] = num_alloc_stmts*num_read_write_data[l]+
//				loop->nstmts+
//				num_ref_count_update_stmts*num_read_write_data[l];
//
//		sep_stmts[l] = (Stmt **)malloc(sep_stmts_size[l]*sizeof(Stmt *));
//		count = 0;
//
//		for (i = 0; i < num_read_write_data[l]; ++i) {
//			for (k = 0; k < num_alloc_stmts; ++k) {
//				pluto_add_given_stmt(prog, data_alloc_stmts[l][i][k]);
//				sep_stmts[l][count++] = data_alloc_stmts[l][i][k];
//			}
//		}
//
//		for (i = 0; i < loop->nstmts; ++i) {
//			sep_stmts[l][count++] = loop->stmts[i];
//		}
//
//		for (i = 0; i < num_read_write_data[l]; ++i) {
//			for (k = 0; k < num_ref_count_update_stmts; ++k) {
//				pluto_add_given_stmt(prog, ref_count_update_stmts[l][i][k]);
//				sep_stmts[l][count++] = ref_count_update_stmts[l][i][k];
//			}
//		}
//
//	}
//
//	for (l=0; l<nloops; l++) {
//		Ploop *loop = loops[l];
//
//		int level = outer_dist_loop_level[l];
//
//		/* Parallel loop is distributed (loop distribution) around these
//		 * statements - except for write-out statements */
//		pluto_separate_stmts(prog, sep_stmts[l],sep_stmts_size[l], level, 0);
//
//		Stmt *first_stmt = loop->stmts[0];
//		for (k = 1; k < loop->nstmts; ++k) {
//			Stmt *orig_stmt = loop->stmts[k];
//
//			orig_stmt->trans->val[level][orig_stmt->trans->ncols-1] =
//					first_stmt->trans->val[level][first_stmt->trans->ncols-1];
//		}
//
//	}

	free(copy_level);

	pluto_loops_free(loops,nloops);

	return (nloops==0);
}

void pluto_dist_generate_copy_back_code(Ploop **loops, int nloops, int *copy_level,
		FILE *headerfp, PlutoProg *prog){

	Stmt ****copy_back_stmts;

	int i, l;
	copy_back_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));

	int *num_write_data = (int *) malloc(nloops*sizeof(int));

	for (l=0; l<nloops; l++) {
		Ploop *loop = loops[l];
		int *num_stmts_per_wacc; // indexed by data variable
		struct stmt_access_pair ***wacc_stmts; // indexed by data variable

		wacc_stmts = get_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_write_data[l], &num_stmts_per_wacc);

		copy_back_stmts[l] = (Stmt ***) malloc(num_write_data[l]*sizeof(Stmt **));

		for (i=0; i<num_write_data[l]; i++)    {
			copy_back_stmts[l][i] =
					gen_copy_back_code(wacc_stmts[i],num_stmts_per_wacc[i],prog,copy_level,l);

			gen_data_tile_copy_back_cloog_code(prog, l, wacc_stmts[i],
					num_stmts_per_wacc[i], copy_level[l], headerfp);
		}

		for (i=0; i<num_write_data[l]; i++) {
			free(wacc_stmts[i]);
		}
	}

	int num_copy_back_stmts = 1;

	gen_copy_back_function(copy_back_stmts, num_write_data, num_copy_back_stmts,
			nloops, prog,  headerfp);

}

/* Distributed memory parallelization */
int pluto_distmem_parallelize(PlutoProg *prog, FILE *sigmafp, FILE *headerfp, FILE *pifp)
{
	int i, j, l, k;
	int nstmts;
	Stmt ****copy_comm_stmts, ****write_out_stmts , ****data_alloc_stmts, 
    ****ref_count_update_stmts, ****copy_back_stmts;

	Stmt ***sep_stmts ;
	int *pi_mappings, *num_write_data, *num_read_data, *num_read_write_data, *copy_level;
	Stmt ***sep_writeout_stmts;
	int  *num_data, *outer_dist_loop_level;

	nstmts = prog->nstmts;

	// print_hyperplane_properties(prog);

	FILE *pidefs = NULL;
	pidefs = fopen("pi_defs.h", "w");
	fclose(pidefs);

	/* Create a new pi.c file */

	fprintf(pifp, "\
#include <math.h>\n\
#include <stdio.h>\n\
#include <assert.h>\n\
#include \"pi_defs.h\"\n\
#include \"polyrt.h\"\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");

	fprintf(pifp, "#ifndef __BLOCK_CYCLIC_BLOCK_SIZE\n");
	fprintf(pifp, "#define __BLOCK_CYCLIC_BLOCK_SIZE 1\n");
	fprintf(pifp, "#endif\n\n");

	if (options->fop_unicast_runtime) {
		fprintf(sigmafp, "#ifndef __FOP_UNICAST_RECV_LIMIT\n");
		fprintf(sigmafp, "#define __FOP_UNICAST_RECV_LIMIT 1\n");
		fprintf(sigmafp, "#endif\n\n");
	}

	fprintf(sigmafp, "\n#ifdef __DYNSCHEDULER_DEBUG_PRINT\n");
	fprintf(sigmafp, "#include <stdio.h>\n");
	fprintf(sigmafp, "#endif\n\n");

	fprintf(sigmafp, "\
#include <math.h>\n\
#include <assert.h>\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");

	pluto_add_distmem_decls(prog, headerfp);


    int nbands = 0;
    int *comm_place_levels;
    pluto_get_dom_parallel_bands(prog, &nbands, &comm_place_levels);

	int nloops=0;
	Ploop **loops = pluto_get_dom_parallel_loops(prog, &nloops);
	/* Loops should be in the increasing order of their depths */
	qsort(loops, nloops, sizeof(Ploop *), pluto_loop_compar);

	IF_DEBUG(printf("distmem: parallelizing loops\n"););
	IF_DEBUG(pluto_loops_print(loops,nloops););

	pi_mappings = malloc(nstmts*sizeof(int));
	init_pi_mappings(prog, loops, nloops, pi_mappings);

	/* Stripmine parallel loop to achieve block-cyclic partitioning */
	if (options->blockcyclic) {
		for (l=0; l<nloops; l++) {
			int tsizes[1];
			Ploop *bloop = pluto_loop_dup(loops[l]);
			bloop->depth += l;
			Band band = {bloop, 1};
			tsizes[0]=options->cyclesize;
			pluto_tile_band(prog, &band, tsizes);
		}
	}

	copy_comm_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	write_out_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	copy_back_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	data_alloc_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	ref_count_update_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
	sep_stmts = (Stmt ***) malloc(nloops*sizeof(Stmt **));
	sep_writeout_stmts = (Stmt ***) malloc(nloops*sizeof(Stmt **));

	num_write_data = (int *) malloc(nloops*sizeof(int));
	num_read_data = (int *) malloc(nloops*sizeof(int));
	num_data = (int *) malloc(nloops*sizeof(int));
	num_read_write_data = (int *) malloc(nloops*sizeof(int));
	copy_level = (int *) malloc(nloops*sizeof(int));
	outer_dist_loop_level = (int *) malloc(nloops*sizeof(int));

	int num_comm_stmts = 0;

	init_packfile();

	init_copy_level(prog, loops, nloops, copy_level, outer_dist_loop_level);

	if (options->dynschedule) {
		int data_dist = options->data_dist;
		options->data_dist = 0;
		pluto_dynschedule_common_parallelize(prog, sigmafp, headerfp, loops, nloops, pi_mappings, copy_level);
		options->data_dist = data_dist;
	}

	pluto_detect_hyperplane_types(prog);


	pluto_dist_gen_intial_code(nloops, loops, copy_level, prog, headerfp);

	pluto_dist_gen_compute_pi(nloops, loops, copy_level,
			prog, pi_mappings, pifp, headerfp);


	for (l=0; l<nloops; l++) {
		Ploop *loop = loops[l];

		for (i=0; i<loop->nstmts; i++) {
			loop->stmts[i]->ploop_id = l;
		}

		int *num_stmts_per_wacc; // indexed by data variable
		int *num_stmts_per_racc; // indexed by data variable
		int *num_stmts_per_acc; // indexed by data variable

		struct stmt_access_pair ***wacc_stmts; // indexed by data variable
		struct stmt_access_pair ***racc_stmts; // indexed by data variable
		struct stmt_access_pair ***acc_stmts; // indexed by data variable

		acc_stmts = get_read_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_read_write_data[l], &num_stmts_per_acc);
		racc_stmts = get_read_access_with_stmts(loop->stmts,
				loop->nstmts, &num_read_data[l], &num_stmts_per_racc);
		wacc_stmts = get_write_access_with_stmts(loop->stmts,
				loop->nstmts, &num_write_data[l], &num_stmts_per_wacc);

		num_data[l] = num_write_data[l];

		if(options->data_dist){
			copy_comm_stmts[l] = (Stmt ***) malloc(num_read_write_data[l]*sizeof(Stmt **));
		} else {

			copy_comm_stmts[l] = (Stmt ***) malloc(num_write_data[l]*sizeof(Stmt **));
		}

		write_out_stmts[l] = (Stmt ***) malloc(num_write_data[l]*sizeof(Stmt **));
		copy_back_stmts[l] = (Stmt ***) malloc(num_write_data[l]*sizeof(Stmt **));

		Stmt *anchor_stmt = get_new_anchor_stmt(loop->stmts, loop->nstmts);

		if (options->data_dist && !options->dynschedule) {

			data_alloc_stmts[l] = (Stmt ***) malloc(num_read_write_data[l]*sizeof(Stmt **));
			ref_count_update_stmts[l] = (Stmt ***) malloc(num_read_write_data[l]*sizeof(Stmt **));

			for (i = 0; i < num_read_write_data[l]; ++i) {
				data_alloc_stmts[l][i] =
						gen_data_tile_alloc_code(acc_stmts[i],num_stmts_per_acc[i],prog,copy_level,l);
			}

			for (i = 0; i < num_read_write_data[l]; ++i) {
				ref_count_update_stmts[l][i] =
						gen_data_tile_ref_count_update_code(acc_stmts[i],num_stmts_per_acc[i],prog,copy_level,l);
			}

		}

		/* if data is distributed then comm needs to be done for RAR and RAW deps
		 * hence we need to iterate over both read and write access
		 */
		if (options->data_dist) {
			for (i=0; i<num_read_write_data[l]; i++)    {

				generate_sigma(acc_stmts[i], num_stmts_per_acc[i],
						copy_level, prog, l, pi_mappings, sigmafp, headerfp);

				if (options->commopt_fop) {
					copy_comm_stmts[l][i] =
							gen_comm_code_opt_fop(i,acc_stmts[i],num_stmts_per_acc[i],nloops, num_read_write_data[l],prog,
									anchor_stmt, copy_level,outer_dist_loop_level[l],l,
									pi_mappings,&num_comm_stmts, sigmafp, headerfp);

				}else if (options->commopt_foifi) {
					copy_comm_stmts[l][i] =
							gen_comm_code_opt_foifi(i,acc_stmts[i],num_stmts_per_acc[i],nloops, num_read_write_data[l], prog,
									anchor_stmt, copy_level,outer_dist_loop_level[l],l,pi_mappings,&num_comm_stmts, headerfp);
				}else {
					copy_comm_stmts[l][i] =
							gen_comm_code_opt(i,acc_stmts[i],num_stmts_per_acc[i], nloops, num_read_write_data[l], prog,
									anchor_stmt, copy_level,outer_dist_loop_level[l],l,pi_mappings,&num_comm_stmts, headerfp);
				}
			}
			if(options->verify_output){
				for (i=0; i<num_write_data[l]; i++)    {
					copy_back_stmts[l][i] =
							gen_copy_back_code(wacc_stmts[i],num_stmts_per_wacc[i],prog,copy_level,l);

					gen_data_tile_copy_back_cloog_code(prog, l, wacc_stmts[i],
							num_stmts_per_wacc[i], copy_level[l], headerfp);
				}
				for (i=0; i<num_write_data[l]; i++)    {
					write_out_stmts[l][i] =
							gen_write_out_code(wacc_stmts[i],num_stmts_per_wacc[i],prog,
									anchor_stmt, copy_level,outer_dist_loop_level[l],l, headerfp);
				}
			}
		}else {
			for (i=0; i<num_write_data[l]; i++)    {
				generate_sigma(wacc_stmts[i], num_stmts_per_wacc[i],
						copy_level, prog, l, pi_mappings, sigmafp, headerfp);
				if (options->commopt_fop) {
					copy_comm_stmts[l][i] =
							gen_comm_code_opt_fop(i, wacc_stmts[i],num_stmts_per_wacc[i],nloops,num_data[l],prog,
									anchor_stmt,copy_level,outer_dist_loop_level[l],l,
									pi_mappings,&num_comm_stmts, sigmafp, headerfp);
				}else if (options->commopt_foifi) {
					copy_comm_stmts[l][i] =
							gen_comm_code_opt_foifi(i, wacc_stmts[i],num_stmts_per_wacc[i],nloops,num_data[l],prog,
									anchor_stmt,copy_level,outer_dist_loop_level[l],l,
									pi_mappings,&num_comm_stmts, headerfp);
				}else {
					copy_comm_stmts[l][i] =
							gen_comm_code_opt(i, wacc_stmts[i],num_stmts_per_wacc[i],nloops,num_data[l],prog,
									anchor_stmt,copy_level,outer_dist_loop_level[l],l,
									pi_mappings,&num_comm_stmts, headerfp);
				}
				write_out_stmts[l][i] =
						gen_write_out_code(wacc_stmts[i],num_stmts_per_wacc[i],prog,
								anchor_stmt,copy_level,outer_dist_loop_level[l],l, headerfp);
			}
		}

		if (options->dynschedule) {

			if(options->data_dist){
				gen_pack_send_text_code(prog, copy_comm_stmts[l], acc_stmts,
						l, num_read_write_data[l], num_comm_stmts, copy_level[l], headerfp, loop);
				gen_unpack_text_code(prog, copy_comm_stmts[l], acc_stmts,
						l, num_read_write_data[l], num_comm_stmts, copy_level[l], headerfp);
			}
			else {
				gen_pack_send_text_code(prog, copy_comm_stmts[l], wacc_stmts,
						l, num_data[l], num_comm_stmts, copy_level[l], headerfp, loop);
				gen_unpack_text_code(prog, copy_comm_stmts[l], wacc_stmts,
						l, num_data[l], num_comm_stmts, copy_level[l], headerfp);
			}
		}


		for (i=0; i<num_write_data[l]; i++) {
			free(wacc_stmts[i]);
		}
		free(wacc_stmts);
		free(num_stmts_per_wacc);

		for (i=0; i<num_read_data[l]; i++) {
			free(racc_stmts[i]);
		}
		free(racc_stmts);
		free(num_stmts_per_racc);

		for (i=0; i<num_read_write_data[l]; i++) {
			free(acc_stmts[i]);
		}
		free(acc_stmts);
		free(num_stmts_per_acc);
	}

	int num_copy_back_stmts = 1;
	int num_ref_count_update_stmts = 1;
	int num_alloc_stmts = 1;
	int count = 0;
	int *sep_stmts_size = (int *)malloc(nloops * sizeof(int));
	int *sep_writeout_stmts_size = (int *)malloc(nloops * sizeof(int));


	PlutoProg *write_out_prog = pluto_prog_alloc();
	for (i=0; i<prog->npar; i++) {
		pluto_prog_add_param(write_out_prog, prog->params[i], write_out_prog->npar);
	}
	if (options->verify_output) {

		gen_copy_back_function(copy_back_stmts, num_write_data, num_copy_back_stmts,
				nloops, prog,  headerfp);

		for (l=0; l<nloops; l++) {
			if(options->data_dist){
				//				sep_writeout_stmts_size[l] = (NUM_WRITE_OUT_STMTS + num_copy_back_stmts)*
				//						num_write_data[l];
				sep_writeout_stmts_size[l] = (NUM_WRITE_OUT_STMTS )* num_data[l];
				sep_writeout_stmts[l] = (Stmt **)malloc(sep_writeout_stmts_size[l] *sizeof(Stmt *));

				//				for (i=0; i<num_write_data[l]; i++) {
				//					for (k=0; k<num_copy_back_stmts; k++) {
				//						pluto_add_given_stmt(write_out_prog, copy_back_stmts[l][i][k]);
				//						sep_writeout_stmts[l][count++] = copy_back_stmts[l][i][k];
				//					}
				//				}

				count = 0;
				for (i=0; i<num_write_data[l]; i++) {
					for (j=0; j<NUM_WRITE_OUT_STMTS; j++) {
						pluto_add_given_stmt(write_out_prog, write_out_stmts[l][i][j]);
						sep_writeout_stmts[l][count++] = write_out_stmts[l][i][j];
					}
				}
			}
			else {
				sep_writeout_stmts_size[l] = (NUM_WRITE_OUT_STMTS )* num_data[l];
				sep_writeout_stmts[l] = (Stmt **)malloc(sep_writeout_stmts_size[l] *sizeof(Stmt *));
				count = 0;
				for (i=0; i<num_data[l]; i++) {
					for (j=0; j<NUM_WRITE_OUT_STMTS; j++) {
						pluto_add_given_stmt(write_out_prog, write_out_stmts[l][i][j]);
						sep_writeout_stmts[l][count++] = write_out_stmts[l][i][j];
					}
				}
			}
		}
	}

	if (!options->dynschedule) {
		for (l=0; l<nloops; l++) {
			Ploop *loop = loops[l];
			IF_DEBUG(printf("separating for loop: \n"););
			IF_DEBUG(pluto_loop_print(loop););

			if(options->data_dist){
				//num_data_alloc_stmts = num_alloc_stmts*num_read_write_data[l];

				/* Add statements to PlutoProg */
				sep_stmts_size[l] = num_alloc_stmts*num_read_write_data[l]+
						loop->nstmts+
						num_comm_stmts*num_read_write_data[l] +
						num_ref_count_update_stmts*num_read_write_data[l];

				sep_stmts[l] = (Stmt **)malloc(sep_stmts_size[l]*sizeof(Stmt *));

			}
			else {
				sep_stmts_size[l] = num_comm_stmts*num_write_data[l];
				sep_stmts[l] = (Stmt **)malloc(sep_stmts_size[l]*sizeof(Stmt *));
			}

			count = 0;

			if(options->data_dist){

				for (i = 0; i < num_read_write_data[l]; ++i) {
					for (k = 0; k < num_alloc_stmts; ++k) {
						pluto_add_given_stmt(prog, data_alloc_stmts[l][i][k]);
						sep_stmts[l][count++] = data_alloc_stmts[l][i][k];
					}
				}

				for (i = 0; i < loop->nstmts; ++i) {
					sep_stmts[l][count++] = loop->stmts[i];
				}

				for (i = 0; i < num_read_write_data[l]; ++i) {
					for (k = 0; k < num_ref_count_update_stmts; ++k) {
						pluto_add_given_stmt(prog, ref_count_update_stmts[l][i][k]);
						sep_stmts[l][count++] = ref_count_update_stmts[l][i][k];
					}
				}
				for (i=0; i<num_read_write_data[l]; i++)    {
					for (k=0; k<num_comm_stmts; k++) {
						pluto_add_given_stmt(prog, copy_comm_stmts[l][i][k]);
						sep_stmts[l][count++] = copy_comm_stmts[l][i][k];
					}
				}
			}
			else {
				for (i=0; i<num_write_data[l]; i++)    {
					for (k=0; k<num_comm_stmts; k++) {
						pluto_add_given_stmt(prog, copy_comm_stmts[l][i][k]);
						sep_stmts[l][count++] = copy_comm_stmts[l][i][k];
					}
				}
			}

		}
	}

	int scalar_dimensions_added[nloops*2];
	int num_scalar_dimensions_added = 0;
	int level = -1;
	int dist_parallel_loop;

	for (l=0; l<nloops; l++) {
		Ploop *loop = loops[l];
		// depth of the loop before any scalar dimensions were added
		if (options->blockcyclic) {
			dist_parallel_loop = loop->depth + l;
		}else{
			dist_parallel_loop = loop->depth;
		}

		// maintain the scalar dimensions added in sorted order
		for (i=0; i<num_scalar_dimensions_added; i++) {
			int orig_dimension = scalar_dimensions_added[i]-(i+1);
			// compare original dimensions before any scalar dimensinos were added
			if (orig_dimension>=outer_dist_loop_level[l]) {
				// update the depth of the loops after scalar dimensions were added by the previous loops
				// i represents the number of dimensions added before original depth of the loop
				outer_dist_loop_level[l] += i;
				for (j=num_scalar_dimensions_added; j>i; j--) {
					scalar_dimensions_added[j] = scalar_dimensions_added[j-1] + 1;
				}
				scalar_dimensions_added[i++] = outer_dist_loop_level[l];
				num_scalar_dimensions_added++;
				break;
			}
		}
		if (i==num_scalar_dimensions_added) {
			// update the depth of the loops after scalar dimensions were added by the previous loops
			// i represents the number of dimensions added before original depth of the loop
			outer_dist_loop_level[l] += i;
			scalar_dimensions_added[i++] = outer_dist_loop_level[l];
			num_scalar_dimensions_added++;
		}

		if(options->verify_output)
			pluto_separate_stmts(write_out_prog, sep_writeout_stmts[l], NUM_WRITE_OUT_STMTS*num_data[l], outer_dist_loop_level[l], 0);


		if (!options->dynschedule) {

            /* Parallel loop is distributed (loop distribution) around these
             * statements - except for write-out statements */
            //                 pluto_separate_stmts(prog, sep_stmts[l], num_comm_stmts*num_data[l], outer_dist_loop_level[l], 0);

            level = outer_dist_loop_level[l];

            /* Parallel loop is distributed (loop distribution) around these
             * statements - except for write-out statements */
            pluto_separate_stmts(prog, sep_stmts[l],sep_stmts_size[l], level, l);

            if(options->data_dist){
                /* orignal stmts should not be seperated
                 * assumes that orig stmts are at the beginning of prog
                 */
                Stmt *first_stmt = loop->stmts[0];
                for (k = 1; k < loop->nstmts; ++k) {
                    Stmt *orig_stmt = loop->stmts[k];

                    orig_stmt->trans->val[level][orig_stmt->trans->ncols-1] =
                            first_stmt->trans->val[level][first_stmt->trans->ncols-1];
                }
            }

            /* Add scalar dimension to separate the fused statements out inside (immediately inside); since they are
             * to be fused for dist_parallel_loop - used only for write-out pack and unpack */
            //       		   pluto_separate_stmts(prog, NULL, 0, level+2, 0);

            /* dist_parallel_loop's position now changes */
            dist_parallel_loop++;

            FILE *outfp = fopen(".distmem", "w");
            if (outfp)  {
                fprintf(outfp, "t%d", level+1);
                fclose(outfp);
            }

		}

	}


	// recompute dependence satisfaction levels
	// after scalar dimensions have been added to the program:
	// should be done before generate_pi since that may require it

	if(!options->dynschedule)
		pluto_compute_dep_satisfaction_complex(prog);

	int inner_dist_loop_level;
	for (l=0; l<nloops; l++) {
		// depth of the loop before any scalar dimensions were added
		inner_dist_loop_level = copy_level[l] - 1;
		if (options->blockcyclic) {
			inner_dist_loop_level -= l;
			dist_parallel_loop = loops[l]->depth + l;
		}else{
			dist_parallel_loop = loops[l]->depth;
		}

		// find the depth of the loop after scalar dimensions were added
		int num_scalar_dimensions_added_before = 0;
		// find the depth of the loop after scalar dimensions were added
		for (i=0; i<num_scalar_dimensions_added; i++) {
			int orig_dimension = scalar_dimensions_added[i]-(i+1);
			// compare original dimensions before any scalar dimensinos were added
			if (orig_dimension>=inner_dist_loop_level) {
				num_scalar_dimensions_added_before = i;
				inner_dist_loop_level += i;
				break;
			}
		}
		if (i==num_scalar_dimensions_added) {
			num_scalar_dimensions_added_before = i;
			inner_dist_loop_level += i;
		}

		/* Generate pi */
		generate_pi(pifp, headerfp, outer_dist_loop_level[l], inner_dist_loop_level, prog, loops[l], l,
				scalar_dimensions_added, num_scalar_dimensions_added_before);
	}

	if (options->verify_output) {
		gen_write_out_cloog_code(prog, write_out_prog, headerfp);
	}

	for (l=0; l<nloops; l++) {
		for (i=0; i<num_data[l]; i++)    {
			free(copy_comm_stmts[l][i]);
		}
		free(copy_comm_stmts[l]);

		if (options->verify_output) {
			for (i=0; i<num_data[l]; i++)    {
				free(write_out_stmts[l][i]);
			}
			free(write_out_stmts[l]);
			free(sep_writeout_stmts[l]);
		}

		if (!options->dynschedule) free(sep_stmts[l]);
	}

	free(copy_comm_stmts);
	free(sep_stmts);

	if(options->verify_output){
		free(write_out_stmts);
		free(sep_writeout_stmts);
	}

	free(num_data);
	free(copy_level);
	free(outer_dist_loop_level);
	free(pi_mappings);

	pluto_mark_statements(prog);

	pluto_loops_free(loops,nloops);

	return (nloops==0);

}


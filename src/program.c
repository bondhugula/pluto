/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution. 
 *
 * program.c
 *
 * This file contains functions that do the job interfacing the PLUTO 
 * core to the frontend and related matters
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "pluto.h"
#include "math_support.h"
#include "program.h"

#include "osl/macros.h"
#include "osl/scop.h"
#include "osl/body.h"
#include "osl/relation_list.h"
#include "osl/extensions/arrays.h"
#include "osl/extensions/dependence.h"
#include "osl/extensions/loop.h"
#include "osl/body.h"
#include "osl/extensions/scatnames.h"

#include "cloog/cloog.h"

#include "candl/candl.h"
#include "candl/scop.h"
#include "candl/options.h"
#include "candl/dependence.h"

#include <isl/map.h>
#include <isl/mat.h>
#include <isl/set.h>
#include <isl/flow.h>
#include <isl/union_map.h>

void pluto_add_dep(PlutoProg *prog, Dep *dep)
{
    dep->id = prog->ndeps;
    prog->ndeps++;
    prog->deps = (Dep **) realloc(prog->deps, sizeof(Dep *)*prog->ndeps);
    prog->deps[prog->ndeps-1] = dep;
}


/*
 * In an [eq -I A c] relation, rows can be ordered any way.
 * Returns the index for the row for the nth output dimension.
 */
int osl_relation_get_row_id_for_nth_dimension(osl_relation_p relation,
                                              int ndim){
  int nb_ndims_found = 0;
  int row_id = -1;
  int i = 0;

  if (relation == NULL)
    return OSL_UNDEFINED;
  
  if ((relation->nb_rows < ndim) || (0 > ndim)) {
    fprintf(stderr, "error: dimension out of bounds");
    exit(1);
  }

  nb_ndims_found = 0;
  for (i = 0; i < relation->nb_rows; i++) {
    if (!osl_int_zero(relation->precision, relation->m[i][ndim])) {
      nb_ndims_found ++;
      row_id = i;
    }
  }
  if (nb_ndims_found == 0) {
    fprintf(stderr, "error: specified dimension not found");
    exit(1);
  }
  if (nb_ndims_found > 1) {
    fprintf(stderr, "error: specified dimension occurs multiple times");
    exit(1);
  }

  return row_id;
}


/*
 * Converts a [eq A c] relation to [A c] Pluto constraints
 */
PlutoConstraints *osl_relation_to_pluto_constraints (osl_relation_p rln){
  
  int i, j = 0;
  PlutoConstraints *cst;

  if(rln==NULL)
    return NULL;

  if(rln->nb_local_dims){
    fprintf(stderr, "Cannot handle Local Dimensions in a relation.\n");
    exit(1);
  }

  cst = pluto_constraints_alloc(rln->nb_rows, rln->nb_columns-1);
  cst->nrows = rln->nb_rows;

  //copy matrix values
  for(i=0; i < rln->nb_rows; i++){
    cst->is_eq[i] = osl_int_zero(rln->precision, rln->m[i][0]);
    for(j=0; j < cst->ncols; j++){
      cst->val[i][j] = osl_int_get_si(rln->precision, rln->m[i][j+1]);
    }
  }

  return cst;
}

/*
 * Converts [A c] PLuto constraints to a [eq A c] domain relation
 */
osl_relation_p pluto_constraints_to_osl_domain(PlutoConstraints *cst, int npar){
  
  int i, j = 0;
  osl_relation_p rln;

  if(cst==NULL)
    return NULL;

  rln = osl_relation_pmalloc(PLUTO_OSL_PRECISION, cst->nrows, cst->ncols+1);

  //copy matrix values
  for(i=0; i < rln->nb_rows; i++){
    osl_int_set_si(rln->precision, &rln->m[i][0], cst->is_eq[i]?0:1);
    for(j=0; j < cst->ncols; j++){
      osl_int_set_si(rln->precision, &rln->m[i][j+1], cst->val[i][j]);
    }
  }

  rln->type = OSL_TYPE_DOMAIN;
  rln->nb_parameters = npar;
  rln->nb_output_dims = rln->nb_columns - rln->nb_parameters - 2;
  rln->nb_input_dims = 0;
  rln->nb_local_dims  = 0;

  return rln;
}


/*
 * Converts a [eq -I A c] osl access relation to [A c] pluto matrix
 * Note: a[c] and a, having two and one output dimensions respectively
 * in osl, are converted to a one-dimensional pluto matrix.
 */
PlutoMatrix *osl_access_relation_to_pluto_matrix(osl_relation_p smat)
{
    int i, j;

    PlutoMatrix *mat;

    if(smat==NULL)
      return NULL;

    if(smat->nb_local_dims){
      fprintf(stderr, "Cannot handle Local Dimensions in a relation.\n");
      exit(1);
    }

    int nrows = smat->nb_rows==1?smat->nb_rows: smat->nb_rows-1; //skp id line
    int ncols = smat->nb_columns - smat->nb_output_dims - 1; //-1: skip 1st col
    mat = pluto_matrix_alloc(nrows, ncols);

    if(smat->nb_rows==1){  //special case for scalars
          for (j=smat->nb_output_dims+1; j<smat->nb_columns; j++)  {
              mat->val[0][j-(smat->nb_output_dims+1)] = 0;
          }
    }
    else{
      //fill in the rest of the information
      for (i=1; i<smat->nb_rows; i++)  {
          int row = osl_relation_get_row_id_for_nth_dimension(smat, i+1);
          for (j=smat->nb_output_dims+1; j<smat->nb_columns; j++)  {
              mat->val[i-1][j-(smat->nb_output_dims+1)] = 
                                osl_int_get_si(smat->precision, smat->m[row][j]);
          }
      }
    }

    return mat;
}


/* Return the number of lines until the next non-zero element
 * in the first column of "access" or until the end of the matrix.
 */
int access_len(PlutoMatrix *access, int first)
{
    int i;

    for (i = first + 1; i < access->nrows; ++i)
        if (access->val[i][0]!=0)
            break;

    return i - first;
}

int is_array(int id, PlutoMatrix* pmat){
  int i=0;
  int j=0;
  int is_array=0;

  for(i=0; i< pmat->nrows; i++){
    if(pmat->val[i][0]==id){
      if(access_len(pmat,i) > 1)
        is_array = 1;
      else{
        for(j=1; j< pmat->ncols; j++)
          if(pmat->val[i][j] != 0)
            is_array = 1;
      }
    }
  }

  return is_array;
}


/*
 * Converts a [A c] pluto matrix to [eq -I A c] osl access relation
 * Note: a[c] and a, having two and one output dimensions respectively
 * in osl, are converted to a one-dimensional pluto matrix.
 */
osl_relation_p pluto_matrix_to_osl_access_relation(PlutoMatrix *pmat)
{
  int i=0;
  int j=0;

  if(pmat==NULL)
    return NULL;

  //check if it's a scalar
  int nrows = pmat->nrows + is_array(pmat->val[0][0], pmat);
  int ncols = 1 + nrows + pmat->ncols;

  osl_relation_p rl = osl_relation_malloc(nrows, ncols);
  //set the dims outside
  rl->nb_output_dims = nrows;
  //set the type outside

  //first row with array_id
  osl_int_set_si(rl->precision, &rl->m[0][1], -1);
  //osl_int_set_si(rl->precision, &rl->m[0][rl->nb_columns-1], pmat->val[0][0]);

  //rest of the rows
  for(i=1; i< rl->nb_rows; i++){
    for(j=0; j< rl->nb_columns; j++){
      if(j==i+1)
        osl_int_set_si(rl->precision, &rl->m[i][j], -1);
      else if(j<=rl->nb_output_dims)
        osl_int_set_si(rl->precision, &rl->m[i][j], 0);
      else if(j<rl->nb_columns)
        osl_int_set_si(rl->precision, &rl->m[i][j], pmat->val[i-1][j-rl->nb_output_dims-1]);
    }
  }

  return rl;
}



/*
 * Converts a [eq -I A c] osl scattering to [A c] pluto transformations
 */
PlutoMatrix *osl_scattering_to_pluto_trans(osl_relation_p smat)
{
    int i, j;
    PlutoMatrix *mat;

    if(!smat)
      return NULL;

    if(smat->nb_local_dims){
      fprintf(stderr, "Cannot handle Local Dimensions in a relation.\n");
      exit(1);
    }

    mat = pluto_matrix_alloc(smat->nb_rows, smat->nb_columns-smat->nb_output_dims-1);
    for (i=0; i<smat->nb_rows; i++)  {
        /* Only equalities in schedule expected */
        assert(osl_int_get_si(smat->precision, smat->m[i][0]) == 0);

        int row = osl_relation_get_row_id_for_nth_dimension(smat, i+1);
        for (j=smat->nb_output_dims+1; j<smat->nb_columns; j++)  {
            mat->val[i][j-smat->nb_output_dims-1] = osl_int_get_si(smat->precision, smat->m[row][j]);
        }
    }

    return mat;
}

/*
 * Converts a [A c] pluto transformations to a [eq -I A c] osl scattering
 */
osl_relation_p pluto_trans_to_osl_scattering(PlutoMatrix *mat, int npar)
{
    int i, j;
    osl_relation_p smat;

    if(!mat)
      return NULL;

    smat = osl_relation_pmalloc(PLUTO_OSL_PRECISION, mat->nrows, mat->nrows+mat->ncols+1);
    smat->type = OSL_TYPE_SCATTERING;
    smat->nb_parameters  = npar;
    smat->nb_output_dims = mat->nrows;
    smat->nb_input_dims  = mat->ncols - npar - 1;
    smat->nb_local_dims  = 0;

    for (i=0; i<smat->nb_rows; i++)  {
      for (j=1; j<smat->nb_columns; j++)  {

        /* Only equalities in schedule expected */
        if(j==0)  // eq/neq (first) column
          osl_int_set_si(smat->precision, &smat->m[i][j], 0);

        //fill out the output dims
        else if(j==i+1)
          osl_int_set_si(smat->precision, &smat->m[i][j], -1);
        else if(j<=smat->nb_output_dims) // non diagonal zeros
          osl_int_set_si(smat->precision, &smat->m[i][j], 0);

        //fill out the intput_dims+params+const
        else
          osl_int_set_si(smat->precision, &smat->m[i][j], mat->val[i][j-smat->nb_output_dims-1]);

      }
    }

    return smat;
}


/*
 * get a list of to-be-vectorized loops from PlutoProg
 */
osl_loop_p pluto_get_vector_loop_list( const PlutoProg *prog)
{
    int i, j, nploops;
    osl_loop_p ret_loop = NULL;

    Ploop **ploops = pluto_get_parallel_loops(prog, &nploops);

    for (i=0; i<nploops; i++) {
        /* Only the innermost ones */
        if (!pluto_is_loop_innermost(ploops[i], prog)) continue;

        IF_DEBUG(printf("[pluto_get_vector_loop_list] marking loop\n"););
        IF_DEBUG(pluto_loop_print(ploops[i]););

        osl_loop_p newloop = osl_loop_malloc();

        char iter[5];
        sprintf(iter, "t%d", ploops[i]->depth+1);
        newloop->iter =  strdup(iter);


        newloop->nb_stmts = ploops[i]->nstmts;
        newloop->stmt_ids = malloc(ploops[i]->nstmts*sizeof(int));
        for (j=0; j<ploops[i]->nstmts; j++) {
            newloop->stmt_ids[j] = ploops[i]->stmts[j]->id+1;
        }

        newloop->directive   += CLAST_PARALLEL_VEC;

        //add new loop to looplist
        osl_loop_add(newloop, &ret_loop);
    }

    pluto_loops_free(ploops, nploops);

    return ret_loop;
}


/*
 * get a list of to-be-parallelized loops frop PlutoProg
 */
osl_loop_p pluto_get_parallel_loop_list(const PlutoProg *prog)
{
    int i, j, nploops;
    osl_loop_p ret_loop = NULL;

    Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);

    IF_DEBUG(printf("[pluto_parallel_loop_list] parallelizable loops\n"););
    IF_DEBUG(pluto_loops_print(ploops, nploops););


    for (i=0; i<nploops; i++) {

        osl_loop_p newloop = osl_loop_malloc();

        char iter[5];
        sprintf(iter, "t%d", ploops[i]->depth+1);
        newloop->iter = strdup(iter);

        newloop->nb_stmts = ploops[i]->nstmts;
        newloop->stmt_ids = malloc(ploops[i]->nstmts*sizeof(int));
        for (j=0; j<ploops[i]->nstmts; j++) {
            newloop->stmt_ids[j] = ploops[i]->stmts[j]->id+1;
        }

        newloop->private_vars = strdup("lbv,ubv");
        newloop->directive   += CLAST_PARALLEL_OMP;

        //add new loop to looplist
        osl_loop_add(newloop, &ret_loop);
    }

    pluto_loops_free(ploops, nploops);

    return ret_loop;
}

/*
 * Replace the original scop's statements' domains and scatterings
 * by those generated by Pluto
 */
void pluto_populate_scop (osl_scop_p scop, PlutoProg *prog,
                           PlutoOptions *options){
  
    int i;
    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    int npar = prog->npar;

    osl_statement_p stm = scop->statement;
    /* Fill domains (may have been changed for tiling purposes). */
    for (i=0; i<nstmts; i++)    {

      //replace domain
      if(stm->domain)
        osl_relation_free(stm->domain);
      stm->domain = pluto_constraints_to_osl_domain(stmts[i]->domain, npar);

      //replace scattering 
      if(stm->scattering)
        osl_relation_free(stm->scattering);
      stm->scattering = pluto_trans_to_osl_scattering(stmts[i]->trans, npar);

      stm = stm->next;
    }

    //update iterators
    //if domains(iterators) chanaged due to optimizations (tiling, etc.)
    for (stm = scop->statement; stm; stm = stm->next)
    {
      int niter = stm->domain->nb_columns - scop->context->nb_columns;
      int nb_orig_it = -1;
      if(stm->body){
        osl_body_p stmt_body = (osl_body_p)(stm->body->data);
        nb_orig_it = osl_strings_size(stmt_body->iterators);
        if (nb_orig_it != niter)
          {//update iterators.

            char** iters =(char**)malloc(sizeof(char*)*(niter+1));//+1 for NULL
            for (i = 0; i < niter - nb_orig_it; ++i)
              {
                iters[i] = (char*) malloc( sizeof(char)*16 );
                sprintf (iters[i], "fk%d", i);

                //update accesses
                osl_relation_list_p rll = stm->access;
                while(rll){
                  osl_relation_insert_blank_column(rll->elt, rll->elt->nb_output_dims);
                  rll->elt->nb_input_dims++;
                  rll = rll->next;
                }

              }
            for (; i < niter; ++i)
              iters[i] = stmt_body->iterators->string[i - niter + nb_orig_it];

            iters[i]=(char*)NULL;
            
            free(stmt_body->iterators->string);
            stmt_body->iterators->string = iters;
          }
       }
    }

    //update scatnames
    //get max scat dims
    int nb_scatt = 0;
    for (stm = scop->statement; stm; stm = stm->next)
    {
        int cur_scatt = stm->scattering->nb_output_dims;
        nb_scatt = nb_scatt > cur_scatt ? nb_scatt : cur_scatt;
    }

    //generate scatt names
    osl_strings_p newnames = osl_strings_generate("t", nb_scatt);
    osl_scatnames_p scatt = osl_scatnames_malloc();
    scatt->names = newnames;
  
    //replace the old scatnames with new one
    osl_generic_remove(&scop->extension, OSL_URI_SCATNAMES);
    osl_generic_p gen = osl_generic_shell(scatt, osl_scatnames_interface());
    osl_generic_add(&scop->extension, gen);
    

    //update loop information
    //get loops to be marked for parallization and vectorization
    osl_loop_p pll = NULL;
    if(options->parallel){
      pll = pluto_get_parallel_loop_list(prog);
    }
    osl_loop_p vll = NULL;
    if (options->prevector ){
      vll = pluto_get_vector_loop_list(prog);
    }
    //concatenate the two lists
    osl_loop_add(vll, &pll);

    if(pll){
      osl_generic_p loopgen = osl_generic_shell(pll, osl_loop_interface());
      osl_generic_add(&scop->extension, loopgen);
    }

}

static int get_osl_write_access_position(osl_relation_list_p rl,
                                   osl_relation_p access)
{
    int num;

    num = -1;

    osl_relation_list_p tmp = rl;
    for (; tmp; tmp = tmp->next)  {

        if ( (tmp->elt->type == OSL_TYPE_WRITE) ||
             (tmp->elt->type == OSL_TYPE_MAY_WRITE) )
            num++;

        if(tmp->elt == access)
          break;
    }
    assert(num >= 0);
    return num;
}

static int get_osl_read_access_position(osl_relation_list_p rl,
                                   osl_relation_p access)
{
    int num;

    num = -1;

    osl_relation_list_p tmp = rl;
    for (; tmp; tmp = tmp->next)  {

        if ( (tmp->elt->type == OSL_TYPE_READ) )
            num++;

        if(tmp->elt == access)
          break;
    }
    assert(num >= 0);
    return num;
}

/*
 * Returns a list of write or may_write access relations in a list
 */
osl_relation_list_p osl_access_list_filter_write(osl_relation_list_p list) {

  osl_relation_list_p copy = osl_relation_list_clone(list);
  osl_relation_list_p filtered = NULL;
  osl_relation_list_p previous = NULL;
  osl_relation_list_p trash;
  int first = 1;

  while (copy != NULL) {
    if ((copy->elt != NULL) &&
        ( (copy->elt->type == OSL_TYPE_WRITE) ||
          (copy->elt->type == OSL_TYPE_MAY_WRITE))) {
      if (first) {
        filtered = copy;
        first = 0;
      }
      
      previous = copy;
      copy = copy->next;
    }
    else {
      trash = copy;
      if (!first)
        previous->next = copy->next;
      copy = copy->next;
      trash->next = NULL;
      osl_relation_list_free(trash);
    }
  }

  return filtered;
}

/*
 * Returns a list of read access relations in a list
 */
osl_relation_list_p osl_access_list_filter_read(osl_relation_list_p list) {

  osl_relation_list_p copy = osl_relation_list_clone(list);
  osl_relation_list_p filtered = NULL;
  osl_relation_list_p previous = NULL;
  osl_relation_list_p trash;
  int first = 1;

  while (copy != NULL) {
    if ((copy->elt != NULL) &&
        (copy->elt->type == OSL_TYPE_READ)) {
      if (first) {
        filtered = copy;
        first = 0;
      }
      
      previous = copy;
      copy = copy->next;
    }
    else {
      trash = copy;
      if (!first)
        previous->next = copy->next;
      copy = copy->next;
      trash->next = NULL;
      osl_relation_list_free(trash);
    }
  }

  return filtered;
}

/*
* Converts an osl dependence domain to Pluto constraints
* See osl/extensions/dependence.h for the osl dependence domain matrix format
*/
PlutoConstraints* osl_dep_domain_to_pluto_constraints(osl_dependence_p in_dep){

 int s_dom_output_dims = in_dep->source_nb_output_dims_domain;
 int t_dom_output_dims = in_dep->target_nb_output_dims_domain;


 int nb_output_dims = in_dep->source_nb_output_dims_domain +
                      in_dep->source_nb_output_dims_access;
 int nb_input_dims  = in_dep->target_nb_output_dims_domain +
                      in_dep->target_nb_output_dims_access;

 /* Compute osl domain indexes */
 int osl_ind_source_local_domain = 1 + nb_output_dims + nb_input_dims;
 int osl_ind_source_local_access = osl_ind_source_local_domain +
                                   in_dep->source_nb_local_dims_domain;
 int osl_ind_target_local_domain = osl_ind_source_local_access +
                                   in_dep->source_nb_local_dims_access;
 int osl_ind_target_local_access = osl_ind_target_local_domain +
                                   in_dep->target_nb_local_dims_domain;
 int osl_ind_params              = osl_ind_target_local_access +
                                   in_dep->target_nb_local_dims_access;

 /* Compute pluto constraints domain indexes */
 int pl_ind_target_domain      = 1 + in_dep->source_nb_output_dims_domain;
 int pl_ind_params             = pl_ind_target_domain +
                                   in_dep->target_nb_output_dims_domain;

 int rows, cols = 0;
 
 int nb_pars = in_dep->stmt_source_ptr->domain->nb_parameters;
 int s_dom_rows = in_dep->stmt_source_ptr->domain->nb_rows;
 int t_dom_rows = in_dep->stmt_target_ptr->domain->nb_rows;
 int s_acc_rows = in_dep->ref_source_access_ptr->nb_rows - 1;
 int depth = in_dep->depth;

 //
 rows = s_dom_rows+t_dom_rows+
      (s_acc_rows==0? 1: s_acc_rows)  //special case for 0-dimention array(scalar)
        +depth;
 cols = s_dom_output_dims+t_dom_output_dims+nb_pars+2; //cols: 2 => eq + const

 PlutoConstraints *cst;

 cst = pluto_constraints_alloc(rows, cols-1);
 cst->nrows = rows;
 cst->ncols = cols-1;

 int i=0;
 int j=0;
 int osl_constraint = 0;
 int pl_constraint = 0;
 int osl_index=0;
 int pl_index=0;


 // copy source domain
 osl_relation_p s_domain = in_dep->stmt_source_ptr->domain;
 for(i=0; i< s_domain->nb_rows; i++){

   //copy first column
   if (osl_int_zero(in_dep->domain->precision, in_dep->domain->m[osl_constraint][0])) {
     cst->is_eq[pl_constraint] = 1;
   }else{
     cst->is_eq[pl_constraint] = 0;
   }

   //start of matrix
   osl_index = 1; //start of src_stmt_domain_output_dims
   pl_index = 1-1; // -1 for pluto
   for (j=0;j<s_dom_output_dims; j++)
     cst->val[pl_constraint][pl_index+j] = osl_int_get_si(in_dep->domain->precision,
                                     in_dep->domain->m[osl_constraint][osl_index+j]);
     
   // copy localdims - not supprted by converter
   if(s_domain->nb_local_dims){
     fprintf(stderr, "local dimensions in domain not supported\n");
     exit(1);
   }

   // copy params + constant
   osl_index = osl_ind_params;
   pl_index = pl_ind_params-1;  // -1 for pluto
   for (j=0; j<nb_pars+1; j++)
     cst->val[pl_constraint][pl_index+j] = osl_int_get_si(in_dep->domain->precision,
                               in_dep->domain->m[osl_constraint][osl_index+j]);

   osl_constraint++;
   pl_constraint++;
 }


 // copy target domain
 osl_relation_p t_domain = in_dep->stmt_target_ptr->domain;
 for(i=0; i< t_domain->nb_rows; i++){

   //copy first column
   if (osl_int_zero(in_dep->domain->precision, in_dep->domain->m[osl_constraint][0])) {
     cst->is_eq[pl_constraint] = 1;
   }else{
     cst->is_eq[pl_constraint] = 0;
   }

   //start of matrix
   osl_index = 1 + nb_output_dims;
   pl_index = pl_ind_target_domain-1; // -1 for pluto
   for (j=0;j<t_dom_output_dims; j++)
     cst->val[pl_constraint][pl_index+j] = osl_int_get_si(in_dep->domain->precision,
                               in_dep->domain->m[osl_constraint][osl_index+j]);
     
   // copy local dims - not supported in converter
   if(t_domain->nb_local_dims){
     fprintf(stderr,"local dimensions in domain not supproted\n");
     exit(1);
   }

   // copy params + constant
   osl_index = osl_ind_params;
   pl_index = pl_ind_params-1; // -1 for pluto
   for (j=0; j<nb_pars+1; j++)
     cst->val[pl_constraint][pl_index+j] = osl_int_get_si(in_dep->domain->precision,
                               in_dep->domain->m[osl_constraint][osl_index+j]);

   pl_constraint++;
   osl_constraint++;
 }


 // copy source as well as target access
 int osl_s_index     = 0;
 int osl_t_index     = 0;
 int pl_s_index = 0;
 int pl_t_index = 0;

 osl_relation_p s_access = in_dep->ref_source_access_ptr;
 osl_relation_p t_access = in_dep->ref_target_access_ptr;

 osl_constraint++; //skip the array_id line

 for(i=0; i < s_acc_rows; i++){

   //copy first column
   if (osl_int_zero(in_dep->domain->precision, in_dep->domain->m[osl_constraint][0])) {
     cst->is_eq[pl_constraint] = 1;
   }else{
     cst->is_eq[pl_constraint] = 0;
   }

   osl_s_index     = 1;
   osl_t_index     = 1 + nb_output_dims;
   pl_s_index = 1-1; // -1 for pluto
   pl_t_index = pl_ind_target_domain-1; // -1 for pluto

   for (j=0; j<s_access->nb_input_dims; j++){
     cst->val[pl_constraint][pl_s_index+j] = osl_int_get_si(in_dep->domain->precision,
                              in_dep->domain->m[osl_constraint][osl_s_index+j]);
   }

   for (j=0; j<t_access->nb_input_dims; j++){ //t_acc_dims==s_acc_dims
     cst->val[pl_constraint][pl_t_index+j] = osl_int_get_si(in_dep->domain->precision,
          in_dep->domain->m[osl_constraint+s_access->nb_rows][osl_t_index+j]);
   }

   //copy local dimensions - not supported by converter
   if(s_access->nb_local_dims || t_access->nb_local_dims){
     fprintf(stderr, "local dimensions in Access not supproted\n");
     exit(1);
   }

   // copy params + constant
   osl_index = osl_ind_params;
   pl_index = pl_ind_params-1; // -1 for pluto
   for (j=0; j<nb_pars+1; j++){
     //get src params
     int src_param = osl_int_get_si(in_dep->domain->precision,
                              in_dep->domain->m[osl_constraint][osl_index+j]);
     //get tgt params
     int tgt_param = osl_int_get_si(in_dep->domain->precision,
             in_dep->domain->m[osl_constraint+s_access->nb_rows][osl_index+j]);

     tgt_param = -tgt_param; //oppose

     cst->val[pl_constraint][pl_index+j] = src_param - tgt_param;

   }

   pl_constraint++;
   osl_constraint++;
 }


 // copy access equalities
 // skip min_depth
 int min_depth = OSL_min(s_access->nb_output_dims, 
                         t_access->nb_output_dims);
 osl_constraint += s_access->nb_rows + min_depth; 
                   
 //s_acc_rows calculated by subtracting 1 from acc.nb_rows
 //in case of a scalar this results in 0, still add a constraint for pluto
 if(s_acc_rows==0) pl_constraint++;

 
 // copy depth
 osl_s_index = 1;
 osl_t_index = 1 + nb_output_dims;
 pl_s_index = 1-1; // -1 for pluto
 pl_t_index = pl_ind_target_domain-1; // -1 for pluto
 for(i=0; i< depth; i++){
   // copy first column
   if (osl_int_zero(in_dep->domain->precision, in_dep->domain->m[osl_constraint][0])) {
     cst->is_eq[pl_constraint] = 1;
   }else{
     cst->is_eq[pl_constraint] = 0;
   }

   // copy subscript equalities
   cst->val[pl_constraint][pl_s_index+i] = osl_int_get_si(in_dep->domain->precision,
                             in_dep->domain->m[osl_constraint][osl_s_index+i]);
   cst->val[pl_constraint][pl_t_index+i] = osl_int_get_si(in_dep->domain->precision,
                             in_dep->domain->m[osl_constraint][osl_t_index+i]);

   // copy params -> not applicable here

   // copy const == last column
   cst->val[pl_constraint][cst->ncols-1] = osl_int_get_si(in_dep->domain->precision,
         in_dep->domain->m[osl_constraint][in_dep->domain->nb_columns-1]);

   osl_constraint++;
   pl_constraint++;
 }

 // return new domain
 return cst;
}



/* Get the position of this access given a CandlStmt access matrix
 * (concatenated)
 * ref: starting row for a particular access in concatenated rows of
 * access functions
 * Return the position of this access in the list  */
/*static int get_access_position(CandlMatrix *accesses, int ref)
{
    int num, i;

    num = -1;
    for (i=0; i<=ref; i++)  {
        if (accesses->p[i][0] != 0)   {
            num++;
        }
    }
    assert(num >= 0);
    return num;
}*/


/* Read dependences from candl structures */
static Dep **deps_read(osl_dependence_p candlDeps, PlutoProg *prog)
{
    int i, ndeps;
    int spos, tpos;
    Dep **deps;
    int npar = prog->npar;
    Stmt **stmts = prog->stmts;

    ndeps = osl_nb_dependences(candlDeps);

    deps = (Dep **) malloc(ndeps*sizeof(Dep *));

    for (i=0; i<ndeps; i++) {
        deps[i] = pluto_dep_alloc();
    }

    osl_dependence_p candl_dep = candlDeps;

    candl_dep = candlDeps;

    IF_DEBUG(candl_dependence_pprint(stdout, candl_dep));

    /* Dependence polyhedra information */
    for (i=0; i<ndeps; i++)  {
        Dep *dep = deps[i];
        dep->id = i;
        dep->type = candl_dep->type;
        dep->src = candl_dep->label_source;
        dep->dest = candl_dep->label_target;

        //candl_matrix_print(stdout, candl_dep->domain);
        dep->dpolytope = osl_dep_domain_to_pluto_constraints(candl_dep);

        switch (dep->type) {
            case OSL_DEPENDENCE_RAW: 
                spos = get_osl_write_access_position(
                             candl_dep->stmt_source_ptr->access,
                             candl_dep->ref_source_access_ptr);
                dep->src_acc = stmts[dep->src]->writes[spos];
                tpos = get_osl_read_access_position(
                             candl_dep->stmt_target_ptr->access,
                             candl_dep->ref_target_access_ptr);
                dep->dest_acc = stmts[dep->dest]->reads[tpos];
                    
                break;
            case OSL_DEPENDENCE_WAW: 
                spos = get_osl_write_access_position(
                             candl_dep->stmt_source_ptr->access,
                             candl_dep->ref_source_access_ptr);
                dep->src_acc = stmts[dep->src]->writes[spos];
                tpos = get_osl_write_access_position(
                             candl_dep->stmt_target_ptr->access,
                             candl_dep->ref_target_access_ptr);
                dep->dest_acc = stmts[dep->dest]->writes[tpos];
                break;
            case OSL_DEPENDENCE_WAR: 
                spos = get_osl_read_access_position(
                             candl_dep->stmt_source_ptr->access,
                             candl_dep->ref_source_access_ptr);
                dep->src_acc = stmts[dep->src]->reads[spos];
                tpos = get_osl_write_access_position(
                             candl_dep->stmt_target_ptr->access,
                             candl_dep->ref_target_access_ptr);
                dep->dest_acc = stmts[dep->dest]->writes[tpos];
                break;
            case OSL_DEPENDENCE_RAR: 
                spos = get_osl_read_access_position(
                             candl_dep->stmt_source_ptr->access,
                             candl_dep->ref_source_access_ptr);
                dep->src_acc = stmts[dep->src]->reads[spos];
                tpos = get_osl_read_access_position(
                             candl_dep->stmt_target_ptr->access,
                             candl_dep->ref_target_access_ptr);
                dep->dest_acc = stmts[dep->dest]->reads[tpos];
                break;
            default:
                assert(0);
        }

        /* Get rid of rows that are all zero */
        int r, c;
        bool *remove = (bool *) malloc(sizeof(bool)*dep->dpolytope->nrows);
        for (r=0; r<dep->dpolytope->nrows; r++) {
            for (c=0; c<dep->dpolytope->ncols; c++) {
                if (dep->dpolytope->val[r][c] != 0) {
                    break;
                }
            }
            if (c == dep->dpolytope->ncols) {
                remove[r] = true;
            }else{
                remove[r] = false;
            }
        }
        int orig_nrows = dep->dpolytope->nrows;
        int del_count = 0;
        for (r=0; r<orig_nrows; r++) {
            if (remove[r])  {
                pluto_constraints_remove_row(dep->dpolytope, r-del_count);
                del_count++;
            }
        }
        free(remove);

        int src_dim = stmts[dep->src]->dim;
        int target_dim = stmts[dep->dest]->dim;

        assert(candl_dep->source_nb_output_dims_domain + 
               candl_dep->target_nb_output_dims_domain +
               candl_dep->stmt_source_ptr->domain->nb_parameters + 1
               == src_dim+target_dim+npar+1);

        candl_dep = candl_dep->next;
    }

    return deps;
}

void pluto_dep_print(FILE *fp, Dep *dep)
{
    fprintf(fp, "--- Dep %d from S%d to S%d; satisfied: %d, sat level: %d; Type: ",
            dep->id+1, dep->src+1, dep->dest+1, dep->satisfied, dep->satisfaction_level);

    switch (dep->type) {
        case OSL_UNDEFINED : fprintf(fp, "UNSET"); break;
        case OSL_DEPENDENCE_RAW   : fprintf(fp, "RAW")  ; break;
        case OSL_DEPENDENCE_WAR   : fprintf(fp, "WAR")  ; break;
        case OSL_DEPENDENCE_WAW   : fprintf(fp, "WAW")  ; break;
        case OSL_DEPENDENCE_RAR   : fprintf(fp, "RAR")  ; break;
        default : fprintf(fp, "unknown"); break;
    }

    fprintf(fp, "\n");
    if (dep->src_acc != NULL) {
        fprintf(fp, "Var: %s\n", dep->src_acc->name);
    }

    fprintf(fp, "Dependence polyhedron\n");
    pluto_constraints_print(fp, dep->dpolytope);
    fprintf(fp, "\n");
}


void pluto_deps_print(FILE *fp, PlutoProg *prog)
{
    int i;
    for (i=0; i<prog->ndeps; i++) {
        pluto_dep_print(fp, prog->deps[i]);
    }
}


/* Read statement info from openscop structures (nvar: max domain dim) */
static Stmt **osl_to_pluto_stmts(const osl_scop_p scop)
{
    int i, j;
    Stmt **stmts;
    int npar, nvar, nstmts, max_sched_rows;
    osl_statement_p scop_stmt;

    npar = scop->context->nb_parameters;
    nstmts = osl_statement_number(scop->statement);

    if (nstmts == 0)    return NULL;

    /* Max dom dimensionality */
    nvar = -1;
    max_sched_rows = 0;
    scop_stmt = scop->statement;
    for (i=0; i<nstmts; i++) {
        nvar = PLMAX(nvar, osl_statement_get_nb_iterators(scop_stmt));
        max_sched_rows = PLMAX(max_sched_rows, scop_stmt->scattering->nb_rows);
        scop_stmt = scop_stmt->next;
    }

    /* Allocate more to account for unroll/jamming later on */
    stmts = (Stmt **) malloc(nstmts*sizeof(Stmt *));

    scop_stmt = scop->statement;

    for(i=0; i<nstmts; i++)  {
        PlutoConstraints *domain = 
            osl_relation_to_pluto_constraints(scop_stmt->domain);
        PlutoMatrix *trans = osl_scattering_to_pluto_trans(scop_stmt->scattering);

        int nb_iter = osl_statement_get_nb_iterators(scop_stmt);

        stmts[i] = pluto_stmt_alloc(nb_iter, domain, trans);

        /* Pad with all zero rows */
        int curr_sched_rows = stmts[i]->trans->nrows;
        for (j=curr_sched_rows; j<max_sched_rows; j++) {
            pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, j);
        }

        pluto_constraints_free(domain);
        pluto_matrix_free(trans);

        Stmt *stmt = stmts[i];

        stmt->id = i;
        stmt->type = ORIG;

        assert(scop_stmt->domain->nb_columns-1 == stmt->dim + npar + 1);

        for (j=0; j<stmt->dim; j++)  {
            stmt->is_orig_loop[j] = true;
        }

        /* Tile it if it's tilable unless turned off by .fst/.precut file */
        stmt->tile = 1;

        osl_body_p stmt_body = (osl_body_p)(scop_stmt->body->data);

        for (j=0; j<stmt->dim; j++)    {
            stmt->iterators[j] = strdup(stmt_body->iterators->string[j]);
        }

        /* Statement text */
        stmt->text = osl_strings_sprint(stmt_body->expression); //appends \n
        stmt->text[strlen(stmt->text)-1] = '\0';  //remove the \n from end


        /* Read/write accesses */
        osl_relation_list_p wlist  = osl_access_list_filter_write(scop_stmt->access);
        osl_relation_list_p rlist  = osl_access_list_filter_read(scop_stmt->access);

        osl_relation_list_p rlist_t, wlist_t;
        rlist_t = rlist;
        wlist_t = wlist;

        stmt->nwrites = osl_relation_list_count(wlist);
        stmt->writes = (PlutoAccess **) malloc(stmt->nwrites*sizeof(PlutoAccess *));

        stmt->nreads =  osl_relation_list_count(rlist);
        stmt->reads = (PlutoAccess **) malloc(stmt->nreads*sizeof(PlutoAccess *));

        osl_arrays_p arrays = osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);

        int count = 0;
        while (wlist != NULL)   {
            PlutoMatrix *wmat = osl_access_relation_to_pluto_matrix(wlist->elt);
            stmt->writes[count] = (PlutoAccess *) malloc(sizeof(PlutoAccess));
            stmt->writes[count]->mat = wmat;

            //stmt->writes[count]->symbol = NULL;
            if(arrays){
              int id = osl_relation_get_array_id(wlist->elt);
              stmt->writes[count]->name = strdup(arrays->names[id-1]);
            }
            else{
              stmt->writes[count]->name = NULL;
            }

            count++;
            wlist = wlist->next;
        }

        count = 0;
        while (rlist != NULL)   {
           
            PlutoMatrix *rmat = osl_access_relation_to_pluto_matrix(rlist->elt);
            stmt->reads[count] = (PlutoAccess *) malloc(sizeof(PlutoAccess));
            stmt->reads[count]->mat = rmat;

            //stmt->reads[count]->symbol = NULL;
            if(arrays){
              int id = osl_relation_get_array_id(rlist->elt);
              stmt->reads[count]->name = strdup(arrays->names[id-1]);
            }
            else{
              stmt->reads[count]->name = NULL;
            }

            count++;
            rlist = rlist->next;
        }

        osl_relation_list_free(wlist_t);
        osl_relation_list_free(rlist_t);

        scop_stmt = scop_stmt->next;
    }

    return stmts;
}

void pluto_stmt_print(FILE *fp, const Stmt *stmt)
{
    int i;

    fprintf(fp, "S%d \"%s\"; ndims: %d; orig_depth: %d\n", 
            stmt->id+1, stmt->text, stmt->dim, stmt->dim_orig);
    fprintf(fp, "Domain\n");
    pluto_constraints_print(fp, stmt->domain);
    fprintf(fp, "Transformation\n");
    pluto_matrix_print(fp, stmt->trans);

    if (stmt->nreads==0) {
        fprintf(fp, "No Read accesses\n");
    }else{
        fprintf(fp, "Read accesses\n");
        for (i=0; i<stmt->nreads; i++)  {
            pluto_matrix_print(fp, stmt->reads[i]->mat);
        }
    }

    if (stmt->nwrites==0) {
        fprintf(fp, "No write access\n");
    }else{
        fprintf(fp, "Write accesses\n");
        for (i=0; i<stmt->nwrites; i++)  {
            pluto_matrix_print(fp, stmt->writes[i]->mat);
        }
    }

    for (i=0; i<stmt->dim; i++) {
        printf("Original loop: %d -> %d\n", i, stmt->is_orig_loop[i]);
    }

    fprintf(fp, "\n");
}


void pluto_stmts_print(FILE *fp, Stmt **stmts, int nstmts)
{
    int i;

    for(i=0; i<nstmts; i++)  {
        pluto_stmt_print(fp, stmts[i]);
    }
}


void pluto_prog_print(PlutoProg *prog)
{
    printf("nvar = %d, npar = %d\n", prog->nvar, prog->npar);

    pluto_stmts_print(stdout, prog->stmts, prog->nstmts);
    pluto_deps_print(stdout, prog);
    pluto_transformations_pretty_print(prog);
}


void pluto_dep_free(Dep *dep)
{
    pluto_constraints_free(dep->dpolytope);
    pluto_constraints_free(dep->depsat_poly);
    if (dep->dirvec) {
        free(dep->dirvec);
    }
    if (dep->dirvec) {
        free(dep->satvec);
    }
    free(dep);
}


/* Set the dimension names of type "type" according to the elements
 * in the array "names".
 */
static __isl_give isl_dim *set_names(__isl_take isl_dim *dim,
        enum isl_dim_type type, char **names)
{
    int i;

    for (i = 0; i < isl_dim_size(dim, type); ++i)
        dim = isl_dim_set_name(dim, type, i, names[i]);

    return dim;
}


/* Convert a osl_relation_p containing the constraints of a domain
 * to an isl_set.
 * One shot only; does not take into account the next ptr.
 */
static __isl_give isl_set *osl_relation_to_isl_set(osl_relation_p relation,
        __isl_take isl_dim *dim)
{
    int i, j;
    int n_eq = 0, n_ineq = 0;
    isl_ctx *ctx;
    isl_mat *eq, *ineq;
    isl_int v;
    isl_basic_set *bset;

    isl_int_init(v);

    ctx = isl_dim_get_ctx(dim);

    for (i = 0; i < relation->nb_rows; ++i)
        if (osl_int_zero(relation->precision, relation->m[i][0]))
            n_eq++;
        else
            n_ineq++;

    eq = isl_mat_alloc(ctx, n_eq, relation->nb_columns - 1);
    ineq = isl_mat_alloc(ctx, n_ineq, relation->nb_columns - 1);

    n_eq = n_ineq = 0;
    for (i = 0; i < relation->nb_rows; ++i) {
        isl_mat **m;
        int row;

        if (osl_int_zero(relation->precision, relation->m[i][0])) {
            m = &eq;
            row = n_eq++;
        } else {
            m = &ineq;
            row = n_ineq++;
        }

        for (j = 0; j < relation->nb_columns - 1; ++j) {
            int t = osl_int_get_si(relation->precision, relation->m[i][1 + j]);
            isl_int_set_si(v, t);
            *m = isl_mat_set_element(*m, row, j, v);
        }
    }

    isl_int_clear(v);

    bset = isl_basic_set_from_constraint_matrices(dim, eq, ineq,
            isl_dim_set, isl_dim_div, isl_dim_param, isl_dim_cst);
    return isl_set_from_basic_set(bset);
}

/* Convert a osl_relation_p describing a union of domains
 * to an isl_set.
 */
static __isl_give isl_set *osl_relation_list_to_isl_set(
        osl_relation_p list, __isl_take isl_dim *dim)
{
    isl_set *set;

    set = isl_set_empty(isl_dim_copy(dim));
    for (; list; list = list->next) {
        isl_set *set_i;
        set_i = osl_relation_to_isl_set(list, isl_dim_copy(dim));
        set = isl_set_union(set, set_i);
    }

    isl_dim_free(dim);
    return set;
}

/* Convert an m x ( n + 1) pluto access_matrix_p [d A c]
 * to an m x (m + n + 1) isl_mat [-I A c].
 */
static __isl_give isl_mat *pluto_extract_equalities(isl_ctx *ctx,
        PlutoMatrix *matrix)
{
    int i, j;
    int n_col, n;
    isl_int v;
    isl_mat *eq;

    n_col = matrix->ncols;
    n = matrix->nrows;

    isl_int_init(v);
    eq = isl_mat_alloc(ctx, n, n + n_col);

    for (i = 0; i < n; ++i) {
        isl_int_set_si(v, 0);
        for (j = 0; j < n; ++j)
            eq = isl_mat_set_element(eq, i, j, v);
        isl_int_set_si(v, -1);
        eq = isl_mat_set_element(eq, i, i, v);
        for (j = 0; j < n_col ; ++j) {
            int t = SCOPVAL_get_si(matrix->val[i][j]);
            isl_int_set_si(v, t);
            eq = isl_mat_set_element(eq, i, n + j, v);
        }
    }

    isl_int_clear(v);

    return eq;
}

/* Convert an m x (1 + m + n + 1) osl_relation_p [d -I A c]
 * to an m x (m + n + 1) isl_mat [-I A c].
 */
static __isl_give isl_mat *extract_equalities_osl(isl_ctx *ctx,
        osl_relation_p relation)
{
    int i, j;
    int n_col, n_row;
    isl_int v;
    isl_mat *eq;

    n_col = relation->nb_columns;
    n_row = relation->nb_rows;

    isl_int_init(v);
    eq = isl_mat_alloc(ctx, n_row, n_col - 1);

    for (i = 0; i < n_row; ++i) {
        for (j = 0; j < n_col - 1; ++j) {
            int row = osl_relation_get_row_id_for_nth_dimension(relation, i+1);
            int t = osl_int_get_si(relation->precision, relation->m[row][1 + j]);
            isl_int_set_si(v, t);
            eq = isl_mat_set_element(eq, i, j, v);
        }
    }

    isl_int_clear(v);

    return eq;
}



/* Convert an m x (1 + m + n + 1) osl_relation_p [d -I A c]
 * to an m x (m + n + 1) isl_mat [-I A c].
 */
static __isl_give isl_mat *extract_equalities_osl_access(isl_ctx *ctx,
        osl_relation_p relation)
{
    int i, j;
    int n_col, n_row;
    isl_int v;
    isl_mat *eq;

    n_row = relation->nb_rows==1?1:relation->nb_rows-1;
    n_col = relation->nb_columns - (relation->nb_rows==1?1:2);

    isl_int_init(v);
    eq = isl_mat_alloc(ctx, n_row, n_col);

    if(relation->nb_rows==1){
      isl_int_set_si(v, -1);
      eq = isl_mat_set_element(eq, 0, 0, v);
      for (j = 1; j < n_col; ++j) {
        isl_int_set_si(v, 0);
        eq = isl_mat_set_element(eq, 0, j, v);
      }
    }
    else{
      for (i = 1; i < relation->nb_rows; ++i) {
          for (j = 2; j < relation->nb_columns; ++j) {
              int row = osl_relation_get_row_id_for_nth_dimension(relation, i+1);
              int t = osl_int_get_si(relation->precision, relation->m[row][j]);
              isl_int_set_si(v, t);
              eq = isl_mat_set_element(eq, i-1, j-2, v);
          }
      }
    }

    isl_int_clear(v);

    return eq;
}


/* Convert a osl_relation_p scattering [0 M A c] to
 * the isl_map { i -> A i + c } in the space prescribed by "dim".
 */
static __isl_give isl_map *osl_scattering_to_isl_map(
        osl_relation_p scattering, __isl_take isl_dim *dim)
{
    int n_col;
    isl_ctx *ctx;
    isl_mat *eq, *ineq;
    isl_basic_map *bmap;

    ctx = isl_dim_get_ctx(dim);
    n_col = scattering->nb_columns;

    ineq = isl_mat_alloc(ctx, 0, n_col - 1);
    eq = extract_equalities_osl(ctx, scattering);

    bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq,
            isl_dim_out, isl_dim_in, isl_dim_div, isl_dim_param, isl_dim_cst);

    return isl_map_from_basic_map(bmap);
}



/* Convert a osl_relation_list_p describing a series of accesses [eq -I B c]
 * to an isl_union_map with domain "dom" (in space "D").
 * The -I columns identify the output dimensions of the access, the first
 * of them being the identity of the array being accessed.  The remaining
 * output dimensions identiy the array subscripts.
 *
 * Let "A" be array identified by the first entry.
 * The input dimension columns have the form [B c].
 * Each such access is converted to a map { D[i] -> A[B i + c] } * dom.
 *
 */
static __isl_give isl_union_map *osl_access_list_to_isl_union_map(
        osl_relation_list_p list, __isl_take isl_set *dom, char **arrays)
{
    int len, n_col;
    isl_ctx *ctx;
    isl_dim *dim;
    isl_mat *eq, *ineq;
    isl_union_map *res;

    ctx = isl_set_get_ctx(dom);

    dim = isl_set_get_dim(dom);
    dim = isl_dim_drop(dim, isl_dim_set, 0, isl_dim_size(dim, isl_dim_set));
    res = isl_union_map_empty(dim);

    for ( ; list; list = list->next) {

        n_col = list->elt->nb_columns - (list->elt->nb_rows==1?1:2);
        len   = list->elt->nb_rows==1?1:list->elt->nb_rows-1;

        isl_basic_map *bmap;
        isl_map *map;
        int arr = osl_relation_get_array_id(list->elt) - 1;

        dim = isl_set_get_dim(dom);
        dim = isl_dim_from_domain(dim);
        dim = isl_dim_add(dim, isl_dim_out, len);
        dim = isl_dim_set_tuple_name(dim, isl_dim_out, arrays[arr]);

        ineq = isl_mat_alloc(ctx, 0, n_col);
        eq = extract_equalities_osl_access(ctx, list->elt);

        bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq,
                isl_dim_out, isl_dim_in, isl_dim_div, isl_dim_param, isl_dim_cst);
        map = isl_map_from_basic_map(bmap);
        map = isl_map_intersect_domain(map, isl_set_copy(dom));
        res = isl_union_map_union(res, isl_union_map_from_map(map));
    }

    isl_set_free(dom);

    return res;
}

/*
 * Like osl_access_list_to_isl_union_map, but just for a single osl access
 * (read or write)
 */
static __isl_give isl_map *osl_basic_access_to_isl_union_map(
        osl_relation_p access, __isl_take isl_set *dom, 
        char **arrays)
{
    int len, n_col;
    isl_ctx *ctx;
    isl_dim *dim;
    isl_mat *eq, *ineq;

    ctx = isl_set_get_ctx(dom);

    n_col = access->nb_columns - (access->nb_rows==1?1:2);
    len   = access->nb_rows==1?1:access->nb_rows-1;


    isl_basic_map *bmap;
    isl_map *map;
    int arr = osl_relation_get_array_id(access) - 1;

    dim = isl_set_get_dim(dom);
    dim = isl_dim_from_domain(dim);
    dim = isl_dim_add(dim, isl_dim_out, len);
    dim = isl_dim_set_tuple_name(dim, isl_dim_out, arrays[arr]);

    ineq = isl_mat_alloc(ctx, 0, n_col);
    eq = extract_equalities_osl_access(ctx, access);

    bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq,
            isl_dim_out, isl_dim_in, isl_dim_div, isl_dim_param, isl_dim_cst);
    map = isl_map_from_basic_map(bmap);
    map = isl_map_intersect_domain(map, dom);

    return map;
}


/*
 * Like osl_access_list_to_isl_union_map, but just for a single pluto access
 * (read or write)
 * pos: position (starting row) of the access in 'access'
 */
static __isl_give isl_map *pluto_basic_access_to_isl_union_map(
        PlutoMatrix  *mat, char* access_name,  __isl_take isl_set *dom)
{
    int len, n_col;
    isl_ctx *ctx;
    isl_dim *dim;
    isl_mat *eq, *ineq;

    ctx = isl_set_get_ctx(dom);

    dim = isl_set_get_dim(dom);
    dim = isl_dim_drop(dim, isl_dim_set, 0, isl_dim_size(dim, isl_dim_set));

    n_col = mat->ncols;

    isl_basic_map *bmap;
    isl_map *map;
    //int arr = SCOPVAL_get_si(access->p[pos][0]) - 1;

    len = mat->nrows;

    dim = isl_set_get_dim(dom);
    dim = isl_dim_from_domain(dim);
    dim = isl_dim_add(dim, isl_dim_out, len);
    dim = isl_dim_set_tuple_name(dim, isl_dim_out, access_name);

    ineq = isl_mat_alloc(ctx, 0, len + n_col);
    eq = pluto_extract_equalities(ctx, mat);

    bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq,
            isl_dim_out, isl_dim_in, isl_dim_div, isl_dim_param, isl_dim_cst);
    map = isl_map_from_basic_map(bmap);
    map = isl_map_intersect_domain(map, dom);

    return map;
}


static int basic_map_count(__isl_take isl_basic_map *bmap, void *user)
{
    int *count = user;

    *count += 1;
    isl_basic_map_free(bmap);
    return 0;
}


int isl_map_count(__isl_take isl_map *map, void *user)
{
    int r;

    r = isl_map_foreach_basic_map(map, &basic_map_count, user);
    isl_map_free(map);
    return r;
}


/* Temporary data structure used inside extract_deps.
 *
 * deps points to the array of Deps being constructed
 * type is the type of the next Dep
 * index is the index of the next Dep in the array.
 */
struct pluto_extra_dep_info {
    Dep **deps;
    Stmt **stmts;
    int type;
    int index;
};


/* Convert an isl_basic_map describing part of a dependence to a Dep.
 * The names of the input and output spaces are of the form S_d or S_d_e
 * with d an integer identifying the statement, e identifying the access
 * (relative to the statement). If it's of the form S_d_e and read/write
 * accesses for the statement are available, source and target accesses 
 * are set for the dependence, otherwise not.
 */
static int basic_map_extract(__isl_take isl_basic_map *bmap, void *user)
{
    Stmt **stmts;
    Dep *dep;
    struct pluto_extra_dep_info *info;
    info = (struct pluto_extra_dep_info *)user;

    stmts = info->stmts;

    bmap = isl_basic_map_remove_divs(bmap);

    dep = info->deps[info->index];

    dep->id = info->index;
    dep->dpolytope = isl_basic_map_to_pluto_constraints(bmap);
    dep->dirvec = NULL;
    dep->type = info->type;
    dep->src = atoi(isl_basic_map_get_tuple_name(bmap, isl_dim_in) + 2);
    dep->dest = atoi(isl_basic_map_get_tuple_name(bmap, isl_dim_out) + 2);

    // pluto_stmt_print(stdout, stmts[dep->src]);
    // pluto_stmt_print(stdout, stmts[dep->dest]);
    // printf("Src acc: %d dest acc: %d\n", src_acc_num, dest_acc_num);

    if (stmts[dep->src]->reads != NULL && stmts[dep->dest]->reads != NULL) {
        /* Extract access function information */
        int src_acc_num, dest_acc_num;
        const char *name;
        name = isl_basic_map_get_tuple_name(bmap, isl_dim_in) + 2;
        while (*name != '\0' && *(name++) != '_');
        if (*name != '\0') src_acc_num = atoi(name+1);
        else assert(0); // access function num not encoded in dependence

        name = isl_basic_map_get_tuple_name(bmap, isl_dim_out) + 2;
        while (*name != '\0' && *(name++) != '_');
        if (*name != '\0') dest_acc_num = atoi(name+1);
        else assert(0); // access function num not encoded in dependence

        switch (info->type) {
            case OSL_DEPENDENCE_RAW: 
                dep->src_acc = stmts[dep->src]->writes[src_acc_num];
                dep->dest_acc = stmts[dep->dest]->reads[dest_acc_num];
                break;
            case OSL_DEPENDENCE_WAW: 
                dep->src_acc = stmts[dep->src]->writes[src_acc_num];
                dep->dest_acc = stmts[dep->dest]->writes[dest_acc_num];
                break;
            case OSL_DEPENDENCE_WAR: 
                dep->src_acc = stmts[dep->src]->reads[src_acc_num];
                dep->dest_acc = stmts[dep->dest]->writes[dest_acc_num];
                break;
            case OSL_DEPENDENCE_RAR: 
                dep->src_acc = stmts[dep->src]->reads[src_acc_num];
                dep->dest_acc = stmts[dep->dest]->reads[dest_acc_num];
                break;
            default:
                assert(0);
        }
    }else{
        dep->src_acc = NULL;
        dep->dest_acc = NULL;
    }

    info->index++;
    isl_basic_map_free(bmap);
    return 0;
}


static int map_extract(__isl_take isl_map *map, void *user)
{
    int r;

    r = isl_map_foreach_basic_map(map, &basic_map_extract, user);
    isl_map_free(map);
    return r;
}


int extract_deps(Dep **deps, int first, Stmt **stmts,
        __isl_keep isl_union_map *umap, int type)
{
    struct pluto_extra_dep_info info = { deps, stmts, type, first };

    isl_union_map_foreach_map(umap, &map_extract, &info);

    return info.index - first;
}


osl_names_p get_scop_names(osl_scop_p scop){

  //generate temp names
  osl_names_p names = osl_scop_names(scop);

  //if scop has names substitute them for temp names
  if(scop->context->nb_parameters){
    osl_strings_free(names->parameters);
    names->parameters = osl_strings_clone((osl_strings_p)scop->parameters->data);
  }

  osl_arrays_p arrays = osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);
  if(arrays){
    osl_strings_free(names->arrays);
    names->arrays = osl_arrays_to_strings(arrays);
  }

  return names;
}

/* Compute dependences based on the iteration domain and access
 * information in "scop" and put the result in "prog".
 *
 * If options->lastwriter is false, then
 *      RAW deps are those from any earlier write to a read
 *      WAW deps are those from any earlier write to a write
 *      WAR deps are those from any earlier read to a write
 *      RAR deps are those from any earlier read to a read
 * If options->lastwriter is true, then
 *      RAW deps are those from the last write to a read
 *      WAW deps are those from the last write to a write
 *      WAR deps are those from any earlier read not masked by an intermediate
 *      write to a write
 *      RAR deps are those from the last read to a read
 *
 * The RAR deps are only computed if options->rar is set.
 */
static void compute_deps(osl_scop_p scop, PlutoProg *prog,
        PlutoOptions *options)
{
    int i, racc_num, wacc_num;
    int nstmts = osl_statement_number(scop->statement);
    isl_ctx *ctx;
    isl_dim *dim;
    isl_space *param_space;
    isl_set *context;
    isl_union_map *empty;
    isl_union_map *write;
    isl_union_map *read;
    isl_union_map *schedule;
    isl_union_map *dep_raw, *dep_war, *dep_waw, *dep_rar, *trans_dep_war;
    isl_union_map *trans_dep_waw;
    osl_statement_p stmt;
    osl_strings_p scop_params = NULL;

    ctx = isl_ctx_alloc();
    assert(ctx);

    osl_names_p names = get_scop_names(scop);


    dim = isl_dim_set_alloc(ctx, scop->context->nb_parameters, 0);
    if(scop->context->nb_parameters){
      scop_params = (osl_strings_p)scop->parameters->data;
      dim = set_names(dim, isl_dim_param, scop_params->string);
    }
    param_space = isl_space_params(isl_space_copy(dim));
    context = osl_relation_to_isl_set(scop->context, param_space);

    if (!options->rar) dep_rar = isl_union_map_empty(isl_dim_copy(dim));
    empty = isl_union_map_empty(isl_dim_copy(dim));
    write = isl_union_map_empty(isl_dim_copy(dim));
    read = isl_union_map_empty(isl_dim_copy(dim));
    schedule = isl_union_map_empty(dim);

    if (options->isldepcompact) {
        /* Leads to fewer dependences. Each dependence may not have a unique
         * source/target access relating to it, since a union is taken
         * across all reads for a statement (and writes) for a particualr
         * array. Relationship between a dependence and associated dependent
         * data / array elements is lost, and some analyses may not work with
         * such a representation
         */
        for (i = 0, stmt = scop->statement; i < nstmts; ++i, stmt = stmt->next) {
            isl_set *dom;
            isl_map *schedule_i;
            isl_union_map *read_i;
            isl_union_map *write_i;
            char name[20];

            snprintf(name, sizeof(name), "S_%d", i);

            int niter = osl_statement_get_nb_iterators(stmt);
            dim = isl_dim_set_alloc(ctx, scop->context->nb_parameters, niter);
            if(scop->context->nb_parameters){
              scop_params = (osl_strings_p)scop->parameters->data;
              dim = set_names(dim, isl_dim_param, scop_params->string);
            }
            if(niter){
              osl_body_p stmt_body = (osl_body_p)(stmt->body->data);
              dim = set_names(dim, isl_dim_set, stmt_body->iterators->string);
            }
            dim = isl_dim_set_tuple_name(dim, isl_dim_set, name);
            dom = osl_relation_list_to_isl_set(stmt->domain, dim);
            dom = isl_set_intersect_params(dom, isl_set_copy(context));

            dim = isl_dim_alloc(ctx, scop->context->nb_parameters, niter,
                    2 * niter + 1);
            if(scop->context->nb_parameters){
              scop_params = (osl_strings_p)scop->parameters->data;
              dim = set_names(dim, isl_dim_param, scop_params->string);
            }
            if(niter){
              osl_body_p stmt_body = (osl_body_p)(stmt->body->data);
              dim = set_names(dim, isl_dim_in, stmt_body->iterators->string);
            }
            dim = isl_dim_set_tuple_name(dim, isl_dim_in, name);
            schedule_i = osl_scattering_to_isl_map(stmt->scattering, dim);

            osl_relation_list_p rlist  = osl_access_list_filter_read(stmt->access);
            osl_relation_list_p wlist  = osl_access_list_filter_write(stmt->access);

            //osl_arrays_p arrays = osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);

            read_i = osl_access_list_to_isl_union_map(rlist, isl_set_copy(dom),
                     names->arrays->string);
            write_i = osl_access_list_to_isl_union_map(wlist, isl_set_copy(dom),
                     names->arrays->string);

            read = isl_union_map_union(read, read_i);
            write = isl_union_map_union(write, write_i);
            schedule = isl_union_map_union(schedule,
                    isl_union_map_from_map(schedule_i));

            osl_relation_list_free(rlist);
            osl_relation_list_free(wlist);
        }
    }else{
        /* Each dependence is for a particular source and target access. Use
         * <stmt, access> pair while relating to accessed data so each
         * dependence can be associated to a unique source and target access
         */

        for (i = 0, stmt = scop->statement; i < nstmts; ++i, stmt = stmt->next) {
            isl_set *dom;

            racc_num = 0;
            wacc_num = 0;

            osl_relation_list_p access = stmt->access;
            for( ; access; access = access->next) {
                isl_map *read_pos;
                isl_map *write_pos;
                isl_map *schedule_i;

                char name[20];

                if (access->elt->type == OSL_TYPE_READ) {
                    snprintf(name, sizeof(name), "S_%d_r%d", i, racc_num);
                }else{
                    snprintf(name, sizeof(name), "S_%d_w%d", i, wacc_num);
                }

                int niter = osl_statement_get_nb_iterators(stmt);
                dim = isl_dim_set_alloc(ctx, scop->context->nb_parameters, niter);
                if(scop->context->nb_parameters){
                  scop_params = (osl_strings_p)scop->parameters->data;
                  dim = set_names(dim, isl_dim_param, scop_params->string);

                  osl_strings_free(names->parameters);
                  names->parameters = osl_strings_clone(scop_params);
                }
                if(niter){
                  osl_body_p stmt_body = (osl_body_p)(stmt->body->data);
                  dim = set_names(dim, isl_dim_set, stmt_body->iterators->string);

                  osl_strings_free(names->iterators);
                  names->iterators = osl_strings_clone(stmt_body->iterators);
                }
                dim = isl_dim_set_tuple_name(dim, isl_dim_set, name);
                dom = osl_relation_list_to_isl_set(stmt->domain, dim);
                dom = isl_set_intersect_params(dom, isl_set_copy(context));

                dim = isl_dim_alloc(ctx, scop->context->nb_parameters, niter,
                        2 * niter + 1);
                if(scop->context->nb_parameters){
                  scop_params = (osl_strings_p)scop->parameters->data;
                  dim = set_names(dim, isl_dim_param, scop_params->string);
                }
                if(niter){
                  osl_body_p stmt_body = (osl_body_p)(stmt->body->data);
                  dim = set_names(dim, isl_dim_in, stmt_body->iterators->string);
                }
                dim = isl_dim_set_tuple_name(dim, isl_dim_in, name);

                schedule_i = osl_scattering_to_isl_map(stmt->scattering, dim);

                //osl_arrays_p arrays = osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);

                if (access->elt->type == OSL_TYPE_READ) {
                    read_pos = osl_basic_access_to_isl_union_map(access->elt, 
                            dom, names->arrays->string);
                    read = isl_union_map_union(read, isl_union_map_from_map(read_pos));
                }else{
                    write_pos = osl_basic_access_to_isl_union_map(access->elt, 
                            dom, names->arrays->string);
                    write = isl_union_map_union(write, isl_union_map_from_map(write_pos));
                }

                schedule = isl_union_map_union(schedule,
                        isl_union_map_from_map(schedule_i));

                if (access->elt->type == OSL_TYPE_READ) {
                    racc_num++;
                }else{
                    wacc_num++;
                }
            }
        }
    }

    if (options->lastwriter) {
        // compute RAW dependences which do not contain transitive dependences
        isl_union_map_compute_flow(isl_union_map_copy(read),
                isl_union_map_copy(write),
                isl_union_map_copy(empty),
                isl_union_map_copy(schedule),
                &dep_raw, NULL, NULL, NULL);
        // compute WAW and WAR dependences which do not contain transitive dependences
        isl_union_map_compute_flow(isl_union_map_copy(write),
                isl_union_map_copy(write),
                isl_union_map_copy(read),
                isl_union_map_copy(schedule),
                &dep_waw, &dep_war, NULL, NULL);
        // compute WAR dependences which may contain transitive dependences
        isl_union_map_compute_flow(isl_union_map_copy(write),
                isl_union_map_copy(empty),
                isl_union_map_copy(read),
                isl_union_map_copy(schedule),
                NULL, &trans_dep_war, NULL, NULL);
        isl_union_map_compute_flow(isl_union_map_copy(write),
                isl_union_map_copy(empty),
                isl_union_map_copy(write),
                isl_union_map_copy(schedule),
                NULL, &trans_dep_waw, NULL, NULL);
        if (options->rar) {
            // compute RAR dependences which do not contain transitive dependences
            isl_union_map_compute_flow(isl_union_map_copy(read),
                    isl_union_map_copy(read),
                    isl_union_map_copy(empty),
                    isl_union_map_copy(schedule),
                    &dep_rar, NULL, NULL, NULL);
        }
    }else{
        // compute RAW dependences which may contain transitive dependences
        isl_union_map_compute_flow(isl_union_map_copy(read),
                isl_union_map_copy(empty),
                isl_union_map_copy(write),
                isl_union_map_copy(schedule),
                NULL, &dep_raw, NULL, NULL);
        // compute WAR dependences which may contain transitive dependences
        isl_union_map_compute_flow(isl_union_map_copy(write),
                isl_union_map_copy(empty),
                isl_union_map_copy(read),
                isl_union_map_copy(schedule),
                NULL, &dep_war, NULL, NULL);
        // compute WAW dependences which may contain transitive dependences
        isl_union_map_compute_flow(isl_union_map_copy(write),
                isl_union_map_copy(empty),
                isl_union_map_copy(write),
                isl_union_map_copy(schedule),
                NULL, &dep_waw, NULL, NULL);
        if (options->rar) {
            // compute RAR dependences which may contain transitive dependences
            isl_union_map_compute_flow(isl_union_map_copy(read),
                    isl_union_map_copy(empty),
                    isl_union_map_copy(read),
                    isl_union_map_copy(schedule),
                    NULL, &dep_rar, NULL, NULL);
        }
    }

    dep_raw = isl_union_map_coalesce(dep_raw);
    dep_war = isl_union_map_coalesce(dep_war);
    dep_waw = isl_union_map_coalesce(dep_waw);
    dep_rar = isl_union_map_coalesce(dep_rar);

    prog->ndeps = 0;
    isl_union_map_foreach_map(dep_raw, &isl_map_count, &prog->ndeps);
    isl_union_map_foreach_map(dep_war, &isl_map_count, &prog->ndeps);
    isl_union_map_foreach_map(dep_waw, &isl_map_count, &prog->ndeps);
    isl_union_map_foreach_map(dep_rar, &isl_map_count, &prog->ndeps);

    prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i] = pluto_dep_alloc();
    }
    prog->ndeps = 0;
    prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_raw, OSL_DEPENDENCE_RAW);
    prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_war, OSL_DEPENDENCE_WAR);
    prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_waw, OSL_DEPENDENCE_WAW);
    prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_rar, OSL_DEPENDENCE_RAR);

    if (options->lastwriter) {
        trans_dep_war = isl_union_map_coalesce(trans_dep_war);
        trans_dep_waw = isl_union_map_coalesce(trans_dep_waw);

        prog->ntransdeps = 0;
        isl_union_map_foreach_map(dep_raw, &isl_map_count, &prog->ntransdeps);
        isl_union_map_foreach_map(trans_dep_war, &isl_map_count, &prog->ntransdeps);
        isl_union_map_foreach_map(trans_dep_waw, &isl_map_count, &prog->ntransdeps);
        isl_union_map_foreach_map(dep_rar, &isl_map_count, &prog->ntransdeps);

        if (prog->ntransdeps >= 1) {
            prog->transdeps = (Dep **)malloc(prog->ntransdeps * sizeof(Dep *));
            for (i=0; i<prog->ntransdeps; i++) {
                prog->transdeps[i] = pluto_dep_alloc();
            }
            prog->ntransdeps = 0;
            prog->ntransdeps += extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts, dep_raw, OSL_DEPENDENCE_RAW);
            prog->ntransdeps += extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts, trans_dep_war, OSL_DEPENDENCE_WAR);
            prog->ntransdeps += extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts, trans_dep_waw, OSL_DEPENDENCE_WAW);
            prog->ntransdeps += extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts, dep_rar, OSL_DEPENDENCE_RAR);
        }

        isl_union_map_free(trans_dep_war);
        isl_union_map_free(trans_dep_waw);
    }

    isl_union_map_free(dep_raw);
    isl_union_map_free(dep_war);
    isl_union_map_free(dep_waw);
    isl_union_map_free(dep_rar);

    isl_union_map_free(empty);
    isl_union_map_free(write);
    isl_union_map_free(read);
    isl_union_map_free(schedule);
    isl_set_free(context);

    if(names) osl_names_free(names);

    isl_ctx_free(ctx);
}

/*scoplib_matrix_p get_identity_schedule(int dim, int npar){
    scoplib_matrix_p smat = scoplib_matrix_malloc(2*dim+1, dim+npar+1+1);

    int i, j;
    for(i =0; i<2*dim+1; i++)
        for(j=0; j<dim+1+npar+1; j++)
            smat->p[i][j] = 0;

    for(i=1; i<dim; i++)
        smat->p[2*i-1][i] = 1;

    return smat;

}*/

/* 
 * Extract necessary information from clan_scop to create PlutoProg - a
 * representation of the program sufficient to be used throughout Pluto. 
 * PlutoProg also includes dependences; so candl is run here.
 */
PlutoProg *scop_to_pluto_prog(osl_scop_p scop, PlutoOptions *options)
{
    int i, max_sched_rows;

    PlutoProg *prog = pluto_prog_alloc();

    prog->nstmts = osl_statement_number(scop->statement);
    prog->options = options;

    /* Data variables in the program */
    osl_arrays_p arrays = osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);
    if(arrays==NULL){
      prog->num_data = 0;
      fprintf(stderr, "warning: arrays extension not found\n");
    }
    else{
      prog->num_data = arrays->nb_names;
      prog->data_names = (char **) malloc (prog->num_data * sizeof(char *));
      for(i=0; i< prog->num_data; i++) {
          prog->data_names[i] = strdup(arrays->names[i]);
      }
    }


    /* Program parameters */
    prog->npar = scop->context->nb_parameters;

    osl_strings_p osl_scop_params = NULL;
    if (prog->npar >= 1)    {
        prog->params = (char **) malloc(sizeof(char *)*prog->npar);
        osl_scop_params = (osl_strings_p)scop->parameters->data;
    }
    for (i=0; i<prog->npar; i++)  {
        prog->params[i] = strdup(osl_scop_params->string[i]);
    }

    pluto_constraints_free(prog->context);
    prog->context = osl_relation_to_pluto_constraints(scop->context);

    if (options->context != -1)	{
        for (i=0; i<prog->npar; i++)  {
            pluto_constraints_add_inequality(prog->context);
            prog->context->val[i][i] = 1;
            prog->context->val[i][prog->context->ncols-1] = -options->context;
        }
    }

    osl_statement_p scop_stmt = scop->statement;

    prog->nvar = osl_statement_get_nb_iterators(scop_stmt);
    max_sched_rows = 0;
    for (i=0; i<prog->nstmts; i++) {
        int stmt_num_iter = osl_statement_get_nb_iterators(scop_stmt); 
        prog->nvar = PLMAX(prog->nvar, stmt_num_iter);
        max_sched_rows = PLMAX(max_sched_rows, scop_stmt->scattering->nb_rows);
        scop_stmt = scop_stmt->next;
    }

    prog->stmts = osl_to_pluto_stmts(scop);
    prog->scop = scop;

    /* Compute dependences */
    if (options->isldep) {
        compute_deps(scop, prog, options);
    }else{
        /*  Using Candl */
        candl_options_p candlOptions = candl_options_malloc();
        if (options->rar)   {
            candlOptions->rar = 1;
        }
        candlOptions->lastwriter = options->lastwriter;
        candlOptions->scalar_privatization = options->scalpriv;
        // candlOptions->verbose = 1;

        /* Add more infos (depth, label, ...) */
        /* Needed by Candl */
        candl_scop_usr_init(scop);

        osl_dependence_p candl_deps = candl_dependence(scop, candlOptions);
        prog->deps = deps_read(candl_deps, prog);
        prog->ndeps = osl_nb_dependences(candl_deps);
        candl_options_free(candlOptions);
        osl_dependence_free(candl_deps);

        candl_scop_usr_cleanup(scop); //undo candl_scop_user_init

        prog->transdeps = NULL;
        prog->ntransdeps = 0;
    }

    /* Add hyperplanes */
    if (prog->nstmts >= 1) {
        for (i=0; i<max_sched_rows; i++) {
            pluto_prog_add_hyperplane(prog,prog->num_hyperplanes,H_UNKNOWN);
            prog->hProps[prog->num_hyperplanes-1].type = 
                (i%2)? H_LOOP: H_SCALAR;
        }
    }

    /* Hack for linearized accesses */
    FILE *lfp = fopen(".linearized", "r");
    FILE *nlfp = fopen(".nonlinearized", "r");
    char tmpstr[256];
    char linearized[256];
    if (lfp && nlfp) {
        for (i=0; i<prog->nstmts; i++)    {
            rewind(lfp);
            rewind(nlfp);
            while (!feof(lfp) && !feof(nlfp))      {
                fgets(tmpstr, 256, nlfp);
                fgets(linearized, 256, lfp);
                if (strstr(tmpstr, prog->stmts[i]->text))        {
                    prog->stmts[i]->text = (char *) realloc(prog->stmts[i]->text, sizeof(char)*(strlen(linearized)+1));
                    strcpy(prog->stmts[i]->text, linearized);
                }
            }
        }
        fclose(lfp);
        fclose(nlfp);
    }

    return prog;
}

/* Get an upper bound for transformation coefficients to prevent spurious
 * transformations that represent shifts or skews proportional to trip counts:
 * this happens when loop bounds are constants
 */
int get_coeff_upper_bound(PlutoProg *prog)
{
    int max, i, r;

    max = 0;
    for (i=0; i<prog->nstmts; i++)  {
        Stmt *stmt = prog->stmts[i];
        for (r=0; r<stmt->domain->nrows; r++) {
            max  = PLMAX(max,stmt->domain->val[r][stmt->domain->ncols-1]);
        }
    }

    return max-1;
}


PlutoProg *pluto_prog_alloc()
{
    PlutoProg *prog = (PlutoProg *) malloc(sizeof(PlutoProg));

    prog->nstmts = 0;
    prog->stmts = NULL;
    prog->npar = 0;
    prog->nvar = 0;
    prog->params = NULL;
    prog->context = pluto_constraints_alloc(1, prog->npar+1);
    prog->deps = NULL;
    prog->ndeps = 0;
    prog->transdeps = NULL;
    prog->ntransdeps = 0;
    prog->ddg = NULL;
    prog->hProps = NULL;
    prog->num_hyperplanes = 0;

    prog->globcst = NULL;
    prog->depcst = NULL;

    return prog;
}



void pluto_prog_free(PlutoProg *prog)
{
    int i;

    /* Free dependences */
    for (i=0; i<prog->ndeps; i++) {
        pluto_dep_free(prog->deps[i]);
    }
    if (prog->deps != NULL) {
        free(prog->deps);
    }

    /* Free DDG */
    if (prog->ddg != NULL)  {
        graph_free(prog->ddg);
    }

    if (prog->hProps != NULL)   {
        free(prog->hProps);
    }

    for (i=0; i<prog->npar; i++)  {
        free(prog->params[i]);
    }
    if (prog->npar >= 1)    {
        free(prog->params);
    }

    /* Statements */
    for (i=0; i<prog->nstmts; i++) {
        pluto_stmt_free(prog->stmts[i]);
    }
    if (prog->nstmts >= 1)  {
        free(prog->stmts);
    }

    pluto_constraints_free(prog->context);
    if (prog->depcst != NULL) {
        pluto_constraints_free(*prog->depcst);
    }
    free(prog->depcst);
    pluto_constraints_free(prog->globcst);

    free(prog);
}


PlutoOptions *pluto_options_alloc()
{
    PlutoOptions *options;

    options  = (PlutoOptions *) malloc(sizeof(PlutoOptions));

    /* Initialize to default */
    options->tile = 0;
    options->intratileopt = 1;
    options->debug = 0;
    options->moredebug = 0;
    options->scancount = 0;
    options->parallel = 0;
    options->innerpar = 0;
    options->identity = 0;
    options->unroll = 0;

    /* Unroll/jam factor */
    options->ufactor = 8;

    /* Ignore input deps */
    options->rar = 0;

    /* Override for first and last levels to tile */
    options->ft = -1;
    options->lt = -1;

    /* Override for first and last cloog options */
    options->cloogf = -1;
    options->cloogl = -1;

    options->cloogsh = 0;

    options->cloogbacktrack = 1;

    options->multipipe = 0;
    options->l2tile = 0;
    options->prevector = 1;
    options->fuse = SMART_FUSE;

    /* Experimental */
    options->polyunroll = 0;

    /* Default context is no context */
    options->context = -1;

    options->forceparallel = -42;

    options->bee = 0;

    options->isldep = 0;
    options->isldepcompact = 0;

    options->islsolve = 0;

    options->readscop = 0;

    options->lastwriter = 0;

    options->nobound = 0;

    options->scalpriv = 0;

    options->silent = 0;

    options->out_file = NULL;

    return options;
}


/* Add global/program parameter at position 'pos' */
void pluto_prog_add_param(PlutoProg *prog, const char *param, int pos)
{
    int i, j;

    for (i=0; i<prog->nstmts; i++) {
        Stmt *stmt = prog->stmts[i];
        pluto_constraints_add_dim(stmt->domain, stmt->domain->ncols-1-prog->npar+pos);
        pluto_matrix_add_col(stmt->trans, stmt->trans->ncols-1-prog->npar+pos);

        for (j=0; j<stmt->nwrites; j++)  {
            pluto_matrix_add_col(stmt->writes[j]->mat, stmt->dim+pos);
        }
        for (j=0; j<stmt->nreads; j++)  {
            pluto_matrix_add_col(stmt->reads[j]->mat, stmt->dim+pos);
        }
    }
    for (i=0; i<prog->ndeps; i++)   {
        pluto_constraints_add_dim(prog->deps[i]->dpolytope, 
                prog->deps[i]->dpolytope->ncols-1-prog->npar+pos);
    }
    pluto_constraints_add_dim(prog->context, prog->context->ncols-1-prog->npar+pos);

    prog->params = (char **) realloc(prog->params, sizeof(char *)*(prog->npar+1));

    for (i=prog->npar-1; i>=pos; i--)    {
        prog->params[i+1] = prog->params[i];
    }

    prog->params[pos] = strdup(param);
    prog->npar++;
}


void pluto_options_free(PlutoOptions *options)
{
    if (options->out_file != NULL)  {
        free(options->out_file);
    }
    free(options);
}


/* pos: position of domain iterator 
 * time_pos: position of time iterator; iter: domain iterator; supply -1
 * if you don't want a scattering function row added for it */
void pluto_stmt_add_dim(Stmt *stmt, int pos, int time_pos, const char *iter,
        PlutoHypType hyp_type, PlutoProg *prog)
{
    int i, npar;

    npar = stmt->domain->ncols - stmt->dim - 1;

    assert(pos <= stmt->dim);
    assert(time_pos <= stmt->trans->nrows);
    assert(stmt->dim + npar + 1 == stmt->domain->ncols);

    pluto_constraints_add_dim(stmt->domain, pos);
    stmt->dim++;
    stmt->iterators = (char **) realloc(stmt->iterators, stmt->dim*sizeof(char *));
    for (i=stmt->dim-2; i>=pos; i--) {
        stmt->iterators[i+1] = stmt->iterators[i];
    }
    stmt->iterators[pos] = strdup(iter);

    /* Stmt should always have a transformation */
    assert(stmt->trans != NULL);
    pluto_matrix_add_col(stmt->trans, pos);

    if (time_pos != -1) {
        pluto_matrix_add_row(stmt->trans, time_pos);
        stmt->trans->val[time_pos][pos] = 1;


        stmt->hyp_types = realloc(stmt->hyp_types, 
                sizeof(int)*stmt->trans->nrows);
        for (i=stmt->trans->nrows-2; i>=time_pos; i--) {
            stmt->hyp_types[i+1] = stmt->hyp_types[i];
        }
        stmt->hyp_types[time_pos] = hyp_type;
    }

    /* Update is_orig_loop */
    stmt->is_orig_loop = realloc(stmt->is_orig_loop, sizeof(bool)*stmt->dim);
    for (i=stmt->dim-2; i>=pos; i--) {
        stmt->is_orig_loop[i+1] = stmt->is_orig_loop[i];
    }
    stmt->is_orig_loop[pos] = true;

    for (i=0; i<stmt->nwrites; i++)   {
        pluto_matrix_add_col(stmt->writes[i]->mat, pos);
    }

    for (i=0; i<stmt->nreads; i++)   {
        pluto_matrix_add_col(stmt->reads[i]->mat, pos);
    }

    for (i=0; i<prog->ndeps; i++) {
        if (prog->deps[i]->src == stmt->id) {
            pluto_constraints_add_dim(prog->deps[i]->dpolytope, pos);
        }
        if (prog->deps[i]->dest == stmt->id) {
            pluto_constraints_add_dim(prog->deps[i]->dpolytope, 
                    prog->stmts[prog->deps[i]->src]->dim+pos);
        }
    }

    for (i=0; i<prog->ntransdeps; i++) {
        assert(prog->transdeps[i] != NULL);
        if (prog->transdeps[i]->src == stmt->id) {
            pluto_constraints_add_dim(prog->transdeps[i]->dpolytope, pos);
        }
        if (prog->transdeps[i]->dest == stmt->id) {
            pluto_constraints_add_dim(prog->transdeps[i]->dpolytope, 
                    prog->stmts[prog->transdeps[i]->src]->dim+pos);
        }
    }
}

/* Warning: use it only to knock off a dummy dimension (unrelated to 
 * anything else */
void pluto_stmt_remove_dim(Stmt *stmt, int pos, PlutoProg *prog)
{
    int i, npar;

    npar = stmt->domain->ncols - stmt->dim - 1;

    assert(pos <= stmt->dim);
    assert(stmt->dim + npar + 1 == stmt->domain->ncols);

    pluto_constraints_remove_dim(stmt->domain, pos);
    stmt->dim--;

    if (stmt->iterators != NULL) {
        free(stmt->iterators[pos]);
        for (i=pos; i<=stmt->dim-1; i++) {
            stmt->iterators[i] = stmt->iterators[i+1];
        }
        stmt->iterators = (char **) realloc(stmt->iterators, stmt->dim*sizeof(char *));
    }

    pluto_matrix_remove_col(stmt->trans, pos);

    /* Update is_orig_loop */
    for (i=pos; i<=stmt->dim-1; i++) {
        stmt->is_orig_loop[i] = stmt->is_orig_loop[i+1];
    }
    stmt->is_orig_loop = realloc(stmt->is_orig_loop, sizeof(bool)*stmt->dim);

    for (i=0; i<stmt->nwrites; i++)   {
        pluto_matrix_remove_col(stmt->writes[i]->mat, pos);
    }

    for (i=0; i<stmt->nreads; i++)   {
        pluto_matrix_remove_col(stmt->reads[i]->mat, pos);
    }

    for (i=0; i<prog->ndeps; i++) {
        if (prog->deps[i]->src == stmt->id) {
            pluto_constraints_remove_dim(prog->deps[i]->dpolytope, pos);
        }
        if (prog->deps[i]->dest == stmt->id) {
            // if (i==0)  printf("removing dim\n");
            pluto_constraints_remove_dim(prog->deps[i]->dpolytope, 
                    prog->stmts[prog->deps[i]->src]->dim+pos);
        }
    }

    for (i=0; i<prog->ntransdeps; i++) {
        assert(prog->transdeps[i] != NULL);
        if (prog->transdeps[i]->src == stmt->id) {
            pluto_constraints_remove_dim(prog->transdeps[i]->dpolytope, pos);
        }
        if (prog->transdeps[i]->dest == stmt->id) {
            // if (i==0)  printf("removing dim\n");
            pluto_constraints_remove_dim(prog->transdeps[i]->dpolytope,
                    prog->stmts[prog->transdeps[i]->src]->dim+pos);
        }
    }
}

void pluto_stmt_add_hyperplane(Stmt *stmt, PlutoHypType type, int pos)
{
    int i;

    assert(pos <= stmt->trans->nrows);

    pluto_matrix_add_row(stmt->trans, pos);

    stmt->hyp_types = realloc(stmt->hyp_types, 
            sizeof(int)*stmt->trans->nrows);
    for (i=stmt->trans->nrows-2; i>=pos; i--) {
        stmt->hyp_types[i+1] = stmt->hyp_types[i];
    }
    stmt->hyp_types[pos] = type;
}


void pluto_prog_add_hyperplane(PlutoProg *prog, int pos, PlutoHypType hyp_type)
{
    int i;

    prog->num_hyperplanes++;
    prog->hProps = (HyperplaneProperties *) realloc(prog->hProps, 
            prog->num_hyperplanes*sizeof(HyperplaneProperties));

    for (i=prog->num_hyperplanes-2; i>=pos; i--) {
        prog->hProps[i+1] = prog->hProps[i];
    }
    /* Initialize some */
    prog->hProps[pos].unroll = NO_UNROLL;
    prog->hProps[pos].prevec = 0;
    prog->hProps[pos].band_num = -1;
    prog->hProps[pos].dep_prop = UNKNOWN;
    prog->hProps[pos].type = hyp_type;
}


/* Pad statement transformations so that they all equal number
 * of rows */
void pluto_pad_stmt_transformations(PlutoProg *prog)
{
    int max_nrows, i, j, nstmts;

    nstmts = prog->nstmts;
    Stmt **stmts = prog->stmts;

    /* Pad all trans if necessary with zeros */
    max_nrows = 0;
    for (i=0; i<nstmts; i++)    {
        if (stmts[i]->trans != NULL)    {
            max_nrows = PLMAX(max_nrows, stmts[i]->trans->nrows);
        }
    }

    if (max_nrows >= 1) {
        for (i=0; i<nstmts; i++)    {
            if (stmts[i]->trans == NULL)    {
                stmts[i]->trans = pluto_matrix_alloc(max_nrows, 
                        stmts[i]->dim+prog->npar+1);
                stmts[i]->trans->nrows = 0;
            }

            int curr_rows = stmts[i]->trans->nrows;

            /* Add all zero rows */
            for (j=curr_rows; j<max_nrows; j++)    {
                pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
            }
        }

        int old_hyp_num = prog->num_hyperplanes;
        for (i=old_hyp_num; i<max_nrows; i++) {
            /* This is not really H_SCALAR, but this is the best we can do */
            pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);
        }
    }
}


/* Add statement to program; can't reuse arg stmt pointer any more */
void pluto_add_given_stmt(PlutoProg *prog, Stmt *stmt)
{
    prog->stmts = (Stmt **) realloc(prog->stmts, ((prog->nstmts+1)*sizeof(Stmt *)));

    stmt->id = prog->nstmts;

    prog->nvar = PLMAX(prog->nvar, stmt->dim);
    prog->stmts[prog->nstmts] = stmt;
    prog->nstmts++;

    pluto_pad_stmt_transformations(prog);

}



/* Create a statement and add it to the program
 * iterators: domain iterators
 * trans: schedule/transformation
 * domain: domain
 * text: statement text
 */
void pluto_add_stmt(PlutoProg *prog, 
        const PlutoConstraints *domain,
        const PlutoMatrix *trans,
        char ** iterators,
        const char *text,
        PlutoStmtType type)
{
    int i, nstmts;

    assert(trans != NULL);
    assert(trans->ncols == domain->ncols);

    nstmts = prog->nstmts;

    prog->stmts = (Stmt **) realloc(prog->stmts, ((nstmts+1)*sizeof(Stmt *)));

    Stmt *stmt = pluto_stmt_alloc(domain->ncols-prog->npar-1, domain, trans);

    stmt->id = nstmts;
    stmt->type = type;

    stmt->text = strdup(text);
    prog->nvar = PLMAX(prog->nvar, stmt->dim);

    for (i=0; i<stmt->dim; i++) {
        stmt->iterators[i] = strdup(iterators[i]);
    }

    prog->stmts[nstmts] = stmt;
    prog->nstmts++;

    pluto_pad_stmt_transformations(prog);
}


Dep *pluto_dep_alloc()
{
    Dep *dep = malloc(sizeof(Dep));

    dep->id = -1;
    dep->satvec = NULL;
    dep->depsat_poly = NULL;
    dep->satisfied = false;
    dep->satisfaction_level = -1;
    dep->dirvec = NULL;

    return dep;
}


Stmt *pluto_stmt_alloc(int dim, const PlutoConstraints *domain, 
        const PlutoMatrix *trans)
{
    int i;

    /* Have to provide a transformation */
    assert(trans != NULL);

    Stmt *stmt = (Stmt *) malloc(sizeof(Stmt));

    /* id will be assigned when added to PlutoProg */
    stmt->id = -1;
    stmt->dim = dim;
    stmt->dim_orig = dim;
    if (domain != NULL) {
        stmt->domain = pluto_constraints_dup(domain);
    }else{
        stmt->domain = NULL;
    }

    stmt->trans = pluto_matrix_dup(trans);

    stmt->hyp_types = malloc(stmt->trans->nrows*sizeof(int));
    for (i=0; i<stmt->trans->nrows; i++) {
        stmt->hyp_types[i] = H_LOOP;
    }

    stmt->text = NULL;
    stmt->tile =  1;
    stmt->num_tiled_loops = 0;
    stmt->reads = NULL;
    stmt->writes = NULL;
    stmt->nreads = 0;
    stmt->nwrites = 0;

    stmt->first_tile_dim = 0;
    stmt->last_tile_dim = -1;

    stmt->type = STMT_UNKNOWN;
    stmt->parent_compute_stmt = NULL;

    if (dim >= 1)   {
        stmt->is_orig_loop = (bool *) malloc(dim*sizeof(bool));
        stmt->iterators = (char **) malloc(sizeof(char *)*dim);
        for (i=0; i<stmt->dim; i++) {
            stmt->iterators[i] = NULL;
        }
    }else{
        stmt->is_orig_loop = NULL;
        stmt->iterators = NULL;
    }

    return stmt;
}


void pluto_access_free(PlutoAccess *acc)
{
    pluto_matrix_free(acc->mat);
    free(acc->name);
    free(acc);
}

void pluto_stmt_free(Stmt *stmt)
{
    int i, j;

    pluto_constraints_free(stmt->domain);

    pluto_matrix_free(stmt->trans);

    free(stmt->hyp_types);

    if (stmt->text != NULL) {
        free(stmt->text);
    }

    for (j=0; j<stmt->dim; j++)    {
        if (stmt->iterators[j] != NULL) {
            free(stmt->iterators[j]);
        }
    }

    /* If dim is zero, iterators, is_orig_loop are NULL */
    if (stmt->iterators != NULL)    {
        free(stmt->iterators);
        free(stmt->is_orig_loop);
    }

    PlutoAccess **writes = stmt->writes;
    PlutoAccess **reads = stmt->reads;

    if (writes != NULL) {
        for (i=0; i<stmt->nwrites; i++)   {
            pluto_access_free(writes[i]);
        }
        free(writes);
    }
    if (reads != NULL) {
        for (i=0; i<stmt->nreads; i++)   {
            pluto_access_free(reads[i]);
        }
        free(reads);
    }

    free(stmt);
}



/*char *get_data_extent(PlutoAccess *acc, char **params, int npars, int dim)
{
    return scoplib_symbol_table_get_bound(acc->symbol, dim, params, npars);
}*/

/* Get Alpha matrix (A matrix - INRIA transformation representation */
PlutoMatrix *get_alpha(const Stmt *stmt, const PlutoProg *prog)
{
    int r, c, i;

    PlutoMatrix *a;
    a = pluto_matrix_alloc(stmt->dim, stmt->dim);

    r=0;
    for (i=0; i<stmt->trans->nrows; i++)    {
        if (stmt->hyp_types[i] == H_LOOP || 
                stmt->hyp_types[i] == H_TILE_SPACE_LOOP) {
            for (c=0; c<stmt->dim; c++) {
                a->val[r][c] = stmt->trans->val[i][c];
            }
            r++;
            if (r==stmt->dim)   break;
        }
    }

    assert(r==stmt->dim);

    return a;
}


int pluto_is_hyperplane_scalar(const Stmt *stmt, int level)
{
    int j;

    assert(level <= stmt->trans->nrows-1);

    for (j=0; j<stmt->dim; j++) {
        if (stmt->trans->val[level][j] != 0) return 0;
    }

    return 1;
}


int pluto_is_hyperplane_loop(const Stmt *stmt, int level)
{
    return !pluto_is_hyperplane_scalar(stmt, level);
}

/* Get the remapping matrix: maps time iterators back to the domain 
 * iterators; divs: divisors for the rows */
PlutoMatrix *pluto_stmt_get_remapping(const Stmt *stmt, int **divs)
{
    int i, j, k, _lcm, factor1, npar;

    PlutoMatrix *remap, *trans;

    trans = stmt->trans;
    remap = pluto_matrix_dup(trans);

    npar = stmt->domain->ncols - stmt->dim - 1;

    *divs = malloc(sizeof(int)*(stmt->dim+npar+1));

    for (i=0; i<remap->nrows; i++)  {
        pluto_matrix_negate_row(remap, remap->nrows-1-i);
        pluto_matrix_add_col(remap, 0);
        remap->val[trans->nrows-1-i][0] = 1;
    }

    /* Bring the stmt iterators to the left */
    for (i=0; i<stmt->dim; i++)  {
        pluto_matrix_move_col(remap, remap->nrows+i, i);
    }

    assert(stmt->dim <= remap->nrows);

    for (i=0; i<stmt->dim; i++)  {
        // pluto_matrix_print(stdout, remap);
        if (remap->val[i][i] == 0) {
            for (k=i+1; k<remap->nrows; k++) {
                if (remap->val[k][i] != 0) break;
            }
            if (k<remap->nrows)    {
                pluto_matrix_interchange_rows(remap, i, k);
            }else{
                /* Can't associate domain iterator with time iterator */
                /* Shouldn't happen with a full-ranked transformation */
                printf("Can't associate domain iterator #%d with time iterators\n", i+1);
                pluto_matrix_print(stdout, remap);
                assert(0);
            }
        }
        //printf("after interchange %d\n", i); 
        //pluto_matrix_print(stdout, remap);
        assert(remap->val[i][i] != 0);
        for (k=i+1; k<remap->nrows; k++) {
            if (remap->val[k][i] == 0) continue;
            _lcm = lcm(remap->val[k][i], remap->val[i][i]);
            factor1 = _lcm/remap->val[k][i];
            for (j=i; j<remap->ncols; j++) {
                remap->val[k][j] = remap->val[k][j]*factor1
                    - remap->val[i][j]*(_lcm/remap->val[i][i]);
            }

        }
        //printf("after iteration %d\n", i); 
        //pluto_matrix_print(stdout, remap);
    }

    //pluto_matrix_print(stdout, remap);

    /* Solve upper triangular system now */
    for (i=stmt->dim-1; i>=0; i--)  {
        assert(remap->val[i][i] != 0);
        for (k=i-1; k>=0; k--) {
            if (remap->val[k][i] == 0) continue;
            _lcm = lcm(remap->val[k][i], remap->val[i][i]);
            factor1 = _lcm/remap->val[k][i];
            for (j=0; j<remap->ncols; j++) {
                remap->val[k][j] = remap->val[k][j]*(factor1) 
                    - remap->val[i][j]*(_lcm/remap->val[i][i]);
            }
        }
    }

    assert(remap->nrows >= stmt->dim);
    for (i=remap->nrows-1; i>=stmt->dim; i--) {
        pluto_matrix_remove_row(remap, remap->nrows-1);
    }
    // pluto_matrix_print(stdout, remap);

    for (i=0; i<stmt->dim; i++) {
        assert(remap->val[i][i] != 0);
        if (remap->val[i][i] <= -1) {
            pluto_matrix_negate_row(remap, i);
        }
        (*divs)[i] = abs(remap->val[i][i]);
    }
    // pluto_matrix_print(stdout, remap);

    for (i=0; i<stmt->dim; i++) {
        pluto_matrix_remove_col(remap, 0);
    }

    for (i=0; i<stmt->dim; i++) {
        pluto_matrix_negate_row(remap, i);
    }

    /* Identity for the parameter and constant part */
    for (i=0; i<npar+1; i++) {
        pluto_matrix_add_row(remap, remap->nrows);
        remap->val[remap->nrows-1][remap->ncols-npar-1+i] = 1;
        (*divs)[stmt->dim+i] = 1;
    }

    // printf("Remapping using new technique is\n");
    // pluto_matrix_print(stdout, remap);

    return remap;
}


void pluto_prog_params_print(const PlutoProg *prog)
{
    int i;
    for (i=0; i<prog->npar; i++) {
        printf("%s\n", prog->params[i]);
    }
}


/* Get new access function */
PlutoMatrix *pluto_get_new_access_func(const Stmt *stmt, 
        const PlutoMatrix *acc, int **divs) 
{
    PlutoMatrix *remap, *newacc;
    int r, c, npar, *remap_divs;

    npar = stmt->domain->ncols - stmt->dim - 1;
    *divs = malloc(sizeof(int)*acc->nrows);

    // printf("Old access function is \n");;
    // pluto_matrix_print(stdout, acc);;

    // printf("Stmt trans\n");
    // pluto_matrix_print(stdout, stmt->trans);

    remap = pluto_stmt_get_remapping(stmt, &remap_divs);
    // printf("Remapping matrix\n");
    // pluto_matrix_print(stdout, remap);
    //

    int _lcm = 1;
    for (r=0; r<remap->nrows; r++) {
        assert(remap_divs[r] != 0);
        _lcm = lcm(_lcm,remap_divs[r]);
    }
    for (r=0; r<remap->nrows; r++) {
        for (c=0; c<remap->ncols; c++) {
            remap->val[r][c] = (remap->val[r][c]*_lcm)/remap_divs[r];
        }
    }

    newacc = pluto_matrix_product(acc, remap);

    for (r=0; r<newacc->nrows; r++) {
        (*divs)[r] = _lcm;
    }

    // printf("New access function is \n");
    // pluto_matrix_print(stdout, newacc);

    assert(newacc->ncols = stmt->trans->nrows+npar+1);

    pluto_matrix_free(remap);
    free(remap_divs);

    return newacc;
}


/* Separates a list of statements at level 'level' */
void pluto_separate_stmts(PlutoProg *prog, Stmt **stmts, int num, int level)
{
    int i, nstmts, k;

    nstmts = prog->nstmts;

    // pluto_matrix_print(stdout, stmt->trans);
    for (i=0; i<nstmts; i++)    {
        pluto_stmt_add_hyperplane(prog->stmts[i], H_SCALAR, level);
    }
    // pluto_matrix_print(stdout, stmt->trans);
    for (k=0; k<num; k++)   {
        stmts[k]->trans->val[level][stmts[k]->trans->ncols-1] = 1+k;
    }

    pluto_prog_add_hyperplane(prog, level, H_SCALAR);
    prog->hProps[level].dep_prop = SEQ;
}


/* Separates a statement from the rest (places it later) at level 'level';
 * this is done by inserting a scalar dimension separating them */
void pluto_separate_stmt(PlutoProg *prog, const Stmt *stmt, int level)
{
    int i, nstmts;

    nstmts = prog->nstmts;

    // pluto_matrix_print(stdout, stmt->trans);
    for (i=0; i<nstmts; i++)    {
        pluto_stmt_add_hyperplane(prog->stmts[i], H_SCALAR, level);
    }
    // pluto_matrix_print(stdout, stmt->trans);
    stmt->trans->val[level][stmt->trans->ncols-1] = 1;

    pluto_prog_add_hyperplane(prog, level, H_SCALAR);
    prog->hProps[level].dep_prop = SEQ;
}

int pluto_stmt_is_member_of(const Stmt *s, Stmt **slist, int len)
{
    int i;
    for (i=0; i<len; i++) {
        if (s->id == slist[i]->id) return 1;
    }
    return 0;
}


int pluto_stmt_is_subset_of(Stmt **s1, int n1, Stmt **s2, int n2)
{
    int i;

    for (i=0; i<n1; i++) {
        if (!pluto_stmt_is_member_of(s1[i], s2, n2)) return 0;
    }

    return 1;
}

void add_if_new(PlutoAccess ***accs, int *num, PlutoAccess *new)
{
    int i;

    for (i=0; i<*num; i++) {
        if (!strcmp((*accs)[i]->name, new->name)) {
            break;
        }
    }

    if (i==*num) {
        *accs = realloc(*accs, (*num+1)*sizeof(PlutoAccess *));
        (*accs)[*num] = new;
        (*num)++;
    }
}


/* Get all write accesses in the program */
PlutoAccess **pluto_get_all_waccs(PlutoProg *prog, int *num)
{
    int i;

    PlutoAccess **accs = NULL;
    *num = 0;

    for (i=0; i<prog->nstmts; i++) {
        assert(prog->stmts[i]->nwrites == 1);
        add_if_new(&accs, num, prog->stmts[i]->writes[0]);
    }
    return accs;
}

/* Temporary data structure used inside extra_stmt_domains
 *
 * stmts points to the array of Stmts being constructed
 * index is the index of the next stmt in the array
 */
struct pluto_extra_stmt_info {
    Stmt **stmts;
    int index;
};

static int extract_basic_set(__isl_take isl_basic_set *bset, void *user)
{
    Stmt **stmts;
    Stmt *stmt;
    PlutoConstraints *bcst;
    struct pluto_extra_stmt_info *info;

    info = (struct pluto_extra_stmt_info *)user;

    stmts = info->stmts;
    stmt = stmts[info->index];

    bcst = isl_basic_set_to_pluto_constraints(bset);
    if (stmt->domain) {
        stmt->domain = pluto_constraints_unionize_simple(stmt->domain, bcst);
        pluto_constraints_free(bcst);
    }else{
        stmt->domain = bcst;
    }

    isl_basic_set_free(bset);
    return 0;
}

static int extract_stmt(__isl_take isl_set *set, void *user)
{
    int r;
    Stmt **stmts;
    int id;

    stmts = (Stmt **) user;

    int dim = isl_set_dim(set, isl_dim_all);
    int npar = isl_set_dim(set, isl_dim_param);
    PlutoMatrix *trans = pluto_matrix_alloc(dim-npar, dim+1);
    pluto_matrix_initialize(trans, 0);
    trans->nrows = 0;

    id = atoi(isl_set_get_tuple_name(set)+2);

    stmts[id] = pluto_stmt_alloc(dim-npar, NULL, trans);

    Stmt *stmt = stmts[id];
    stmt->type = ORIG;
    stmt->id = id;

    struct pluto_extra_stmt_info info = {stmts, id};
    r = isl_set_foreach_basic_set(set, &extract_basic_set, &info);

    pluto_matrix_free(trans);

    int j;
    for (j=0; j<stmt->dim; j++)  {
        stmt->is_orig_loop[j] = true;
    }

    isl_set_free(set);

    return r;
}

int extract_stmts(__isl_keep isl_union_set *domains, Stmt **stmts)
{
    isl_union_set_foreach_set(domains, &extract_stmt, stmts);

    return 0;
}

int pluto_get_max_ind_hyps_non_scalar(const PlutoProg *prog)
{
    int max, i;

    max = 0;

    for (i=0; i<prog->nstmts; i++) {
        max = PLMAX(max, pluto_stmt_get_num_ind_hyps_non_scalar(prog->stmts[i]));
    }

    return max;
}

int pluto_get_max_ind_hyps(const PlutoProg *prog)
{
    int max, i;

    max = 0;

    for (i=0; i<prog->nstmts; i++) {
        max = PLMAX(max, pluto_stmt_get_num_ind_hyps(prog->stmts[i]));
    }

    return max;
}

int pluto_stmt_get_num_ind_hyps_non_scalar(const Stmt *stmt)
{
    int isols, i,j=0;

    PlutoMatrix *tprime = pluto_matrix_dup(stmt->trans);

    /* Ignore padding dimensions, params, and constant part */
    for (i=stmt->dim_orig; i<stmt->trans->ncols; i++) {
        pluto_matrix_remove_col(tprime, stmt->dim_orig);
    }
    for (i=0; i<stmt->trans->nrows; i++) {
        if (stmt->hyp_types[i]==H_SCALAR) {   
            pluto_matrix_remove_row(tprime, i-j); 
            j++; 
        }
    }

    isols = pluto_matrix_get_rank(tprime);
    pluto_matrix_free(tprime);

    return isols;
}

int pluto_stmt_get_num_ind_hyps(const Stmt *stmt)
{
    int isols, i;

    PlutoMatrix *tprime = pluto_matrix_dup(stmt->trans);

    /* Ignore padding dimensions, params, and constant part */
    for (i=stmt->dim_orig; i<stmt->trans->ncols; i++) {
        pluto_matrix_remove_col(tprime, stmt->dim_orig);
    }

    isols = pluto_matrix_get_rank(tprime);
    pluto_matrix_free(tprime);

    return isols;
}

int pluto_transformations_full_ranked(PlutoProg *prog)
{
    int i;

    for (i=0; i<prog->nstmts; i++) {
        if (pluto_stmt_get_num_ind_hyps(prog->stmts[i]) < prog->stmts[i]->dim_orig) {
            return 0;
        }
    }

    return 1;
}

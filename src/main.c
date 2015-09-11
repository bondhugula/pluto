/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007--2008 Uday Kumar Bondhugula
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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pluto.h"

#include "clan/clan.h"
#include "candl/candl.h"

#include "math_support.h"
#include "post_transform.h"
#include "ddg.h"
#include "program.h"
#include "version.h"

PlutoOptions *options;

void usage_message();

int main(int argc, char *argv[])
{
    int i,count,j,k;

    FILE *src_fp;

    int option;
    int option_index = 0;

    char *srcFileName;

    FILE *cloogfp, *outfp;

    if (argc <= 1)  {
        usage_message();
        return 1;
    }

    options = pluto_options_alloc();

    const struct option pluto_options[] =
    {
        {"tile", no_argument, &options->tile, 1},
        {"notile", no_argument, &options->tile, 0},
        {"debug", no_argument, &options->debug, true},
        {"moredebug", no_argument, &options->moredebug, true},
        {"rar", no_argument, &options->rar, 1},
        {"identity", no_argument, &options->identity, 1},
        {"nofuse", no_argument, &options->fuse, NO_FUSE},
        {"maxfuse", no_argument, &options->fuse, MAXIMAL_FUSE},
        {"smartfuse", no_argument, &options->fuse, SMART_FUSE},
        {"parallel", no_argument, &options->parallel, 1},
        {"parallelize", no_argument, &options->parallel, 1},
        {"unroll", no_argument, &options->unroll, 1},
        {"nounroll", no_argument, &options->unroll, 0},
        {"polyunroll", no_argument, &options->polyunroll, 1},
        {"bee", no_argument, &options->bee, 1},
        {"ufactor", required_argument, 0, 'u'},
        {"prevector", no_argument, &options->prevector, 1},
        {"noprevector", no_argument, &options->prevector, 0},
        {"context", required_argument, 0, 'c'},
        {"cloogf", required_argument, 0, 'F'},
        {"cloogl", required_argument, 0, 'L'},
        {"ft", required_argument, 0, 'f'},
        {"lt", required_argument, 0, 'l'},
        {"multipipe", no_argument, &options->multipipe, 1},
        {"l2tile", no_argument, &options->l2tile, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"indent", no_argument, 0, 'i'},
        {"silent", no_argument, &options->silent, 1},
        {"lastwriter", no_argument, &options->lastwriter, 1},
        {"nobound", no_argument, &options->nobound, 1},
        {"scalpriv", no_argument, &options->scalpriv, 1},
        {"isldep", no_argument, &options->isldep, 1},
        {"isldepcompact", no_argument, &options->isldepcompact, 1},
        {"readscoplib", no_argument, &options->readscoplib, 1},
        {"islsolve", no_argument, &options->islsolve, 1},
        {0, 0, 0, 0}
    };


    /* Read command-line options */
    while (1) {
        option = getopt_long(argc, argv, "bhiqvf:l:F:L:c:o:", pluto_options,
                &option_index);

        if (option == -1)   {
            break;
        }

        switch (option) {
            case 0:
                break;
            case 'F':
                options->cloogf = atoi(optarg);
                break;
            case 'L':
                options->cloogl = atoi(optarg);
                break;
            case 'b':
                options->bee = 1;
                break;
            case 'c':
                options->context = atoi(optarg);
                break;
            case 'd':
                break;
            case 'f':
                options->ft = atoi(optarg);
                break;
            case 'g':
                break;
            case 'h':
                usage_message();
                return 2;
            case 'i':
                /* Handled in polycc */
                break;
            case 'l':
                options->lt = atoi(optarg);
                break;
            case 'm':
                break;
            case 'n':
                break;
            case 'o':
                options->out_file = strdup(optarg);
                break;
            case 'p':
                break;
            case 'q':
                options->silent = 1;
                break;
            case 's':
                break;
            case 'u':
                options->ufactor = atoi(optarg);
                break;
            case 'v':
                printf("PLUTO %s - An automatic parallelizer and locality optimizer\n\
Copyright (C) 2007--2008  Uday Kumar Bondhugula\n\
This is free software; see the source for copying conditions.  There is NO\n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n", PLUTO_VERSION);
                pluto_options_free(options);
                return 3;
            default:
                usage_message();
                pluto_options_free(options);
                return 4;
        }
    }


    if (optind <= argc-1)   {
        srcFileName = alloca(strlen(argv[optind])+1);
        strcpy(srcFileName, argv[optind]);
    }else{
        /* No non-option argument was specified */
        usage_message();
        pluto_options_free(options);
        return 5;
    }

    src_fp  = fopen(srcFileName, "r");

    if (!src_fp)   {
        fprintf(stderr, "pluto: error opening source file: '%s'\n", srcFileName);
        pluto_options_free(options);
        return 6;
    }

    /* Extract polyhedral representation from input program */
    scoplib_scop_p scop;

    clan_options_p clanOptions = clan_options_malloc();

    if (options->readscoplib) scop = scoplib_scop_read(src_fp);
    else scop = clan_scop_extract(src_fp, clanOptions);

    if (!scop || !scop->statement)   {
        fprintf(stderr, "Error extracting polyhedra from source file: \'%s'\n",
                srcFileName);
        pluto_options_free(options);
        return 7;
    }
    FILE *srcfp = fopen(".srcfilename", "w");
    if (srcfp)    {
        fprintf(srcfp, "%s\n", srcFileName);
        fclose(srcfp);
    }
        
    /* IF_DEBUG(clan_scop_print_dot_scop(stdout, scop, clanOptions)); */

double t_start,t_end;

bug1("before 'start'");
t_start = rtclock();
    /* Convert clan scop to Pluto program */
    PlutoProg *prog = scop_to_pluto_prog(scop, options);

cutAllSccs=0;

prog->stmts_dim = (int*)malloc(prog->nstmts*sizeof(int));
max_dim=0;
for(i=0;i<prog->nstmts;i++)
{
max_dim=max(max_dim,prog->stmts[i]->dim);
prog->stmts_dim[i]=prog->stmts[i]->dim;
}

    clan_options_free(clanOptions);

    /* Backup irregular program portion in .scop. */
    char* irroption = scoplib_scop_tag_content(scop, "<irregular>",
            "</irregular>");

    scoplib_scop_free(scop);

    IF_DEBUG2(pluto_deps_print(stdout, prog->deps, prog->ndeps));
    IF_DEBUG2(pluto_stmts_print(stdout, prog->stmts, prog->nstmts));

    /* Create the data dependence graph */
    prog->ddg = ddg_create(prog);
    ddg_compute_scc(prog);
//pluto_deps_print(stdout, prog->deps, prog->ndeps);


/* Code for part (c) starts ... */

for(i=0;i<prog->nstmts;i++)
prog->stmts[i]->orig_scc_id = prog->stmts[i]->scc_id;

int nrdeps=0;

for(i=prog->ndeps-1;i>=0;i--)
{
count = 0;
for(j=0;j<prog->deps[i]->src_acc->mat->nrows;j++)
{
for(k=0;k<prog->nvar;k++)
{
if(prog->deps[i]->src_acc->mat->val[j][k])
count++;
}
}
if(i==prog->ndeps-1 || !(pluto_domain_equality1(prog->deps[i]->src_acc->mat,prog->deps[i+1]->src_acc->mat) && pluto_domain_equality1(prog->deps[i]->dest_acc->mat,prog->deps[i+1]->dest_acc->mat) && (prog->deps[i]->src == prog->deps[i+1]->src) && (prog->deps[i]->dest == prog->deps[i+1]->dest) && (count < prog->stmts[prog->deps[i]->src]->dim) /*&& !strcmp(prog->deps[i]->src_acc->name,prog->deps[i+1]->src_acc->name)*/))
{
//if(prog->deps[i]->src <= prog->deps[i]->dest)
nrdeps++;
}
}

prog->nrdeps = nrdeps;
prog->rdeps = (int*)malloc(sizeof(int)*nrdeps); 

bug("%d",nrdeps);

for(i=prog->ndeps-1;i>=0;i--)
{
count=0;
for(j=0;j<prog->deps[i]->src_acc->mat->nrows;j++)
{
for(k=0;k<prog->nvar;k++)
{
if(prog->deps[i]->src_acc->mat->val[j][k])
count++;
}
}
//if(count < prog->stmts[prog->deps[i]->src]->dim) // that is if the source of the dep is an incomplete array
//{
if(i==prog->ndeps-1 || !(pluto_domain_equality1(prog->deps[i]->src_acc->mat,prog->deps[i+1]->src_acc->mat) && pluto_domain_equality1(prog->deps[i]->dest_acc->mat,prog->deps[i+1]->dest_acc->mat) && (prog->deps[i]->src == prog->deps[i+1]->src) && (prog->deps[i]->dest == prog->deps[i+1]->dest) && (count < prog->stmts[prog->deps[i]->src]->dim) /*&& !strcmp(prog->deps[i]->src_acc->name,prog->deps[i+1]->src_acc->name)*/))
{
//if(prog->deps[i]->src <= prog->deps[i]->dest)
prog->rdeps[--nrdeps]=i;
}
//}
}

/* Code for part (c) ends ...*/


count=1;
for(i=0;i<prog->nstmts-1;i++)
{
//pluto_constraints_print(stdout,prog->stmts[i]->domain);

/* counting the number of loops */
if(!pluto_domain_equality(prog->stmts[i]->domain,prog->stmts[i+1]->domain))
count++;

}


prog->nloops = count;

prog->loops = (int*) malloc(sizeof(int)*count); // loops[i] contains the ending location of each loop
count=0;

for(i=0;i<prog->nstmts-1;i++)
{

if(!pluto_domain_equality(prog->stmts[i]->domain,prog->stmts[i+1]->domain))
prog->loops[count++]=i;

}

prog->loops[count]=i;

/* printing loops */
for(i=0;i<prog->nloops;i++)
{
bug("loop %d: Statements: %d - %d",i,(i==0?0:prog->loops[i-1]+1),prog->loops[i]);
}

prog->distr = (int*) malloc(prog->nloops*sizeof(int));
for(i=0;i<prog->nloops;i++)
{
prog->distr[i]=0;
}



struct dist* d;
int count1,count2,num,num1=-1;
dist_ = (struct dist***) malloc(sizeof(struct dist**)*prog->nloops);
for(i=0;i<prog->nloops;i++)
{
dist_[i]=(struct dist**)malloc(sizeof(struct dist*)*prog->nloops);
for(j=0;j<prog->nloops;j++)
dist_[i][j]=NULL;
}

for(i=0;i<prog->ndeps;i++)
{
if(prog->deps[i]->src==prog->deps[i]->dest && prog->deps[i]->type==CANDL_WAW)
continue;
d = dist_[which_loop(prog,prog->deps[i]->src)][which_loop(prog,prog->deps[i]->dest)];

if(d==NULL)
{

dist_[which_loop(prog,prog->deps[i]->src)][which_loop(prog,prog->deps[i]->dest)] = (struct dist*) malloc(sizeof(struct dist));
d=dist_[which_loop(prog,prog->deps[i]->src)][which_loop(prog,prog->deps[i]->dest)];
d->dep=i;
d->next=NULL;

count2=0;
for(j=0;j<prog->deps[i]->dpolytope->nrows;j++)
{
count1=0;
for(k=0;k<prog->deps[i]->dpolytope->ncols-prog->npar-1;k++)
{
if(prog->deps[i]->dpolytope->val[j][k])
count1++;
}
if((count1>1) || prog->deps[i]->dpolytope->is_eq[j])
count2++;
}

d->value = pluto_constraints_alloc(count2,prog->deps[i]->dpolytope->ncols);
d->value->nrows=count2;

count2=0;
for(j=0;j<prog->deps[i]->dpolytope->nrows;j++)
{
count1=0;
for(k=0;k<prog->deps[i]->dpolytope->ncols-prog->npar-1;k++)
{
if(prog->deps[i]->dpolytope->val[j][k])
count1++;
}
if((count1>1) || prog->deps[i]->dpolytope->is_eq[j])
{
for(k=0;k<prog->deps[i]->dpolytope->ncols;k++)
d->value->val[count2][k]=prog->deps[i]->dpolytope->val[j][k];
count2++;
}
}

bug("dep: %d (%d, %d)",i,prog->deps[i]->src,prog->deps[i]->dest);
//pluto_constraints_print(stdout, d->value);
} // end if

else
{

while(d->next!=NULL)
d=d->next;

d->next = (struct dist*) malloc(sizeof(struct dist));
d->next->dep = i;
d->next->next = NULL;

count2=0;
for(j=0;j<prog->deps[i]->dpolytope->nrows;j++)
{
count1=0;
for(k=0;k<prog->deps[i]->dpolytope->ncols-prog->npar-1;k++)
{
if(prog->deps[i]->dpolytope->val[j][k])
count1++;
}
if((count1>1) || prog->deps[i]->dpolytope->is_eq[j])
count2++;
}

d->next->value = pluto_constraints_alloc(count2,prog->deps[i]->dpolytope->ncols);
d->next->value->nrows=count2;

count2=0;
for(j=0;j<prog->deps[i]->dpolytope->nrows;j++)
{
count1=0;
for(k=0;k<prog->deps[i]->dpolytope->ncols-prog->npar-1;k++)
{
if(prog->deps[i]->dpolytope->val[j][k])
count1++;
}
if((count1>1) || prog->deps[i]->dpolytope->is_eq[j])
{
for(k=0;k<prog->deps[i]->dpolytope->ncols;k++)
d->next->value->val[count2][k]=prog->deps[i]->dpolytope->val[j][k];
count2++;
}
}

struct dist* d1 = dist_[which_loop(prog,prog->deps[i]->src)][which_loop(prog,prog->deps[i]->dest)];

while(d1!=d->next)
{
if(pluto_domain_equality(d1->value,d->next->value) /*&& pluto_domain_equality1(prog->stmts[prog->deps[d1->dep]->src]->writes[0]->mat,prog->stmts[prog->deps[d->next->dep]->src]->writes[0]->mat)*/)
break;
d1=d1->next;
}


if(d1!=d->next)
{
free(d->next);  
d->next=NULL;
}

} // else

}

//for(i=0;i<prog->nstmts;i++)
//pluto_constraints_print(stdout, prog->stmts[i]->domain);
bug("Original");
//pluto_deps_print(stdout, prog->deps,prog->ndeps);

marker = fopen("marker","w");
keys = fopen("keys","w");
skipdeps = fopen("skipdeps.txt","w");
FILE* marker2;
prog->nrstmts=0;

for(i=0;i<prog->nloops;i++)
{
for(j=0;j<prog->nloops;j++)
{
d = dist_[i][j];
while(d!=NULL)
{
fprintf(marker,"%d\n",d->dep);

//// removing inter-scc deps on scalars and 1Ds
//if(prog->stmts[prog->deps[d->dep]->src]->scc_id!=prog->stmts[prog->deps[d->dep]->dest]->scc_id && prog->deps[d->dep]->src_acc->mat->nrows<=(max_dim-2))
//{
//if(pluto_domain_equality2(prog->stmts[prog->deps[d->dep]->src]->domain,prog->stmts[prog->deps[d->dep]->dest]->domain,2,prog->stmts[prog->deps[d->dep]->src]->dim,prog->stmts[prog->deps[d->dep]->dest]->dim))
//{
//bug("cutting inter-scc dep: %d",d->dep);
//pluto_constraints_add_equality(prog->deps[d->dep]->dpolytope,prog->deps[d->dep]->dpolytope->nrows);
//prog->deps[d->dep]->dpolytope->val[prog->deps[d->dep]->dpolytope->nrows-1][0]=1;
//prog->deps[d->dep]->dpolytope->val[prog->deps[d->dep]->dpolytope->nrows-1][prog->stmts[prog->deps[d->dep]->src]->dim]=-1;
//}
//else
//{
//bug("cutting between sccs containing same scalar: %d",d->dep);
//cut_between_sccs(prog,prog->ddg,prog->stmts[prog->deps[d->dep]->src]->scc_id,prog->stmts[prog->deps[d->dep]->dest]->scc_id);
//}
//}

bug("dep-stmts: %d - %d",prog->deps[d->dep]->src,prog->deps[d->dep]->dest);
d=d->next;
}
}
}
fclose(marker);


//marker = fopen("marker","r");
//
//fscanf(marker, "%d",&num);
//marker2 = fopen("rstmts","w");
//fprintf(marker2,"%d\n",prog->deps[num]->src);
//prog->nrstmts++;
//fclose(marker2);
//
//while(!feof(marker)) {
//fscanf(marker, "%d", &num);
//
//marker2 = fopen("rstmts","r");
//
//while(!feof(marker2)) {
//fscanf(marker2, "%d", &num1);
//if(prog->deps[num]->src==num1)
//break;
//}
//
//fclose(marker2);
//
//if(prog->deps[num]->src!=num1)
//{
//marker2 = fopen("rstmts","a");
//fprintf(marker2, "%d\n",prog->deps[num]->src);
//fclose(marker2);
//prog->nrstmts++;
//}
//
//marker2 = fopen("rstmts","r");
//
//while(!feof(marker2)) {
//fscanf(marker2, "%d", &num1);
//if(prog->deps[num]->dest==num1)
//break;
//}
//
//fclose(marker2);
//
//if(prog->deps[num]->dest!=num1)
//{
//marker2 = fopen("rstmts","a");
//fprintf(marker2, "%d\n",prog->deps[num]->dest);
//fclose(marker2);
//prog->nrstmts++;
//}
//
//}
//fclose(marker);

//bug("nrstmts: %d",prog->nrstmts);

prog->keys_ = (int* ) malloc(prog->nloops*sizeof(int));

for(i=0;i<prog->nloops;i++)
{
marker = fopen("marker","r");
while(!feof(marker)) {
fscanf(marker, "%d", &num);
if(which_loop(prog,prog->deps[num]->src)==i)
{
bug("%d",num);
prog->keys_[i]=prog->deps[num]->src;
fprintf(keys,"%d\n",prog->deps[num]->src);
break;
}
}
fclose(marker);
}

fclose(keys);

CST_WIDTH= prog->npar+1+prog->nloops*(prog->nvar+1)+1;

for(i=0;i<prog->ndeps;i++)
{
marker = fopen("marker","r");

while(!feof(marker)) {
fscanf(marker, "%d", &num);
if(num==i)
break;
}

if(feof(marker))
fprintf(skipdeps,"%d\n",i);
fclose(marker);
}

i=0;j=0;

while(i<prog->ndeps)
{
while(i<prog->rdeps[j])
{
fprintf(skipdeps,"%d\n",i);
i++;
}
j++;
i++;
}


fclose(skipdeps);

    int dim_sum=0;
    for (i=0; i<prog->nstmts; i++) {
        dim_sum += prog->stmts[i]->dim;
    }

    /* Make options consistent */
    if (options->multipipe == 1 && options->parallel == 0)    {
        fprintf(stdout, "Warning: multipipe needs parallel to be on; turning on parallel\n");
        options->parallel = 1;
    }

    /* Disable pre-vectorization if tile is not on */
    if (options->tile == 0 && options->prevector == 1) {
        /* If code will not be tiled, pre-vectorization does not make
         * sense */
        if (!options->silent)   {
            fprintf(stdout, "[Pluto] Warning: pre-vectorization does not fit (--tile is off)\n");
        }
        options->prevector = 0;
    }

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Number of statements: %d\n", prog->nstmts);
        fprintf(stdout, "[Pluto] Total number of loops: %d\n", dim_sum);
        fprintf(stdout, "[Pluto] Number of deps: %d\n", prog->ndeps);
        fprintf(stdout, "[Pluto] Maximum domain dimensionality: %d\n", prog->nvar);
        fprintf(stdout, "[Pluto] Number of parameters: %d\n", prog->npar);
    }

    /* Auto transformation */
    if (!options->identity) {
noOuterLoop = 0;
        pluto_auto_transform(prog);
t_end = rtclock();
bug1("after 'end';\n Time taken: %0.6lfs\n", t_end - t_start);
bug("Max_rows_fme: %d",max_rows_fme);
pluto_constraints_set_var(prog->context,0,102);
pluto_constraints_set_var(prog->context,1,102);
pluto_constraints_set_var(prog->context,2,102);

count = 0;

for(i=0;i<prog->ndeps;i++)
{
bug("%d",i);
skipdeps = fopen("skipdeps.txt","r");
while(!feof(skipdeps))
{
fscanf(skipdeps,"%d",&num);
if(num==i)
break;
}
if(feof(skipdeps))
count++;
fclose(skipdeps);
}

bug1("The number of dependences analyzed, molecules: %d %d",count, prog->nloops);
    }
    pluto_detect_transformation_properties(prog);

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Affine transformations [<iter coeff's> <const>]\n\n");
    }

    /* Print out transformations */
    if (!options->silent)   {
        pluto_transformations_pretty_print(prog);
        pluto_print_hyperplane_properties(prog);

    }
    if (options->moredebug) {
        pluto_print_dep_directions(prog->deps, prog->ndeps, prog->num_hyperplanes);
    }

    if (options->tile)   {
        pluto_tile(prog);
    }

    if (options->parallel)   {
        int outermostBandStart, outermostBandEnd;
        getOutermostTilableBand(prog, &outermostBandStart, &outermostBandEnd);

        /* Obtain pipelined parallelization by skewing the tile space */
        bool retval = create_tile_schedule(prog, outermostBandStart, outermostBandEnd);

        /* Even if the user hasn't supplied --tile and there is only pipelined
         * parallelism, we will warn the user, but anyway do fine-grained 
         * parallelization
         */
        if (retval && options->tile == 0)   {
            printf("WARNING: pipelined parallelism exists and --tile is not used.\n");
            printf("use --tile option for better parallelization \n");
            pluto_print_hyperplane_properties(prog);
        }
    }

    if (options->prevector) {
        pre_vectorize(prog);
    }else{
        /* Create an empty .vectorize file */
        fopen(".vectorize", "w");
    }

    if (options->tile && !options->silent)  {
        fprintf(stdout, "[Pluto] After tiling:\n");
        pluto_transformations_pretty_print(prog);
        pluto_print_hyperplane_properties(prog);
    }

    if (options->parallel)  {
        /* Generate meta info for insertion of OpenMP pragmas */
        pluto_omp_parallelize(prog);
    }

    if (options->moredebug) {
        pluto_detect_transformation_properties(prog);
        pluto_print_dep_directions(prog->deps, prog->ndeps, prog->num_hyperplanes);
    }

    if (options->unroll || options->polyunroll)    {
        /* Will generate a .unroll file */
        /* plann/plorc needs a .params */
        FILE *paramsFP = fopen(".params", "w");
        if (paramsFP)   {
            int i;
            for (i=0; i<prog->npar; i++)  {
                fprintf(paramsFP, "%s\n", prog->params[i]);
            }
            fclose(paramsFP);
        }
        detect_mark_unrollable_loops(prog);
    }

    if (options->polyunroll)    {
        /* Experimental */
        for (i=0; i<prog->num_hyperplanes; i++)   {
            if (prog->hProps[i].unroll)  {
                unroll_phis(prog, i, options->ufactor);
            }
        }
    }


    /* NO MORE TRANSFORMATIONS BEYOND THIS POINT */
    /* Since meta info about loops
     * is printed to be processed by scripts - if transformations are
     * performed, changed loop order/iterator names will be missed  */
    gen_unroll_file(prog);
    gen_vecloop_file(prog);

    char *outFileName;
    char *cloogFileName;
    if (options->out_file == NULL)  {
        /* Get basename, remove .c extension and append a new one */
        char *basec, *bname;
        basec = strdup(srcFileName);
        bname = basename(basec);

        /* max size when tiled.* */
        outFileName = alloca(strlen(bname)+strlen(".pluto.c")+1);
        cloogFileName = alloca(strlen(bname)+strlen(".pluto.cloog")+1);

        if (strlen(bname) >= 2 && !strcmp(bname+strlen(bname)-2, ".c")) {
            strncpy(outFileName, bname, strlen(bname)-2);
            strncpy(cloogFileName, bname, strlen(bname)-2);
            outFileName[strlen(bname)-2] = '\0';
            cloogFileName[strlen(bname)-2] = '\0';
        }else{
            strcpy(outFileName, bname);
            strcpy(cloogFileName, bname);
        }
        strcat(outFileName, ".pluto.c");
        free(basec);
    }else{
        outFileName = options->out_file;
        cloogFileName = alloca(strlen(options->out_file)+1);
        strcpy(cloogFileName, options->out_file);
    }

    strcat(cloogFileName, ".pluto.cloog");

    outfp = fopen(outFileName, "w");
    cloogfp = fopen(cloogFileName, "w+");

    if (!cloogfp)   {
        fprintf(stderr, "Can't open .cloog file: %s\n", cloogFileName);
        pluto_options_free(options);
        pluto_prog_free(prog);
        return 8;
    }

    /* Generate .cloog file */
    pluto_gen_cloog_file(cloogfp, prog);
    /* Add the <irregular> tag from clan, if any */
    if (irroption != NULL) {
        fprintf(cloogfp, "<irregular>\n%s\n</irregular>\n\n", irroption);
        free(irroption);
    }

    /* Generate code using Cloog and add necessary stuff before/after code */
    rewind(cloogfp);

    if (!outfp) {
        fprintf(stderr, "Can't open file %s for writing\n", outFileName);
        pluto_options_free(options);
        pluto_prog_free(prog);
        return 9;
    }

    /* Generate code using Cloog and add necessary stuff before/after code */
    pluto_multicore_codegen(cloogfp, outfp, prog);

    FILE *tmpfp = fopen(".outfilename", "w");
    if (tmpfp)    {
        fprintf(tmpfp, "%s\n", outFileName);
        fclose(tmpfp);
    }

    fclose(cloogfp);
    fclose(outfp);

    pluto_options_free(options);

    pluto_prog_free(prog);

    return 0;
}


void usage_message(void)
{
    fprintf(stdout, "Usage: polycc <input.c> [options] [-o output]\n");
    fprintf(stdout, "\nOptions:\n");
    fprintf(stdout, "       --tile                 Tile for locality\n");
    fprintf(stdout, "       --parallel             Automatically parallelize using OpenMP pragmas\n");
    fprintf(stdout, "       | --parallelize\n");
    fprintf(stdout, "       --l2tile               Tile a second time (typically for L2 cache) - disabled by default \n");
    fprintf(stdout, "       --multipipe            Extract two degrees of pipelined parallelism if possible;\n");
    fprintf(stdout, "                                 by default one degree is extracted (if it exists)\n");
    fprintf(stdout, "       --rar                  Consider RAR dependences too (disabled by default)\n");
    fprintf(stdout, "       --[no]unroll           Unroll-jam (disabled by default)\n");
    fprintf(stdout, "       --ufactor=<factor>     Unroll-jam factor (default is 8)\n");
    fprintf(stdout, "       --[no]prevector        Make code amenable to compiler auto-vectorization (with ICC) - enabled by default\n");
    fprintf(stdout, "       --context=<context>    Parameters are at least as much as <context>\n");
    fprintf(stdout, "       --isldep               Use ISL-based dependence tester\n");
    fprintf(stdout, "       --islsolve             Use ISL as ilp solver\n");
    fprintf(stdout, "       --readscoplib          Read input from a scoplib file\n");
    fprintf(stdout, "       --lastwriter           Work with refined dependences (last conflicting access is computed for RAW/WAW)\n");
    fprintf(stdout, "       --bee                  Generate pragmas for Bee+Cl@k\n\n");
    fprintf(stdout, "       --indent  | -i         Indent generated code (disabled by default)\n");
    fprintf(stdout, "       --silent  | -q         Silent mode; no output as long as everything goes fine (disabled by default)\n");
    fprintf(stdout, "       --help    | -h         Print this help menu\n");
    fprintf(stdout, "       --version | -v         Display version number\n");
    fprintf(stdout, "\n   Fusion                Options to control fusion heuristic\n");
    fprintf(stdout, "       --nofuse               Do not fuse across SCCs of data dependence graph\n");
    fprintf(stdout, "       --maxfuse              Maximal fusion\n");
    fprintf(stdout, "       --smartfuse [default]  Heuristic (in between nofuse and maxfuse)\n");
    fprintf(stdout, "\n   Debugging\n");
    fprintf(stdout, "       --debug        Verbose output\n");
    fprintf(stdout, "       --moredebug    More verbose output\n");
    fprintf(stdout, "\nTo report bugs, please send an email to <pluto-development@googlegroups.com>\n\n");
}

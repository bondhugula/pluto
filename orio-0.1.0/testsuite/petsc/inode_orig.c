#define PETSCMAT_DLL

/*
  This file provides high performance routines for the Inode format (compressed sparse row)
  by taking advantage of rows with identical nonzero structure (I-nodes).
*/
#include "src/mat/impls/aij/seq/aij.h"

#undef __FUNCT__
#define __FUNCT__ "Mat_CreateColInode"
static PetscErrorCode Mat_CreateColInode(Mat A,PetscInt* size,PetscInt ** ns) {
    Mat_SeqAIJ      *a = (Mat_SeqAIJ*)A->data;
    PetscErrorCode ierr;
    PetscInt       i,count,m,n,min_mn,*ns_row,*ns_col;

    PetscFunctionBegin;
    n      = A->cmap.n;
    m      = A->rmap.n;
    ns_row = a->inode.size;

    min_mn = (m < n) ? m : n;
    if (!ns) {
        for (count=0,i=0; count<min_mn; count+=ns_row[i],i++);
        for(; count+1 < n; count++,i++);
        if (count < n)  {
            i++;
        }
        *size = i;
        PetscFunctionReturn(0);
    }
    ierr = PetscMalloc((n+1)*sizeof(PetscInt),&ns_col);
    CHKERRQ(ierr);

    /* Use the same row structure wherever feasible. */
    for (count=0,i=0; count<min_mn; count+=ns_row[i],i++) {
        ns_col[i] = ns_row[i];
    }

    /* if m < n; pad up the remainder with inode_limit */
    for(; count+1 < n; count++,i++) {
        ns_col[i] = 1;
    }
    /* The last node is the odd ball. padd it up with the remaining rows; */
    if (count < n)  {
        ns_col[i] = n - count;
        i++;
    } else if (count > n) {
        /* Adjust for the over estimation */
        ns_col[i-1] += n - count;
    }
    *size = i;
    *ns   = ns_col;
    PetscFunctionReturn(0);
}


/*
      This builds symmetric version of nonzero structure,
*/
#undef __FUNCT__
#define __FUNCT__ "MatGetRowIJ_Inode_Symmetric"
static PetscErrorCode MatGetRowIJ_Inode_Symmetric(Mat A,PetscInt *iia[],PetscInt *jja[],PetscInt ishift,PetscInt oshift) {
    Mat_SeqAIJ      *a = (Mat_SeqAIJ*)A->data;
    PetscErrorCode ierr;
    PetscInt       *work,*ia,*ja,*j,nz,nslim_row,nslim_col,m,row,col,*jmax,n;
    PetscInt       *tns,*tvc,*ns_row = a->inode.size,*ns_col,nsz,i1,i2,*ai= a->i,*aj = a->j;

    PetscFunctionBegin;
    nslim_row = a->inode.node_count;
    m         = A->rmap.n;
    n         = A->cmap.n;
    if (m != n) SETERRQ(PETSC_ERR_SUP,"MatGetRowIJ_Inode_Symmetric: Matrix should be square");

    /* Use the row_inode as column_inode */
    nslim_col = nslim_row;
    ns_col    = ns_row;

    /* allocate space for reformated inode structure */
    ierr = PetscMalloc((nslim_col+1)*sizeof(PetscInt),&tns);
    CHKERRQ(ierr);
    ierr = PetscMalloc((n+1)*sizeof(PetscInt),&tvc);
    CHKERRQ(ierr);
    for (i1=0,tns[0]=0; i1<nslim_col; ++i1) tns[i1+1] = tns[i1]+ ns_row[i1];

    for (i1=0,col=0; i1<nslim_col; ++i1) {
        nsz = ns_col[i1];
        for (i2=0; i2<nsz; ++i2,++col)
            tvc[col] = i1;
    }
    /* allocate space for row pointers */
    ierr = PetscMalloc((nslim_row+1)*sizeof(PetscInt),&ia);
    CHKERRQ(ierr);
    *iia = ia;
    ierr = PetscMemzero(ia,(nslim_row+1)*sizeof(PetscInt));
    CHKERRQ(ierr);
    ierr = PetscMalloc((nslim_row+1)*sizeof(PetscInt),&work);
    CHKERRQ(ierr);

    /* determine the number of columns in each row */
    ia[0] = oshift;
    for (i1=0,row=0 ; i1<nslim_row; row+=ns_row[i1],i1++) {

        j    = aj + ai[row] + ishift;
        jmax = aj + ai[row+1] + ishift;
        i2   = 0;
        col  = *j++ + ishift;
        i2   = tvc[col];
        while (i2<i1 && j<jmax) { /* 1.[-xx-d-xx--] 2.[-xx-------],off-diagonal elemets */
            ia[i1+1]++;
            ia[i2+1]++;
            i2++;                     /* Start col of next node */
            while(((col=*j+ishift)<tns[i2]) && (j<jmax)) ++j;
            i2 = tvc[col];
        }
        if(i2 == i1) ia[i2+1]++;    /* now the diagonal element */
    }

    /* shift ia[i] to point to next row */
    for (i1=1; i1<nslim_row+1; i1++) {
        row        = ia[i1-1];
        ia[i1]    += row;
        work[i1-1] = row - oshift;
    }

    /* allocate space for column pointers */
    nz   = ia[nslim_row] + (!ishift);
    ierr = PetscMalloc(nz*sizeof(PetscInt),&ja);
    CHKERRQ(ierr);
    *jja = ja;

    /* loop over lower triangular part putting into ja */
    for (i1=0,row=0; i1<nslim_row; row += ns_row[i1],i1++) {
        j    = aj + ai[row] + ishift;
        jmax = aj + ai[row+1] + ishift;
        i2   = 0;                     /* Col inode index */
        col  = *j++ + ishift;
        i2   = tvc[col];
        while (i2<i1 && j<jmax) {
            ja[work[i2]++] = i1 + oshift;
            ja[work[i1]++] = i2 + oshift;
            ++i2;
            while(((col=*j+ishift)< tns[i2])&&(j<jmax)) ++j; /* Skip rest col indices in this node */
            i2 = tvc[col];
        }
        if (i2 == i1) ja[work[i1]++] = i2 + oshift;

    }
    ierr = PetscFree(work);
    CHKERRQ(ierr);
    ierr = PetscFree(tns);
    CHKERRQ(ierr);
    ierr = PetscFree(tvc);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/*
      This builds nonsymmetric version of nonzero structure,
*/
#undef __FUNCT__
#define __FUNCT__ "MatGetRowIJ_Inode_Nonsymmetric"
static PetscErrorCode MatGetRowIJ_Inode_Nonsymmetric(Mat A,PetscInt *iia[],PetscInt *jja[],PetscInt ishift,PetscInt oshift) {
    Mat_SeqAIJ      *a = (Mat_SeqAIJ*)A->data;
    PetscErrorCode ierr;
    PetscInt       *work,*ia,*ja,*j,nz,nslim_row,n,row,col,*ns_col,nslim_col;
    PetscInt       *tns,*tvc,*ns_row = a->inode.size,nsz,i1,i2,*ai= a->i,*aj = a->j;

    PetscFunctionBegin;
    nslim_row = a->inode.node_count;
    n         = A->cmap.n;

    /* Create The column_inode for this matrix */
    ierr = Mat_CreateColInode(A,&nslim_col,&ns_col);
    CHKERRQ(ierr);

    /* allocate space for reformated column_inode structure */
    ierr = PetscMalloc((nslim_col +1)*sizeof(PetscInt),&tns);
    CHKERRQ(ierr);
    ierr = PetscMalloc((n +1)*sizeof(PetscInt),&tvc);
    CHKERRQ(ierr);
    for (i1=0,tns[0]=0; i1<nslim_col; ++i1) tns[i1+1] = tns[i1] + ns_col[i1];

    for (i1=0,col=0; i1<nslim_col; ++i1) {
        nsz = ns_col[i1];
        for (i2=0; i2<nsz; ++i2,++col)
            tvc[col] = i1;
    }
    /* allocate space for row pointers */
    ierr = PetscMalloc((nslim_row+1)*sizeof(PetscInt),&ia);
    CHKERRQ(ierr);
    *iia = ia;
    ierr = PetscMemzero(ia,(nslim_row+1)*sizeof(PetscInt));
    CHKERRQ(ierr);
    ierr = PetscMalloc((nslim_row+1)*sizeof(PetscInt),&work);
    CHKERRQ(ierr);

    /* determine the number of columns in each row */
    ia[0] = oshift;
    for (i1=0,row=0; i1<nslim_row; row+=ns_row[i1],i1++) {
        j   = aj + ai[row] + ishift;
        col = *j++ + ishift;
        i2  = tvc[col];
        nz  = ai[row+1] - ai[row];
        while (nz-- > 0) {           /* off-diagonal elemets */
            ia[i1+1]++;
            i2++;                     /* Start col of next node */
            while (((col = *j++ + ishift) < tns[i2]) && nz > 0) {
                nz--;
            }
            if (nz > 0) i2 = tvc[col];
        }
    }

    /* shift ia[i] to point to next row */
    for (i1=1; i1<nslim_row+1; i1++) {
        row        = ia[i1-1];
        ia[i1]    += row;
        work[i1-1] = row - oshift;
    }

    /* allocate space for column pointers */
    nz   = ia[nslim_row] + (!ishift);
    ierr = PetscMalloc(nz*sizeof(PetscInt),&ja);
    CHKERRQ(ierr);
    *jja = ja;

    /* loop over matrix putting into ja */
    for (i1=0,row=0; i1<nslim_row; row+=ns_row[i1],i1++) {
        j   = aj + ai[row] + ishift;
        i2  = 0;                     /* Col inode index */
        col = *j++ + ishift;
        i2  = tvc[col];
        nz  = ai[row+1] - ai[row];
        while (nz-- > 0) {
            ja[work[i1]++] = i2 + oshift;
            ++i2;
            while(((col = *j++ + ishift) < tns[i2]) && nz > 0) {
                nz--;
            }
            if (nz > 0) i2 = tvc[col];
        }
    }
    ierr = PetscFree(ns_col);
    CHKERRQ(ierr);
    ierr = PetscFree(work);
    CHKERRQ(ierr);
    ierr = PetscFree(tns);
    CHKERRQ(ierr);
    ierr = PetscFree(tvc);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatGetRowIJ_Inode"
static PetscErrorCode MatGetRowIJ_Inode(Mat A,PetscInt oshift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt *ja[],PetscTruth *done) {
    Mat_SeqAIJ      *a = (Mat_SeqAIJ*)A->data;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    *n     = a->inode.node_count;
    if (!ia) PetscFunctionReturn(0);
    if (!blockcompressed) {
        ierr = MatGetRowIJ_SeqAIJ(A,oshift,symmetric,blockcompressed,n,ia,ja,done);
        CHKERRQ(ierr);;
    } else if (symmetric) {
        ierr = MatGetRowIJ_Inode_Symmetric(A,ia,ja,0,oshift);
        CHKERRQ(ierr);
    } else {
        ierr = MatGetRowIJ_Inode_Nonsymmetric(A,ia,ja,0,oshift);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatRestoreRowIJ_Inode"
static PetscErrorCode MatRestoreRowIJ_Inode(Mat A,PetscInt oshift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt *ja[],PetscTruth *done) {
    PetscErrorCode ierr;

    PetscFunctionBegin;
    if (!ia) PetscFunctionReturn(0);

    if (!blockcompressed) {
        ierr = MatRestoreRowIJ_SeqAIJ(A,oshift,symmetric,blockcompressed,n,ia,ja,done);
        CHKERRQ(ierr);;
    } else {
        ierr = PetscFree(*ia);
        CHKERRQ(ierr);
        ierr = PetscFree(*ja);
        CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

/* ----------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "MatGetColumnIJ_Inode_Nonsymmetric"
static PetscErrorCode MatGetColumnIJ_Inode_Nonsymmetric(Mat A,PetscInt *iia[],PetscInt *jja[],PetscInt ishift,PetscInt oshift) {
    Mat_SeqAIJ      *a = (Mat_SeqAIJ*)A->data;
    PetscErrorCode ierr;
    PetscInt       *work,*ia,*ja,*j,nz,nslim_row, n,row,col,*ns_col,nslim_col;
    PetscInt       *tns,*tvc,*ns_row = a->inode.size,nsz,i1,i2,*ai= a->i,*aj = a->j;

    PetscFunctionBegin;
    nslim_row = a->inode.node_count;
    n         = A->cmap.n;

    /* Create The column_inode for this matrix */
    ierr = Mat_CreateColInode(A,&nslim_col,&ns_col);
    CHKERRQ(ierr);

    /* allocate space for reformated column_inode structure */
    ierr = PetscMalloc((nslim_col + 1)*sizeof(PetscInt),&tns);
    CHKERRQ(ierr);
    ierr = PetscMalloc((n + 1)*sizeof(PetscInt),&tvc);
    CHKERRQ(ierr);
    for (i1=0,tns[0]=0; i1<nslim_col; ++i1) tns[i1+1] = tns[i1] + ns_col[i1];

    for (i1=0,col=0; i1<nslim_col; ++i1) {
        nsz = ns_col[i1];
        for (i2=0; i2<nsz; ++i2,++col)
            tvc[col] = i1;
    }
    /* allocate space for column pointers */
    ierr = PetscMalloc((nslim_col+1)*sizeof(PetscInt),&ia);
    CHKERRQ(ierr);
    *iia = ia;
    ierr = PetscMemzero(ia,(nslim_col+1)*sizeof(PetscInt));
    CHKERRQ(ierr);
    ierr = PetscMalloc((nslim_col+1)*sizeof(PetscInt),&work);
    CHKERRQ(ierr);

    /* determine the number of columns in each row */
    ia[0] = oshift;
    for (i1=0,row=0; i1<nslim_row; row+=ns_row[i1],i1++) {
        j   = aj + ai[row] + ishift;
        col = *j++ + ishift;
        i2  = tvc[col];
        nz  = ai[row+1] - ai[row];
        while (nz-- > 0) {           /* off-diagonal elemets */
            /* ia[i1+1]++; */
            ia[i2+1]++;
            i2++;
            while (((col = *j++ + ishift) < tns[i2]) && nz > 0) {
                nz--;
            }
            if (nz > 0) i2 = tvc[col];
        }
    }

    /* shift ia[i] to point to next col */
    for (i1=1; i1<nslim_col+1; i1++) {
        col        = ia[i1-1];
        ia[i1]    += col;
        work[i1-1] = col - oshift;
    }

    /* allocate space for column pointers */
    nz   = ia[nslim_col] + (!ishift);
    ierr = PetscMalloc(nz*sizeof(PetscInt),&ja);
    CHKERRQ(ierr);
    *jja = ja;

    /* loop over matrix putting into ja */
    for (i1=0,row=0; i1<nslim_row; row+=ns_row[i1],i1++) {
        j   = aj + ai[row] + ishift;
        i2  = 0;                     /* Col inode index */
        col = *j++ + ishift;
        i2  = tvc[col];
        nz  = ai[row+1] - ai[row];
        while (nz-- > 0) {
            /* ja[work[i1]++] = i2 + oshift; */
            ja[work[i2]++] = i1 + oshift;
            i2++;
            while(((col = *j++ + ishift) < tns[i2]) && nz > 0) {
                nz--;
            }
            if (nz > 0) i2 = tvc[col];
        }
    }
    ierr = PetscFree(ns_col);
    CHKERRQ(ierr);
    ierr = PetscFree(work);
    CHKERRQ(ierr);
    ierr = PetscFree(tns);
    CHKERRQ(ierr);
    ierr = PetscFree(tvc);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatGetColumnIJ_Inode"
static PetscErrorCode MatGetColumnIJ_Inode(Mat A,PetscInt oshift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt *ja[],PetscTruth *done) {
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = Mat_CreateColInode(A,n,PETSC_NULL);
    CHKERRQ(ierr);
    if (!ia) PetscFunctionReturn(0);

    if (!blockcompressed) {
        ierr = MatGetColumnIJ_SeqAIJ(A,oshift,symmetric,blockcompressed,n,ia,ja,done);
        CHKERRQ(ierr);;
    } else if (symmetric) {
        /* Since the indices are symmetric it does'nt matter */
        ierr = MatGetRowIJ_Inode_Symmetric(A,ia,ja,0,oshift);
        CHKERRQ(ierr);
    } else {
        ierr = MatGetColumnIJ_Inode_Nonsymmetric(A,ia,ja,0,oshift);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatRestoreColumnIJ_Inode"
static PetscErrorCode MatRestoreColumnIJ_Inode(Mat A,PetscInt oshift,PetscTruth symmetric,PetscTruth blockcompressed,PetscInt *n,PetscInt *ia[],PetscInt *ja[],PetscTruth *done) {
    PetscErrorCode ierr;

    PetscFunctionBegin;
    if (!ia) PetscFunctionReturn(0);
    if (!blockcompressed) {
        ierr = MatRestoreColumnIJ_SeqAIJ(A,oshift,symmetric,blockcompressed,n,ia,ja,done);
        CHKERRQ(ierr);;
    } else {
        ierr = PetscFree(*ia);
        CHKERRQ(ierr);
        ierr = PetscFree(*ja);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

/* ----------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "MatMult_Inode"
static PetscErrorCode MatMult_Inode(Mat A,Vec xx,Vec yy) {
    Mat_SeqAIJ        *a = (Mat_SeqAIJ*)A->data;
    PetscScalar       sum1,sum2,sum3,sum4,sum5,tmp0,tmp1;
    PetscScalar       *y;
    const PetscScalar *x;
    const MatScalar   *v1,*v2,*v3,*v4,*v5;
    PetscErrorCode    ierr;
    PetscInt          *idx,i1,i2,n,i,row,node_max,*ns,*ii,nsz,sz,nonzerorow=0;

#if defined(PETSC_HAVE_PRAGMA_DISJOINT)
#pragma disjoint(*x,*y,*v1,*v2,*v3,*v4,*v5)
#endif

    PetscFunctionBegin;
    if (!a->inode.size) SETERRQ(PETSC_ERR_COR,"Missing Inode Structure");
    node_max = a->inode.node_count;
    ns       = a->inode.size;     /* Node Size array */
    ierr = VecGetArray(xx,(PetscScalar**)&x);
    CHKERRQ(ierr);
    ierr = VecGetArray(yy,&y);
    CHKERRQ(ierr);
    idx  = a->j;
    v1   = a->a;
    ii   = a->i;

    for (i = 0,row = 0; i< node_max; ++i) {
        nsz  = ns[i];
        n    = ii[1] - ii[0];
        nonzerorow += (n>0)*nsz;
        ii  += nsz;
        sz   = n;                   /* No of non zeros in this row */
        /* Switch on the size of Node */
        switch (nsz) {              /* Each loop in 'case' is unrolled */
        case 1 :
            sum1  = 0;

            for(n = 0; n< sz-1; n+=2) {
                i1   = idx[0];          /* The instructions are ordered to */
                i2   = idx[1];          /* make the compiler's job easy */
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
            }

            if (n == sz-1) {         /* Take care of the last nonzero  */
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
            }
            y[row++]=sum1;
            break;
        case 2:
            sum1  = 0;
            sum2  = 0;
            v2    = v1 + n;

            for (n = 0; n< sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
            }
            if (n == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            v1      =v2;              /* Since the next block to be processed starts there*/
            idx    +=sz;
            break;
        case 3:
            sum1  = 0;
            sum2  = 0;
            sum3  = 0;
            v2    = v1 + n;
            v3    = v2 + n;

            for (n = 0; n< sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
                sum3 += v3[0] * tmp0 + v3[1] * tmp1;
                v3 += 2;
            }
            if (n == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
                sum3 += *v3++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            y[row++]=sum3;
            v1       =v3;             /* Since the next block to be processed starts there*/
            idx     +=2*sz;
            break;
        case 4:
            sum1  = 0;
            sum2  = 0;
            sum3  = 0;
            sum4  = 0;
            v2    = v1 + n;
            v3    = v2 + n;
            v4    = v3 + n;

            for (n = 0; n< sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] *tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] *tmp1;
                v2 += 2;
                sum3 += v3[0] * tmp0 + v3[1] *tmp1;
                v3 += 2;
                sum4 += v4[0] * tmp0 + v4[1] *tmp1;
                v4 += 2;
            }
            if (n == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
                sum3 += *v3++ * tmp0;
                sum4 += *v4++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            y[row++]=sum3;
            y[row++]=sum4;
            v1      =v4;              /* Since the next block to be processed starts there*/
            idx    +=3*sz;
            break;
        case 5:
            sum1  = 0;
            sum2  = 0;
            sum3  = 0;
            sum4  = 0;
            sum5  = 0;
            v2    = v1 + n;
            v3    = v2 + n;
            v4    = v3 + n;
            v5    = v4 + n;

            for (n = 0; n<sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] *tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] *tmp1;
                v2 += 2;
                sum3 += v3[0] * tmp0 + v3[1] *tmp1;
                v3 += 2;
                sum4 += v4[0] * tmp0 + v4[1] *tmp1;
                v4 += 2;
                sum5 += v5[0] * tmp0 + v5[1] *tmp1;
                v5 += 2;
            }
            if (n == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
                sum3 += *v3++ * tmp0;
                sum4 += *v4++ * tmp0;
                sum5 += *v5++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            y[row++]=sum3;
            y[row++]=sum4;
            y[row++]=sum5;
            v1      =v5;       /* Since the next block to be processed starts there */
            idx    +=4*sz;
            break;
        default :
            SETERRQ(PETSC_ERR_COR,"Node size not yet supported");
        }
    }
    ierr = VecRestoreArray(xx,(PetscScalar**)&x);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(yy,&y);
    CHKERRQ(ierr);
    ierr = PetscLogFlops(2*a->nz - nonzerorow);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
/* ----------------------------------------------------------- */
/* Almost same code as the MatMult_Inode() */
#undef __FUNCT__
#define __FUNCT__ "MatMultAdd_Inode"
static PetscErrorCode MatMultAdd_Inode(Mat A,Vec xx,Vec zz,Vec yy) {
    Mat_SeqAIJ      *a = (Mat_SeqAIJ*)A->data;
    PetscScalar    sum1,sum2,sum3,sum4,sum5,tmp0,tmp1;
    MatScalar      *v1,*v2,*v3,*v4,*v5;
    PetscScalar    *x,*y,*z,*zt;
    PetscErrorCode ierr;
    PetscInt       *idx,i1,i2,n,i,row,node_max,*ns,*ii,nsz,sz;

    PetscFunctionBegin;
    if (!a->inode.size) SETERRQ(PETSC_ERR_COR,"Missing Inode Structure");
    node_max = a->inode.node_count;
    ns       = a->inode.size;     /* Node Size array */
    ierr = VecGetArray(xx,&x);
    CHKERRQ(ierr);
    ierr = VecGetArray(yy,&y);
    CHKERRQ(ierr);
    if (zz != yy) {
        ierr = VecGetArray(zz,&z);
        CHKERRQ(ierr);
    } else {
        z = y;
    }
    zt = z;

    idx  = a->j;
    v1   = a->a;
    ii   = a->i;

    for (i = 0,row = 0; i< node_max; ++i) {
        nsz  = ns[i];
        n    = ii[1] - ii[0];
        ii  += nsz;
        sz   = n;                   /* No of non zeros in this row */
        /* Switch on the size of Node */
        switch (nsz) {              /* Each loop in 'case' is unrolled */
        case 1 :
            sum1  = *zt++;

            for(n = 0; n< sz-1; n+=2) {
                i1   = idx[0];          /* The instructions are ordered to */
                i2   = idx[1];          /* make the compiler's job easy */
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
            }

            if(n   == sz-1) {         /* Take care of the last nonzero  */
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
            }
            y[row++]=sum1;
            break;
        case 2:
            sum1  = *zt++;
            sum2  = *zt++;
            v2    = v1 + n;

            for(n = 0; n< sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
            }
            if(n   == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            v1      =v2;              /* Since the next block to be processed starts there*/
            idx    +=sz;
            break;
        case 3:
            sum1  = *zt++;
            sum2  = *zt++;
            sum3  = *zt++;
            v2    = v1 + n;
            v3    = v2 + n;

            for (n = 0; n< sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
                sum3 += v3[0] * tmp0 + v3[1] * tmp1;
                v3 += 2;
            }
            if (n == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
                sum3 += *v3++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            y[row++]=sum3;
            v1       =v3;             /* Since the next block to be processed starts there*/
            idx     +=2*sz;
            break;
        case 4:
            sum1  = *zt++;
            sum2  = *zt++;
            sum3  = *zt++;
            sum4  = *zt++;
            v2    = v1 + n;
            v3    = v2 + n;
            v4    = v3 + n;

            for (n = 0; n< sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] *tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] *tmp1;
                v2 += 2;
                sum3 += v3[0] * tmp0 + v3[1] *tmp1;
                v3 += 2;
                sum4 += v4[0] * tmp0 + v4[1] *tmp1;
                v4 += 2;
            }
            if (n == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
                sum3 += *v3++ * tmp0;
                sum4 += *v4++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            y[row++]=sum3;
            y[row++]=sum4;
            v1      =v4;              /* Since the next block to be processed starts there*/
            idx    +=3*sz;
            break;
        case 5:
            sum1  = *zt++;
            sum2  = *zt++;
            sum3  = *zt++;
            sum4  = *zt++;
            sum5  = *zt++;
            v2    = v1 + n;
            v3    = v2 + n;
            v4    = v3 + n;
            v5    = v4 + n;

            for (n = 0; n<sz-1; n+=2) {
                i1   = idx[0];
                i2   = idx[1];
                idx += 2;
                tmp0 = x[i1];
                tmp1 = x[i2];
                sum1 += v1[0] * tmp0 + v1[1] *tmp1;
                v1 += 2;
                sum2 += v2[0] * tmp0 + v2[1] *tmp1;
                v2 += 2;
                sum3 += v3[0] * tmp0 + v3[1] *tmp1;
                v3 += 2;
                sum4 += v4[0] * tmp0 + v4[1] *tmp1;
                v4 += 2;
                sum5 += v5[0] * tmp0 + v5[1] *tmp1;
                v5 += 2;
            }
            if(n   == sz-1) {
                tmp0  = x[*idx++];
                sum1 += *v1++ * tmp0;
                sum2 += *v2++ * tmp0;
                sum3 += *v3++ * tmp0;
                sum4 += *v4++ * tmp0;
                sum5 += *v5++ * tmp0;
            }
            y[row++]=sum1;
            y[row++]=sum2;
            y[row++]=sum3;
            y[row++]=sum4;
            y[row++]=sum5;
            v1      =v5;       /* Since the next block to be processed starts there */
            idx    +=4*sz;
            break;
        default :
            SETERRQ(PETSC_ERR_COR,"Node size not yet supported");
        }
    }
    ierr = VecRestoreArray(xx,&x);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(yy,&y);
    CHKERRQ(ierr);
    if (zz != yy) {
        ierr = VecRestoreArray(zz,&z);
        CHKERRQ(ierr);
    }
    ierr = PetscLogFlops(2*a->nz);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/* ----------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MatSolve_Inode"
PetscErrorCode MatSolve_Inode(Mat A,Vec bb,Vec xx) {
    Mat_SeqAIJ        *a = (Mat_SeqAIJ*)A->data;
    IS                iscol = a->col,isrow = a->row;
    PetscErrorCode    ierr;
    PetscInt          *r,*c,i,j,n = A->rmap.n,*ai = a->i,nz,*a_j = a->j;
    PetscInt          node_max,*ns,row,nsz,aii,*vi,*ad,*aj,i0,i1,*rout,*cout;
    PetscScalar       *x,*tmp,*tmps,tmp0,tmp1;
    PetscScalar       sum1,sum2,sum3,sum4,sum5;
    const MatScalar   *v1,*v2,*v3,*v4,*v5,*a_a = a->a,*aa;
    const PetscScalar *b;

    PetscFunctionBegin;
    if (A->factor!=FACTOR_LU) SETERRQ(PETSC_ERR_ARG_WRONGSTATE,"Not for unfactored matrix");
    if (!a->inode.size) SETERRQ(PETSC_ERR_COR,"Missing Inode Structure");
    node_max = a->inode.node_count;
    ns       = a->inode.size;     /* Node Size array */

    ierr = VecGetArray(bb,(PetscScalar**)&b);
    CHKERRQ(ierr);
    ierr = VecGetArray(xx,&x);
    CHKERRQ(ierr);
    tmp  = a->solve_work;

    ierr = ISGetIndices(isrow,&rout);
    CHKERRQ(ierr);
    r = rout;
    ierr = ISGetIndices(iscol,&cout);
    CHKERRQ(ierr);
    c = cout + (n-1);

    /* forward solve the lower triangular */
    tmps = tmp ;
    aa   = a_a ;
    aj   = a_j ;
    ad   = a->diag;

    for (i = 0,row = 0; i< node_max; ++i) {
        nsz = ns[i];
        aii = ai[row];
        v1  = aa + aii;
        vi  = aj + aii;
        nz  = ad[row]- aii;

        switch (nsz) {              /* Each loop in 'case' is unrolled */
        case 1 :
            sum1 = b[*r++];
            /*      while (nz--) sum1 -= *v1++ *tmps[*vi++];*/
            for(j=0; j<nz-1; j+=2) {
                i0   = vi[0];
                i1   = vi[1];
                vi  +=2;
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
            }
            if(j == nz-1) {
                tmp0 = tmps[*vi++];
                sum1 -= *v1++ *tmp0;
            }
            tmp[row ++]=sum1;
            break;
        case 2:
            sum1 = b[*r++];
            sum2 = b[*r++];
            v2   = aa + ai[row+1];

            for(j=0; j<nz-1; j+=2) {
                i0   = vi[0];
                i1   = vi[1];
                vi  +=2;
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
            }
            if(j == nz-1) {
                tmp0 = tmps[*vi++];
                sum1 -= *v1++ *tmp0;
                sum2 -= *v2++ *tmp0;
            }
            sum2 -= *v2++ * sum1;
            tmp[row ++]=sum1;
            tmp[row ++]=sum2;
            break;
        case 3:
            sum1 = b[*r++];
            sum2 = b[*r++];
            sum3 = b[*r++];
            v2   = aa + ai[row+1];
            v3   = aa + ai[row+2];

            for (j=0; j<nz-1; j+=2) {
                i0   = vi[0];
                i1   = vi[1];
                vi  +=2;
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
                sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                v3 += 2;
            }
            if (j == nz-1) {
                tmp0 = tmps[*vi++];
                sum1 -= *v1++ *tmp0;
                sum2 -= *v2++ *tmp0;
                sum3 -= *v3++ *tmp0;
            }
            sum2 -= *v2++ * sum1;
            sum3 -= *v3++ * sum1;
            sum3 -= *v3++ * sum2;
            tmp[row ++]=sum1;
            tmp[row ++]=sum2;
            tmp[row ++]=sum3;
            break;

        case 4:
            sum1 = b[*r++];
            sum2 = b[*r++];
            sum3 = b[*r++];
            sum4 = b[*r++];
            v2   = aa + ai[row+1];
            v3   = aa + ai[row+2];
            v4   = aa + ai[row+3];

            for (j=0; j<nz-1; j+=2) {
                i0   = vi[0];
                i1   = vi[1];
                vi  +=2;
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
                sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                v3 += 2;
                sum4 -= v4[0] * tmp0 + v4[1] * tmp1;
                v4 += 2;
            }
            if (j == nz-1) {
                tmp0 = tmps[*vi++];
                sum1 -= *v1++ *tmp0;
                sum2 -= *v2++ *tmp0;
                sum3 -= *v3++ *tmp0;
                sum4 -= *v4++ *tmp0;
            }
            sum2 -= *v2++ * sum1;
            sum3 -= *v3++ * sum1;
            sum4 -= *v4++ * sum1;
            sum3 -= *v3++ * sum2;
            sum4 -= *v4++ * sum2;
            sum4 -= *v4++ * sum3;

            tmp[row ++]=sum1;
            tmp[row ++]=sum2;
            tmp[row ++]=sum3;
            tmp[row ++]=sum4;
            break;
        case 5:
            sum1 = b[*r++];
            sum2 = b[*r++];
            sum3 = b[*r++];
            sum4 = b[*r++];
            sum5 = b[*r++];
            v2   = aa + ai[row+1];
            v3   = aa + ai[row+2];
            v4   = aa + ai[row+3];
            v5   = aa + ai[row+4];

            for (j=0; j<nz-1; j+=2) {
                i0   = vi[0];
                i1   = vi[1];
                vi  +=2;
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                v1 += 2;
                sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                v2 += 2;
                sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                v3 += 2;
                sum4 -= v4[0] * tmp0 + v4[1] * tmp1;
                v4 += 2;
                sum5 -= v5[0] * tmp0 + v5[1] * tmp1;
                v5 += 2;
            }
            if (j == nz-1) {
                tmp0 = tmps[*vi++];
                sum1 -= *v1++ *tmp0;
                sum2 -= *v2++ *tmp0;
                sum3 -= *v3++ *tmp0;
                sum4 -= *v4++ *tmp0;
                sum5 -= *v5++ *tmp0;
            }

            sum2 -= *v2++ * sum1;
            sum3 -= *v3++ * sum1;
            sum4 -= *v4++ * sum1;
            sum5 -= *v5++ * sum1;
            sum3 -= *v3++ * sum2;
            sum4 -= *v4++ * sum2;
            sum5 -= *v5++ * sum2;
            sum4 -= *v4++ * sum3;
            sum5 -= *v5++ * sum3;
            sum5 -= *v5++ * sum4;

            tmp[row ++]=sum1;
            tmp[row ++]=sum2;
            tmp[row ++]=sum3;
            tmp[row ++]=sum4;
            tmp[row ++]=sum5;
            break;
        default:
            SETERRQ(PETSC_ERR_COR,"Node size not yet supported \n");
        }
    }
    /* backward solve the upper triangular */
    for (i=node_max -1,row = n-1 ; i>=0; i--) {
        nsz = ns[i];
        aii = ai[row+1] -1;
        v1  = aa + aii;
        vi  = aj + aii;
        nz  = aii- ad[row];
        switch (nsz) {              /* Each loop in 'case' is unrolled */
        case 1 :
            sum1 = tmp[row];

            for(j=nz ; j>1; j-=2) {
                vi  -=2;
                i0   = vi[2];
                i1   = vi[1];
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                v1   -= 2;
                sum1 -= v1[2] * tmp0 + v1[1] * tmp1;
            }
            if (j==1) {
                tmp0  = tmps[*vi--];
                sum1 -= *v1-- * tmp0;
            }
            x[*c--] = tmp[row] = sum1*a_a[ad[row]];
            row--;
            break;
        case 2 :
            sum1 = tmp[row];
            sum2 = tmp[row -1];
            v2   = aa + ai[row]-1;
            for (j=nz ; j>1; j-=2) {
                vi  -=2;
                i0   = vi[2];
                i1   = vi[1];
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                v1   -= 2;
                v2   -= 2;
                sum1 -= v1[2] * tmp0 + v1[1] * tmp1;
                sum2 -= v2[2] * tmp0 + v2[1] * tmp1;
            }
            if (j==1) {
                tmp0  = tmps[*vi--];
                sum1 -= *v1-- * tmp0;
                sum2 -= *v2-- * tmp0;
            }

            tmp0    = x[*c--] = tmp[row] = sum1*a_a[ad[row]];
            row--;
            sum2   -= *v2-- * tmp0;
            x[*c--] = tmp[row] = sum2*a_a[ad[row]];
            row--;
            break;
        case 3 :
            sum1 = tmp[row];
            sum2 = tmp[row -1];
            sum3 = tmp[row -2];
            v2   = aa + ai[row]-1;
            v3   = aa + ai[row -1]-1;
            for (j=nz ; j>1; j-=2) {
                vi  -=2;
                i0   = vi[2];
                i1   = vi[1];
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                v1   -= 2;
                v2   -= 2;
                v3   -= 2;
                sum1 -= v1[2] * tmp0 + v1[1] * tmp1;
                sum2 -= v2[2] * tmp0 + v2[1] * tmp1;
                sum3 -= v3[2] * tmp0 + v3[1] * tmp1;
            }
            if (j==1) {
                tmp0  = tmps[*vi--];
                sum1 -= *v1-- * tmp0;
                sum2 -= *v2-- * tmp0;
                sum3 -= *v3-- * tmp0;
            }
            tmp0    = x[*c--] = tmp[row] = sum1*a_a[ad[row]];
            row--;
            sum2   -= *v2-- * tmp0;
            sum3   -= *v3-- * tmp0;
            tmp0    = x[*c--] = tmp[row] = sum2*a_a[ad[row]];
            row--;
            sum3   -= *v3-- * tmp0;
            x[*c--] = tmp[row] = sum3*a_a[ad[row]];
            row--;

            break;
        case 4 :
            sum1 = tmp[row];
            sum2 = tmp[row -1];
            sum3 = tmp[row -2];
            sum4 = tmp[row -3];
            v2   = aa + ai[row]-1;
            v3   = aa + ai[row -1]-1;
            v4   = aa + ai[row -2]-1;

            for (j=nz ; j>1; j-=2) {
                vi  -=2;
                i0   = vi[2];
                i1   = vi[1];
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                v1  -= 2;
                v2  -= 2;
                v3  -= 2;
                v4  -= 2;
                sum1 -= v1[2] * tmp0 + v1[1] * tmp1;
                sum2 -= v2[2] * tmp0 + v2[1] * tmp1;
                sum3 -= v3[2] * tmp0 + v3[1] * tmp1;
                sum4 -= v4[2] * tmp0 + v4[1] * tmp1;
            }
            if (j==1) {
                tmp0  = tmps[*vi--];
                sum1 -= *v1-- * tmp0;
                sum2 -= *v2-- * tmp0;
                sum3 -= *v3-- * tmp0;
                sum4 -= *v4-- * tmp0;
            }

            tmp0    = x[*c--] = tmp[row] = sum1*a_a[ad[row]];
            row--;
            sum2   -= *v2-- * tmp0;
            sum3   -= *v3-- * tmp0;
            sum4   -= *v4-- * tmp0;
            tmp0    = x[*c--] = tmp[row] = sum2*a_a[ad[row]];
            row--;
            sum3   -= *v3-- * tmp0;
            sum4   -= *v4-- * tmp0;
            tmp0    = x[*c--] = tmp[row] = sum3*a_a[ad[row]];
            row--;
            sum4   -= *v4-- * tmp0;
            x[*c--] = tmp[row] = sum4*a_a[ad[row]];
            row--;
            break;
        case 5 :
            sum1 = tmp[row];
            sum2 = tmp[row -1];
            sum3 = tmp[row -2];
            sum4 = tmp[row -3];
            sum5 = tmp[row -4];
            v2   = aa + ai[row]-1;
            v3   = aa + ai[row -1]-1;
            v4   = aa + ai[row -2]-1;
            v5   = aa + ai[row -3]-1;
            for (j=nz ; j>1; j-=2) {
                vi  -= 2;
                i0   = vi[2];
                i1   = vi[1];
                tmp0 = tmps[i0];
                tmp1 = tmps[i1];
                v1   -= 2;
                v2   -= 2;
                v3   -= 2;
                v4   -= 2;
                v5   -= 2;
                sum1 -= v1[2] * tmp0 + v1[1] * tmp1;
                sum2 -= v2[2] * tmp0 + v2[1] * tmp1;
                sum3 -= v3[2] * tmp0 + v3[1] * tmp1;
                sum4 -= v4[2] * tmp0 + v4[1] * tmp1;
                sum5 -= v5[2] * tmp0 + v5[1] * tmp1;
            }
            if (j==1) {
                tmp0  = tmps[*vi--];
                sum1 -= *v1-- * tmp0;
                sum2 -= *v2-- * tmp0;
                sum3 -= *v3-- * tmp0;
                sum4 -= *v4-- * tmp0;
                sum5 -= *v5-- * tmp0;
            }

            tmp0    = x[*c--] = tmp[row] = sum1*a_a[ad[row]];
            row--;
            sum2   -= *v2-- * tmp0;
            sum3   -= *v3-- * tmp0;
            sum4   -= *v4-- * tmp0;
            sum5   -= *v5-- * tmp0;
            tmp0    = x[*c--] = tmp[row] = sum2*a_a[ad[row]];
            row--;
            sum3   -= *v3-- * tmp0;
            sum4   -= *v4-- * tmp0;
            sum5   -= *v5-- * tmp0;
            tmp0    = x[*c--] = tmp[row] = sum3*a_a[ad[row]];
            row--;
            sum4   -= *v4-- * tmp0;
            sum5   -= *v5-- * tmp0;
            tmp0    = x[*c--] = tmp[row] = sum4*a_a[ad[row]];
            row--;
            sum5   -= *v5-- * tmp0;
            x[*c--] = tmp[row] = sum5*a_a[ad[row]];
            row--;
            break;
        default:
            SETERRQ(PETSC_ERR_COR,"Node size not yet supported \n");
        }
    }
    ierr = ISRestoreIndices(isrow,&rout);
    CHKERRQ(ierr);
    ierr = ISRestoreIndices(iscol,&cout);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(bb,(PetscScalar**)&b);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(xx,&x);
    CHKERRQ(ierr);
    ierr = PetscLogFlops(2*a->nz - A->cmap.n);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatLUFactorNumeric_Inode"
PetscErrorCode MatLUFactorNumeric_Inode(Mat A,MatFactorInfo *info,Mat *B) {
    Mat               C = *B;
    Mat_SeqAIJ        *a = (Mat_SeqAIJ*)A->data,*b = (Mat_SeqAIJ*)C->data;
    IS                iscol = b->col,isrow = b->row,isicol = b->icol;
    PetscErrorCode    ierr;
    PetscInt          *r,*ic,*c,n = A->rmap.n,*bi = b->i;
    PetscInt          *bj = b->j,*nbj=b->j +1,*ajtmp,*bjtmp,nz,nz_tmp,row,prow;
    PetscInt          *ics,i,j,idx,*ai = a->i,*aj = a->j,*bd = b->diag,node_max,nodesz;
    PetscInt          *ns,*tmp_vec1,*tmp_vec2,*nsmap,*pj;
    PetscScalar       mul1,mul2,mul3,tmp;
    MatScalar         *pc1,*pc2,*pc3,*ba = b->a,*pv,*rtmp11,*rtmp22,*rtmp33;
    const MatScalar   *v1,*v2,*v3,*aa = a->a,*rtmp1;
    PetscReal         rs=0.0;
    LUShift_Ctx       sctx;
    PetscInt          newshift;

    PetscFunctionBegin;
    sctx.shift_top  = 0;
    sctx.nshift_max = 0;
    sctx.shift_lo   = 0;
    sctx.shift_hi   = 0;

    /* if both shift schemes are chosen by user, only use info->shiftpd */
    if (info->shiftpd && info->shiftnz) info->shiftnz = 0.0;
    if (info->shiftpd) { /* set sctx.shift_top=max{rs} */
        sctx.shift_top = 0;
        for (i=0; i<n; i++) {
            /* calculate rs = sum(|aij|)-RealPart(aii), amt of shift needed for this row */
            rs    = 0.0;
            ajtmp = aj + ai[i];
            rtmp1 = aa + ai[i];
            nz = ai[i+1] - ai[i];
            for (j=0; j<nz; j++) {
                if (*ajtmp != i) {
                    rs += PetscAbsScalar(*rtmp1++);
                } else {
                    rs -= PetscRealPart(*rtmp1++);
                }
                ajtmp++;
            }
            if (rs>sctx.shift_top) sctx.shift_top = rs;
        }
        if (sctx.shift_top == 0.0) sctx.shift_top += 1.e-12;
        sctx.shift_top *= 1.1;
        sctx.nshift_max = 5;
        sctx.shift_lo   = 0.;
        sctx.shift_hi   = 1.;
    }
    sctx.shift_amount = 0;
    sctx.nshift       = 0;

    ierr  = ISGetIndices(isrow,&r);
    CHKERRQ(ierr);
    ierr  = ISGetIndices(iscol,&c);
    CHKERRQ(ierr);
    ierr  = ISGetIndices(isicol,&ic);
    CHKERRQ(ierr);
    ierr  = PetscMalloc((3*n+1)*sizeof(PetscScalar),&rtmp11);
    CHKERRQ(ierr);
    ierr  = PetscMemzero(rtmp11,(3*n+1)*sizeof(PetscScalar));
    CHKERRQ(ierr);
    ics   = ic ;
    rtmp22 = rtmp11 + n;
    rtmp33 = rtmp22 + n;

    node_max = a->inode.node_count;
    ns       = a->inode.size ;
    if (!ns) {
        SETERRQ(PETSC_ERR_PLIB,"Matrix without inode information");
    }

    /* If max inode size > 3, split it into two inodes.*/
    /* also map the inode sizes according to the ordering */
    ierr = PetscMalloc((n+1)* sizeof(PetscInt),&tmp_vec1);
    CHKERRQ(ierr);
    for (i=0,j=0; i<node_max; ++i,++j) {
        if (ns[i]>3) {
            tmp_vec1[j] = ns[i]/2; /* Assuming ns[i] < =5  */
            ++j;
            tmp_vec1[j] = ns[i] - tmp_vec1[j-1];
        } else {
            tmp_vec1[j] = ns[i];
        }
    }
    /* Use the correct node_max */
    node_max = j;

    /* Now reorder the inode info based on mat re-ordering info */
    /* First create a row -> inode_size_array_index map */
    ierr = PetscMalloc(n*sizeof(PetscInt)+1,&nsmap);
    CHKERRQ(ierr);
    ierr = PetscMalloc(node_max*sizeof(PetscInt)+1,&tmp_vec2);
    CHKERRQ(ierr);
    for (i=0,row=0; i<node_max; i++) {
        nodesz = tmp_vec1[i];
        for (j=0; j<nodesz; j++,row++) {
            nsmap[row] = i;
        }
    }
    /* Using nsmap, create a reordered ns structure */
    for (i=0,j=0; i< node_max; i++) {
        nodesz       = tmp_vec1[nsmap[r[j]]];    /* here the reordered row_no is in r[] */
        tmp_vec2[i]  = nodesz;
        j           += nodesz;
    }
    ierr = PetscFree(nsmap);
    CHKERRQ(ierr);
    ierr = PetscFree(tmp_vec1);
    CHKERRQ(ierr);
    /* Now use the correct ns */
    ns = tmp_vec2;

    do {
        sctx.lushift = PETSC_FALSE;
        /* Now loop over each block-row, and do the factorization */
        for (i=0,row=0; i<node_max; i++) {
            nodesz = ns[i];
            nz     = bi[row+1] - bi[row];
            bjtmp  = bj + bi[row];

            switch (nodesz) {
            case 1:
                for  (j=0; j<nz; j++) {
                    idx        = bjtmp[j];
                    rtmp11[idx] = 0.0;
                }

                /* load in initial (unfactored row) */
                idx    = r[row];
                nz_tmp = ai[idx+1] - ai[idx];
                ajtmp  = aj + ai[idx];
                v1     = aa + ai[idx];

                for (j=0; j<nz_tmp; j++) {
                    idx        = ics[ajtmp[j]];
                    rtmp11[idx] = v1[j];
                }
                rtmp11[ics[r[row]]] += sctx.shift_amount;

                prow = *bjtmp++ ;
                while (prow < row) {
                    pc1 = rtmp11 + prow;
                    if (*pc1 != 0.0) {
                        pv   = ba + bd[prow];
                        pj   = nbj + bd[prow];
                        mul1 = *pc1 * *pv++;
                        *pc1 = mul1;
                        nz_tmp = bi[prow+1] - bd[prow] - 1;
                        ierr = PetscLogFlops(2*nz_tmp);
                        CHKERRQ(ierr);
                        for (j=0; j<nz_tmp; j++) {
                            tmp = pv[j];
                            idx = pj[j];
                            rtmp11[idx] -= mul1 * tmp;
                        }
                    }
                    prow = *bjtmp++ ;
                }
                pj  = bj + bi[row];
                pc1 = ba + bi[row];

                sctx.pv    = rtmp11[row];
                rtmp11[row] = 1.0/rtmp11[row]; /* invert diag */
                rs         = 0.0;
                for (j=0; j<nz; j++) {
                    idx    = pj[j];
                    pc1[j] = rtmp11[idx]; /* rtmp11 -> ba */
                    if (idx != row) rs += PetscAbsScalar(pc1[j]);
                }
                sctx.rs  = rs;
                ierr = MatLUCheckShift_inline(info,sctx,row,newshift);
                CHKERRQ(ierr);
                if (newshift == 1) goto endofwhile;
                break;

            case 2:
                for (j=0; j<nz; j++) {
                    idx        = bjtmp[j];
                    rtmp11[idx] = 0.0;
                    rtmp22[idx] = 0.0;
                }

                /* load in initial (unfactored row) */
                idx    = r[row];
                nz_tmp = ai[idx+1] - ai[idx];
                ajtmp  = aj + ai[idx];
                v1     = aa + ai[idx];
                v2     = aa + ai[idx+1];
                for (j=0; j<nz_tmp; j++) {
                    idx        = ics[ajtmp[j]];
                    rtmp11[idx] = v1[j];
                    rtmp22[idx] = v2[j];
                }
                rtmp11[ics[r[row]]]   += sctx.shift_amount;
                rtmp22[ics[r[row+1]]] += sctx.shift_amount;

                prow = *bjtmp++ ;
                while (prow < row) {
                    pc1 = rtmp11 + prow;
                    pc2 = rtmp22 + prow;
                    if (*pc1 != 0.0 || *pc2 != 0.0) {
                        pv   = ba + bd[prow];
                        pj   = nbj + bd[prow];
                        mul1 = *pc1 * *pv;
                        mul2 = *pc2 * *pv;
                        ++pv;
                        *pc1 = mul1;
                        *pc2 = mul2;

                        nz_tmp = bi[prow+1] - bd[prow] - 1;
                        for (j=0; j<nz_tmp; j++) {
                            tmp = pv[j];
                            idx = pj[j];
                            rtmp11[idx] -= mul1 * tmp;
                            rtmp22[idx] -= mul2 * tmp;
                        }
                        ierr = PetscLogFlops(4*nz_tmp);
                        CHKERRQ(ierr);
                    }
                    prow = *bjtmp++ ;
                }

                /* Now take care of diagonal 2x2 block. Note: prow = row here */
                pc1 = rtmp11 + prow;
                pc2 = rtmp22 + prow;

                sctx.pv = *pc1;
                pj      = bj + bi[prow];
                rs      = 0.0;
                for (j=0; j<nz; j++) {
                    idx = pj[j];
                    if (idx != prow) rs += PetscAbsScalar(rtmp11[idx]);
                }
                sctx.rs = rs;
                ierr = MatLUCheckShift_inline(info,sctx,row,newshift);
                CHKERRQ(ierr);
                if (newshift == 1) goto endofwhile;

                if (*pc2 != 0.0) {
                    pj     = nbj + bd[prow];
                    mul2   = (*pc2)/(*pc1); /* since diag is not yet inverted.*/
                    *pc2   = mul2;
                    nz_tmp = bi[prow+1] - bd[prow] - 1;
                    for (j=0; j<nz_tmp; j++) {
                        idx = pj[j] ;
                        tmp = rtmp11[idx];
                        rtmp22[idx] -= mul2 * tmp;
                    }
                    ierr = PetscLogFlops(2*nz_tmp);
                    CHKERRQ(ierr);
                }

                pj  = bj + bi[row];
                pc1 = ba + bi[row];
                pc2 = ba + bi[row+1];

                sctx.pv = rtmp22[row+1];
                rs = 0.0;
                rtmp11[row]   = 1.0/rtmp11[row];
                rtmp22[row+1] = 1.0/rtmp22[row+1];
                /* copy row entries from dense representation to sparse */
                for (j=0; j<nz; j++) {
                    idx    = pj[j];
                    pc1[j] = rtmp11[idx];
                    pc2[j] = rtmp22[idx];
                    if (idx != row+1) rs += PetscAbsScalar(pc2[j]);
                }
                sctx.rs = rs;
                ierr = MatLUCheckShift_inline(info,sctx,row+1,newshift);
                CHKERRQ(ierr);
                if (newshift == 1) goto endofwhile;
                break;

            case 3:
                for  (j=0; j<nz; j++) {
                    idx        = bjtmp[j];
                    rtmp11[idx] = 0.0;
                    rtmp22[idx] = 0.0;
                    rtmp33[idx] = 0.0;
                }
                /* copy the nonzeros for the 3 rows from sparse representation to dense in rtmp*[] */
                idx    = r[row];
                nz_tmp = ai[idx+1] - ai[idx];
                ajtmp = aj + ai[idx];
                v1    = aa + ai[idx];
                v2    = aa + ai[idx+1];
                v3    = aa + ai[idx+2];
                for (j=0; j<nz_tmp; j++) {
                    idx        = ics[ajtmp[j]];
                    rtmp11[idx] = v1[j];
                    rtmp22[idx] = v2[j];
                    rtmp33[idx] = v3[j];
                }
                rtmp11[ics[r[row]]]   += sctx.shift_amount;
                rtmp22[ics[r[row+1]]] += sctx.shift_amount;
                rtmp33[ics[r[row+2]]] += sctx.shift_amount;

                /* loop over all pivot row blocks above this row block */
                prow = *bjtmp++ ;
                while (prow < row) {
                    pc1 = rtmp11 + prow;
                    pc2 = rtmp22 + prow;
                    pc3 = rtmp33 + prow;
                    if (*pc1 != 0.0 || *pc2 != 0.0 || *pc3 !=0.0) {
                        pv   = ba  + bd[prow];
                        pj   = nbj + bd[prow];
                        mul1 = *pc1 * *pv;
                        mul2 = *pc2 * *pv;
                        mul3 = *pc3 * *pv;
                        ++pv;
                        *pc1 = mul1;
                        *pc2 = mul2;
                        *pc3 = mul3;

                        nz_tmp = bi[prow+1] - bd[prow] - 1;
                        /* update this row based on pivot row */
                        for (j=0; j<nz_tmp; j++) {
                            tmp = pv[j];
                            idx = pj[j];
                            rtmp11[idx] -= mul1 * tmp;
                            rtmp22[idx] -= mul2 * tmp;
                            rtmp33[idx] -= mul3 * tmp;
                        }
                        ierr = PetscLogFlops(6*nz_tmp);
                        CHKERRQ(ierr);
                    }
                    prow = *bjtmp++ ;
                }

                /* Now take care of diagonal 3x3 block in this set of rows */
                /* note: prow = row here */
                pc1 = rtmp11 + prow;
                pc2 = rtmp22 + prow;
                pc3 = rtmp33 + prow;

                sctx.pv = *pc1;
                pj      = bj + bi[prow];
                rs      = 0.0;
                for (j=0; j<nz; j++) {
                    idx = pj[j];
                    if (idx != row) rs += PetscAbsScalar(rtmp11[idx]);
                }
                sctx.rs = rs;
                ierr = MatLUCheckShift_inline(info,sctx,row,newshift);
                CHKERRQ(ierr);
                if (newshift == 1) goto endofwhile;

                if (*pc2 != 0.0 || *pc3 != 0.0) {
                    mul2 = (*pc2)/(*pc1);
                    mul3 = (*pc3)/(*pc1);
                    *pc2 = mul2;
                    *pc3 = mul3;
                    nz_tmp = bi[prow+1] - bd[prow] - 1;
                    pj     = nbj + bd[prow];
                    for (j=0; j<nz_tmp; j++) {
                        idx = pj[j] ;
                        tmp = rtmp11[idx];
                        rtmp22[idx] -= mul2 * tmp;
                        rtmp33[idx] -= mul3 * tmp;
                    }
                    ierr = PetscLogFlops(4*nz_tmp);
                    CHKERRQ(ierr);
                }
                ++prow;

                pc2 = rtmp22 + prow;
                pc3 = rtmp33 + prow;
                sctx.pv = *pc2;
                pj      = bj + bi[prow];
                rs      = 0.0;
                for (j=0; j<nz; j++) {
                    idx = pj[j];
                    if (idx != prow) rs += PetscAbsScalar(rtmp22[idx]);
                }
                sctx.rs = rs;
                ierr = MatLUCheckShift_inline(info,sctx,row+1,newshift);
                CHKERRQ(ierr);
                if (newshift == 1) goto endofwhile;

                if (*pc3 != 0.0) {
                    mul3   = (*pc3)/(*pc2);
                    *pc3   = mul3;
                    pj     = nbj + bd[prow];
                    nz_tmp = bi[prow+1] - bd[prow] - 1;
                    for (j=0; j<nz_tmp; j++) {
                        idx = pj[j] ;
                        tmp = rtmp22[idx];
                        rtmp33[idx] -= mul3 * tmp;
                    }
                    ierr = PetscLogFlops(4*nz_tmp);
                    CHKERRQ(ierr);
                }

                pj  = bj + bi[row];
                pc1 = ba + bi[row];
                pc2 = ba + bi[row+1];
                pc3 = ba + bi[row+2];

                sctx.pv = rtmp33[row+2];
                rs = 0.0;
                rtmp11[row]   = 1.0/rtmp11[row];
                rtmp22[row+1] = 1.0/rtmp22[row+1];
                rtmp33[row+2] = 1.0/rtmp33[row+2];
                /* copy row entries from dense representation to sparse */
                for (j=0; j<nz; j++) {
                    idx    = pj[j];
                    pc1[j] = rtmp11[idx];
                    pc2[j] = rtmp22[idx];
                    pc3[j] = rtmp33[idx];
                    if (idx != row+2) rs += PetscAbsScalar(pc3[j]);
                }

                sctx.rs = rs;
                ierr = MatLUCheckShift_inline(info,sctx,row+2,newshift);
                CHKERRQ(ierr);
                if (newshift == 1) goto endofwhile;
                break;

            default:
                SETERRQ(PETSC_ERR_COR,"Node size not yet supported \n");
            }
            row += nodesz;                 /* Update the row */
        }
endofwhile:
        ;
    } while (sctx.lushift);
    ierr = PetscFree(rtmp11);
    CHKERRQ(ierr);
    ierr = PetscFree(tmp_vec2);
    CHKERRQ(ierr);
    ierr = ISRestoreIndices(isicol,&ic);
    CHKERRQ(ierr);
    ierr = ISRestoreIndices(isrow,&r);
    CHKERRQ(ierr);
    ierr = ISRestoreIndices(iscol,&c);
    CHKERRQ(ierr);
    C->factor      = FACTOR_LU;
    C->assembled   = PETSC_TRUE;
    if (sctx.nshift) {
        if (info->shiftnz) {
            ierr = PetscInfo2(A,"number of shift_nz tries %D, shift_amount %G\n",sctx.nshift,sctx.shift_amount);
            CHKERRQ(ierr);
        } else if (info->shiftpd) {
            ierr = PetscInfo4(A,"number of shift_pd tries %D, shift_amount %G, diagonal shifted up by %e fraction top_value %e\n",sctx.nshift,sctx.shift_amount,info->shift_fraction,sctx.shift_top);
            CHKERRQ(ierr);
        }
    }
    ierr = PetscLogFlops(C->cmap.n);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/*
     Makes a longer coloring[] array and calls the usual code with that
*/
#undef __FUNCT__
#define __FUNCT__ "MatColoringPatch_Inode"
PetscErrorCode MatColoringPatch_Inode(Mat mat,PetscInt ncolors,PetscInt nin,ISColoringValue coloring[],ISColoring *iscoloring) {
    Mat_SeqAIJ       *a = (Mat_SeqAIJ*)mat->data;
    PetscErrorCode  ierr;
    PetscInt        n = mat->cmap.n,m = a->inode.node_count,j,*ns = a->inode.size,row;
    PetscInt        *colorused,i;
    ISColoringValue *newcolor;

    PetscFunctionBegin;
    ierr = PetscMalloc((n+1)*sizeof(PetscInt),&newcolor);
    CHKERRQ(ierr);
    /* loop over inodes, marking a color for each column*/
    row = 0;
    for (i=0; i<m; i++) {
        for (j=0; j<ns[i]; j++) {
            newcolor[row++] = coloring[i] + j*ncolors;
        }
    }

    /* eliminate unneeded colors */
    ierr = PetscMalloc(5*ncolors*sizeof(PetscInt),&colorused);
    CHKERRQ(ierr);
    ierr = PetscMemzero(colorused,5*ncolors*sizeof(PetscInt));
    CHKERRQ(ierr);
    for (i=0; i<n; i++) {
        colorused[newcolor[i]] = 1;
    }

    for (i=1; i<5*ncolors; i++) {
        colorused[i] += colorused[i-1];
    }
    ncolors = colorused[5*ncolors-1];
    for (i=0; i<n; i++) {
        newcolor[i] = colorused[newcolor[i]]-1;
    }
    ierr = PetscFree(colorused);
    CHKERRQ(ierr);
    ierr = ISColoringCreate(((PetscObject)mat)->comm,ncolors,n,newcolor,iscoloring);
    CHKERRQ(ierr);
    ierr = PetscFree(coloring);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#include "src/inline/ilu.h"

#undef __FUNCT__
#define __FUNCT__ "MatRelax_Inode"
PetscErrorCode MatRelax_Inode(Mat A,Vec bb,PetscReal omega,MatSORType flag,PetscReal fshift,PetscInt its,PetscInt lits,Vec xx) {
    Mat_SeqAIJ         *a = (Mat_SeqAIJ*)A->data;
    PetscScalar        *x,*xs,sum1,sum2,sum3,sum4,sum5,tmp0,tmp1,tmp2,tmp3;
    MatScalar          *ibdiag,*bdiag;
    PetscScalar        *b,*xb,tmp4,tmp5,x1,x2,x3,x4,x5;
    const MatScalar    *v = a->a,*v1,*v2,*v3,*v4,*v5;
    PetscReal          zeropivot = 1.0e-15, shift = 0.0;
    PetscErrorCode     ierr;
    PetscInt           n,m = a->inode.node_count,*sizes = a->inode.size,cnt = 0,i,j,row,i1,i2;
    PetscInt           *idx,*diag = a->diag,*ii = a->i,sz,k;

    PetscFunctionBegin;
    if (omega != 1.0) SETERRQ(PETSC_ERR_SUP,"No support for omega != 1.0; use -mat_no_inode");
    if (fshift != 0.0) SETERRQ(PETSC_ERR_SUP,"No support for fshift != 0.0; use -mat_no_inode");
    if (flag & SOR_EISENSTAT) SETERRQ(PETSC_ERR_SUP,"No support for Eisenstat trick; use -mat_no_inode");

    if (!a->inode.ibdiagvalid) {
        if (!a->inode.ibdiag) {
            /* calculate space needed for diagonal blocks */
            for (i=0; i<m; i++) {
                cnt += sizes[i]*sizes[i];
            }
            a->inode.bdiagsize = cnt;
            ierr   = PetscMalloc2(cnt,MatScalar,&a->inode.ibdiag,cnt,MatScalar,&a->inode.bdiag);
            CHKERRQ(ierr);
        }

        /* copy over the diagonal blocks and invert them */
        ibdiag = a->inode.ibdiag;
        bdiag  = a->inode.bdiag;
        cnt = 0;
        for (i=0, row = 0; i<m; i++) {
            for (j=0; j<sizes[i]; j++) {
                for (k=0; k<sizes[i]; k++) {
                    bdiag[cnt+k*sizes[i]+j] = v[diag[row+j] - j + k];
                }
            }
            ierr = PetscMemcpy(ibdiag+cnt,bdiag+cnt,sizes[i]*sizes[i]*sizeof(MatScalar));
            CHKERRQ(ierr);

            switch(sizes[i]) {
            case 1:
                /* Create matrix data structure */
                if (PetscAbsScalar(ibdiag[cnt]) < zeropivot) SETERRQ1(PETSC_ERR_MAT_LU_ZRPVT,"Zero pivot on row %D",row);
                ibdiag[cnt] = 1.0/ibdiag[cnt];
                break;
            case 2:
                ierr = Kernel_A_gets_inverse_A_2(ibdiag+cnt,shift);
                CHKERRQ(ierr);
                break;
            case 3:
                ierr = Kernel_A_gets_inverse_A_3(ibdiag+cnt,shift);
                CHKERRQ(ierr);
                break;
            case 4:
                ierr = Kernel_A_gets_inverse_A_4(ibdiag+cnt,shift);
                CHKERRQ(ierr);
                break;
            case 5:
                ierr = Kernel_A_gets_inverse_A_5(ibdiag+cnt,shift);
                CHKERRQ(ierr);
                break;
            default:
                SETERRQ1(PETSC_ERR_SUP,"Inode size %D not supported",sizes[i]);
            }
            cnt += sizes[i]*sizes[i];
            row += sizes[i];
        }
        a->inode.ibdiagvalid = PETSC_TRUE;
    }
    ibdiag = a->inode.ibdiag;
    bdiag  = a->inode.bdiag;

    ierr = VecGetArray(xx,&x);
    CHKERRQ(ierr);
    if (xx != bb) {
        ierr = VecGetArray(bb,(PetscScalar**)&b);
        CHKERRQ(ierr);
    } else {
        b = x;
    }

    /* We count flops by assuming the upper triangular and lower triangular parts have the same number of nonzeros */
    xs   = x;
    if (flag & SOR_ZERO_INITIAL_GUESS) {
        if (flag & SOR_FORWARD_SWEEP || flag & SOR_LOCAL_FORWARD_SWEEP) {

            for (i=0, row=0; i<m; i++) {
                sz  = diag[row] - ii[row];
                v1  = a->a + ii[row];
                idx = a->j + ii[row];

                /* see comments for MatMult_Inode() for how this is coded */
                switch (sizes[i]) {
                case 1:

                    sum1  = b[row];
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= *v1 * tmp0;
                    }
                    x[row++] = sum1*(*ibdiag++);
                    break;
                case 2:
                    v2    = a->a + ii[row+1];
                    sum1  = b[row];
                    sum2  = b[row+1];
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= v1[0] * tmp0;
                        sum2 -= v2[0] * tmp0;
                    }
                    x[row++] = sum1*ibdiag[0] + sum2*ibdiag[2];
                    x[row++] = sum1*ibdiag[1] + sum2*ibdiag[3];
                    ibdiag  += 4;
                    break;
                case 3:
                    v2    = a->a + ii[row+1];
                    v3    = a->a + ii[row+2];
                    sum1  = b[row];
                    sum2  = b[row+1];
                    sum3  = b[row+2];
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                        sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                        v3 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= v1[0] * tmp0;
                        sum2 -= v2[0] * tmp0;
                        sum3 -= v3[0] * tmp0;
                    }
                    x[row++] = sum1*ibdiag[0] + sum2*ibdiag[3] + sum3*ibdiag[6];
                    x[row++] = sum1*ibdiag[1] + sum2*ibdiag[4] + sum3*ibdiag[7];
                    x[row++] = sum1*ibdiag[2] + sum2*ibdiag[5] + sum3*ibdiag[8];
                    ibdiag  += 9;
                    break;
                case 4:
                    v2    = a->a + ii[row+1];
                    v3    = a->a + ii[row+2];
                    v4    = a->a + ii[row+3];
                    sum1  = b[row];
                    sum2  = b[row+1];
                    sum3  = b[row+2];
                    sum4  = b[row+3];
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                        sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                        v3 += 2;
                        sum4 -= v4[0] * tmp0 + v4[1] * tmp1;
                        v4 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= v1[0] * tmp0;
                        sum2 -= v2[0] * tmp0;
                        sum3 -= v3[0] * tmp0;
                        sum4 -= v4[0] * tmp0;
                    }
                    x[row++] = sum1*ibdiag[0] + sum2*ibdiag[4] + sum3*ibdiag[8] + sum4*ibdiag[12];
                    x[row++] = sum1*ibdiag[1] + sum2*ibdiag[5] + sum3*ibdiag[9] + sum4*ibdiag[13];
                    x[row++] = sum1*ibdiag[2] + sum2*ibdiag[6] + sum3*ibdiag[10] + sum4*ibdiag[14];
                    x[row++] = sum1*ibdiag[3] + sum2*ibdiag[7] + sum3*ibdiag[11] + sum4*ibdiag[15];
                    ibdiag  += 16;
                    break;
                case 5:
                    v2    = a->a + ii[row+1];
                    v3    = a->a + ii[row+2];
                    v4    = a->a + ii[row+3];
                    v5    = a->a + ii[row+4];
                    sum1  = b[row];
                    sum2  = b[row+1];
                    sum3  = b[row+2];
                    sum4  = b[row+3];
                    sum5  = b[row+4];
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                        sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                        v3 += 2;
                        sum4 -= v4[0] * tmp0 + v4[1] * tmp1;
                        v4 += 2;
                        sum5 -= v5[0] * tmp0 + v5[1] * tmp1;
                        v5 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= v1[0] * tmp0;
                        sum2 -= v2[0] * tmp0;
                        sum3 -= v3[0] * tmp0;
                        sum4 -= v4[0] * tmp0;
                        sum5 -= v5[0] * tmp0;
                    }
                    x[row++] = sum1*ibdiag[0] + sum2*ibdiag[5] + sum3*ibdiag[10] + sum4*ibdiag[15] + sum5*ibdiag[20];
                    x[row++] = sum1*ibdiag[1] + sum2*ibdiag[6] + sum3*ibdiag[11] + sum4*ibdiag[16] + sum5*ibdiag[21];
                    x[row++] = sum1*ibdiag[2] + sum2*ibdiag[7] + sum3*ibdiag[12] + sum4*ibdiag[17] + sum5*ibdiag[22];
                    x[row++] = sum1*ibdiag[3] + sum2*ibdiag[8] + sum3*ibdiag[13] + sum4*ibdiag[18] + sum5*ibdiag[23];
                    x[row++] = sum1*ibdiag[4] + sum2*ibdiag[9] + sum3*ibdiag[14] + sum4*ibdiag[19] + sum5*ibdiag[24];
                    ibdiag  += 25;
                    break;
                default:
                    SETERRQ1(PETSC_ERR_SUP,"Inode size %D not supported",sizes[i]);
                }
            }

            xb = x;
            ierr = PetscLogFlops(a->nz);
            CHKERRQ(ierr);
        } else xb = b;
        if ((flag & SOR_FORWARD_SWEEP || flag & SOR_LOCAL_FORWARD_SWEEP) &&
                (flag & SOR_BACKWARD_SWEEP || flag & SOR_LOCAL_BACKWARD_SWEEP)) {
            cnt = 0;
            for (i=0, row=0; i<m; i++) {

                switch (sizes[i]) {
                case 1:
                    x[row++] *= bdiag[cnt++];
                    break;
                case 2:
                    x1   = x[row];
                    x2 = x[row+1];
                    tmp1 = x1*bdiag[cnt] + x2*bdiag[cnt+2];
                    tmp2 = x1*bdiag[cnt+1] + x2*bdiag[cnt+3];
                    x[row++] = tmp1;
                    x[row++] = tmp2;
                    cnt += 4;
                    break;
                case 3:
                    x1   = x[row];
                    x2 = x[row+1];
                    x3 = x[row+2];
                    tmp1 = x1*bdiag[cnt] + x2*bdiag[cnt+3] + x3*bdiag[cnt+6];
                    tmp2 = x1*bdiag[cnt+1] + x2*bdiag[cnt+4] + x3*bdiag[cnt+7];
                    tmp3 = x1*bdiag[cnt+2] + x2*bdiag[cnt+5] + x3*bdiag[cnt+8];
                    x[row++] = tmp1;
                    x[row++] = tmp2;
                    x[row++] = tmp3;
                    cnt += 9;
                    break;
                case 4:
                    x1   = x[row];
                    x2 = x[row+1];
                    x3 = x[row+2];
                    x4 = x[row+3];
                    tmp1 = x1*bdiag[cnt] + x2*bdiag[cnt+4] + x3*bdiag[cnt+8] + x4*bdiag[cnt+12];
                    tmp2 = x1*bdiag[cnt+1] + x2*bdiag[cnt+5] + x3*bdiag[cnt+9] + x4*bdiag[cnt+13];
                    tmp3 = x1*bdiag[cnt+2] + x2*bdiag[cnt+6] + x3*bdiag[cnt+10] + x4*bdiag[cnt+14];
                    tmp4 = x1*bdiag[cnt+3] + x2*bdiag[cnt+7] + x3*bdiag[cnt+11] + x4*bdiag[cnt+15];
                    x[row++] = tmp1;
                    x[row++] = tmp2;
                    x[row++] = tmp3;
                    x[row++] = tmp4;
                    cnt += 16;
                    break;
                case 5:
                    x1   = x[row];
                    x2 = x[row+1];
                    x3 = x[row+2];
                    x4 = x[row+3];
                    x5 = x[row+4];
                    tmp1 = x1*bdiag[cnt] + x2*bdiag[cnt+5] + x3*bdiag[cnt+10] + x4*bdiag[cnt+15] + x5*bdiag[cnt+20];
                    tmp2 = x1*bdiag[cnt+1] + x2*bdiag[cnt+6] + x3*bdiag[cnt+11] + x4*bdiag[cnt+16] + x5*bdiag[cnt+21];
                    tmp3 = x1*bdiag[cnt+2] + x2*bdiag[cnt+7] + x3*bdiag[cnt+12] + x4*bdiag[cnt+17] + x5*bdiag[cnt+22];
                    tmp4 = x1*bdiag[cnt+3] + x2*bdiag[cnt+8] + x3*bdiag[cnt+13] + x4*bdiag[cnt+18] + x5*bdiag[cnt+23];
                    tmp5 = x1*bdiag[cnt+4] + x2*bdiag[cnt+9] + x3*bdiag[cnt+14] + x4*bdiag[cnt+19] + x5*bdiag[cnt+24];
                    x[row++] = tmp1;
                    x[row++] = tmp2;
                    x[row++] = tmp3;
                    x[row++] = tmp4;
                    x[row++] = tmp5;
                    cnt += 25;
                    break;
                default:
                    SETERRQ1(PETSC_ERR_SUP,"Inode size %D not supported",sizes[i]);
                }
            }
            ierr = PetscLogFlops(m);
            CHKERRQ(ierr);
        }
        if (flag & SOR_BACKWARD_SWEEP || flag & SOR_LOCAL_BACKWARD_SWEEP) {

            ibdiag = a->inode.ibdiag+a->inode.bdiagsize;
            for (i=m-1, row=A->rmap.n-1; i>=0; i--) {
                ibdiag -= sizes[i]*sizes[i];
                sz      = ii[row+1] - diag[row] - 1;
                v1      = a->a + diag[row] + 1;
                idx     = a->j + diag[row] + 1;

                /* see comments for MatMult_Inode() for how this is coded */
                switch (sizes[i]) {
                case 1:

                    sum1  = xb[row];
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= *v1*tmp0;
                    }
                    x[row--] = sum1*(*ibdiag);
                    break;

                case 2:

                    sum1  = xb[row];
                    sum2  = xb[row-1];
                    /* note that sum1 is associated with the second of the two rows */
                    v2    = a->a + diag[row-1] + 2;
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= *v1*tmp0;
                        sum2 -= *v2*tmp0;
                    }
                    x[row--] = sum2*ibdiag[1] + sum1*ibdiag[3];
                    x[row--] = sum2*ibdiag[0] + sum1*ibdiag[2];
                    break;
                case 3:

                    sum1  = xb[row];
                    sum2  = xb[row-1];
                    sum3  = xb[row-2];
                    v2    = a->a + diag[row-1] + 2;
                    v3    = a->a + diag[row-2] + 3;
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                        sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                        v3 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= *v1*tmp0;
                        sum2 -= *v2*tmp0;
                        sum3 -= *v3*tmp0;
                    }
                    x[row--] = sum3*ibdiag[2] + sum2*ibdiag[5] + sum1*ibdiag[8];
                    x[row--] = sum3*ibdiag[1] + sum2*ibdiag[4] + sum1*ibdiag[7];
                    x[row--] = sum3*ibdiag[0] + sum2*ibdiag[3] + sum1*ibdiag[6];
                    break;
                case 4:

                    sum1  = xb[row];
                    sum2  = xb[row-1];
                    sum3  = xb[row-2];
                    sum4  = xb[row-3];
                    v2    = a->a + diag[row-1] + 2;
                    v3    = a->a + diag[row-2] + 3;
                    v4    = a->a + diag[row-3] + 4;
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                        sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                        v3 += 2;
                        sum4 -= v4[0] * tmp0 + v4[1] * tmp1;
                        v4 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= *v1*tmp0;
                        sum2 -= *v2*tmp0;
                        sum3 -= *v3*tmp0;
                        sum4 -= *v4*tmp0;
                    }
                    x[row--] = sum4*ibdiag[3] + sum3*ibdiag[7] + sum2*ibdiag[11] + sum1*ibdiag[15];
                    x[row--] = sum4*ibdiag[2] + sum3*ibdiag[6] + sum2*ibdiag[10] + sum1*ibdiag[14];
                    x[row--] = sum4*ibdiag[1] + sum3*ibdiag[5] + sum2*ibdiag[9] + sum1*ibdiag[13];
                    x[row--] = sum4*ibdiag[0] + sum3*ibdiag[4] + sum2*ibdiag[8] + sum1*ibdiag[12];
                    break;
                case 5:

                    sum1  = xb[row];
                    sum2  = xb[row-1];
                    sum3  = xb[row-2];
                    sum4  = xb[row-3];
                    sum5  = xb[row-4];
                    v2    = a->a + diag[row-1] + 2;
                    v3    = a->a + diag[row-2] + 3;
                    v4    = a->a + diag[row-3] + 4;
                    v5    = a->a + diag[row-4] + 5;
                    for(n = 0; n<sz-1; n+=2) {
                        i1   = idx[0];
                        i2   = idx[1];
                        idx += 2;
                        tmp0 = x[i1];
                        tmp1 = x[i2];
                        sum1 -= v1[0] * tmp0 + v1[1] * tmp1;
                        v1 += 2;
                        sum2 -= v2[0] * tmp0 + v2[1] * tmp1;
                        v2 += 2;
                        sum3 -= v3[0] * tmp0 + v3[1] * tmp1;
                        v3 += 2;
                        sum4 -= v4[0] * tmp0 + v4[1] * tmp1;
                        v4 += 2;
                        sum5 -= v5[0] * tmp0 + v5[1] * tmp1;
                        v5 += 2;
                    }

                    if (n == sz-1) {
                        tmp0  = x[*idx];
                        sum1 -= *v1*tmp0;
                        sum2 -= *v2*tmp0;
                        sum3 -= *v3*tmp0;
                        sum4 -= *v4*tmp0;
                        sum5 -= *v5*tmp0;
                    }
                    x[row--] = sum5*ibdiag[4] + sum4*ibdiag[9] + sum3*ibdiag[14] + sum2*ibdiag[19] + sum1*ibdiag[24];
                    x[row--] = sum5*ibdiag[3] + sum4*ibdiag[8] + sum3*ibdiag[13] + sum2*ibdiag[18] + sum1*ibdiag[23];
                    x[row--] = sum5*ibdiag[2] + sum4*ibdiag[7] + sum3*ibdiag[12] + sum2*ibdiag[17] + sum1*ibdiag[22];
                    x[row--] = sum5*ibdiag[1] + sum4*ibdiag[6] + sum3*ibdiag[11] + sum2*ibdiag[16] + sum1*ibdiag[21];
                    x[row--] = sum5*ibdiag[0] + sum4*ibdiag[5] + sum3*ibdiag[10] + sum2*ibdiag[15] + sum1*ibdiag[20];
                    break;
                default:
                    SETERRQ1(PETSC_ERR_SUP,"Inode size %D not supported",sizes[i]);
                }
            }

            ierr = PetscLogFlops(a->nz);
            CHKERRQ(ierr);
        }
        its--;
    }
    if (its) SETERRQ(PETSC_ERR_SUP,"Currently no support for multiply SOR sweeps using inode version of AIJ matrix format;\n run with the option -mat_no_inode");
    ierr = VecRestoreArray(xx,&x);
    CHKERRQ(ierr);
    if (bb != xx) {
        ierr = VecRestoreArray(bb,(PetscScalar**)&b);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}


/*
    samestructure indicates that the matrix has not changed its nonzero structure so we
    do not need to recompute the inodes
*/
#undef __FUNCT__
#define __FUNCT__ "Mat_CheckInode"
PetscErrorCode Mat_CheckInode(Mat A,PetscTruth samestructure) {
    Mat_SeqAIJ     *a = (Mat_SeqAIJ*)A->data;
    PetscErrorCode ierr;
    PetscInt       i,j,m,nzx,nzy,*idx,*idy,*ns,*ii,node_count,blk_size;
    PetscTruth     flag;

    PetscFunctionBegin;
    if (!a->inode.use)                     PetscFunctionReturn(0);
    if (a->inode.checked && samestructure) PetscFunctionReturn(0);


    m = A->rmap.n;
    if (a->inode.size) {
        ns = a->inode.size;
    } else {
        ierr = PetscMalloc((m+1)*sizeof(PetscInt),&ns);
        CHKERRQ(ierr);
    }

    i          = 0;
    node_count = 0;
    idx        = a->j;
    ii         = a->i;
    while (i < m) {               /* For each row */
        nzx = ii[i+1] - ii[i];       /* Number of nonzeros */
        /* Limits the number of elements in a node to 'a->inode.limit' */
        for (j=i+1,idy=idx,blk_size=1; j<m && blk_size <a->inode.limit; ++j,++blk_size) {
            nzy     = ii[j+1] - ii[j]; /* Same number of nonzeros */
            if (nzy != nzx) break;
            idy  += nzx;             /* Same nonzero pattern */
            ierr = PetscMemcmp(idx,idy,nzx*sizeof(PetscInt),&flag);
            CHKERRQ(ierr);
            if (!flag) break;
        }
        ns[node_count++] = blk_size;
        idx += blk_size*nzx;
        i    = j;
    }
    /* If not enough inodes found,, do not use inode version of the routines */
    if (!a->inode.size && m && node_count > .9*m) {
        ierr = PetscFree(ns);
        CHKERRQ(ierr);
        a->inode.node_count     = 0;
        a->inode.size           = PETSC_NULL;
        a->inode.use            = PETSC_FALSE;
        ierr = PetscInfo2(A,"Found %D nodes out of %D rows. Not using Inode routines\n",node_count,m);
        CHKERRQ(ierr);
    } else {
        A->ops->mult            = MatMult_Inode;
        A->ops->relax           = MatRelax_Inode;
        A->ops->multadd         = MatMultAdd_Inode;
        A->ops->solve           = MatSolve_Inode;
        A->ops->lufactornumeric = MatLUFactorNumeric_Inode;
        A->ops->getrowij        = MatGetRowIJ_Inode;
        A->ops->restorerowij    = MatRestoreRowIJ_Inode;
        A->ops->getcolumnij     = MatGetColumnIJ_Inode;
        A->ops->restorecolumnij = MatRestoreColumnIJ_Inode;
        A->ops->coloringpatch   = MatColoringPatch_Inode;
        a->inode.node_count     = node_count;
        a->inode.size           = ns;
        ierr = PetscInfo3(A,"Found %D nodes of %D. Limit used: %D. Using Inode routines\n",node_count,m,a->inode.limit);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

/*
     This is really ugly. if inodes are used this replaces the
  permutations with ones that correspond to rows/cols of the matrix
  rather then inode blocks
*/
#undef __FUNCT__
#define __FUNCT__ "MatInodeAdjustForInodes"
PetscErrorCode PETSCMAT_DLLEXPORT MatInodeAdjustForInodes(Mat A,IS *rperm,IS *cperm) {
    PetscErrorCode ierr,(*f)(Mat,IS*,IS*);

    PetscFunctionBegin;
    ierr = PetscObjectQueryFunction((PetscObject)A,"MatInodeAdjustForInodes_C",(void (**)(void))&f);
    CHKERRQ(ierr);
    if (f) {
        ierr = (*f)(A,rperm,cperm);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatAdjustForInodes_Inode"
PetscErrorCode PETSCMAT_DLLEXPORT MatInodeAdjustForInodes_Inode(Mat A,IS *rperm,IS *cperm) {
    Mat_SeqAIJ      *a=(Mat_SeqAIJ*)A->data;
    PetscErrorCode ierr;
    PetscInt       m = A->rmap.n,n = A->cmap.n,i,j,*ridx,*cidx,nslim_row = a->inode.node_count;
    PetscInt       row,col,*permr,*permc,*ns_row =  a->inode.size,*tns,start_val,end_val,indx;
    PetscInt       nslim_col,*ns_col;
    IS             ris = *rperm,cis = *cperm;

    PetscFunctionBegin;
    if (!a->inode.size) PetscFunctionReturn(0); /* no inodes so return */
    if (a->inode.node_count == m) PetscFunctionReturn(0); /* all inodes are of size 1 */

    ierr  = Mat_CreateColInode(A,&nslim_col,&ns_col);
    CHKERRQ(ierr);
    ierr  = PetscMalloc((((nslim_row>nslim_col)?nslim_row:nslim_col)+1)*sizeof(PetscInt),&tns);
    CHKERRQ(ierr);
    ierr  = PetscMalloc((m+n+1)*sizeof(PetscInt),&permr);
    CHKERRQ(ierr);
    permc = permr + m;

    ierr  = ISGetIndices(ris,&ridx);
    CHKERRQ(ierr);
    ierr  = ISGetIndices(cis,&cidx);
    CHKERRQ(ierr);

    /* Form the inode structure for the rows of permuted matric using inv perm*/
    for (i=0,tns[0]=0; i<nslim_row; ++i) tns[i+1] = tns[i] + ns_row[i];

    /* Construct the permutations for rows*/
    for (i=0,row = 0; i<nslim_row; ++i) {
        indx      = ridx[i];
        start_val = tns[indx];
        end_val   = tns[indx + 1];
        for (j=start_val; j<end_val; ++j,++row) permr[row]= j;
    }

    /* Form the inode structure for the columns of permuted matrix using inv perm*/
    for (i=0,tns[0]=0; i<nslim_col; ++i) tns[i+1] = tns[i] + ns_col[i];

    /* Construct permutations for columns */
    for (i=0,col=0; i<nslim_col; ++i) {
        indx      = cidx[i];
        start_val = tns[indx];
        end_val   = tns[indx + 1];
        for (j = start_val; j<end_val; ++j,++col) permc[col]= j;
    }

    ierr = ISCreateGeneral(PETSC_COMM_SELF,n,permr,rperm);
    CHKERRQ(ierr);
    ISSetPermutation(*rperm);
    ierr = ISCreateGeneral(PETSC_COMM_SELF,n,permc,cperm);
    CHKERRQ(ierr);
    ISSetPermutation(*cperm);

    ierr  = ISRestoreIndices(ris,&ridx);
    CHKERRQ(ierr);
    ierr  = ISRestoreIndices(cis,&cidx);
    CHKERRQ(ierr);

    ierr = PetscFree(ns_col);
    CHKERRQ(ierr);
    ierr = PetscFree(permr);
    CHKERRQ(ierr);
    ierr = ISDestroy(cis);
    CHKERRQ(ierr);
    ierr = ISDestroy(ris);
    CHKERRQ(ierr);
    ierr = PetscFree(tns);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "MatInodeGetInodeSizes"
/*@C
   MatInodeGetInodeSizes - Returns the inode information of the Inode matrix.

   Collective on Mat

   Input Parameter:
.  A - the Inode matrix or matrix derived from the Inode class -- e.g., SeqAIJ

   Output Parameter:
+  node_count - no of inodes present in the matrix.
.  sizes      - an array of size node_count,with sizes of each inode.
-  limit      - the max size used to generate the inodes.

   Level: advanced

   Notes: This routine returns some internal storage information
   of the matrix, it is intended to be used by advanced users.
   It should be called after the matrix is assembled.
   The contents of the sizes[] array should not be changed.
   PETSC_NULL may be passed for information not requested.

.keywords: matrix, seqaij, get, inode

.seealso: MatGetInfo()
@*/
PetscErrorCode PETSCMAT_DLLEXPORT MatInodeGetInodeSizes(Mat A,PetscInt *node_count,PetscInt *sizes[],PetscInt *limit) {
    PetscErrorCode ierr,(*f)(Mat,PetscInt*,PetscInt*[],PetscInt*);

    PetscFunctionBegin;
    if (!A->assembled) SETERRQ(PETSC_ERR_ARG_WRONGSTATE,"Not for unassembled matrix");
    ierr = PetscObjectQueryFunction((PetscObject)A,"MatInodeGetInodeSizes_C",(void (**)(void))&f);
    CHKERRQ(ierr);
    if (f) {
        ierr = (*f)(A,node_count,sizes,limit);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatInodeGetInodeSizes_Inode"
PetscErrorCode PETSCMAT_DLLEXPORT MatInodeGetInodeSizes_Inode(Mat A,PetscInt *node_count,PetscInt *sizes[],PetscInt *limit) {
    Mat_SeqAIJ *a = (Mat_SeqAIJ*)A->data;

    PetscFunctionBegin;
    if (node_count) *node_count = a->inode.node_count;
    if (sizes)      *sizes      = a->inode.size;
    if (limit)      *limit      = a->inode.limit;
    PetscFunctionReturn(0);
}
EXTERN_C_END

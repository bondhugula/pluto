#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <polylib/polylib.h>


void Union_Read( Polyhedron **P, Polyhedron **C, char ***param_name )
{
	Matrix *pm;
	Polyhedron *ptmp;
	unsigned NbRows, NbColumns;
	char s[1025], param[1025];
	int i, j, c, f;

	*P = NULL;
	pm = Matrix_Read();
	f=1;
	while( f )
	{
		do
		{
			if( fgets(s, 1024, stdin) == 0  )
				f=0;
		}
		while ( (*s=='#' || *s=='\n') && f );

		if( f && sscanf(s, "%d %d", &NbRows, &NbColumns)==2 )
		{
			/* gets old pm and add it to the union */
			if( *P )
				if( pm->NbColumns != ((*P)->Dimension)+2 )
				{
					fprintf( stderr,
						"Polyhedra must be in the same dimension space !\n");
					exit(0);
				}
			ptmp = Constraints2Polyhedron(pm, 200);
			ptmp->next = *P;
			*P = ptmp;
			Matrix_Free(pm);

			/* reads the new pm */
			pm = Matrix_Alloc(NbRows, NbColumns);
			Matrix_Read_Input( pm );
		}
		else
			break;
	}

	/* Context : last read pm */
	*C = Constraints2Polyhedron(pm, 200);
	Matrix_Free(pm);


	if( f )
	{
		/* read the parameter names */
		*param_name = (char **)malloc( (*C)->Dimension*sizeof(char *) );
		c = 0;
		for( i=0 ; i<(*C)->Dimension ; ++i )
		{
			j=0;
			for( ; ; ++c )
			{
				if( s[c]==' ' || s[c]=='\n' || s[c]==0 ) {
					if( j==0 )
						continue;
					else
						break;
				}
				param[j++] = s[c];
			}

			/* else, no parameters (use default) */
			if( j==0 )
				break;
			param[j] = 0;
			(*param_name)[i] = (char *)malloc( j );
			strcpy( (*param_name)[i], param );
		}
		if( i != (*C)->Dimension )
		{
			free( *param_name );
			*param_name = Read_ParamNames(NULL,(*C)->Dimension);
		}
	}
	else
		*param_name = Read_ParamNames(NULL,(*C)->Dimension);

}

void recurse ( 	Polyhedron *C, char **param_name, Enumeration *e, 
                  Value *pmin, Value *pmax, Value *p, int l )
{
	Value z, *tmp; int k;
	value_init( z );
	if( l == C->Dimension )
	{
		fprintf(stdout,"EP( ");
		value_print(stdout,VALUE_FMT,p[0]);
		for(k=1;k<C->Dimension;++k) {
		  fprintf(stdout,",");
		  value_print(stdout,VALUE_FMT,p[k]);
		}  
		fprintf(stdout," ) = ");
		value_print(stdout,VALUE_FMT,*(tmp=compute_poly(e,p)));
		value_clear( *tmp );
		free(tmp);
		fprintf(stdout,"\n");  
	}
	else
	{
		for( value_assign( z, pmin[l]) ; value_le(z,pmax[l]) ; value_increment(z,z) )
		{
			value_assign( p[l], z );
			recurse ( 	C, param_name, e, pmin, pmax, p, l+1 );
		}
	}

}



int main( int argc, char **argv)
{
	Polyhedron *P, *C;
	char **param_name;
	Enumeration *e, *en;
	Value *pmin, *pmax, *p; int i, k; char str[256], *s;

	if( argc != 1 )
	{
		fprintf( stderr, " Usage : %s [< file]\n", argv[0] );
		return( -1 );
	}

	Union_Read( &P, &C, &param_name );

	e = Domain_Enumerate( P, C, 200, param_name );

	for( en=e ; en ; en=en->next )
	{
	  Print_Domain(stdout,en->ValidityDomain, param_name);
	  print_evalue(stdout,&en->EP, param_name);
	  printf( "\n-----------------------------------\n" );
	}

    if( isatty(0) && C->Dimension != 0)
        {  /* no tty input or no polyhedron -> no evaluation. */
            printf("Evaluation of the Ehrhart polynomial :\n");
            pmin = (Value *)malloc(sizeof(Value) * (C->Dimension));
            pmax = (Value *)malloc(sizeof(Value) * (C->Dimension));
            p = (Value *)malloc(sizeof(Value) * (C->Dimension));
            for(i=0;i<C->Dimension;i++) 
            {
               value_init(pmin[i]);
               value_init(pmax[i]);
               value_init(p[i]);
            }
            FOREVER {
                fflush(stdin);
                printf("Enter %d parameters (or intervals, comma separated) : ",C->Dimension);
                for(k=0;k<C->Dimension;++k)
                {
                    scanf("%s",str);
                    if( (s=strpbrk(str,",")) )
                    {  *s = 0;
                       value_read(pmin[k],str);
                       value_read(pmax[k],(s+1));
                    }
                    else
                    {  value_read(pmin[k],str);
                       value_assign(pmax[k],pmin[k]);
                    }
                }

                recurse( C, param_name, e, pmin, pmax, p, 0 );

            }
        }

	return( 0 );
}


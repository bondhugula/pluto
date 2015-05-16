
#define N 256
int main(){
	int i,j;
	double A[N+1][N+1];
#pragma scop
	for (i=0;i<N-1;i++){
		for(j=0;j<N-1;j++){
			A[i+1][j+1]=0.6*A[i+1][j+1]+0.4*A[i][j];
		}
	}
#pragma endscop
	return 0;
}



/*@ begin Loop(
 transform RegTile(loops=['i','j','k'], ufactors=[2,2,2])
 for (i=0; i<=M-1; i++)
 {
   for (j=0; j<=N-1; j++)
   {
     for (k=i; k<=O-1; k++)
     {
       S(i,j,k);
     }
   }
 }
) @*/{
    for (it=0; it<=M-2; it=it+2) {
        for (jt=0; jt<=N-2; jt=jt+2) {
            newlb_k=-2147483648;
            newub_k=O-1;
            for (i=it; i<=it+1; i=i+1)
                newlb_k=max(newlb_k,i);
            for (i=it; i<=it+1; i=i+1)
                for (k=i; k<=newlb_k-1; k=k+1) {
                    S(i,jt,k);
                    S(i,(jt+1),k);
                }
            for (kt=newlb_k; kt<=newub_k-1; kt=kt+2) {
                S(it,jt,kt);
                S(it,jt,(kt+1));
                S(it,(jt+1),kt);
                S(it,(jt+1),(kt+1));
                S((it+1),jt,kt);
                S((it+1),jt,(kt+1));
                S((it+1),(jt+1),kt);
                S((it+1),(jt+1),(kt+1));
            }
            for (k=kt; k<=newub_k; k=k+1) {
                S(it,jt,k);
                S(it,(jt+1),k);
                S((it+1),jt,k);
                S((it+1),(jt+1),k);
            }
            for (i=it; i<=it+1; i=i+1)
                for (k=newub_k+1; k<=O-1; k=k+1) {
                    S(i,jt,k);
                    S(i,(jt+1),k);
                }
        }
        for (j=jt; j<=N-1; j=j+1) {
            newlb_k=-2147483648;
            newub_k=O-1;
            for (i=it; i<=it+1; i=i+1)
                newlb_k=max(newlb_k,i);
            for (i=it; i<=it+1; i=i+1)
                for (k=i; k<=newlb_k-1; k=k+1) {
                    S(i,j,k);
                }
            for (kt=newlb_k; kt<=newub_k-1; kt=kt+2) {
                S(it,j,kt);
                S(it,j,(kt+1));
                S((it+1),j,kt);
                S((it+1),j,(kt+1));
            }
            for (k=kt; k<=newub_k; k=k+1) {
                S(it,j,k);
                S((it+1),j,k);
            }
            for (i=it; i<=it+1; i=i+1)
                for (k=newub_k+1; k<=O-1; k=k+1) {
                    S(i,j,k);
                }
        }
    }
    for (i=it; i<=M-1; i=i+1) {
        for (jt=0; jt<=N-2; jt=jt+2) {
            for (kt=i; kt<=O-2; kt=kt+2) {
                S(i,jt,kt);
                S(i,jt,(kt+1));
                S(i,(jt+1),kt);
                S(i,(jt+1),(kt+1));
            }
            for (k=kt; k<=O-1; k=k+1) {
                S(i,jt,k);
                S(i,(jt+1),k);
            }
        }
        for (j=jt; j<=N-1; j=j+1) {
            for (kt=i; kt<=O-2; kt=kt+2) {
                S(i,j,kt);
                S(i,j,(kt+1));
            }
            for (k=kt; k<=O-1; k=k+1) {
                S(i,j,k);
            }
        }
    }
}
/*@ end @*/



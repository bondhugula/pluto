for (v1=0; v1<=V-1; v1=v1+1)
    for (v2=0; v2<=V-1; v2=v2+1)
        for (o1=0; o1<=O-1; o1=o1+1)
            for (o2=0; o2<=O-1; o2=o2+1)
                for (ox=0; ox<=O-1; ox=ox+1)
                    R[v1][v2][o1][o2]=R[v1][v2][o1][o2]+T[v1][ox][o1][o2]*A2[v2][ox];

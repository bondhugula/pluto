constant n, sqrt2, w;

    for(i=0;i<n;i++)
        vecp[i] = 0;

    while(w>1)
    {
        /* w = w / 2; */
        for(i=0;i<w;i++)
        {
            vecp[i] = (vec[2*i] + vec[2*i+1])/sqrt2;
            vecp[i+w] = (vec[2*i] - vec[2*i+1])/sqrt2;
        }

        for(i=0;i<(w*2);i++)
            vec[i] = vecp[i]; 
    }

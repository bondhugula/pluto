/* pluto start (Cz,Cym,Cxm) */
for(iz=0; iz<Cz; iz++) {
    for(iy=0; iy<Cym; iy++) {
        for(ix=0; ix<Cxm; ix++) {
            clf[iz][iy]=Ex[iz][iy][ix]-Ex[iz][iy+1][ix]+Ey[iz][iy][ix+1]-Ey[iz][iy][ix];
            tmp[iz][iy]=(cymh[iy]/cyph[iy])*Bza[iz][iy][ix]-(ch/cyph[iy])*clf[iz][iy];
            Hz[iz][iy][ix]=(cxmh[ix]/cxph[ix])*Hz[iz][iy][ix]
                           +(mui*czp[iz]/cxph[ix])*tmp[iz][iy]
                           -(mui*czm[iz]/cxph[ix])*Bza[iz][iy][ix];
            Bza[iz][iy][ix]=tmp[iz][iy];
        }
        clf[iz][iy]=Ex[iz][iy][Cxm]-Ex[iz][iy+1][Cxm]+Ry[iz][iy]-Ey[iz][iy][Cxm];
        tmp[iz][iy]=(cymh[iy]/cyph[iy])*Bza[iz][iy][Cxm]-(ch/cyph[iy])*clf[iz][iy];
        Hz[iz][iy][Cxm]=(cxmh[Cxm]/cxph[Cxm])*Hz[iz][iy][Cxm]
                        +(mui*czp[iz]/cxph[Cxm])*tmp[iz][iy]
                        -(mui*czm[iz]/cxph[Cxm])*Bza[iz][iy][Cxm];
        Bza[iz][iy][Cxm]=tmp[iz][iy];
    }
    for(ix=0; ix<Cxm; ix++) {
        clf[iz][iy]=Ex[iz][Cym][ix]-Ax[iz][ix]+Ey[iz][Cym][ix+1]-Ey[iz][Cym][ix];
        tmp[iz][iy]=(cymh[Cym]/cyph[iy])*Bza[iz][iy][ix]-(ch/cyph[iy])*clf[iz][iy];
        Hz[iz][Cym][ix]=(cxmh[ix]/cxph[ix])*Hz[iz][Cym][ix]
                        +(mui*czp[iz]/cxph[ix])*tmp[iz][iy]
                        -(mui*czm[iz]/cxph[ix])*Bza[iz][Cym][ix];
        Bza[iz][Cym][ix]=tmp[iz][iy];
    }
    clf[iz][iy]=Ex[iz][Cym][Cxm]-Ax[iz][Cxm]+Ry[iz][Cym]-Ey[iz][Cym][Cxm];
    tmp[iz][iy]=(cymh[Cym]/cyph[Cym])*Bza[iz][Cym][Cxm]-(ch/cyph[Cym])*clf[iz][iy];
    Hz[iz][Cym][Cxm]=(cxmh[Cxm]/cxph[Cxm])*Hz[iz][Cym][Cxm]
                     +(mui*czp[iz]/cxph[Cxm])*tmp[iz][iy]
                     -(mui*czm[iz]/cxph[Cxm])*Bza[iz][Cym][Cxm];
    Bza[iz][Cym][Cxm]=tmp[iz][iy];
}
/* pluto end */

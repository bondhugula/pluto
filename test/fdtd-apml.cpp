#include "pfdtd.h"

#ifdef APML
 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <time.h>
  #ifdef MPIRUN
    #ifdef TESTMPI 
      #include "mpi.h"
    #else
      #include <mpi.h>
  #endif
#endif

/* Required Microsoft include files */
#ifdef WIN95
  #include <io.h>
#endif
#include "prototype.h"


/*****************************************************************************/

void setapmlid(Global_Block *gblock)
{
  int ii,ix,iy,iz;
 
  // Set up apmlid
  ii=0;
  iz=0;
  for (iy=0;iy<gblock->BBy;iy++){
    for (ix=0;ix<gblock->BBx;ix++){
      if(iy==0){
	if(ix==0){
	  gblock->iblock[ii].apmlid=0;
	}
	else if(ix==gblock->BBx-1){
	  gblock->iblock[ii].apmlid=2;
	}
	else {
	  gblock->iblock[ii].apmlid=1;
	}
      }
      else if(iy==gblock->BBy-1){
	if(ix==0){
	  gblock->iblock[ii].apmlid=6;
	}
	else if(ix==gblock->BBx-1){
	  gblock->iblock[ii].apmlid=8;
	}
	else {
	  gblock->iblock[ii].apmlid=7;
	}
      }
      else {
	if(ix==0){
	  gblock->iblock[ii].apmlid=3;
	}
	else if(ix==gblock->BBx-1){
	  gblock->iblock[ii].apmlid=5;
	}
	else {
	  gblock->iblock[ii].apmlid=4;
	}
      }
      ii++;
    }
  }

  for(iz=1;iz<gblock->BBz-1;iz++){
    for (iy=0;iy<gblock->BBy;iy++){
      for (ix=0;ix<gblock->BBx;ix++){
	if(iy==0){
	  if(ix==0){
	    gblock->iblock[ii].apmlid=9;
	  }
	  else if(ix==gblock->BBx-1){
	    gblock->iblock[ii].apmlid=11;
	  }
	  else {
	    gblock->iblock[ii].apmlid=10;
	  }
	}
	else if(iy==gblock->BBy-1){
	  if(ix==0){
	    gblock->iblock[ii].apmlid=15;
	  }
	  else if(ix==gblock->BBx-1){
	    gblock->iblock[ii].apmlid=17;
	  }
	  else {
	    gblock->iblock[ii].apmlid=16;
	  }
	}
	else {
	  if(ix==0){
	    gblock->iblock[ii].apmlid=12;
	  }
	  else if(ix==gblock->BBx-1){
	    gblock->iblock[ii].apmlid=14;
	  }
	  else {
	    gblock->iblock[ii].apmlid=13;
	  }
	}
	ii++;
      }
    }
  }

  iz=gblock->BBz-1;
  for (iy=0;iy<gblock->BBy;iy++){
    for (ix=0;ix<gblock->BBx;ix++){
      if(iy==0){
	if(ix==0){
	  gblock->iblock[ii].apmlid=18;
	}
	else if(ix==gblock->BBx-1){
	  gblock->iblock[ii].apmlid=20;
	}
	else {
	  gblock->iblock[ii].apmlid=19;
	}
      }
      else if(iy==gblock->BBy-1){
	if(ix==0){
	  gblock->iblock[ii].apmlid=24;
	}
	else if(ix==gblock->BBx-1){
	  gblock->iblock[ii].apmlid=26;
	}
	else {
	  gblock->iblock[ii].apmlid=25;
	}
      }
      else {
	if(ix==0){
	  gblock->iblock[ii].apmlid=21;
	}
	else if(ix==gblock->BBx-1){
	  gblock->iblock[ii].apmlid=23;
	}
	else {
	  gblock->iblock[ii].apmlid=22;
	}
      }
      ii++;
    }
  }
  
}


/*****************************************************************************/

void define_apml_space(I_Block *iblock)
{
  REALNO ***Ex, ***Ey, ***Ez;
  REALNO ***Hx, ***Hy, ***Hz;
  REALNO ***Px, ***Py, ***Pz;
  REALNO ***Pxp, ***Pyp, ***Pzp;
  REALNO ***Bxa, ***Bya, ***Bza;
  REALNO *cxm,*cxp,*cym,*cyp,*czm,*czp;
  REALNO *cxmh,*cxph,*cymh,*cyph,*czmh,*czph;
  //REALNO ce;
  REALNO **Tx, **Ty, **Bx, **By;
  REALNO **Ry, **Rz, **Ly, **Lz;
  REALNO **Fx, **Fz, **Ax, **Az;
  int ix, iy, iz;
  //REALNO tmp1,tmp2,tmp3;
  REALNO tmp1,tmp3;

  /* shortcut variable names */
  int Cx=iblock->length.x;
  int Cy=iblock->length.y;
  int Cz=iblock->length.z;
  REALNO dt=iblock->gblock->dt;
  //REALNO dx=iblock->gblock->dx;
  int ID=iblock->gblock->ID[iblock->id];
  
  iblock->apml.cxm=(REALNO *)(calloc(Cx,sizeof(REALNO )));
  iblock->apml.cxp=(REALNO *)(calloc(Cx,sizeof(REALNO )));
  iblock->apml.cym=(REALNO *)(calloc(Cy,sizeof(REALNO )));
  iblock->apml.cyp=(REALNO *)(calloc(Cy,sizeof(REALNO )));
  iblock->apml.czm=(REALNO *)(calloc(Cz,sizeof(REALNO )));
  iblock->apml.czp=(REALNO *)(calloc(Cz,sizeof(REALNO )));
  iblock->apml.cxmh=(REALNO *)(calloc(Cx,sizeof(REALNO )));
  iblock->apml.cxph=(REALNO *)(calloc(Cx,sizeof(REALNO )));
  iblock->apml.cymh=(REALNO *)(calloc(Cy,sizeof(REALNO )));
  iblock->apml.cyph=(REALNO *)(calloc(Cy,sizeof(REALNO )));
  iblock->apml.czmh=(REALNO *)(calloc(Cz,sizeof(REALNO )));
  iblock->apml.czph=(REALNO *)(calloc(Cz,sizeof(REALNO )));

  cxm=iblock->apml.cxm; cxp=iblock->apml.cxp;
  cym=iblock->apml.cym; cyp=iblock->apml.cyp;
  czm=iblock->apml.czm; czp=iblock->apml.czp;
  cxmh=iblock->apml.cxmh; cxph=iblock->apml.cxph;
  cymh=iblock->apml.cymh; cyph=iblock->apml.cyph;
  czmh=iblock->apml.czmh; czph=iblock->apml.czph;
  // Setting X profile
  for(ix=0;ix<Cx;ix++){
    cxm[ix]=2*CNST_ep; cxp[ix]=2*CNST_ep;
    cxmh[ix]=2*CNST_ep; cxph[ix]=2*CNST_ep;
  }
  if (ID & 2){ // left PML
    for (ix=0; ix<Cx; ++ix){
      tmp3=(REALNO)(Cx-ix)/Cx;
      tmp1=tmp3*tmp3*iblock->gblock->SGx;
      cxm[ix]=2*CNST_ep-tmp1*dt;
      cxp[ix]=2*CNST_ep+tmp1*dt;
      tmp3=(REALNO)(Cx-ix-0.5)/Cx;
      tmp1=tmp3*tmp3*iblock->gblock->SGx;
      cxmh[ix]=2*CNST_ep-tmp1*dt;
      cxph[ix]=2*CNST_ep+tmp1*dt;
    }
  }
  if (ID & 4){ // right PML
    for (ix=0; ix<Cx; ++ix){
      tmp3=(REALNO)(ix+0.5)/Cx;
      tmp1=tmp3*tmp3*iblock->gblock->SGx;
      cxmh[ix]=2*CNST_ep-tmp1*dt;
      cxph[ix]=2*CNST_ep+tmp1*dt;
      tmp3=(REALNO)(ix)/Cx;
      tmp1=tmp3*tmp3*iblock->gblock->SGx;
      cxm[ix]=2*CNST_ep-tmp1*dt;
      cxp[ix]=2*CNST_ep+tmp1*dt;
    }
  }
  // Setting Y profile
  for(iy=0;iy<Cy;iy++){
    cym[iy]=2*CNST_ep; cyp[iy]=2*CNST_ep;
    cymh[iy]=2*CNST_ep; cyph[iy]=2*CNST_ep;
  }
  if (ID & 8){ // front PML
    for (iy=0; iy<Cy; ++iy){
      tmp3=(REALNO)(Cy-iy-0.5)/Cy;
      tmp1=tmp3*tmp3*iblock->gblock->SGy;
      cymh[iy]=2*CNST_ep-tmp1*dt;
      cyph[iy]=2*CNST_ep+tmp1*dt;
      tmp3=(REALNO)(Cy-iy)/Cy;
      tmp1=tmp3*tmp3*iblock->gblock->SGy;
      cym[iy]=2*CNST_ep-tmp1*dt;
      cyp[iy]=2*CNST_ep+tmp1*dt;
    }
  }
  if (ID & 16){ // back PML
    for (iy=0; iy<Cy; ++iy){
      tmp3=(REALNO)(iy+0.5)/Cy;
      tmp1=tmp3*tmp3*iblock->gblock->SGy;
      cymh[iy]=2*CNST_ep-tmp1*dt;
      cyph[iy]=2*CNST_ep+tmp1*dt;
      tmp3=(REALNO)(iy)/Cy;
      tmp1=tmp3*tmp3*iblock->gblock->SGy;
      cym[iy]=2*CNST_ep-tmp1*dt;
      cyp[iy]=2*CNST_ep+tmp1*dt;
    }
  }
  // Setting Z profile
  for(iz=0;iz<Cz;iz++){
    czm[iz]=2*CNST_ep; czp[iz]=2*CNST_ep;
    czmh[iz]=2*CNST_ep; czph[iz]=2*CNST_ep;
  }
  if (ID & 32){ // bottom PML
    for (iz=0; iz<Cz; ++iz){
      tmp3=(REALNO)(Cz-iz-0.5)/Cz;
      tmp1=tmp3*tmp3*iblock->gblock->SGz;
      czmh[iz]=2*CNST_ep-tmp1*dt;
      czph[iz]=2*CNST_ep+tmp1*dt;
      tmp3=(REALNO)(Cz-iz)/Cz;
      tmp1=tmp3*tmp3*iblock->gblock->SGz;
      czm[iz]=2*CNST_ep-tmp1*dt;
      czp[iz]=2*CNST_ep+tmp1*dt;
    }
  }
  if (ID & 64){ // top PML
    for (iz=0; iz<Cz; ++iz){
      tmp3=(REALNO)(iz+0.5)/Cz;
      tmp1=tmp3*tmp3*iblock->gblock->SGz;
      czmh[iz]=2*CNST_ep-tmp1*dt;
      czph[iz]=2*CNST_ep+tmp1*dt;
      tmp3=(REALNO)(iz)/Cz;
      tmp1=tmp3*tmp3*iblock->gblock->SGz;
      czm[iz]=2*CNST_ep-tmp1*dt;
      czp[iz]=2*CNST_ep+tmp1*dt;
    }
  }
  
//   printf(" bkpt : testing \n");
//   printf(" %e %e %e %e %e \n",CNST_ep,dt,iblock->gblock->SGx,iblock->gblock->SGy,iblock->gblock->SGz);
//   printf(" apmlid : %d \n",iblock->apmlid); getchar();
//   for(iz=0;iz<Cz;iz++){
//     printf(" %d   %e  %e \n",iz,czm[iz],czp[iz]);
//   }
//   getchar();
//   for(iy=0;iy<Cy;iy++){
//     printf(" %d   %e  %e \n",iy,cym[iy],cyp[iy]);
//   }
//   getchar();
//   for(ix=0;ix<Cx;ix++){
//     printf(" %d   %e  %e \n",ix,cxm[ix],cxp[ix]);
//   }
//   printf("done \n ");getchar();
  

  // setup coefficients for P's

  switch (iblock->apmlid) {
  case 0 :
    define_apml_case00(iblock);
    break;
  case 2 :
    define_apml_case02(iblock);
    break;
  case 6 :
    define_apml_case06(iblock);
    break;
  case 8 :
    define_apml_case08(iblock);
    break;
  case 18 :
    define_apml_case18(iblock);
    break;
  case 20 :
    define_apml_case20(iblock);
    break;
  case 24 :
    define_apml_case24(iblock);
    break;
  case 26 :
    define_apml_case26(iblock);
    break;
  case 1 :
    define_apml_case01(iblock);
    break;
  case 7 :
    define_apml_case07(iblock);
    break;
  case 19 :
    define_apml_case19(iblock);
    break;
  case 25 :
    define_apml_case25(iblock);
    break;
  case 3 :
    define_apml_case03(iblock);
    break;
  case 5 :
    define_apml_case05(iblock);
    break;
  case 21 :
    define_apml_case21(iblock);
    break;
  case 23 :
    define_apml_case23(iblock);
    break;
  case 9 :
    define_apml_case09(iblock);
    break;
  case 11 :
    define_apml_case11(iblock);
    break;
  case 15 :
    define_apml_case15(iblock);
    break;
  case 17 :
    define_apml_case17(iblock);
    break;
  case 4 :
    define_apml_case04(iblock);
    break;
  case 22 :
    define_apml_case22(iblock);
    break;
  case 10 :
    define_apml_case10(iblock);
    break;
  case 16 :
    define_apml_case16(iblock);
    break;
  case 12 :
    define_apml_case12(iblock);
    break;
  case 14 :
    define_apml_case14(iblock);
    break;
  default :
    printf(" Error in apmlid %d %d \n",iblock->id,iblock->apmlid);
  }

  /*
  printf(" bkpt : testing \n");
  printf(" ep:%e dt:%e \n",ep,dt); getchar(); 
  switch(iblock->apmlid){
  case 0 :
  case 2 :
  case 6 :
  case 8 :
  case 18 :
  case 20 :
  case 24 :
  case 26 :
    printf(" apmlid : %d \n",iblock->apmlid);
    printf(" csgm : %e ",iblock->apmlc.csgm);
    printf(" csgp : %e \n",iblock->apmlc.csgp);
    getchar();
    break;
  case 1 :
  case 7 :
  case 19 :
  case 25 :
    printf(" apmlid : %d \n",iblock->apmlid);
    for(ix=0;ix<iblock->length.x;ix++){
      printf(" csgm[%d] : %e ",ix,iblock->apmle.csgm[ix]);
      printf(" csgp[%d] : %e \n",ix,iblock->apmle.csgp[ix]);     
    }
    getchar();
    break;
  case 3 :
  case 5 :
  case 21 :
  case 23 :
    printf(" apmlid : %d \n",iblock->apmlid);
    for(iy=0;iy<iblock->length.y;iy++){
      printf(" csgm[%d] : %e ",iy,iblock->apmle.csgm[iy]);
      printf(" csgp[%d] : %e \n",iy,iblock->apmle.csgp[iy]);
    }
    getchar();
    break;
  case 9 :
  case 11 :
  case 15 :
  case 17 :
    printf(" apmlid : %d \n",iblock->apmlid);
    for(iz=0;iz<iblock->length.z;iz++){
      printf(" csgm[%d] : %e ",iz,iblock->apmle.csgm[iz]);
      printf(" csgp[%d] : %e \n",iz,iblock->apmle.csgp[iz]);
    }
    getchar();
    break;
  case 12 :
  case 14 :
    printf(" apmlid : %d \n",iblock->apmlid);
    for(iz=0;iz<iblock->length.z;iz++){
      printf(" \n iz = %d \n",iz);
      for(iy=0;iy<iblock->length.y;iy++){
	printf(" csgm[%d][%d] : %e ",iz,iy,iblock->apmlf.csgm[iz][iy]);
	printf(" csgp[%d][%d] : %e \n",iz,iy,iblock->apmlf.csgp[iz][iy]);
      }
    }
    getchar();
    break;
  case 10 :
  case 16 :
    printf(" apmlid : %d \n",iblock->apmlid);
    for(iz=0;iz<iblock->length.z;iz++){
      printf(" \n iz = %d \n",iz);
      for(ix=0;ix<iblock->length.x;ix++){
	printf(" csgm[%d][%d] : %e \n",iz,ix,iblock->apmlf.csgm[iz][ix]);
	printf(" csgp[%d][%d] : %e \n",iz,ix,iblock->apmlf.csgp[iz][ix]);
      }
    }
    getchar();
    break;
  case 4 :
  case 22 :
    printf(" apmlid : %d \n",iblock->apmlid);
    for(iy=0;iy<iblock->length.y;iy++){
      printf(" \n iy = %d \n",iy);
      for(ix=0;ix<iblock->length.x;ix++){
	printf(" csgm[%d][%d] : %e \n",iy,ix,iblock->apmlf.csgm[iy][ix]);
	printf(" csgp[%d][%d] : %e \n",iy,ix,iblock->apmlf.csgp[iy][ix]);
      }
    }
    getchar();
    break;
  default :
    printf(" ERROR \n");
  }
  printf(" Done P coeff \n"); getchar();
  */

  iblock->ib_type.pml=2;  // Anisotropic PML

  iblock->max.Ex=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->max.Ey=(REALNO ***)(calloc(Cz,sizeof(REALNO **)));
  iblock->max.Ez=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->max.Hx=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->max.Hy=(REALNO ***)(calloc(Cz,sizeof(REALNO **)));
  iblock->max.Hz=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->apml.Px=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->apml.Py=(REALNO ***)(calloc(Cz,sizeof(REALNO **)));
  iblock->apml.Pz=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->apml.Pxp=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->apml.Pyp=(REALNO ***)(calloc(Cz,sizeof(REALNO **)));
  iblock->apml.Pzp=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->apml.Bxa=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 
  iblock->apml.Bya=(REALNO ***)(calloc(Cz,sizeof(REALNO **)));
  iblock->apml.Bza=(REALNO ***)(calloc(Cz,sizeof(REALNO **))); 

  iblock->ib_type.pass=1;
  iblock->pass.Tx=(REALNO **)(calloc(Cy,sizeof(REALNO *))); 
  if(iblock->pass.Tx==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  iblock->pass.Ty=(REALNO **)(calloc(Cy,sizeof(REALNO *))); 
  if(iblock->pass.Ty==NULL){
    printf("OOM\n"); fflush(NULL);  exit(1);
  }
  iblock->pass.Bx=(REALNO **)(calloc(Cy,sizeof(REALNO *))); 
  if(iblock->pass.Bx==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
    }
  iblock->pass.By=(REALNO **)(calloc(Cy,sizeof(REALNO *))); 
  if(iblock->pass.By==NULL){
    printf("OOM\n");  fflush(NULL);  exit(1);
  }
  iblock->pass.Ry=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Ry==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  iblock->pass.Rz=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Rz==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  iblock->pass.Ly=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Ly==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  iblock->pass.Lz=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Lz==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  iblock->pass.Fx=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Fx==NULL){
    printf("OOM\n"); fflush(NULL);  exit(1);
  }
  iblock->pass.Fz=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Fz==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  iblock->pass.Ax=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Ax==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  iblock->pass.Az=(REALNO **)(calloc(Cz,sizeof(REALNO *))); 
  if(iblock->pass.Az==NULL){
    printf("OOM\n"); fflush(NULL); exit(1);
  }
  
  Ex=iblock->max.Ex;     Ey=iblock->max.Ey;     Ez=iblock->max.Ez;
  Hx=iblock->max.Hx;     Hy=iblock->max.Hy;     Hz=iblock->max.Hz;
  Px=iblock->apml.Px;    Py=iblock->apml.Py;    Pz=iblock->apml.Pz;
  Pxp=iblock->apml.Pxp;  Pyp=iblock->apml.Pyp;  Pzp=iblock->apml.Pzp;
  Bxa=iblock->apml.Bxa;  Bya=iblock->apml.Bya;  Bza=iblock->apml.Bza;

  Tx=iblock->pass.Tx;  Ty=iblock->pass.Ty;
  Bx=iblock->pass.Bx;  By=iblock->pass.By;
  Ry=iblock->pass.Ry;  Rz=iblock->pass.Rz;
  Ly=iblock->pass.Ly;  Lz=iblock->pass.Lz;
  Fx=iblock->pass.Fx;  Fz=iblock->pass.Fz;
  Ax=iblock->pass.Ax;  Az=iblock->pass.Az;

  for (iz=0; iz<Cz; ++iz){
    Ry[iz]=(REALNO *)(calloc(Cy,sizeof(REALNO )));
    Rz[iz]=(REALNO *)(calloc(Cy,sizeof(REALNO )));
    Ly[iz]=(REALNO *)(calloc(Cy,sizeof(REALNO )));
    Lz[iz]=(REALNO *)(calloc(Cy,sizeof(REALNO )));
    Fx[iz]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    Fz[iz]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    Ax[iz]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    Az[iz]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    if(!Az[iz]) {
      printf("[CPU %d] define_apml_space  Allocation error tag=2\n",iblock->gblock->rank);
      fflush(NULL);
      exit(1);
    }
  }
  for (iy=0; iy<Cy; ++iy){
    Tx[iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    Ty[iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    Bx[iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    By[iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
    if(!By[iy]) {
      printf("[CPU %d] define_apml_space  Allocation error  tag=3\n",iblock->gblock->rank);
      fflush(NULL);
      exit(1);
    }
  }

  /* Set to zero */
  for (iz=0; iz<Cz; ++iz){
    for (iy=0; iy<Cy; ++iy){
      Ry[iz][iy]=INITIAL;
      Rz[iz][iy]=INITIAL;
      Ly[iz][iy]=INITIAL;
      Lz[iz][iy]=INITIAL;
    }
    for (ix=0; ix<Cx; ++ix){
      Fx[iz][ix]=INITIAL;
      Fz[iz][ix]=INITIAL;
      Ax[iz][ix]=INITIAL;
      Az[iz][ix]=INITIAL;
    }
  }
  for (iy=0; iy<Cy; ++iy){
    for (ix=0; ix<Cx; ++ix){
      Tx[iy][ix]=INITIAL;
      Ty[iy][ix]=INITIAL;
      Bx[iy][ix]=INITIAL;
      By[iy][ix]=INITIAL;
    }
  }

  for (iz=0; iz<Cz; ++iz){
    Ex[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Ey[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Ez[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Hx[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Hy[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Hz[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Px[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Py[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Pz[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Pxp[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Pyp[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Pzp[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Bxa[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Bya[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    Bza[iz]=(REALNO **)(calloc(Cy,sizeof(REALNO *)));
    if(!Bza[iz]) {
      printf("[CPU %d] define_apml_space  Allocation error  tag=4\n",iblock->gblock->rank);
      fflush(NULL); exit(1);
    }
    for (iy=0; iy<Cy; ++iy){
      Ex[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Ey[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Ez[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Hx[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Hy[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Hz[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Px[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Py[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Pz[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Pxp[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Pyp[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Pzp[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Bxa[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Bya[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      Bza[iz][iy]=(REALNO *)(calloc(Cx,sizeof(REALNO )));
      if(!Bza[iz][iy]) {
	printf("[CPU %d] define_apml_space  Allocation error  tag=5\n",iblock->gblock->rank);
	fflush(NULL);exit(1);
      }
      for (ix=0; ix<Cx; ++ix){
        Ex[iz][iy][ix]=INITIAL;
	Ey[iz][iy][ix]=INITIAL;
	Ez[iz][iy][ix]=INITIAL;
	Hx[iz][iy][ix]=INITIAL;
	Hy[iz][iy][ix]=INITIAL;
	Hz[iz][iy][ix]=INITIAL;
	Px[iz][iy][ix]=INITIAL;
	Py[iz][iy][ix]=INITIAL;
	Pz[iz][iy][ix]=INITIAL;
	Pxp[iz][iy][ix]=INITIAL;
	Pyp[iz][iy][ix]=INITIAL;
	Pzp[iz][iy][ix]=INITIAL;
	Bxa[iz][iy][ix]=INITIAL;
	Bya[iz][iy][ix]=INITIAL;
	Bza[iz][iy][ix]=INITIAL;
      }
    }
  }

  
  iblock->ib_type.pml_1d=2;

  iblock->physics.inhomogenous_type=0;
  iblock->physics.def_er=1.0;
  iblock->physics.def_sigma=0.0;
  iblock->physics.def_ur=1.0;
  iblock->physics.def_mag=0.0;
  iblock->physics.input_format=0;
  iblock->physics.field_type=1;
  
}



/*****************************************************************************/

void calc_apml_E(I_Block *iblock)
{
  int ix,iy,iz,ixm,iym,izm;
  REALNO ***Ex, ***Ey, ***Ez, ***Hx, ***Hy, ***Hz;
  REALNO ***Px, ***Py, ***Pz, ***Pxp, ***Pyp, ***Pzp;
  REALNO **Bx, **By, **Ly, **Lz, **Fx, **Fz;
  REALNO *cxm, *cxp, *cym, *cyp, *czm, *czp;
  REALNO *cxmh, *cxph, *cymh, *cyph, *czmh, *czph;
  int Cx, Cy, Cz;
  REALNO ch,clf,tmp1,tmp2;
  REALNO csgmc,csgpc;
  REALNO *csgme,*csgpe,*csgmte,*csgpte;
  REALNO **csgmf,**csgpf,**csgmt1f,**csgpt1f,**csgmt2f,**csgpt2f;


  ch=2*CNST_ep*iblock->gblock->dt/iblock->gblock->dx;

  Ex =iblock->max.Ex;    Ey =iblock->max.Ey;   Ez =iblock->max.Ez;
  Hx =iblock->max.Hx;    Hy =iblock->max.Hy;   Hz =iblock->max.Hz;
  Px =iblock->apml.Px;   Py =iblock->apml.Py;  Pz =iblock->apml.Pz;
  Pxp=iblock->apml.Pxp;  Pyp=iblock->apml.Pyp; Pzp=iblock->apml.Pzp;
  cxm  = iblock->apml.cxm;  cxp  = iblock->apml.cxp;
  cym  = iblock->apml.cym;  cyp  = iblock->apml.cyp;
  czm  = iblock->apml.czm;  czp  = iblock->apml.czp;
  cxmh = iblock->apml.cxmh; cxph = iblock->apml.cxph;
  cymh = iblock->apml.cymh; cyph = iblock->apml.cyph;
  czmh = iblock->apml.czmh; czph = iblock->apml.czph;
  Bx=iblock->pass.Bx;  By=iblock->pass.By;
  Ly=iblock->pass.Ly;  Lz=iblock->pass.Lz;
  Fx=iblock->pass.Fx;  Fz=iblock->pass.Fz;
  Cx=iblock->length.x; Cy=iblock->length.y; Cz=iblock->length.z;

  switch (iblock->apmlid) {
  case 0 :
  case 2 :
  case 6 :
  case 8 :
  case 18 :
  case 20 :
  case 24 :
  case 26 :
    // Corner PML
    csgmc=iblock->apmlc.csgm; csgpc=iblock->apmlc.csgp;
    

//     printf(" bkpt : %d csgmc:%e csgpc:%e \n",iblock->apmlid,csgmc,csgpc);
//     getchar();

    // Ex
    iz=0;
    iy=0;
    for(ix=0;ix<Cx;ix++){
      clf=Hz[iz][iy][ix]-Fz[iz][ix]+By[iy][ix]-Hy[iz][iy][ix];
      tmp1=(csgmc/csgpc)*Pxp[iz][iy][ix]+(ch/csgpc)*clf;
      tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
      Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	(cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
      Pxp[iz][iy][ix]=tmp1;
      Px[iz][iy][ix]=tmp2;
    }
    for(iy=1;iy<Cy;iy++){
      iym=iy-1;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+By[iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmc/csgpc)*Pxp[iz][iy][ix]+(ch/csgpc)*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      iy=0;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Fz[iz][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmc/csgpc)*Pxp[iz][iy][ix]+(ch/csgpc)*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	for(ix=0;ix<Cx;ix++){
	  clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	  tmp1=(csgmc/csgpc)*Pxp[iz][iy][ix]+(ch/csgpc)*clf;
	  tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	  Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	    (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	  Pxp[iz][iy][ix]=tmp1;
	  Px[iz][iy][ix]=tmp2;
	}
      }
    }
    
    // Ey
    iz=0;
    for(iy=0;iy<Cy;iy++){
      ix=0;
      clf=Hx[iz][iy][ix]-Bx[iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
      tmp1=(csgmc/csgpc)*Pyp[iz][iy][ix]+(ch/csgpc)*clf;
      tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
      Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	(cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
      Pyp[iz][iy][ix]=tmp1;
      Py[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hx[iz][iy][ix]-Bx[iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	tmp1=(csgmc/csgpc)*Pyp[iz][iy][ix]+(ch/csgpc)*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      for(iy=0;iy<Cy;iy++){
	ix=0;
	clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
	tmp1=(csgmc/csgpc)*Pyp[iz][iy][ix]+(ch/csgpc)*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	  tmp1=(csgmc/csgpc)*Pyp[iz][iy][ix]+(ch/csgpc)*clf;
	  tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	  Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	    (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	  Pyp[iz][iy][ix]=tmp1;
	  Py[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ez
    for(iz=0;iz<Cz;iz++){
      iy=0;
      ix=0;
      clf=Hy[iz][iy][ix]-Ly[iz][iy]+Fx[iz][ix]-Hx[iz][iy][ix];
      tmp1=(csgmc/csgpc)*Pzp[iz][iy][ix]+(ch/csgpc)*clf;
      tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
      Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	(czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
      Pzp[iz][iy][ix]=tmp1;
      Pz[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Fx[iz][ix]-Hx[iz][iy][ix];
	tmp1=(csgmc/csgpc)*Pzp[iz][iy][ix]+(ch/csgpc)*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	ix=0;
	clf=Hy[iz][iy][ix]-Ly[iz][iy]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	tmp1=(csgmc/csgpc)*Pzp[iz][iy][ix]+(ch/csgpc)*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	  tmp1=(csgmc/csgpc)*Pzp[iz][iy][ix]+(ch/csgpc)*clf;
	  tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	  Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	    (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	  Pzp[iz][iy][ix]=tmp1;
	  Pz[iz][iy][ix]=tmp2;
	}
      }
    }

    break;

  case 1 :
  case 7 :
  case 19 :
  case 25 :
    // Edge PML - along x
    csgme=iblock->apmle.csgm; csgpe=iblock->apmle.csgp;
    csgmte=iblock->apmle.csgmt; csgpte=iblock->apmle.csgpt;

//     printf(" bkpt : %d \n",iblock->apmlid);
//     for(ix=0;ix<Cx;ix++){
//       printf(" bkpt : %d csgme:%e csgpe:%e \n",ix,csgme[ix],csgpe[ix]);
//     }
//     getchar();

    // Ex
    iz=0;
    iy=0;
    for(ix=0;ix<Cx;ix++){
      clf=Hz[iz][iy][ix]-Fz[iz][ix]+By[iy][ix]-Hy[iz][iy][ix];
      tmp1=(csgme[ix]/csgpe[ix])*Pxp[iz][iy][ix]+(ch/csgpe[ix])*clf;
      tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
      Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	(cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
      Pxp[iz][iy][ix]=tmp1;
      Px[iz][iy][ix]=tmp2;
    }
    for(iy=1;iy<Cy;iy++){
      iym=iy-1;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+By[iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgme[ix]/csgpe[ix])*Pxp[iz][iy][ix]+(ch/csgpe[ix])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      iy=0;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Fz[iz][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgme[ix]/csgpe[ix])*Pxp[iz][iy][ix]+(ch/csgpe[ix])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	for(ix=0;ix<Cx;ix++){
	  clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	  tmp1=(csgme[ix]/csgpe[ix])*Pxp[iz][iy][ix]+(ch/csgpe[ix])*clf;
	  tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	  Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	    (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	  Pxp[iz][iy][ix]=tmp1;
	  Px[iz][iy][ix]=tmp2;
	}
      }
    }
    
    // Ey
    iz=0;
    for(iy=0;iy<Cy;iy++){
      ix=0;
      clf=Hx[iz][iy][ix]-Bx[iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
      tmp1=(csgmte[ix]/csgpte[ix])*Pyp[iz][iy][ix]+(ch/csgpte[ix])*clf;
      tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
      Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	(cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
      Pyp[iz][iy][ix]=tmp1;
      Py[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hx[iz][iy][ix]-Bx[iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	tmp1=(csgmte[ix]/csgpte[ix])*Pyp[iz][iy][ix]+(ch/csgpte[ix])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      for(iy=0;iy<Cy;iy++){
	ix=0;
	clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
	tmp1=(csgmte[ix]/csgpte[ix])*Pyp[iz][iy][ix]+(ch/csgpte[ix])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	  tmp1=(csgmte[ix]/csgpte[ix])*Pyp[iz][iy][ix]+(ch/csgpte[ix])*clf;
	  tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	  Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	    (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	  Pyp[iz][iy][ix]=tmp1;
	  Py[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ez
    for(iz=0;iz<Cz;iz++){
      iy=0;
      ix=0;
      clf=Hy[iz][iy][ix]-Ly[iz][iy]+Fx[iz][ix]-Hx[iz][iy][ix];
      tmp1=(csgmte[ix]/csgpte[ix])*Pzp[iz][iy][ix]+(ch/csgpte[ix])*clf;
      tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
      Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	(czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
      Pzp[iz][iy][ix]=tmp1;
      Pz[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Fx[iz][ix]-Hx[iz][iy][ix];
	tmp1=(csgmte[ix]/csgpte[ix])*Pzp[iz][iy][ix]+(ch/csgpte[ix])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	ix=0;
	clf=Hy[iz][iy][ix]-Ly[iz][iy]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	tmp1=(csgmte[ix]/csgpte[ix])*Pzp[iz][iy][ix]+(ch/csgpte[ix])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	  tmp1=(csgmte[ix]/csgpte[ix])*Pzp[iz][iy][ix]+(ch/csgpte[ix])*clf;
	  tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	  Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	    (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	  Pzp[iz][iy][ix]=tmp1;
	  Pz[iz][iy][ix]=tmp2;
	}
      }
    }

    break;
 
  case 3 :
  case 5 :
  case 21 :
  case 23 :
    // edge PML - along y
    csgme=iblock->apmle.csgm; csgpe=iblock->apmle.csgp;
    csgmte=iblock->apmle.csgmt; csgpte=iblock->apmle.csgpt;

//     printf(" bkpt : %d \n",iblock->apmlid);
//     for(iy=0;iy<Cy;iy++){
//       printf(" bkpt : %d csgme:%e csgpe:%e \n",iy,csgme[iy],csgpe[iy]);
//     }
//     getchar();
   
    // Ex
    iz=0;
    iy=0;
    for(ix=0;ix<Cx;ix++){
      clf=Hz[iz][iy][ix]-Fz[iz][ix]+By[iy][ix]-Hy[iz][iy][ix];
      tmp1=(csgmte[iy]/csgpte[iy])*Pxp[iz][iy][ix]+(ch/csgpte[iy])*clf;
      tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
      Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	(cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
      Pxp[iz][iy][ix]=tmp1;
      Px[iz][iy][ix]=tmp2;
    }
    for(iy=1;iy<Cy;iy++){
      iym=iy-1;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+By[iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmte[iy]/csgpte[iy])*Pxp[iz][iy][ix]+(ch/csgpte[iy])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      iy=0;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Fz[iz][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmte[iy]/csgpte[iy])*Pxp[iz][iy][ix]+(ch/csgpte[iy])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	for(ix=0;ix<Cx;ix++){
	  clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	  tmp1=(csgmte[iy]/csgpte[iy])*Pxp[iz][iy][ix]+(ch/csgpte[iy])*clf;
	  tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	  Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	    (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	  Pxp[iz][iy][ix]=tmp1;
	  Px[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ey
    iz=0;
    for(iy=0;iy<Cy;iy++){
      ix=0;
      clf=Hx[iz][iy][ix]-Bx[iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
      tmp1=(csgme[iy]/csgpe[iy])*Pyp[iz][iy][ix]+(ch/csgpe[iy])*clf;
      tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
      Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	(cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
      Pyp[iz][iy][ix]=tmp1;
      Py[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hx[iz][iy][ix]-Bx[iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	tmp1=(csgme[iy]/csgpe[iy])*Pyp[iz][iy][ix]+(ch/csgpe[iy])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      for(iy=0;iy<Cy;iy++){
	ix=0;
	clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
	tmp1=(csgme[iy]/csgpe[iy])*Pyp[iz][iy][ix]+(ch/csgpe[iy])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	  tmp1=(csgme[iy]/csgpe[iy])*Pyp[iz][iy][ix]+(ch/csgpe[iy])*clf;
	  tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	  Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	    (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	  Pyp[iz][iy][ix]=tmp1;
	  Py[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ez
    for(iz=0;iz<Cz;iz++){
      iy=0;
      ix=0;
      clf=Hy[iz][iy][ix]-Ly[iz][iy]+Fx[iz][ix]-Hx[iz][iy][ix];
      tmp1=(csgmte[iy]/csgpte[iy])*Pzp[iz][iy][ix]+(ch/csgpte[iy])*clf;
      tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
      Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	(czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
      Pzp[iz][iy][ix]=tmp1;
      Pz[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Fx[iz][ix]-Hx[iz][iy][ix];
	tmp1=(csgmte[iy]/csgpte[iy])*Pzp[iz][iy][ix]+(ch/csgpte[iy])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	ix=0;
	clf=Hy[iz][iy][ix]-Ly[iz][iy]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	tmp1=(csgmte[iy]/csgpte[iy])*Pzp[iz][iy][ix]+(ch/csgpte[iy])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	  tmp1=(csgmte[iy]/csgpte[iy])*Pzp[iz][iy][ix]+(ch/csgpte[iy])*clf;
	  tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	  Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	    (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	  Pzp[iz][iy][ix]=tmp1;
	  Pz[iz][iy][ix]=tmp2;
	}
      }
    }

    break;
    
  case 9 :
  case 11 :
  case 15 :
  case 17 :
    // edge PML - along z
    csgme=iblock->apmle.csgm; csgpe=iblock->apmle.csgp;
    csgmte=iblock->apmle.csgmt; csgpte=iblock->apmle.csgpt;

//     printf(" bkpt : %d \n",iblock->apmlid);
//     for(iz=0;iz<Cz;iz++){
//       printf(" bkpt : %d csgme:%e csgpe:%e \n",iz,csgme[iz],csgpe[iz]);
//     }
//     getchar();
   
    // Ex
    iz=0;
    iy=0;
    for(ix=0;ix<Cx;ix++){
      clf=Hz[iz][iy][ix]-Fz[iz][ix]+By[iy][ix]-Hy[iz][iy][ix];
      tmp1=(csgmte[iz]/csgpte[iz])*Pxp[iz][iy][ix]+(ch/csgpte[iz])*clf;
      tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
      Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	(cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
      Pxp[iz][iy][ix]=tmp1;
      Px[iz][iy][ix]=tmp2;
    }
    for(iy=1;iy<Cy;iy++){
      iym=iy-1;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+By[iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmte[iz]/csgpte[iz])*Pxp[iz][iy][ix]+(ch/csgpte[iz])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      iy=0;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Fz[iz][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmte[iz]/csgpte[iz])*Pxp[iz][iy][ix]+(ch/csgpte[iz])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	for(ix=0;ix<Cx;ix++){
	  clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	  tmp1=(csgmte[iz]/csgpte[iz])*Pxp[iz][iy][ix]+(ch/csgpte[iz])*clf;
	  tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	  Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	    (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	  Pxp[iz][iy][ix]=tmp1;
	  Px[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ey
    iz=0;
    for(iy=0;iy<Cy;iy++){
      ix=0;
      clf=Hx[iz][iy][ix]-Bx[iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
      tmp1=(csgmte[iz]/csgpte[iz])*Pyp[iz][iy][ix]+(ch/csgpte[iz])*clf;
      tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
      Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	(cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
      Pyp[iz][iy][ix]=tmp1;
      Py[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hx[iz][iy][ix]-Bx[iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	tmp1=(csgmte[iz]/csgpte[iz])*Pyp[iz][iy][ix]+(ch/csgpte[iz])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      for(iy=0;iy<Cy;iy++){
	ix=0;
	clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
	tmp1=(csgmte[iz]/csgpte[iz])*Pyp[iz][iy][ix]+(ch/csgpte[iz])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	  tmp1=(csgmte[iz]/csgpte[iz])*Pyp[iz][iy][ix]+(ch/csgpte[iz])*clf;
	  tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	  Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	    (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	  Pyp[iz][iy][ix]=tmp1;
	  Py[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ez
    for(iz=0;iz<Cz;iz++){
      iy=0;
      ix=0;
      clf=Hy[iz][iy][ix]-Ly[iz][iy]+Fx[iz][ix]-Hx[iz][iy][ix];
      tmp1=(csgme[iz]/csgpe[iz])*Pzp[iz][iy][ix]+(ch/csgpe[iz])*clf;
      tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
      Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	(czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
      Pzp[iz][iy][ix]=tmp1;
      Pz[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Fx[iz][ix]-Hx[iz][iy][ix];
	tmp1=(csgme[iz]/csgpe[iz])*Pzp[iz][iy][ix]+(ch/csgpe[iz])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	ix=0;
	clf=Hy[iz][iy][ix]-Ly[iz][iy]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	tmp1=(csgme[iz]/csgpe[iz])*Pzp[iz][iy][ix]+(ch/csgpe[iz])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	  tmp1=(csgme[iz]/csgpe[iz])*Pzp[iz][iy][ix]+(ch/csgpe[iz])*clf;
	  tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	  Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	    (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	  Pzp[iz][iy][ix]=tmp1;
	  Pz[iz][iy][ix]=tmp2;
	}
      }
    }

    break;
 
  case 4 :
  case 22 :
    // Face PML - z
    csgmf=iblock->apmlf.csgm; csgpf=iblock->apmlf.csgp;
    csgmt1f=iblock->apmlf.csgmt1; csgpt1f=iblock->apmlf.csgpt1;
    csgmt2f=iblock->apmlf.csgmt2; csgpt2f=iblock->apmlf.csgpt2;

//     printf(" bkpt : %d \n",iblock->apmlid);
//     for(iy=0;iy<Cy;iy++){
//       for(ix=0;ix<Cx;ix++){
// 	printf(" bkpt : %d %d csgmf:%e csgpf:%e \n",iy,ix,csgmf[iy][ix],csgpf[iy][ix]);
//       }
//       getchar();
//     }
//     getchar();
   
    iz=0;
    iy=0;
    for(ix=0;ix<Cx;ix++){
      clf=Hz[iz][iy][ix]-Fz[iz][ix]+By[iy][ix]-Hy[iz][iy][ix];
      tmp1=(csgmt1f[iy][ix]/csgpt1f[iy][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iy][ix])*clf;
      tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
      Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	(cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
      Pxp[iz][iy][ix]=tmp1;
      Px[iz][iy][ix]=tmp2;
    }
    for(iy=1;iy<Cy;iy++){
      iym=iy-1;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+By[iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmt1f[iy][ix]/csgpt1f[iy][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iy][ix])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      iy=0;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Fz[iz][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmt1f[iy][ix]/csgpt1f[iy][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iy][ix])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	for(ix=0;ix<Cx;ix++){
	  clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	  tmp1=(csgmt1f[iy][ix]/csgpt1f[iy][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iy][ix])*clf;
	  tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	  Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	    (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	  Pxp[iz][iy][ix]=tmp1;
	  Px[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ey
    iz=0;
    for(iy=0;iy<Cy;iy++){
      ix=0;
      clf=Hx[iz][iy][ix]-Bx[iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
      tmp1=(csgmt2f[iy][ix]/csgpt2f[iy][ix])*Pyp[iz][iy][ix]+(ch/csgpt2f[iy][ix])*clf;
      tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
      Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	(cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
      Pyp[iz][iy][ix]=tmp1;
      Py[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hx[iz][iy][ix]-Bx[iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	tmp1=(csgmt2f[iy][ix]/csgpt2f[iy][ix])*Pyp[iz][iy][ix]+(ch/csgpt2f[iy][ix])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      for(iy=0;iy<Cy;iy++){
	ix=0;
	clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
	tmp1=(csgmt2f[iy][ix]/csgpt2f[iy][ix])*Pyp[iz][iy][ix]+(ch/csgpt2f[iy][ix])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	  tmp1=(csgmt2f[iy][ix]/csgpt2f[iy][ix])*Pyp[iz][iy][ix]+(ch/csgpt2f[iy][ix])*clf;
	  tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	  Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	    (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	  Pyp[iz][iy][ix]=tmp1;
	  Py[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ez
    for(iz=0;iz<Cz;iz++){
      iy=0;
      ix=0;
      clf=Hy[iz][iy][ix]-Ly[iz][iy]+Fx[iz][ix]-Hx[iz][iy][ix];
      tmp1=(csgmf[iy][ix]/csgpf[iy][ix])*Pzp[iz][iy][ix]+(ch/csgpf[iy][ix])*clf;
      tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
      Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	(czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
      Pzp[iz][iy][ix]=tmp1;
      Pz[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Fx[iz][ix]-Hx[iz][iy][ix];
	tmp1=(csgmf[iy][ix]/csgpf[iy][ix])*Pzp[iz][iy][ix]+(ch/csgpf[iy][ix])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	ix=0;
	clf=Hy[iz][iy][ix]-Ly[iz][iy]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	tmp1=(csgmf[iy][ix]/csgpf[iy][ix])*Pzp[iz][iy][ix]+(ch/csgpf[iy][ix])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	  tmp1=(csgmf[iy][ix]/csgpf[iy][ix])*Pzp[iz][iy][ix]+(ch/csgpf[iy][ix])*clf;
	  tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	  Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	    (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	  Pzp[iz][iy][ix]=tmp1;
	  Pz[iz][iy][ix]=tmp2;
	}
      }
    }

    break;
 
  case 10 :
  case 16 :
    // Face PML - y
    csgmf=iblock->apmlf.csgm;     csgpf=iblock->apmlf.csgp;
    csgmt1f=iblock->apmlf.csgmt1; csgpt1f=iblock->apmlf.csgpt1;
    csgmt2f=iblock->apmlf.csgmt2; csgpt2f=iblock->apmlf.csgpt2;

//     printf(" bkpt : %d \n",iblock->apmlid);
//     for(iz=0;iz<Cz;iz++){
//       for(ix=0;ix<Cx;ix++){
// 	printf(" bkpt : %d %d csgmf:%e csgpf:%e \n",iz,ix,csgmf[iz][ix],csgpf[iz][ix]);
//       }
//       getchar();
//     }
//     getchar();

   
    // Ex
    iz=0;
    iy=0;
    for(ix=0;ix<Cx;ix++){
      clf=Hz[iz][iy][ix]-Fz[iz][ix]+By[iy][ix]-Hy[iz][iy][ix];
      tmp1=(csgmt1f[iz][ix]/csgpt1f[iz][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iz][ix])*clf;
      tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
      Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	(cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
      Pxp[iz][iy][ix]=tmp1;
      Px[iz][iy][ix]=tmp2;
    }
    for(iy=1;iy<Cy;iy++){
      iym=iy-1;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+By[iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmt1f[iz][ix]/csgpt1f[iz][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iz][ix])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      iy=0;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Fz[iz][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmt1f[iz][ix]/csgpt1f[iz][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iz][ix])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	for(ix=0;ix<Cx;ix++){
	  clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	  tmp1=(csgmt1f[iz][ix]/csgpt1f[iz][ix])*Pxp[iz][iy][ix]+(ch/csgpt1f[iz][ix])*clf;
	  tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	  Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	    (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	  Pxp[iz][iy][ix]=tmp1;
	  Px[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ey
    iz=0;
    for(iy=0;iy<Cy;iy++){
      ix=0;
      clf=Hx[iz][iy][ix]-Bx[iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
      tmp1=(csgmf[iz][ix]/csgpf[iz][ix])*Pyp[iz][iy][ix]+(ch/csgpf[iz][ix])*clf;
      tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
      Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	(cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
      Pyp[iz][iy][ix]=tmp1;
      Py[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hx[iz][iy][ix]-Bx[iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	tmp1=(csgmf[iz][ix]/csgpf[iz][ix])*Pyp[iz][iy][ix]+(ch/csgpf[iz][ix])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      for(iy=0;iy<Cy;iy++){
	ix=0;
	clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
	tmp1=(csgmf[iz][ix]/csgpf[iz][ix])*Pyp[iz][iy][ix]+(ch/csgpf[iz][ix])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	  tmp1=(csgmf[iz][ix]/csgpf[iz][ix])*Pyp[iz][iy][ix]+(ch/csgpf[iz][ix])*clf;
	  tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	  Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	    (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	  Pyp[iz][iy][ix]=tmp1;
	  Py[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ez
    for(iz=0;iz<Cz;iz++){
      iy=0;
      ix=0;
      clf=Hy[iz][iy][ix]-Ly[iz][iy]+Fx[iz][ix]-Hx[iz][iy][ix];
      tmp1=(csgmt2f[iz][ix]/csgpt2f[iz][ix])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][ix])*clf;
      tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
      Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	(czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
      Pzp[iz][iy][ix]=tmp1;
      Pz[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Fx[iz][ix]-Hx[iz][iy][ix];
	tmp1=(csgmt2f[iz][ix]/csgpt2f[iz][ix])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][ix])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	ix=0;
	clf=Hy[iz][iy][ix]-Ly[iz][iy]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	tmp1=(csgmt2f[iz][ix]/csgpt2f[iz][ix])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][ix])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	  tmp1=(csgmt2f[iz][ix]/csgpt2f[iz][ix])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][ix])*clf;
	  tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	  Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	    (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	  Pzp[iz][iy][ix]=tmp1;
	  Pz[iz][iy][ix]=tmp2;
	}
      }
    }

    break;

  case 12 :
  case 14 :
    // Face PML - x
    csgmf=iblock->apmlf.csgm;     csgpf=iblock->apmlf.csgp;
    csgmt1f=iblock->apmlf.csgmt1; csgpt1f=iblock->apmlf.csgpt1;
    csgmt2f=iblock->apmlf.csgmt2; csgpt2f=iblock->apmlf.csgpt2;


//     printf(" bkpt : %d \n",iblock->apmlid);
//     for(iz=0;iz<Cz;iz++){
//       for(iy=0;iy<Cy;iy++){
// 	printf(" bkpt : %d %d csgmf:%e csgpf:%e \n",iz,iy,csgmf[iz][iy],csgpf[iz][iy]);
//       }
//       getchar();
//     }
//     getchar();

   
    // Ex
    iz=0;
    iy=0;
    for(ix=0;ix<Cx;ix++){
      clf=Hz[iz][iy][ix]-Fz[iz][ix]+By[iy][ix]-Hy[iz][iy][ix];
      tmp1=(csgmf[iz][iy]/csgpf[iz][iy])*Pxp[iz][iy][ix]+(ch/csgpf[iz][iy])*clf;
      tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
      Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	(cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
      Pxp[iz][iy][ix]=tmp1;
      Px[iz][iy][ix]=tmp2;
    }
    for(iy=1;iy<Cy;iy++){
      iym=iy-1;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+By[iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmf[iz][iy]/csgpf[iz][iy])*Pxp[iz][iy][ix]+(ch/csgpf[iz][iy])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      iy=0;
      for(ix=0;ix<Cx;ix++){
	clf=Hz[iz][iy][ix]-Fz[iz][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	tmp1=(csgmf[iz][iy]/csgpf[iz][iy])*Pxp[iz][iy][ix]+(ch/csgpf[iz][iy])*clf;
	tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	  (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	Pxp[iz][iy][ix]=tmp1;
	Px[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	for(ix=0;ix<Cx;ix++){
	  clf=Hz[iz][iy][ix]-Hz[iz][iym][ix]+Hy[izm][iy][ix]-Hy[iz][iy][ix];
	  tmp1=(csgmf[iz][iy]/csgpf[iz][iy])*Pxp[iz][iy][ix]+(ch/csgpf[iz][iy])*clf;
	  tmp2=(cym[iy]/cyp[iy])*Px[iz][iy][ix]+(2/cyp[iy])*(tmp1-Pxp[iz][iy][ix]);
	  Ex[iz][iy][ix]=(czm[iz]/czp[iz])*Ex[iz][iy][ix]+
	    (cxph[ix]/czp[iz])*tmp2-(cxmh[ix]/czp[iz])*Px[iz][iy][ix];
	  Pxp[iz][iy][ix]=tmp1;
	  Px[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ey
    iz=0;
    for(iy=0;iy<Cy;iy++){
      ix=0;
      clf=Hx[iz][iy][ix]-Bx[iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
      tmp1=(csgmt1f[iz][iy]/csgpt1f[iz][iy])*Pyp[iz][iy][ix]+(ch/csgpt1f[iz][iy])*clf;
      tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
      Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	(cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
      Pyp[iz][iy][ix]=tmp1;
      Py[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hx[iz][iy][ix]-Bx[iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	tmp1=(csgmt1f[iz][iy]/csgpt1f[iz][iy])*Pyp[iz][iy][ix]+(ch/csgpt1f[iz][iy])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
      }
    }
    for(iz=1;iz<Cz;iz++){
      izm=iz-1;
      for(iy=0;iy<Cy;iy++){
	ix=0;
	clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Lz[iz][iy]-Hz[iz][iy][ix];
	tmp1=(csgmt1f[iz][iy]/csgpt1f[iz][iy])*Pyp[iz][iy][ix]+(ch/csgpt1f[iz][iy])*clf;
	tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	  (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	Pyp[iz][iy][ix]=tmp1;
	Py[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hx[iz][iy][ix]-Hx[izm][iy][ix]+Hz[iz][iy][ixm]-Hz[iz][iy][ix];
	  tmp1=(csgmt1f[iz][iy]/csgpt1f[iz][iy])*Pyp[iz][iy][ix]+(ch/csgpt1f[iz][iy])*clf;
	  tmp2=(czm[iz]/czp[iz])*Py[iz][iy][ix]+(2/czp[iz])*(tmp1-Pyp[iz][iy][ix]);
	  Ey[iz][iy][ix]=(cxm[ix]/cxp[ix])*Ey[iz][iy][ix]+
	    (cyph[iy]/cxp[ix])*tmp2-(cymh[iy]/cxp[ix])*Py[iz][iy][ix];
	  Pyp[iz][iy][ix]=tmp1;
	  Py[iz][iy][ix]=tmp2;
	}
      }
    }

    // Ez
    for(iz=0;iz<Cz;iz++){
      iy=0;
      ix=0;
      clf=Hy[iz][iy][ix]-Ly[iz][iy]+Fx[iz][ix]-Hx[iz][iy][ix];
      tmp1=(csgmt2f[iz][iy]/csgpt2f[iz][iy])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][iy])*clf;
      tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
      Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	(czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
      Pzp[iz][iy][ix]=tmp1;
      Pz[iz][iy][ix]=tmp2;
      for(ix=1;ix<Cx;ix++){
	ixm=ix-1;
	clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Fx[iz][ix]-Hx[iz][iy][ix];
	tmp1=(csgmt2f[iz][iy]/csgpt2f[iz][iy])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][iy])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
      }
      for(iy=1;iy<Cy;iy++){
	iym=iy-1;
	ix=0;
	clf=Hy[iz][iy][ix]-Ly[iz][iy]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	tmp1=(csgmt2f[iz][iy]/csgpt2f[iz][iy])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][iy])*clf;
	tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	  (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	Pzp[iz][iy][ix]=tmp1;
	Pz[iz][iy][ix]=tmp2;
	for(ix=1;ix<Cx;ix++){
	  ixm=ix-1;
	  clf=Hy[iz][iy][ix]-Hy[iz][iy][ixm]+Hx[iz][iym][ix]-Hx[iz][iy][ix];
	  tmp1=(csgmt2f[iz][iy]/csgpt2f[iz][iy])*Pzp[iz][iy][ix]+(ch/csgpt2f[iz][iy])*clf;
	  tmp2=(cxm[ix]/cxp[ix])*Pz[iz][iy][ix]+(2/cxp[ix])*(tmp1-Pzp[iz][iy][ix]);
	  Ez[iz][iy][ix]=(cym[iy]/cyp[iy])*Ez[iz][iy][ix]+
	    (czph[iz]/cyp[iy])*tmp2-(czmh[iz]/cyp[iy])*Pz[iz][iy][ix];
	  Pzp[iz][iy][ix]=tmp1;
	  Pz[iz][iy][ix]=tmp2;
	}
      }
    }

    break;

  default :
    printf("Error in apmlid %d %d \n",iblock->id,iblock->apmlid);

  }
}



/*****************************************************************************/

void calc_apml_H(I_Block *iblock)
{
  int ix, iy, iz, ixp, iyp, izp;
  REALNO ***Ex, ***Ey, ***Ez, ***Hx, ***Hy, ***Hz;
  REALNO ***Bxa, ***Bya, ***Bza;
  REALNO **Tx, **Ty;
  REALNO **Ry, **Rz;
  REALNO **Ax, **Az;
  REALNO *cxm, *cxp, *cym, *cyp, *czm, *czp;
  REALNO *cxmh, *cxph, *cymh, *cyph, *czmh, *czph;
  int Cx, Cy, Cz, Cxm, Cym, Czm;
  REALNO ch,mui;
  REALNO tmp[][], clf[][];
  int i;

  //Short cuts
  Ex=iblock->max.Ex; Ey=iblock->max.Ey; Ez=iblock->max.Ez;
  Hx=iblock->max.Hx; Hy=iblock->max.Hy; Hz=iblock->max.Hz;
  Bxa=iblock->apml.Bxa; Bya=iblock->apml.Bya; Bza=iblock->apml.Bza;
  cxm=iblock->apml.cxm; cxp=iblock->apml.cxp;
  cym=iblock->apml.cym; cyp=iblock->apml.cyp;
  czm=iblock->apml.czm; czp=iblock->apml.czp;
  cxmh=iblock->apml.cxmh; cxph=iblock->apml.cxph;
  cymh=iblock->apml.cymh; cyph=iblock->apml.cyph;
  czmh=iblock->apml.czmh; czph=iblock->apml.czph;
  Tx=iblock->pass.Tx; Ty=iblock->pass.Ty;
  Ry=iblock->pass.Ry; Rz=iblock->pass.Rz;
  Ax=iblock->pass.Ax; Az=iblock->pass.Az;

  Cx=iblock->length.x; Cy=iblock->length.y; Cz=iblock->length.z;
  Czm=Cz-1; Cym=Cy-1; Cxm=Cx-1;

  tmp = (REALNO **) malloc(Cz*sizeof(REALNO *));
  clf = (REALNO **) malloc(Cz*sizeof(REALNO *));

  for (i=0; i<Cz; i++)  {
      tmp[i] = (REALNO *) malloc(Cym*sizeof(REALNO));
      clf[i] = (REALNO *) malloc(Cym*sizeof(REALNO));
  }

  ch=2*CNST_ep*iblock->gblock->dt/iblock->gblock->dx;
  mui=1/CNST_mu;

  // Hx
  for(iz=0;iz<Czm;iz++){
    izp=iz+1;
    for(iy=0;iy<Cym;iy++){
      iyp=iy+1;
      for(ix=0;ix<Cx;ix++){
	clf=Ez[iz][iyp][ix]-Ez[iz][iy][ix]+Ey[iz][iy][ix]-Ey[izp][iy][ix];
	tmp=(czmh[iz]/czph[iz])*Bxa[iz][iy][ix]-(ch/czph[iz])*clf;
	Hx[iz][iy][ix]=(cymh[iy]/cyph[iy])*Hx[iz][iy][ix]
	  +(mui*cxp[ix]/cyph[iy])*tmp
	  -(mui*cxm[ix]/cyph[iy])*Bxa[iz][iy][ix];
	Bxa[iz][iy][ix]=tmp;
      }
    }
    iy=Cym;
    for(ix=0;ix<Cx;ix++){
      clf=Az[iz][ix]-Ez[iz][iy][ix]+Ey[iz][iy][ix]-Ey[izp][iy][ix];
      tmp=(czmh[iz]/czph[iz])*Bxa[iz][iy][ix]-(ch/czph[iz])*clf;
      Hx[iz][iy][ix]=(cymh[iy]/cyph[iy])*Hx[iz][iy][ix]
	+(mui*cxp[ix]/cyph[iy])*tmp
	-(mui*cxm[ix]/cyph[iy])*Bxa[iz][iy][ix];
      Bxa[iz][iy][ix]=tmp;
    }
  }
  iz=Czm;
  for(iy=0;iy<Cym;iy++){
    iyp=iy+1;
    for(ix=0;ix<Cx;ix++){
      clf=Ez[iz][iyp][ix]-Ez[iz][iy][ix]+Ey[iz][iy][ix]-Ty[iy][ix];
      tmp=(czmh[iz]/czph[iz])*Bxa[iz][iy][ix]-(ch/czph[iz])*clf;
      Hx[iz][iy][ix]=(cymh[iy]/cyph[iy])*Hx[iz][iy][ix]
	+(mui*cxp[ix]/cyph[iy])*tmp
	-(mui*cxm[ix]/cyph[iy])*Bxa[iz][iy][ix];
      Bxa[iz][iy][ix]=tmp;
    }
  }
  iy=Cym;
  for(ix=0;ix<Cx;ix++){
    clf=Az[iz][ix]-Ez[iz][iy][ix]+Ey[iz][iy][ix]-Ty[iy][ix];
    tmp=(czmh[iz]/czph[iz])*Bxa[iz][iy][ix]-(ch/czph[iz])*clf;
    Hx[iz][iy][ix]=(cymh[iy]/cyph[iy])*Hx[iz][iy][ix]
      +(mui*cxp[ix]/cyph[iy])*tmp
      -(mui*cxm[ix]/cyph[iy])*Bxa[iz][iy][ix];
    Bxa[iz][iy][ix]=tmp;
  }

  // Hy
  for(iz=0;iz<Czm;iz++){
    izp=iz+1;
    for(iy=0;iy<Cy;iy++){
      for(ix=0;ix<Cxm;ix++){
	ixp=ix+1;
	clf=Ez[iz][iy][ix]-Ez[iz][iy][ixp]+Ex[izp][iy][ix]-Ex[iz][iy][ix];
	tmp=(cxmh[ix]/cxph[ix])*Bya[iz][iy][ix]-(ch/cxph[ix])*clf;
	Hy[iz][iy][ix]=(czmh[iz]/czph[iz])*Hy[iz][iy][ix]
	  +(mui*cyp[iy]/czph[iz])*tmp
	  -(mui*cym[iy]/czph[iz])*Bya[iz][iy][ix];
	Bya[iz][iy][ix]=tmp;
      }
      ix=Cxm;
      clf=Ez[iz][iy][ix]-Rz[iz][iy]+Ex[izp][iy][ix]-Ex[iz][iy][ix];
      tmp=(cxmh[ix]/cxph[ix])*Bya[iz][iy][ix]-(ch/cxph[ix])*clf;
      Hy[iz][iy][ix]=(czmh[iz]/czph[iz])*Hy[iz][iy][ix]
	+(mui*cyp[iy]/czph[iz])*tmp
	-(mui*cym[iy]/czph[iz])*Bya[iz][iy][ix];
      Bya[iz][iy][ix]=tmp;
    }
  }
  iz=Czm;
  for(iy=0;iy<Cy;iy++){
    for(ix=0;ix<Cxm;ix++){
      ixp=ix+1;
      clf=Ez[iz][iy][ix]-Ez[iz][iy][ixp]+Tx[iy][ix]-Ex[iz][iy][ix];
      tmp=(cxmh[ix]/cxph[ix])*Bya[iz][iy][ix]-(ch/cxph[ix])*clf;
      Hy[iz][iy][ix]=(czmh[iz]/czph[iz])*Hy[iz][iy][ix]
	+(mui*cyp[iy]/czph[iz])*tmp
	-(mui*cym[iy]/czph[iz])*Bya[iz][iy][ix];
      Bya[iz][iy][ix]=tmp;
    }
    ix=Cxm;
    clf=Ez[iz][iy][ix]-Rz[iz][iy]+Tx[iy][ix]-Ex[iz][iy][ix];
    tmp=(cxmh[ix]/cxph[ix])*Bya[iz][iy][ix]-(ch/cxph[ix])*clf;
    Hy[iz][iy][ix]=(czmh[iz]/czph[iz])*Hy[iz][iy][ix]
      +(mui*cyp[iy]/czph[iz])*tmp
      -(mui*cym[iy]/czph[iz])*Bya[iz][iy][ix];
    Bya[iz][iy][ix]=tmp;
  }

  /* pluto start (Cz,Cym,Cxm) */
  for(iz=0;iz<Cz;iz++){
      for(iy=0;iy<Cym;iy++){
          for(ix=0;ix<Cxm;ix++){
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
      for(ix=0;ix<Cxm;ix++){
          clf[iz][Cym]=Ex[iz][Cym][ix]-Ax[iz][ix]+Ey[iz][Cym][ix+1]-Ey[iz][Cym][ix];
          tmp[iz][Cym]=(cymh[Cym]/cyph[Cym])*Bza[iz][Cym][ix]-(ch/cyph[Cym])*clf[iz][Cym];
          Hz[iz][Cym][ix]=(cxmh[ix]/cxph[ix])*Hz[iz][Cym][ix]
              +(mui*czp[iz]/cxph[ix])*tmp[iz][Cym]
              -(mui*czm[iz]/cxph[ix])*Bza[iz][Cym][ix];
          Bza[iz][Cym][ix]=tmp[iz][Cym];
      }
      clf[iz][Cym]=Ex[iz][Cym][Cxm]-Ax[iz][Cxm]+Ry[iz][Cym]-Ey[iz][Cym][Cxm];
      tmp[iz][Cym]=(cymh[Cym]/cyph[Cym])*Bza[iz][Cym][Cxm]-(ch/cyph[Cym])*clf[iz][Cym];
      Hz[iz][Cym][Cxm]=(cxmh[Cxm]/cxph[Cxm])*Hz[iz][Cym][Cxm]
          +(mui*czp[iz]/cxph[Cxm])*tmp[iz][Cym]
          -(mui*czm[iz]/cxph[Cxm])*Bza[iz][Cym][Cxm];
      Bza[iz][Cym][Cxm]=tmp[iz][Cym];
  }
  /* pluto end */

}

/*****************************************************************************/

void send_aPML_E_fields(I_Block *iblock, int flag)
{
    int ix, iy, iz, idd;
    REALNO **Tx, **Ty, **Ry, **Rz, **Ax, **Az;

    //short cuts
    REALNO ***Ex, ***Ey, ***Ez, ***Hx, ***Hy, ***Hz;
    int Cx, Cy, Cz, id;
    int BBx,BBy,bx,by,bz;
    Ex=iblock->max.Ex; Ey=iblock->max.Ey; Ez=iblock->max.Ez;
    Hx=iblock->max.Hx; Hy=iblock->max.Hy; Hz=iblock->max.Hz;
    Cx=iblock->length.x; Cy=iblock->length.y; Cz=iblock->length.z;
    id=iblock->id;
    BBx=iblock->gblock->BBx; BBy=iblock->gblock->BBy;
    bz=(int)id/(BBx*BBy);
    by=(int)(id-bz*BBx*BBy)/BBx;
    bx=(int)(id-bz*BBx*BBy-by*BBx);

#ifdef MPIRUN
    int BBz=iblock->gblock->BBz;
    MPI_Status pstatus;
    REALNO *buffer;  
    MPI_Request req;
    //static REALNO buffer[BNDRYbuffy];

    if(flag>0){ /* send mode */

        // if((2*Cx*Cy>BNDRYbuffy) || (2*Cx*Cz>BNDRYbuffy) || (2*Cy*Cz>BNDRYbuffy) ){
        //       printf("Error in static memory in routine send_PML_E_fields \n");
        //       printf("Must increase buffer from %d to max of\n",BNDRYbuffy);
        //       printf("%d %d %d\n",2*Cx*Cy,2*Cx*Cz,2*Cy*Cz);
        //       exit(1);
        //     }

        /* Top - pass down */
        if(bz>0){
            idd=id-BBx*BBy;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]==iblock->gblock->rank){ 
                    // Pass to the Same CPU
                    Tx=iblock->gblock->iblock[idd].pass.Tx; 
                    Ty=iblock->gblock->iblock[idd].pass.Ty;
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            Tx[iy][ix]=Ex[0][iy][ix];
                            Ty[iy][ix]=Ey[0][iy][ix];
                        } 
                    } 
                }else{
                    // Make sure last group was sent before packing
                    MPI_Wait(&(iblock->send_id.request[1]),&pstatus);
                    buffer=iblock->send.face[1];
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[iy*Cx+ix]=Ex[0][iy][ix];
                        }
                    }
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[Cx*Cy+iy*Cx+ix]=Ey[0][iy][ix];
                        }
                    }
                    MPI_Isend(buffer,2*Cx*Cy,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*id+0,MPI_COMM_WORLD,&req);  
                    iblock->send_id.request[1]=req;
                }
        }

        /* bAck - pass front */
        if(by>0){
            idd=id-BBx;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]==iblock->gblock->rank){ 
                    // Pass to the Same CPU
                    Ax=iblock->gblock->iblock[idd].pass.Ax; 
                    Az=iblock->gblock->iblock[idd].pass.Az;
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            Ax[iz][ix]=Ex[iz][0][ix];
                            Az[iz][ix]=Ez[iz][0][ix];
                        }
                    }
                }else{
                    // Make sure last group was sent before packing
                    MPI_Wait(&(iblock->send_id.request[3]),&pstatus);
                    buffer=iblock->send.face[3];
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[iz*Cx+ix]=Ex[iz][0][ix];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[Cx*Cz+iz*Cx+ix]=Ez[iz][0][ix];
                        }
                    }
                    MPI_Isend(buffer,2*Cz*Cx,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*id+2,MPI_COMM_WORLD,&req);
                    iblock->send_id.request[3]=req;
                }
        }

        /* Right - pass left */
        if(bx>0){
            idd=id-1;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]==iblock->gblock->rank){ 
                    // Pass to the Same CPU
                    Ry=iblock->gblock->iblock[idd].pass.Ry; 
                    Rz=iblock->gblock->iblock[idd].pass.Rz;
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            Ry[iz][iy]=Ey[iz][iy][0];
                            Rz[iz][iy]=Ez[iz][iy][0];
                        }
                    }
                }else{
                    // Make sure last group was sent before packing
                    MPI_Wait(&(iblock->send_id.request[5]),&pstatus);
                    buffer=iblock->send.face[5];
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            buffer[iz*Cy+iy]=Ey[iz][iy][0];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            buffer[Cy*Cz+iz*Cy+iy]=Ez[iz][iy][0];
                        }
                    }
                    MPI_Isend(buffer,2*Cz*Cy,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*id+4,MPI_COMM_WORLD,&req);
                    iblock->send_id.request[5]=req;
                }
        }
    }else{ //Recieve Only
        /* Top - pass down */
        if(bz<(BBz-1)){ 
            idd=id+BBx*BBy;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]!=iblock->gblock->rank){ 
                    //Some other CPU sent this
                    Tx=iblock->gblock->iblock[id].pass.Tx; 
                    Ty=iblock->gblock->iblock[id].pass.Ty;
                    // Make sure data was received before unpacking
                    MPI_Wait(&(iblock->receive_id.request[0]),&pstatus);
                    buffer=iblock->receive.face[0];
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            Tx[iy][ix]=buffer[iy*Cx+ix];
                        }
                    }
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            Ty[iy][ix]=buffer[Cx*Cy+iy*Cx+ix];
                        }
                    }
                    // Post the next receive
                    MPI_Irecv(buffer,2*Cx*Cy,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*idd+0,MPI_COMM_WORLD,&req);
                    iblock->receive_id.request[0]=req;
                }
        }

        /* bAck - pass front */
        if(by<(BBy-1)){ //recieve data
            idd=id+BBx;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]!=iblock->gblock->rank){ 
                    //Some other CPU sent this
                    Ax=iblock->gblock->iblock[id].pass.Ax; 
                    Az=iblock->gblock->iblock[id].pass.Az;
                    // Make sure data was received before unpacking
                    MPI_Wait(&(iblock->receive_id.request[2]),&pstatus);
                    buffer=iblock->receive.face[2];	
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            Ax[iz][ix]=buffer[iz*Cx+ix];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            Az[iz][ix]=buffer[Cx*Cz+iz*Cx+ix];
                        }
                    }
                    // Post the next receive
                    MPI_Irecv(buffer,2*Cx*Cz,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*idd+2,MPI_COMM_WORLD,&req);
                    iblock->receive_id.request[2]=req;
                }
        }

        /* Right - pass left */
        if(bx<(BBx-1)){ //recieve data
            idd=id+1;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]!=iblock->gblock->rank){ 
                    //Some other CPU sent this
                    Ry=iblock->gblock->iblock[id].pass.Ry; 
                    Rz=iblock->gblock->iblock[id].pass.Rz;
                    // Make sure data was received before unpacking
                    MPI_Wait(&(iblock->receive_id.request[4]),&pstatus);
                    buffer=iblock->receive.face[4];
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            Ry[iz][iy]=buffer[iz*Cy+iy];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            Rz[iz][iy]=buffer[Cy*Cz+iz*Cy+iy];
                        }
                    }
                    // Post the next receive
                    MPI_Irecv(buffer,2*Cy*Cz,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*idd+4,MPI_COMM_WORLD,&req);
                    iblock->receive_id.request[4]=req;
                }
        }
    }
#else
    /* Top - pass down */
    if(bz>0){
        idd=id-BBx*BBy;
        if(iblock->gblock->iblock[idd].physics.field_type<40){
            // Pass to the Same CPU
            Tx=iblock->gblock->iblock[idd].pass.Tx; 
            Ty=iblock->gblock->iblock[idd].pass.Ty;
            for ( iy=0; iy<Cy; ++iy){
                for ( ix=0; ix<Cx; ++ix){
                    Tx[iy][ix]=Ex[0][iy][ix];
                    Ty[iy][ix]=Ey[0][iy][ix];
                } 
            } 
        }
    }

    /* bAck - pass front */
    if(by>0){
        idd=id-BBx;
        if(iblock->gblock->iblock[idd].physics.field_type<40){
            // Pass to the Same CPU
            Ax=iblock->gblock->iblock[idd].pass.Ax; 
            Az=iblock->gblock->iblock[idd].pass.Az;
            for ( iz=0; iz<Cz; ++iz){
                for ( ix=0; ix<Cx; ++ix){
                    Ax[iz][ix]=Ex[iz][0][ix];
                    Az[iz][ix]=Ez[iz][0][ix];
                }
            }
        }
    }

    /* Right - pass left */
    if(bx>0){
        idd=id-1;
        if(iblock->gblock->iblock[idd].physics.field_type<40){
            // Pass to the Same CPU
            Ry=iblock->gblock->iblock[idd].pass.Ry; 
            Rz=iblock->gblock->iblock[idd].pass.Rz;
            for ( iz=0; iz<Cz; ++iz){
                for ( iy=0; iy<Cy; ++iy){
                    Ry[iz][iy]=Ey[iz][iy][0];
                    Rz[iz][iy]=Ez[iz][iy][0];
                }
            }
        }
    }
#endif
}     

/*****************************************************************************/

void send_aPML_H_fields(I_Block *iblock, int flag)
{
    int ix, iy, iz;
    int Cxm, Cym, Czm, idd;
    REALNO **Bx, **By, **Ly, **Lz, **Fx, **Fz;

    //short cuts
    REALNO ***Hx, ***Hy, ***Hz;
    int Cx, Cy, Cz, id;
    int BBx,BBy,BBz,bx,by,bz;
    Hx=iblock->max.Hx; Hy=iblock->max.Hy; Hz=iblock->max.Hz;
    Cx=iblock->length.x; Cy=iblock->length.y; Cz=iblock->length.z;
    id=iblock->id;
    BBx=iblock->gblock->BBx; BBy=iblock->gblock->BBy; BBz=iblock->gblock->BBz;
    bz=(int)id/(BBx*BBy);
    by=(int)(id-bz*BBx*BBy)/BBx;
    bx=(int)(id-bz*BBx*BBy-by*BBx);
#ifdef MPIRUN
    MPI_Status pstatus;
    //static REALNO buffer[BNDRYbuffy];
    REALNO *buffer;
    MPI_Request req;
#endif

    Czm=Cz-1; Cym=Cy-1; Cxm=Cx-1;
#ifdef MPIRUN
    if(flag>0){ /* send mode */
        if((2*Cx*Cy>BNDRYbuffy) || (2*Cx*Cz>BNDRYbuffy) || (2*Cy*Cz>BNDRYbuffy) ){
            printf("Error in static memory in routine send_PML_H_fields \n");
            printf("Must increase buffer from %d to max of\n",BNDRYbuffy);
            printf("%d %d %d\n",2*Cx*Cy,2*Cx*Cz,2*Cy*Cz);
            exit(1);
        }

        //Tags
        // 10 *(1+rank of cpu sent to) + numbers 0 ->11
        // where 0->11 are defined for
        //  0 =>Tx  1=>Ty
        //  2 =>Bx  3=>By
        //  4 =>Ax  5=>Az
        //  6 =>Fx  7=>Fz
        //  8 =>Ry  9=>Rz
        // 10 =>Ly 11=>Lz

        // define a unique tag as
        // 10*(the sending block id)+id of the face send to
        // where the faces are define as
        // 0 top   1 bottom
        // 2 back  3 front
        // 4 right 5 left

        /* bottom - pass up */
        if(bz<(BBz-1)){
            /*    printf("      valid block to pass up\n"); */
            idd=id+BBx*BBy;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]==iblock->gblock->rank){ 
                    // Pass to same CPU
                    Bx=iblock->gblock->iblock[idd].pass.Bx; 
                    By=iblock->gblock->iblock[idd].pass.By;
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            Bx[iy][ix]=Hx[Czm][iy][ix];
                            By[iy][ix]=Hy[Czm][iy][ix];
                        }
                    }
                }else{
                    // Make sure last group was sent before packing
                    MPI_Wait(&(iblock->send_id.request[0]),&pstatus);
                    buffer=iblock->send.face[0];
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[iy*Cx+ix]=Hx[Czm][iy][ix];
                        }
                    }
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[Cx*Cy+iy*Cx+ix]=Hy[Czm][iy][ix];
                        }
                    }
                    MPI_Isend(buffer,2*Cx*Cy,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*id+1,MPI_COMM_WORLD,&req);  
                    iblock->send_id.request[0]=req;
                }
        }

        /* front - pass back */
        if(by<(BBy-1)){
            idd=id+BBx;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]==iblock->gblock->rank){ 
                    // Pass to same CPU
                    Fx=iblock->gblock->iblock[idd].pass.Fx; 
                    Fz=iblock->gblock->iblock[idd].pass.Fz;
                    for (iz=0; iz<Cz; ++iz){
                        for (ix=0; ix<Cx; ++ix){
                            Fx[iz][ix]=Hx[iz][Cym][ix];
                            Fz[iz][ix]=Hz[iz][Cym][ix];
                        }
                    }
                }else{
                    // Make sure last group was sent before packing
                    MPI_Wait(&(iblock->send_id.request[2]),&pstatus);
                    buffer=iblock->send.face[2];
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[iz*Cx+ix]=Hx[iz][Cym][ix];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            buffer[Cx*Cz+iz*Cx+ix]=Hz[iz][Cym][ix];
                        }
                    }
                    MPI_Isend(buffer,2*Cz*Cx,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*id+3,MPI_COMM_WORLD,&req);
                    iblock->send_id.request[2]=req;
                }
        }

        /* left - pass right */
        if(bx<(BBx-1)){
            idd=id+1;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]==iblock->gblock->rank){ 
                    // Pass to same CPU
                    Ly=iblock->gblock->iblock[idd].pass.Ly; 
                    Lz=iblock->gblock->iblock[idd].pass.Lz;
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            Ly[iz][iy]=Hy[iz][iy][Cxm];
                            Lz[iz][iy]=Hz[iz][iy][Cxm];
                        } 
                    } 
                }else{
                    // Make sure last group was sent before packing
                    MPI_Wait(&(iblock->send_id.request[4]),&pstatus);
                    buffer=iblock->send.face[4];
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            buffer[iz*Cy+iy]=Hy[iz][iy][Cxm];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            buffer[Cy*Cz+iz*Cy+iy]=Hz[iz][iy][Cxm];
                        }
                    }
                    MPI_Isend(buffer,2*Cz*Cy,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*id+5,MPI_COMM_WORLD,&req);
                    iblock->send_id.request[4]=req;
                }
        }
    }else{ /* recieve mode */
        /* bottom - pass up */
        if(bz>0){ //recieve data
            idd=id-BBx*BBy;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]!=iblock->gblock->rank){ 
                    // Some other cpu sent this
                    Bx=iblock->gblock->iblock[id].pass.Bx; 
                    By=iblock->gblock->iblock[id].pass.By;
                    // Make sure data was received before unpacking
                    MPI_Wait(&(iblock->receive_id.request[1]),&pstatus);
                    buffer=iblock->receive.face[1];
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            Bx[iy][ix]=buffer[iy*Cx+ix];
                        }
                    }
                    for ( iy=0; iy<Cy; ++iy){
                        for ( ix=0; ix<Cx; ++ix){
                            By[iy][ix]=buffer[Cx*Cy+iy*Cx+ix];
                        }
                    }
                    // Post the next receive
                    MPI_Irecv(buffer,2*Cx*Cy,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*idd+1,MPI_COMM_WORLD,&req);
                    iblock->receive_id.request[1]=req;
                }
        }

        /* front - pass back */
        if(by>0){ //recieve data
            idd=id-BBx;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]!=iblock->gblock->rank){ 
                    // Some other cpu sent this
                    Fx=iblock->gblock->iblock[id].pass.Fx; 
                    Fz=iblock->gblock->iblock[id].pass.Fz;
                    // Make sure data was received before unpacking
                    MPI_Wait(&(iblock->receive_id.request[3]),&pstatus);
                    buffer=iblock->receive.face[3];
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            Fx[iz][ix]=buffer[iz*Cx+ix];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( ix=0; ix<Cx; ++ix){
                            Fz[iz][ix]=buffer[Cx*Cz+iz*Cx+ix];
                        }
                    }
                    // Post the next receive
                    MPI_Irecv(buffer,2*Cx*Cz,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*idd+3,MPI_COMM_WORLD,&req);
                    iblock->receive_id.request[3]=req;
                }
        }

        /* left - pass right */
        if(bx>0){ //recieve data
            idd=id-1;
            if(iblock->gblock->iblock[idd].physics.field_type<40)
                if(iblock->gblock->useblock[idd]!=iblock->gblock->rank){ 
                    // Some other cpu sent this
                    Ly=iblock->gblock->iblock[id].pass.Ly; 
                    Lz=iblock->gblock->iblock[id].pass.Lz;
                    // Make sure data was received before unpacking
                    MPI_Wait(&(iblock->receive_id.request[5]),&pstatus);
                    buffer=iblock->receive.face[5];
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            Ly[iz][iy]=buffer[iz*Cy+iy];
                        }
                    }
                    for ( iz=0; iz<Cz; ++iz){
                        for ( iy=0; iy<Cy; ++iy){
                            Lz[iz][iy]=buffer[Cy*Cz+iz*Cy+iy];
                        }
                    }
                    // Post the next receive
                    MPI_Irecv(buffer,2*Cy*Cz,MPI_REALNO,iblock->gblock->useblock[idd],
                            10*idd+5,MPI_COMM_WORLD,&req);
                    iblock->receive_id.request[5]=req;
                }
        }    
    }
#else
    /* bottom - pass up */
    if(bz<(BBz-1)){
        idd=id+BBx*BBy;
        if(iblock->gblock->iblock[idd].physics.field_type<40){
            // Pass to same CPU
            Bx=iblock->gblock->iblock[idd].pass.Bx; 
            By=iblock->gblock->iblock[idd].pass.By;
            for ( iy=0; iy<Cy; ++iy){
                for ( ix=0; ix<Cx; ++ix){
                    Bx[iy][ix]=Hx[Czm][iy][ix];
                    By[iy][ix]=Hy[Czm][iy][ix];
                }
            }
        }
    }

    /* front - pass back */
    if(by<(BBy-1)){
        idd=id+BBx;
        if(iblock->gblock->iblock[idd].physics.field_type<40){
            // Pass to same CPU
            Fx=iblock->gblock->iblock[idd].pass.Fx; 
            Fz=iblock->gblock->iblock[idd].pass.Fz;
            for (iz=0; iz<Cz; ++iz){
                for (ix=0; ix<Cx; ++ix){
                    Fx[iz][ix]=Hx[iz][Cym][ix];
                    Fz[iz][ix]=Hz[iz][Cym][ix];
                }
            }
        }
    }

    /* left - pass right */
    if(bx<(BBx-1)){
        idd=id+1;
        if(iblock->gblock->iblock[idd].physics.field_type<40){
            // Pass to same CPU
            Ly=iblock->gblock->iblock[idd].pass.Ly; 
            Lz=iblock->gblock->iblock[idd].pass.Lz;
            for ( iz=0; iz<Cz; ++iz){
                for ( iy=0; iy<Cy; ++iy){
                    Ly[iz][iy]=Hy[iz][iy][Cxm];
                    Lz[iz][iy]=Hz[iz][iy][Cxm];
                } 
            } 
        }
    }

#endif
}     

/*****************************************************************************/
#endif

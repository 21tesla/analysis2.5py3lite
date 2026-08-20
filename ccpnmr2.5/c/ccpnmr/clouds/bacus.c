/*
======================COPYRIGHT/LICENSE START==========================

bacus.c: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Carlos Steren, Alexei Grishaev, Miguel Llinas, Wayne Boucher and Tim Stevens (Carnegie Mellon University and University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
*/

#include "bacus.h"

#ifndef CCPN_FALSE
#define CCPN_FALSE     0
#endif
#ifndef CCPN_TRUE
#define CCPN_TRUE      1
#endif

/*
#define NPK  5000
#define NFR  500
#define NP 200
*/
#define NPK  10000
#define NFR  2500
#define NP 500

#define NPR 100
#define NPRF 50


static float thrhn=0.02;
static float thrhh=0.04;
static float thrnn=0.25;
static float thrcc=0.1;

static float errdiag=0.04;

static float pcosy=0.9998;
static float ptocsy=0.6907;

static float pj11=1.0;
static float pj12=0.412964;
static float pj13=0.277762;
static float pj22=0.32238;
static float pj23=0.244056;
static float pj33=0.201548;

static float pnj11=1.0;
static float pnj12=0.009575;
static float pnj13=0.011237;
static float pnj22=0.01089;
static float pnj23=0.011839;
static float pnj33=0.012398;

static CcpnStatus alloc_memory(int **cylnk, int **tylnk, int **prioras1, int **prioras2,
                int **postas1, int **postas2, float **priorprob,
                float **postprob, float **probmap, float **lkprob)
{
    int i;

    for (i = 0; i < NFR; i++)
    {
        MALLOC(cylnk[i], int, NP);
        MALLOC(tylnk[i], int, NP);

        MALLOC(probmap[i], float, NFR);
    }

    for (i = 0; i < NPK; i++)
    {
        MALLOC(prioras1[i], int, NP);
        MALLOC(prioras2[i], int, NP);
        MALLOC(postas1[i], int, NP);
        MALLOC(postas2[i], int, NP);

        MALLOC(priorprob[i], float, NP);
        MALLOC(postprob[i], float, NP);
        MALLOC(lkprob[i], float, NP);
    }

    return CCPN_OK;
}

static void free_memory(int **cylnk, int **tylnk, int **prioras1, int **prioras2,
                int **postas1, int **postas2, float **priorprob,
                float **postprob, float **probmap, float **lkprob)
{
    int i;

    for (i = 0; i < NFR; i++)
    {
        FREE(cylnk[i], int);
        FREE(tylnk[i], int);

        FREE(probmap[i], float);
    }

    for (i = 0; i < NPK; i++)
    {
        FREE(prioras1[i], int);
        FREE(prioras2[i], int);
        FREE(postas1[i], int);
        FREE(postas2[i], int);

        FREE(priorprob[i], float);
        FREE(postprob[i], float);
        FREE(lkprob[i], float);
    }
}

static void print_msg(char stip[])
{
    printf(" finished step %s\n",stip);
    fflush(stdout);
    return;
}

static CcpnStatus readnoes(char *name, char *exten, int *ndata, float data[][3],
      int newnoe[],float newvol[], char *error_msg)
{
    FILE *fp;
    char infile[20];
    int DONE=CCPN_FALSE,it;
    strcpy(infile,name);  
    strcat(infile,exten);  
    while(!DONE)
    {
        printf("file to open: %s  \n",infile);
        if((fp=fopen(infile,"r"))==NULL)
        {
            sprintf(error_msg, "readnoes: cannot open file '%s'", infile);
            return CCPN_ERROR;
        }
        else
        {
            DONE= CCPN_TRUE;
        }
    }
    it=0;
    while(!feof(fp))
    {
        if (it >= NPK)
        {
            sprintf(error_msg, "readnoes: code compiled for no more than %d peaks", NPK);
            return CCPN_ERROR;
        }
        fscanf(fp," %d %f %f %f %e ",
            &newnoe[it],&data[it][0],&data[it][2],&data[it][1],&newvol[it]);
/*     printf("%d %d %f %f %f %e \n",
       it,newnoe[it],data[it][0],data[it][1],data[it][2],newvol[it]); */

        it++;
    }    
    fclose(fp);
    *ndata=it;

    return CCPN_OK; 
}

static CcpnStatus readfreqlst(char *name,char *exten,int *nfreqs,float freqs[][2],char heavy[][2],
       int ngroup[],char group[][2],char atom[][7], char *error_msg)
{
    FILE *fp;
    char infile[20];
    int DONE,it,tmpint;
    strcpy(infile,name);  
    strcat(infile,exten);
    DONE=CCPN_FALSE;
    while(!DONE)
    {
        printf("file to open: %s  \n",infile);
        if((fp=fopen(infile,"r"))==NULL)
        {
            sprintf(error_msg, "readfreqlst: cannot open file '%s'", infile);
            return CCPN_ERROR;
        }
        else
        {
            DONE= CCPN_TRUE;
        }
    }
    it=0;
    while(!feof(fp))
    {
        if (it >= NFR)
        {
            sprintf(error_msg, "readfreqlst: code compiled for no more than %d freqs", NFR);
            return CCPN_ERROR;
        }
        fscanf(fp,"%d %f %f %s %s %d %s ",
            &tmpint,&freqs[it][1],&freqs[it][0],atom[it],group[it],&ngroup[it],heavy[it]);
/*
     printf("%d %f %f %s %s %d %s \n",
         tmpint,freqs[it][0],freqs[it][0],atom[it],group[it],ngroup[it],heavy[it]);
*/
        it++;
    }    
    fclose(fp);
    *nfreqs=it;
    printf(" %d frequencies read from %s ! \n",*nfreqs,infile);
    printf("read: freq[0]:%f,freq[1]:%f  \n",freqs[0][0],freqs[0][1]);

    return CCPN_OK; 
}

static CcpnStatus readlinkfile(char *name,char *exten,int nfreqs,int numlnk[],int *datlnk[], char *error_msg)
{
    FILE *fp;
    char infile[20];
    int DONE=CCPN_FALSE,it,ij,tmpint;
    float tmpf;
    strcpy(infile,name);  
    strcat(infile,exten);  
    while(!DONE)
    {
        printf("file to open: %s  \n",infile);
        if((fp=fopen(infile,"r"))==NULL)
        {
            sprintf(error_msg, "readlink: cannot open file '%s'", infile);
            return CCPN_ERROR;
        }
        else
        {
            DONE= CCPN_TRUE;
        }
    }

    for(it=0;it<nfreqs;it++)
    {
        if (it >= NFR)
        {
            sprintf(error_msg, "readlink: code compiled for no more than %d freqs", NFR);
            return CCPN_ERROR;
        }

        fscanf(fp,"%d %f %d ",&tmpint,&tmpf,&numlnk[it]);
        if (numlnk[it] > NP)
        {
            sprintf(error_msg, "readlink: code compiled for no more than %d linkages", NP);
            return CCPN_ERROR;
        }

/*
       printf("n=%d,freq=%f,num=%d \n",tmpint,tmpf,numlnk[it]); 
*/
        if(numlnk[it]>0) for(ij=0;ij<numlnk[it];ij++) 
            fscanf(fp,"%d ",&datlnk[it][ij]);
    }    
    
    fclose(fp);
    printf(" %d data from %s read !! \n",nfreqs,infile);

    return CCPN_OK; 
} 


static CcpnStatus readpriors(char *name, char *exten, int *ndata, int numnoe[],int npriors[],float *priors[],
      int *prioas1[],int *prioas2[], char *error_msg)
{
    FILE *fp;
    char infile[20],jk1[10];
    int DONE=CCPN_FALSE,it,ij,ntmp;
    float tmp1,tmp2,tmp3;

    strcpy(infile,name);  
    strcat(infile,exten);  
    while(!DONE)
    {
        printf("file to open: %s  \n",infile);
        if((fp=fopen(infile,"r"))==NULL)
        {
            sprintf(error_msg, "readpriors: cannot open file '%s'", infile);
            return CCPN_ERROR;
        }
        else
        {
            DONE= CCPN_TRUE;
        }
    }
    it=0;
    while(!feof(fp))
    {
        fscanf(fp,"%d %d",&numnoe[it],&npriors[it]);
/*
      printf("%d %d %d \n",it,numnoe[it],npriors[it]);
*/
        if(npriors[it]>0) for(ij=0;ij<npriors[it];ij++)
        {
            fscanf(fp,"%d %d %d %s %s %s %s %s %s ",
                &ntmp,&prioas1[it][ij],&prioas2[it][ij],jk1,jk1,jk1,jk1,jk1,jk1);

/*    printf("%d %d %d %s %s %s %s %s %s  \n",
    ntmp,prioas1[it][ij],prioas2[it][ij],jk1,jk1,jk1,jk1,jk1,jk1);
    getch(); */

            fscanf(fp,"%f %f %f %e",&priors[it][ij],&tmp1,&tmp2,&tmp3);
/*    printf("%f %f %f %e \n",priors[it][ij],tmp1,tmp2,tmp3); 
    getch(); */
        }
        it++;
    }    
    fclose(fp);
    *ndata=it;
    printf("%d data from %s read !! \n",*ndata,infile);

    return CCPN_OK; 
}


static void renorm(int ndata,int nprob[],float *prob[])
{
    int it,ij;
    float sum;

    for(it=0;it<ndata;it++) if(nprob[it]>0)
    {
        sum=0.0;
        for(ij=0;ij<nprob[it];ij++) sum=sum+prob[it][ij];
        for(ij=0;ij<nprob[it];ij++) prob[it][ij]=prob[it][ij]/sum;
    }
}     


static void filterdiag(int newnum,float newdata[][3],int newdiag[])
{
    int it;
    float dif;
    for(it=0;it<newnum;it++) 
    {
        newdiag[it]=0;
        dif=fabs(newdata[it][1]-newdata[it][2])/errdiag;
        if(dif<=2.0) newdiag[it]=1;
    }
}


static CcpnStatus addtolist(int *totnoe,float data[][3],char type[][2],int nunoe[],int rownum[],float vol[],int newnum,
    int newnoe[],float newdata[][3],char newtype[],float newvol[],int newdiag[], char *error_msg)
{
    int it,count, orig_count;
    count=orig_count=*totnoe;
    printf("start  addlist count:%d \n",count);

    for(it=0;it<newnum;it++) if(newdiag[it]==0)
    {
        if (count >= NPK)
        {
            sprintf(error_msg, "addtolist: code compiled for no more than %d peaks", NPK);
            return CCPN_ERROR;
        }

        data[count][0]=newdata[it][0];
        data[count][1]=newdata[it][1];
        data[count][2]=newdata[it][2];
        nunoe[count]=newnoe[it];
        vol[count]=newvol[it];
        strcpy(type[count],newtype);
        rownum[count] = it + orig_count;
        count++;
    }
    *totnoe=count;
    printf("end in addlist count:%d \n",count);

    return CCPN_OK;
}


static void inclupr(float *pro,int idt,int beli[],float priors[],int *num,float *min)
{
        /*
    int it,ij,posit,belpr[101],count,minus,manis;
    float vectpr[101],maxi;
        */
    int it,ij,posit,belpr[NPR+1],count,minus,manis;
    float vectpr[NPR+1],maxi;
    
    count=*num;

    if (count<NPR)
    {
        beli[count]=idt;
        priors[count]=*pro;
        count++;
        *num=count;
        *min=0.0;

        return;  /* wb104 */
    }
    minus=(NPR-1);
    /* O(N^2) descending sort on probability */
    if (count>minus) 
    {
        for(it=0;it<NPR;it++)
        {
            vectpr[it]=priors[it];
            belpr[it]=beli[it];
        }
        vectpr[NPR]=*pro;
        belpr[NPR]=idt;
    
        for(it=0;it<NPR;it++)
        {
            maxi=0.0;
            manis=NPR;
            for(ij=0;ij<manis;ij++)
            {
                if(vectpr[ij]>maxi)
                {
                    maxi=vectpr[ij];
                    posit=ij;
                }
            }
            priors[it]=maxi;
            beli[it]=belpr[posit];
            vectpr[posit]=-999;
        }
        *min=priors[minus];
    }
}

    
static void calcprof2(float noeval[],char *noetype,int nfreq,float freqlst[][2],
        int belin[],char heavy[][2],float priorin[])
{
    float dif1,dif2,pr1,pr2,pro,minu,sigma;
    int it,num;
    minu=0.0;
    num=0;
    if(!strcmp(noetype,"N")) sigma=thrnn;
    if(!strcmp(noetype,"C")) sigma=thrcc;

    for(it=0;it<NPR;it++) priorin[it]=0.0;
    for(it=0;it<nfreq;it++) if(!strcmp(heavy[it],noetype)) 
    {
        dif1=fabs(noeval[0]-freqlst[it][0])/sigma;
        dif2=fabs(noeval[1]-freqlst[it][1])/thrhn;
        if((dif1<2.0) & (dif2<2.0))     
        {
            pr1=exp(-0.5*dif1*dif1);
            pr2=exp(-0.5*dif2*dif2);
            pro=pr1*pr2;
            if(pro>minu) inclupr(&pro,it,belin,priorin,&num,&minu);
        }
    }
}


static void calcprof1(float noeval[],int nfreq,float freqlst[][2],int belin[],float priorin[])
{
    float dif1,pr1,pro,min;
    int it,num;
    min=0.0;
    num=0;
    for(it=0;it<NPR;it++) priorin[it]=0.0;
    for(it=0;it<nfreq;it++) 
    {
        dif1=fabs(noeval[0]-freqlst[it][1])/thrhh;
        if(dif1<2.0)     
        {
            pr1=exp(-0.5*dif1*dif1);
            pro=pr1;
            inclupr(&pro,it,belin,priorin,&num,&min);
        }
    }
}


static CcpnStatus calctot(int idt,float priors1[],float priors2[],int belin1[],int belin2[],
        int *beli1[],int *beli2[],float *priors[],int nupriors[], char *error_msg)
{
    int it,ij,count1,count2,belpr[NFR][2],posi,indi1,indi2;
    float vect[NFR],probi,maxi;

    count1=0;

    for(it=0;it<NPRF;it++) for(ij=0;ij<NPRF;ij++) 
    {
        probi=priors1[it]*priors2[ij];
        if(probi>0.0)
        {
            if (count1 >= NFR)
            {
                sprintf(error_msg, "calctot: code compiled for no more than %d probs", NFR);
                return CCPN_ERROR;
            }
            vect[count1]=probi;
            belpr[count1][0]=belin1[it];
            belpr[count1][1]=belin2[ij];
            count1++;
        }
    }
    
    count2=0;
    if(count1>1)
    {
        indi1=count1;
        indi2=count1;
        if(count1>NPRF) indi1=NPRF;
        for(it=0;it<indi1;it++) 
        {
            maxi=0.0;
            for(ij=0;ij<indi2;ij++)
            {
                if(vect[ij]>maxi) 
                {
                    maxi=vect[ij];
                        posi=ij;
                }
            }
            if (count2 >= NP)
            {
                sprintf(error_msg, "calctot: code compiled for no more than %d probs2", NP);
                return CCPN_ERROR;
            }
            priors[idt][count2]=vect[posi];
            beli1[idt][count2]=belpr[posi][0];
            beli2[idt][count2]=belpr[posi][1];
            vect[posi]=-999;
            count2++;    
        }
    }

    if(count1==1)
    {
        priors[idt][0]=vect[0];
        beli1[idt][0]=belpr[0][0];
        beli2[idt][0]=belpr[0][1];
        count2++;
    }

    if(count1==0)
    {
        priors[idt][0]=0.0;
        beli1[idt][0]=-999;
        beli2[idt][0]=-999;
        if(!(priors1[it]==0.0))     beli1[idt][0]=999;
        if(!(priors2[it]==0.0))     beli2[idt][0]=999;
    }
    nupriors[idt]=count2;
    
    return CCPN_OK;
}


static CcpnStatus calcpriors(int ndata,float data[][3],char type[][2],int nfreqs,
        float freqlst[][2],char heavy[][2],int nupriors[], int *beli1[],
        int *beli2[],float *priors[], char *error_msg)
{
    int it,belin1[NPR],belin2[NPR];
    float temp1[2],temp2[2],priors1[NPR],priors2[NPR];
    char noetype[2];

    for(it=0;it<ndata;it++) nupriors[it]=0;

    for(it=0;it<ndata;it++)
    { 
        temp1[0]=data[it][0];
        temp1[1]=data[it][1];
        strcpy(noetype,type[it]);

        temp2[0]=data[it][2];
        temp2[1]=-999.00;

        calcprof2(temp1,noetype,nfreqs,freqlst,belin1,heavy,priors1);

        calcprof1(temp2,nfreqs,freqlst,belin2,priors2);

/*    for(ij=0;ij<10;ij++) if(!(priors2[ij]==0.0))
    printf("belin2[ij]:%d  ",belin2[ij]);
    printf(" \n");
    getch();  */

        CHECK_STATUS(calctot(it,priors1,priors2,belin1,belin2,beli1,beli2,priors,nupriors,error_msg));

/*    for(ij=0;ij<nupriors[it];ij++) 
    printf("beli:%d %d   ",beli1[it][ij],beli2[it][ij]);
    printf(" \n");
    getch();  */
        /*
    printf(" it:%d number of priors:%d  \n",it,nupriors[it]);
        */

    }

    return CCPN_OK;
}


static void shannon(int ndata, int nuprobs[],float *probs[], float *shan)
{
    int it,ij;
    float sum[NPK],totsum;
    totsum=0.0;
    for(it=0;it<ndata;it++) sum[it]=0.0;
    
    for(it=0;it<ndata;it++)  if(nuprobs[it]>0) 
    for(ij=0;ij<nuprobs[it];ij++) 
    {
    if(probs[it][ij]>0)  
        sum[it]=sum[it]-probs[it][ij]*log(probs[it][ij]);
    }
    for(it=0;it<ndata;it++) sum[it]=exp(sum[it]);
    for(it=0;it<ndata;it++) totsum=totsum+sum[it];
    *shan=totsum;
}


static int islink(int idt1,int idt2,int nlink[],int *link[])
{
    int it,num;
    num=nlink[idt1];
    for(it=0;it<num;it++)     if(link[idt1][it]==idt2)  return 1;
    return 0;
}


static int issamegroup(int idt1,int idt2,int ncylnk[],int *cylnk[],
        int ntylnk[],int *tylnk[],float *likep)
{
    if(islink(idt1,idt2,ncylnk,cylnk))
    {
        *likep=pcosy;
        return 1;
    }
    if(islink(idt1,idt2,ntylnk,tylnk))
    {
        *likep=ptocsy;
        return 1;
    }
    return 0;
}


static void findreporters(int idt,int ncylnk[],int *cylnk[],int ntylnk[],
        int *tylnk[],int *nvect,int vect[],int clvect[])
{
    int it,count;
    count=0;
    vect[count]=idt;
    clvect[count]=1;
    count++;
    if(ncylnk[idt]>0)  for(it=0;it<ncylnk[idt];it++)
    {
        vect[count]=cylnk[idt][it];
        clvect[count]=2;
        count++;
    }

    if(ntylnk[idt]>0) for(it=0;it<ntylnk[idt];it++)
    {
        vect[count]=tylnk[idt][it];
        clvect[count]=3;
        count++;
    }
    *nvect=count;
}     



static void calclkreport(int idt1,int idt2,int nvect1,int vect1[],int clvect1[],int nvect2,int vect2[],
/*
        int clvect2[],float *totsum,float probmap[][NFR])
*/
        int clvect2[],float *totsum,float *probmap[])
{
    int it,ij,vet1,vet2,cl1,cl2;
    float sum,norm,pscore,pj[4][4],pnj[4][4],term1,term2;

    pj[1][1]=pj11;
    pj[1][2]=pj12;
    pj[1][3]=pj13;
    pj[2][2]=pj22;
    pj[2][3]=pj23;
    pj[3][3]=pj33;
    pj[2][1]=pj[1][2];
    pj[3][1]=pj[1][3];
    pj[3][2]=pj[2][3];

    pnj[1][1]=pnj11;
    pnj[1][2]=pnj12;
    pnj[1][3]=pnj13;
    pnj[2][2]=pnj22;
    pnj[2][3]=pnj23;
    pnj[3][3]=pnj33;
    pnj[2][1]=pnj[1][2];
    pnj[3][1]=pnj[1][3];
    pnj[3][2]=pnj[2][3];

    sum=0.0;
    norm=0.0;


    for(it=0;it<nvect1;it++)
    {
        for(ij=0;ij<nvect2;ij++)
        {
            if(!((clvect1[it]==1) && (clvect2[ij]==1)))
            {

/*    if((idt1==250) && (idt2==246)) 
    {
          printf("\n"); 
    printf("initial sum:%f \n",sum); 
    } */

                vet1=vect1[it];
                vet2=vect2[ij];
                cl1=clvect1[it];
                cl2=clvect2[ij];
                pscore=probmap[vet1][vet2];
/*
          calcpscore(vet1,vet2,nvect1,nvect2,bel1,bel2,priors,&pscore);
*/
                term1=pscore*pj[cl1][cl2];
                term2=(1.0-pscore)*pnj[cl1][cl2];
                sum=sum+term1+term2;
/*
                sum=sum+pscore*pj[cl1][cl2]+(1.0-pscore)*pnj[cl1][cl2];
*/

            norm=norm+pj[cl1][cl2];

/*    if((idt1==250) && (idt2==246)) 
    {
     printf("idt1:%d  idt2:%d  sum:%f ,norm:%f pscore:%f \n",idt1,idt2,sum,norm,pscore); 
    printf("cl1:%d,cl2:%d,term1:%f,term2:%f, pnj:%f \n",cl1,cl2,term1,term2,pnj[cl1][cl2]); 
    getch();
    } */

            }
        }
    }

    sum=sum/norm;
    *totsum=sum;
/*    if((idt1==250) && (idt2==246)) 
    {
     for(it=0;it<nvect1;it++) 
            printf("vect1:%d,clvect1:%d \n",vect1[it],clvect1[it]);
     for(it=0;it<nvect2;it++) 
            printf("vect2:%d,clvect2:%d \n",vect2[it],clvect2[it]);

     printf("idt1:%d  idt2:%d  sum:%f ,norm:%f  \n",idt1,idt2,sum,norm); 
    getch();
    }  */
}



/*
int likehood(int totnoe,int ncylnk[],int cylnk[][NP],int ntylnk[],int tylnk[][NP],
        int npost[],float postprob[][NP],int postas1[][NP],int postas2[][NP],
        float probmap[][NFR],float lkprob[][NP])
*/
static void likehood(int totnoe,int ncylnk[],int *cylnk[],int ntylnk[],int *tylnk[],
        int npost[],float *postprob[],int *postas1[],int *postas2[],
        float *probmap[],float *lkprob[])
{
    int it,ij,idt1,idt2,nvect1,nvect2,clvect1[NP],clvect2[NP];
    int vect1[NP],vect2[NP];
    float likep,totsum;

    for(it=0;it<totnoe;it++) if(npost[it]>1) 
    {
        for(ij=0;ij<npost[it];ij++)
        {
            idt1=postas1[it][ij];
            idt2=postas2[it][ij];
            if(issamegroup(idt1,idt2,ncylnk,cylnk,ntylnk,tylnk,&likep)) 
            {
                lkprob[it][ij]=likep;
            }
            else 
            {    
                findreporters(idt1,ncylnk,cylnk,ntylnk,tylnk,&nvect1,vect1,clvect1);
                findreporters(idt2,ncylnk,cylnk,ntylnk,tylnk,&nvect2,vect2,clvect2);
                calclkreport(idt1,idt2,nvect1,vect1,clvect1,nvect2,vect2,clvect2,&totsum,probmap);
                lkprob[it][ij]=totsum;

/*    if(it==8) 
    {
    printf("idt1:%d  idt2:%d  nvect1:%d ,nvect2:%d prob:%5.3f \n",idt1,idt2,nvect1,nvect2,totsum); 
    getch();
    } */

            } 
        }
    }
}   




/*
int calcpscore(int idt1,int idt2,int nv1,int nv2,int bel1[][NP],int bel2[][NP],
        float priors[][NP],float *pscore)
*/
static int calcpscore(int idt1,int idt2,int nv1,int nv2,int *bel1[],int *bel2[],
        float *priors[],float *pscore)
{
    int it,ij;
    for(it=0;it<nv1;it++)
    {
    for(ij=0;ij<nv2;ij++)
    {
    if((bel1[it][ij]==idt1) & (bel2[it][ij]==idt2))
    {
    *pscore=priors[it][ij];
    return 0;
    }
    }    
    }
    *pscore=0;
    return 0;
}


/*
int writeass(char *name,char *exten,int totnoe,int npost[],float postprob[][NP],
      int postas1[][NP],int postas2[][NP],int ngroup[],char group[][2],char atom[][7])
*/
static int writeass(char *name,char *exten,int totnoe,int rownum[], int npost[],float *postprob[],
      int *postas1[],int *postas2[],int ngroup[],char group[][2],char atom[][7])
{
    FILE *fp;
    int it,idt1,idt2,count;
    char outfile[20];
    strcpy(outfile,name);  
    strcat(outfile,exten);  
    printf(" writing the file: %s \n",outfile);
        if((fp=fopen(outfile,"w"))==NULL) {
        puts("cannot open file\n");
        return 0;
        }
        
     count=1;
     for(it=0;it<totnoe;it++)     if(npost[it]==1)
     {
     idt1=postas1[it][0];
     idt2=postas2[it][0];
/*
       fprintf(fp,"%5d %4d %4d    %1s %2d %6s %1s %2d %6s %f %f %f %e\n",count,idt1,idt2,
*/
     fprintf(fp,"%5d %5d %4d %4d    %1s %2d %-6s %1s %2d %-6s %f 0.00 0.000 1.0e03\n",count,rownum[it],idt1+1,idt2+1,
         group[idt1],ngroup[idt1],atom[idt1],group[idt2],ngroup[idt2],atom[idt2],postprob[it][0]);
     count++;
     }
     fclose(fp);
     return 0;
}


/*
int writepos(char *name,char *exten,int totnoe,int npost[],float postprob[][NP],
      int postas1[][NP],int postas2[][NP],int ngroup[],char group[][2],char atom[][7])
*/
static int writepos(char *name,char *exten,int totnoe,int npost[],float *postprob[],
      int *postas1[],int *postas2[],int ngroup[],char group[][2],char atom[][7])
{
    FILE *fp;
    int it,ij,idt1,idt2,count;
    char outfile[20];
    strcpy(outfile,name);  
    strcat(outfile,exten);  
    printf(" writing the file: %s \n",outfile);
        if((fp=fopen(outfile,"w"))==NULL) {
        puts("cannot open file\n");
        return 0;
        }
     printf("in writepos totnoe: %d \n",totnoe);
     count=1;
     for(it=0;it<totnoe;it++)     
     {
/*
       printf("%5d %4d \n",count,npost[it]);
*/
     fprintf(fp,"%5d %4d \n",count,npost[it]);
     if(npost[it]>0) for(ij=0;ij<npost[it];ij++)     
     {
     idt1=postas1[it][ij];
     idt2=postas2[it][ij];
/*
       fprintf(fp,"%5d %4d %4d    %1s %2d %6s %1s %2d %6s %f %f %f %e\n",count,idt1,idt2,
*/
     fprintf(fp,"%5d %4d %4d    %1s %2d %-6s %1s %2d %-6s %f 0.00 0.000 1.0e03\n",ij+1,idt1,idt2,
         group[idt1],ngroup[idt1],atom[idt1],group[idt2],ngroup[idt2],atom[idt2],postprob[it][ij]);
     }
     count++;
     }
     fclose(fp);
     return 0;
}


/*    int filterodds(int totnoe,int npost[],float postprob[][NP])
    {
    int it,ij,ik;
    float odds;
    for(it=0;it<totnoe;it++) if(npost[it]>1)
    {
    for(ij=0;ij<npost[it];ij++)     for(ik=0;ik<npost[it];ik++) 
    {
    odds=postprob[it][ij]/postprob[it][ik];
    if(odds>2.0) 
    {
    postprob[it][ik]=0.0;
    goto fin;
    }
    }
    fin: 
    continue;
    }
    return 0;
    }  */


/*
int filterodds(int totnoe,int npost[],float postprob[][NP])
*/
static int filterodds(int totnoe,int npost[],float *postprob[])
{
    int it,ij,posi;
    float odds,maxi,mini;
    for(it=0;it<totnoe;it++) if(npost[it]>1)
    {
    mini=1.0;
    maxi=0.0;
    for(ij=0;ij<npost[it];ij++) 
    {
    if(postprob[it][ij]>maxi)  maxi=postprob[it][ij];
    if(postprob[it][ij]<mini)  
    {
    mini=postprob[it][ij];
    posi=ij;
    }
    }
    odds=maxi/mini; 
    if(odds>2.0) postprob[it][posi]=0.0;
    }
    return 0;
}


/*
int initpost(int npost[],float postprob[][NP],int posas1[][NP],int posas2[][NP])
*/
static void initpost(int npost[],float *postprob[],int *posas1[],int *posas2[])
{
    int it,ij;
    for(it=0;it<NPK;it++)
    {
        npost[it]=0;
        for(ij=0;ij<NP;ij++)
        {
            postprob[it][ij]=0.0;
            posas1[it][ij]=0;
            posas2[it][ij]=0;
        }
    }
}

/*
int zeroprobmap(float probmap[][NFR])
*/
static void zeroprobmap(float *probmap[])
{
    int it,ij;
    for(it=0;it<NFR;it++)
    {
        for(ij=0;ij<NFR;ij++)
            probmap[it][ij]=0.0;
    }
}

/*
int initprobmap(int totnoe,int npriors[],float priors[][NP],
        int prioas1[][NP],int prioas2[][NP],int npost[],float postprob[][NP],
        int postas1[][NP],int postas2[][NP],float probmap[][NFR])
*/
static void initprobmap(int totnoe,int npriors[],float *priors[],
        int *prioas1[],int *prioas2[],int npost[],float *postprob[],
        int *postas1[],int *postas2[],float *probmap[])
{
    int it,ij,idt1,idt2;
    for(it=0;it<totnoe;it++) if(npriors[it]>0)
    {
        for(ij=0;ij<npriors[it];ij++)
        {
            postas1[it][ij]=prioas1[it][ij];
            postas2[it][ij]=prioas2[it][ij];    
            postprob[it][ij]=priors[it][ij];
            idt1=prioas1[it][ij];
            idt2=prioas2[it][ij];
            probmap[idt1][idt2]=MAX(probmap[idt1][idt2],priors[it][ij]);
            probmap[idt2][idt1]=MAX(probmap[idt2][idt1],priors[it][ij]);
        }
        npost[it]=npriors[it];
    }
}


static int totalcount(int totnoe,int npriors[])
{
    int it,coun1;
    coun1=0;
    for(it=0;it<totnoe;it++) if(npriors[it]>0) coun1=npriors[it]+coun1;
/*
      printf(" count : %d  \n", *count);
*/
    return coun1;
}


/*
int includefake(int totnoe,int ntylnk[],int tylnk[][NP],int npost[],
        int postas1[][NP],int postas2[][NP])
*/
static int includefake(int totnoe,int ntylnk[],int *tylnk[],int npost[],
        int *postas1[],int *postas2[])
{
    int it,ij,idt1,idt2,flag1,flag2,count;
    for(it=0;it<totnoe;it++)  if(npost[it]==1)
    {
    idt1=postas1[it][0];
    idt2=postas2[it][0];
    flag1=0;
    flag2=0;
    for(ij=0;ij<ntylnk[idt1];ij++) if(tylnk[idt1][ij]==idt2) flag1=1;
    if(flag1==0) 
    {
    count=ntylnk[idt1];
    tylnk[idt1][count]=idt2;
    ntylnk[idt1]=ntylnk[idt1]+1;
    }

    for(ij=0;ij<ntylnk[idt2];ij++) if(tylnk[idt2][ij]==idt1) flag2=1;
    if(flag2==0) 
    {
    count=ntylnk[idt2];
    tylnk[idt2][count]=idt1;
    ntylnk[idt2]=ntylnk[idt2]+1;
    }
    }
    return 0;
}

/*
int shiftass(int totnoe,int npost[],float postprob[][NP],int postas1[][NP],int postas2[][NP],
        int nprior[],float priorprob[][NP],int prioras1[][NP],int prioras2[][NP])
*/
static int shiftass(int totnoe,int npost[], float *postprob[],int *postas1[],int *postas2[],
        int nprior[],float *priorprob[],int *prioras1[],int *prioras2[])
{
    int it,ik,count;
    for(it=0;it<totnoe;it++)  if(npost[it]>1)
    {
/*    for(ij=0;ij<npost[it];ij++) if(postprob[it][ij]==0.0) count++;
    if(count>1) 
    {    
    printf("there is an error in the shiftass \n");
    } */

    count=0;
    for(ik=0;ik<npost[it];ik++)    if(!(postprob[it][ik]==0.0))
    {
    postprob[it][count]=postprob[it][ik];
    postas1[it][count]=postas1[it][ik];
    postas2[it][count]=postas2[it][ik];

    priorprob[it][count]=priorprob[it][ik];
    prioras1[it][count]=prioras1[it][ik];
    prioras2[it][count]=prioras2[it][ik];
    count++;
    }
    npost[it]=count;
    nprior[it]=count;
    }
    return 0;
}



static void calcpost(int totnoe,int npost[],float *postprob[],float *likeprob[],
                     float *priorprob[])
{
    int it,ij;
    for(it=0;it<totnoe;it++)  if(npost[it]>1)
    {
        for(ij=0;ij<npost[it];ij++)
            postprob[it][ij]=priorprob[it][ij]*likeprob[it][ij];
    }
}


static CcpnStatus writenotfound(char *name,char *exten,int totnoe,int npost[],
        int *postas1[],int *postas2[],int noe[],float data[][3],char *error_msg)
{
    FILE *fp;
    int it;
    char outfile[20];
    strcpy(outfile,name);  
    strcat(outfile,exten);  
    printf(" writing the file: %s \n",outfile);
    if((fp=fopen(outfile,"w"))==NULL)
    {
        sprintf(error_msg, "writenotfound: cannot open file '%s' for writing", outfile);
        return CCPN_ERROR;
    }

    fprintf(fp,"frequencies not found \n");
    for(it=0;it<totnoe;it++) if(npost[it]==0) 
    {
        if(postas1[it][0]< 0) 
            fprintf(fp,"noe:%5d: frequencies not found:%f , %f \n ",noe[it],data[it][0],data[it][1]);
        if(postas2[it][0]< 0) 
            fprintf(fp,"noe:%5d: frequency:%f \n ",noe[it],data[it][2]);
    }
    fclose(fp);

    return CCPN_OK;
}     


static CcpnStatus do_bacus(char *name, int **cylnk, int **tylnk, int **prioras1, int **prioras2,
                int **postas1, int **postas2, float **priorprob,
                float **postprob, float **probmap, float **lkprob, char *error_msg)
{
/*************************** variables ****************************************/
    char atom[NFR][7],group[NFR][2],heavy[NFR][2],newtype[2],type[NPK][2];
    int DONE;
    int nfreqs,ngroup[NFR];
    int ncylnk[NFR],ntylnk[NFR];
    float freqs[NFR][2];
    int totnoe,newnum,newnoe[NPK],noe[NPK],newdiag[NPK],rownum[NPK];
    float newdata[NPK][3],data[NPK][3];
    float  newvol[NPK],vol[NPK];
    int nprior[NPK];
    int npost[NPK],count;
    float dif,shan1,shan2;

    printf(" \n");
    printf(" program BACUS (in c) Version 0 \n");
    printf(" written by Carlos A. Steren 2004 \n");
    printf(" Copyright: Miguel Llinas \n");
    printf(" Dept. of Chemistry \n Carnegie Mellon University \n Pittsburgh, PA, USA \n");
    printf(" \n");

/************************ reading data ****************************************/    
    
/*
    filename(name);
*/
    
    CHECK_STATUS(readnoes(name,"ncp.inp",&newnum,newdata,newnoe,newvol,error_msg));
    print_msg("readnoes");
    filterdiag(newnum,newdata,newdiag);
    print_msg("filterdiag");
    strcpy(newtype,"N");
    totnoe=0;
    CHECK_STATUS(addtolist(&totnoe,data,type,noe,rownum,vol,newnum,newnoe,newdata,newtype,newvol,newdiag,error_msg));
    print_msg("addtolist");

    CHECK_STATUS(readnoes(name,"ccp.inp",&newnum,newdata,newnoe,newvol,error_msg));
    print_msg("readnoes");
    filterdiag(newnum,newdata,newdiag);
    print_msg("filterdiag");
    strcpy(newtype,"C");
    printf("newtype: %s \n",newtype);
    CHECK_STATUS(addtolist(&totnoe,data,type,noe,rownum,vol,newnum,newnoe,newdata,newtype,newvol,newdiag,error_msg));
    print_msg("addtolist"); 
    printf("total number of noes = %d\n", totnoe);

    CHECK_STATUS(readfreqlst(name,"shi.inp",&nfreqs,freqs,heavy,ngroup,group,atom,error_msg));
    print_msg("readlst");

    CHECK_STATUS(readlinkfile(name,"clk.inp",nfreqs,ncylnk,cylnk, error_msg));
    print_msg("readlink");
    
    CHECK_STATUS(readlinkfile(name,"tlk.inp",nfreqs,ntylnk,tylnk, error_msg));
    print_msg("readlink"); 

/************************ first block ****************************************/

    CHECK_STATUS(calcpriors(totnoe,data,type,nfreqs,freqs,heavy,nprior,prioras1,prioras2,priorprob,error_msg));
    print_msg("calcpriors");

    CHECK_STATUS(writenotfound(name,"noo.out",totnoe,nprior,prioras1,prioras2,noe,data,error_msg));
    renorm(totnoe,nprior,priorprob);    
    print_msg("renorm");

    writepos(name,"prs.out",totnoe,nprior,priorprob,prioras1,prioras2,ngroup,group,atom);
    print_msg("writepriors");
    /*getch();*/
/************************ second block ****************************************/
        
     

    CHECK_STATUS(readpriors(name,"prs.out",&totnoe,noe,nprior,priorprob,prioras1,prioras2,error_msg));
    print_msg("readpriors");
    renorm(totnoe,nprior,priorprob);
    print_msg("renorm priors");

    DONE=CCPN_FALSE;
    while(!DONE)
    { 
        initpost(npost,postprob,postas1,postas2);
        print_msg("initpost");

        zeroprobmap(probmap);
        print_msg("zeroprobmap");

        initprobmap(totnoe,nprior,priorprob,prioras1,prioras2,npost,postprob,postas1,postas2,probmap);
        print_msg("initprobmap");

        count = totalcount(totnoe,nprior);
        print_msg("totalcount");
        printf("total assigments:%d \n",count);

        likehood(totnoe,ncylnk,cylnk,ntylnk,tylnk,nprior,priorprob,postas1,postas2,probmap,lkprob);
        print_msg("likehood");
        renorm(totnoe,npost,lkprob);        
        print_msg("renorm likehood");

        calcpost(totnoe,npost,postprob,lkprob,priorprob);
        print_msg("calcpost");
        renorm(totnoe,npost,postprob);
        print_msg("renorm posterios");

        shannon(totnoe,npost,postprob,&shan1);
        print_msg("shannon 1");
        filterodds(totnoe,npost,postprob);
        print_msg("filterodds");

        shiftass(totnoe,npost,postprob,postas1,postas2,nprior,priorprob,prioras1,prioras2);
        print_msg("shiftass");
        renorm(totnoe,npost,postprob);
        print_msg("renorm posteriors, second time");
        renorm(totnoe,nprior,priorprob);
        print_msg("renorm posteriors, second time");

        includefake(totnoe,ntylnk,tylnk,npost,postas1,postas2);
        print_msg("include fake tocsys");

        shannon(totnoe,npost,postprob,&shan2);
        print_msg("shannon 2");

        printf("shan1:%f, shan2:%f  \n",shan1,shan2);
        dif=fabs(shan1-shan2);
        if(dif<=1.0e-3)  DONE=CCPN_TRUE;
    }
    writeass(name,"ass.out",totnoe,rownum,npost,postprob,postas1,postas2,ngroup,group,atom);
    print_msg("writereass");
        /*
    writepos(name,"pos.out",totnoe,npost,postprob,postas1,postas2,ngroup,group,atom);
    print_msg("writerepost");
        */
    /*getch();*/

    return CCPN_OK;
}


CcpnStatus bacus(char *name, char *error_msg)
{
    int *cylnk[NFR],*tylnk[NFR];
    int *prioras1[NPK],*prioras2[NPK];
    int *postas1[NPK],*postas2[NPK];
    float *priorprob[NPK],*postprob[NPK],*probmap[NFR],*lkprob[NPK];
    CcpnStatus status;

    if (strlen(name) > 10)
    {
        sprintf(error_msg, "name cannot be more than length 10");
        return CCPN_ERROR;
    }

    CHECK_STATUS(alloc_memory(cylnk, tylnk, prioras1, prioras2, postas1, postas2, priorprob, postprob, probmap, lkprob));
 
    status = do_bacus(name, cylnk, tylnk, prioras1, prioras2, postas1, postas2, priorprob, postprob, probmap, lkprob, error_msg);

    free_memory(cylnk, tylnk, prioras1, prioras2, postas1, postas2, priorprob, postprob, probmap, lkprob);

    return status;
}

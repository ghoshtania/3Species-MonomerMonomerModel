//
//  MonomerMonomer.c
//  
//
//  Created by Tania Ghosh on 10/25/22.
//

#include <stdio.h>
#include <stdlib.h>      // standard library with lots of functions
#include <math.h>      // standard math library
#define NRANSI          // needed for NR
#include "my_nrutil.h"    // include NR header files

long L;
int *Lattice;
int *Lcopy;
int *VacList;   //pointer to integer for vacant sites
int VacNum;     //store the number for number of Vacant sites
double TotalVacNum; //count the totoal vac number per run
int *Size;
int *StillAlive;
double *LocalSlope;
double *Localeta;
double *Localz;
double *eta;
double *z;
int *iteration;
long seed;
double P_AB =0.5;
double P_C =0.121;
int tmax;
int MaxRun = 500000;

int TotalVacNum1;
float ran2(long *idum);      // typecast ran2
void initL();
void plotL();
void plotLists();
void initLists();
//int MCupdate(int);
void copy(int *Lattice, int *Lcopy);
void InterfaceSize(int *Lattice, int *Lcopy);
int Dynamics(int *Lcopy);
int CalSize(int *Lattice, int *Lcopy, int VacNum);
int main(int argc, char *argv[])
{
    FILE *fout;    // pointer to filename
    if (argc == 5)
    {
        L = atol(argv[1]);
        seed = atol(argv[2]);
        tmax = atol(argv[3]);
        fout = fopen(argv[4],"w");
    }
    else
    {
        fprintf(stderr,"\n Initialization error\n");
        fprintf(stderr,"Usage: XIEex.x L seed tmax Filename \n");  // correct input syntax
        return 1;
    }
    Lattice = ivector(-L,L);
    Lcopy = ivector(-L,L);
    Size = ivector(0,1);
    VacList = ivector(-L,L);
    StillAlive = ivector(0,tmax);
    LocalSlope = dvector(0,tmax);
    Localeta = dvector(0,tmax);
    Localz = dvector(0,tmax);
    eta = dvector(0,tmax);
    z = dvector(0,tmax);
    iteration = ivector(0,tmax);
    fprintf(stderr,"Initialization:\n");
    initL();
    //plotL();
    copy(Lattice, Lcopy);//copy list
    initLists();
    //plotL();
    //plotLists();
    double step;
    //tmax = 10;
    
    for (int i=0; i<tmax; i++)
    {
        StillAlive[i]=0; //initialize StillAlive vector with zero
        LocalSlope[i] = 0.0;
        Localeta[i] = 0.0;
        Localz[i] = 0.0;
        eta[i] = 0.0;
        iteration[i] = 0;
        z[i]=0.0;
    }
    //int MaxRun = 10;
    int LeftSize, RightSize, change, Size1, SqSize;
    for (int i=0; i<MaxRun; i++)
    {
        TotalVacNum =0.0;
        fprintf(stderr,"Runnumber %d\n",i);
        InterfaceSize(Lattice, Lcopy); //will note down where the changes happened - first element and last element
        initLists();
        double realtime = 0.0;
        int tt=0;
        LeftSize=0;
        RightSize=0;
        change= 0;
        while (realtime<tmax)
        {
            //fprintf(stderr,"start %d\n");
            change = Dynamics(Lcopy);       //target will return where the cange occured
//            fprintf(stderr,"change %d\n", change);
//            plotL();
//            plotLists();
            //fprintf(stderr,"Size %d SqSize %d\n", Size, SqSize);
            
            if (realtime>= tt)
            {
                StillAlive[tt] +=1;
                
                if (change<=LeftSize) LeftSize = change;
                if (change>=RightSize) RightSize = change;
                //fprintf(stderr,"LeftSize %d RightSize %d\n", LeftSize, RightSize);
                Size1 = ((RightSize-LeftSize)+2);
                SqSize = Size1*Size1;
                //eta[tt]+= (double)VacNum;
                z[tt]+= (double)SqSize;
                tt +=1;
            }
            if (VacNum!=0)
            {
                step = 1.0/((double)VacNum);
                realtime = realtime + step;
            }
           
            else
            {
                break;
            }
//
        }
    }

    double LocalScope = 0;
    for (int i=tmax; i>=5; i--)
        if (i % 5 ==0)
        {
            //fprintf(stderr,"StillAlive %d %d %f \n", i,StillAlive[i], (double)i/5.0);
            //z[i] = (z[i]/(double)StillAlive[i]);
            Localz[i] = log((double)z[i]/(double)z[i/5])/log(5);
            //fprintf(stderr,"Time %d SurvivalProb[t] %d SurvivalProb[t/5] %d LocalSlope %f\n",i, StillAlive[i], StillAlive[i/5], LocalSlope[i]);
            //fprintf(stderr,"actual t %d t/5 %d 1000/t %f\n", i,i/5, 1000.0/(double)i);
            fprintf(fout,"%f %f\n", (1000.0/(double)i), Localz[i]);
        }
//    for (int i=0; i<tmax; i++)
//        fprintf(stderr,"t %d eta %f \n", i, eta[i]);
    
    fclose(fout);
    
    free_ivector(Lattice,-L,L);
    free_ivector(Lcopy,-L,L);
    free_ivector(Size,0,1);
    free_ivector(VacList,-L,L);
    free_ivector(StillAlive,0,tmax);
    free_dvector(LocalSlope,0,tmax);
    free_dvector(Localeta,0,tmax);
    free_dvector(Localz,0,tmax);
    free_dvector(eta,0,tmax);
    free_dvector(z,0,tmax);
    free_ivector(iteration,0,tmax);
    return 0;
}


void initL()
{
   long i;
    for (i=-L; i<=L; ++i)   //initialize the lattice with C except the center which
    {Lattice[i] = 1;  }         //is vacant site
    Lattice[0] = 0;
    

    return;
}

void initLists()
{
    long i,j;
    VacNum = 0;
    for (i=-L; i<=L; ++i){
        if (Lcopy[i]==0)
        {
            VacList[VacNum] = i; //note the vacant site through VacList
            ++VacNum;        //increase the vacancy number by 1
        }
    }
    
    return;
}

void plotL()
{
    long i;
    for (i=-L; i<=L; ++i)
    {
        if (Lcopy[i] ==0) fprintf(stderr, "V");
        else if (Lcopy[i] ==1) fprintf(stderr,"A");
        else if (Lcopy[i] ==2) fprintf(stderr, "B");
        else fprintf (stderr,"C");
    }
    fprintf(stderr,"\n");
    return;
}

void plotLists()
{
    long i;
//    fprintf(stderr, "Total Vacancy %d\n", VacNum);
//    fprintf(stderr,"VacList: \n");
    for (i=0; i<VacNum;++i)
    {
        fprintf (stderr, "%d", VacList[i]);
    }
    fprintf(stderr,"\n");
    return;
}

void copy(int *Lattice, int *Lcopy)
{
    for (int i=-L; i<=L; i++) Lcopy[i]=Lattice[i];
}

void InterfaceSize(int *Lattice, int *Lcopy)
{
//    Size[0]=0;
//    Size[1]=0;
//    fprintf(stderr,"original\n");
//    for (int i=-L; i<=L; i++) fprintf(stderr,"%d",Lattice[i]);
//    fprintf(stderr,"\n");
//    fprintf(stderr,"beforecopy\n");
//    for (int i=-L; i<=L; i++) fprintf(stderr,"%d",Lcopy[i]);
//    fprintf(stderr,"\n");
    for (int i=-L; i<=L; i++)
    {
        if (Lcopy[i]!=Lattice[i])
        {
            Size[0]=i;
            break;
        }
    }
    for (int i=L; i>=-L; i--)
    {
        if (Lcopy[i]!=Lattice[i])
        {
            Size[1]=i;
            break;
        }
    }
      //fprintf(stderr,"change is %d %d\n", Size[0], Size[1]);
      for (int i=Size[0]; i<=Size[1]; i++) Lcopy[i]=Lattice[i];
//    fprintf(stderr,"aftercopy\n");
//    for (int i=-L; i<=L; i++) fprintf(stderr,"%d",Lcopy[i]);
//    fprintf(stderr,"\n");
    return ;
}

int CalSize(int *Lattice, int *Lcopy, int VacNum)
{while(VacNum>0)
    {for (int i=-L; i<=L; i++)
    {
        if (Lcopy[i]!=3)
        {
            Size[0]=i;
            break;
        }
    }
    for (int i=L; i>-L; i--)
    {
        if (Lcopy[i]!=3)
        {
            Size[1]=i;
            break;
        }
    }
//    fprintf(stderr,"Lcopy");
//    fprintf(stderr,"\n");
//    for (int i=-L; i<=L; i++)
//    {
//        fprintf(stderr,"%d", Lcopy[i]);
//    }
//    fprintf(stderr,"\n");
//    fprintf(stderr,"Lattice");
//    fprintf(stderr,"\n");
//    for (int i=-L; i<=L; i++)
//    {
//        fprintf(stderr,"%d", Lattice[i]);
//    }
//    fprintf(stderr,"\n");
//    fprintf(stderr,"Size[0] %d Size[1] %d\n", Size[0],Size[1]);
        int Result = Size[1]-Size[0]+1;
        return Result;
        }
}

int Dynamics(int *Lcopy)
{
    
    int i,j,targetSite,l,SurvivalNum, change;
    //fprintf(stderr,"Print the vac num %d\n", VacNum);
    if (VacNum==0)
    {
        fprintf(stderr,"It is in sarurated phase - no vaccancy\n");
    
    }else
    {
        j=VacNum*ran2(&seed);
        //fprintf(stderr,"j is %d\n" ,j);
        targetSite = VacList[j];
        //fprintf(stderr,"TargetSite is %d\n", targetSite);
        if (j<VacNum)
        {
            float roll1,roll2;
            roll1=ran2(&seed);
            //roll1 = 0.35;
            //fprintf(stderr,"first roll %f\n", roll1);
            if (roll1<P_C) //C will come
            {   //fprintf(stderr,"C will come\n");
                if ((Lcopy[targetSite-1]==2 || Lcopy[targetSite-1] ==1) && (Lcopy[targetSite+1] ==2 ||Lcopy[targetSite+1] ==1)) //AVB
                    {
                        //fprintf(stderr,"First-both are different :site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                        int side = 2*ran2(&seed);
                        //fprintf(stderr,"which side %d\n", side);
                        if (side==0)
                        {
                            Lcopy[targetSite-1] = 0;
                            ++VacNum;
                            VacList[VacNum-1] = targetSite-1;
                            change = targetSite-1;
                        }
                        else
                        {
                            Lcopy[targetSite+1] = 0;
                            ++VacNum;
                            VacList[VacNum-1] = targetSite+1;
                            change = targetSite+1;
                        }
                    }
                else if ((Lcopy[targetSite-1]==2 || Lcopy[targetSite-1]==1) && (Lcopy[targetSite+1]!=2 || Lcopy[targetSite+1]!=1))
                    {
                        //fprintf(stderr,"Second-left is different only:site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                        Lcopy[targetSite-1] = 0;
                        ++VacNum;
                        VacList[VacNum-1]=targetSite-1;
                        change = targetSite-1;
                    }
                else if ((Lcopy[targetSite+1]==2||Lcopy[targetSite+1]==1) && (Lcopy[targetSite-1]!=2||Lcopy[targetSite-1]!=1))
                    {
                        //fprintf(stderr,"Third-right is different only :site after and before target site is %d %d\n",Lattice[targetSite+1], Lattice[targetSite-1]);
                        Lcopy[targetSite+1] = 0;
                        ++VacNum;
                        VacList[VacNum-1]=targetSite+1;
                        change = targetSite+1;
                    }
                
                else
                    {
                        //fprintf(stderr,"Fourth-both are same or vacant:site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                        Lcopy[targetSite] = 3;
                        l = VacList[VacNum-1];  //last vacant site is l
                        VacList[j] = l; //update the filled site with the last vacant site
                        change = targetSite;
                        --VacNum;
                    }
            }
            else
            {
                roll2=ran2(&seed); //C will not come
                //fprintf(stderr,"second roll %f\n", roll2);
                if (roll2<P_AB)  //A will come
                {   //fprintf(stderr,"A will come\n");
                    if ((Lcopy[targetSite-1]==2 || Lcopy[targetSite-1] ==3) && (Lcopy[targetSite+1] ==2 ||Lcopy[targetSite+1] ==3)) //CVB
                        {
                            //fprintf(stderr,"First-both are different :site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                            int side = 2*ran2(&seed);
                            //fprintf(stderr,"which side %d\n", side);
                            if (side==0)
                            {
                                Lcopy[targetSite-1] = 0;
                                ++VacNum;
                                VacList[VacNum-1] = targetSite-1;
                                change = targetSite-1;
                            }
                            else
                            {
                                Lcopy[targetSite+1] = 0;
                                ++VacNum;
                                VacList[VacNum-1] = targetSite+1;
                                change = targetSite+1;
                            }
                        }
                    else if ((Lcopy[targetSite-1]==2 || Lcopy[targetSite-1]==3) && (Lcopy[targetSite+1]!=2 || Lcopy[targetSite+1]!=3))
                        {
                            //fprintf(stderr,"Second-left is different only:site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                            Lcopy[targetSite-1] = 0;
                            ++VacNum;
                            VacList[VacNum-1]=targetSite-1;
                            change = targetSite-1;
                        }
                    else if ((Lcopy[targetSite+1]==2||Lcopy[targetSite+1]==3) && (Lcopy[targetSite-1]!=2||Lcopy[targetSite-1]!=3))
                        {
                            //fprintf(stderr,"Third-right is different only :site after and before target site is %d %d\n",Lattice[targetSite+1], Lattice[targetSite-1]);
                            Lcopy[targetSite+1] = 0;
                            ++VacNum;
                            VacList[VacNum-1]=targetSite+1;
                            change = targetSite+1;
                        }
                    
                    else
                        {
                            //fprintf(stderr,"Fourth-both are same or vacant:site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                            Lcopy[targetSite] = 1;
                            l = VacList[VacNum-1];  //last vacant site is l
                            VacList[j] = l; //update the filled site with the last vacant site
                            --VacNum;
                            change = targetSite;
                        }
                }
                else   //B will come
                {   //fprintf(stderr,"B will come\n");
                    if ((Lcopy[targetSite-1]==3 || Lcopy[targetSite-1] ==1) && (Lcopy[targetSite+1] ==3 ||Lcopy[targetSite+1] ==1))
                        {
                            //fprintf(stderr,"First-both are different :site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                            int side = 2*ran2(&seed);
                            //fprintf(stderr,"which side %d\n", side);
                            if (side==0)
                            {
                                Lcopy[targetSite-1] = 0;
                                ++VacNum;
                                VacList[VacNum-1] = targetSite-1;
                                change = targetSite-1;
                            }
                            else
                            {
                                Lcopy[targetSite+1] = 0;
                                ++VacNum;
                                VacList[VacNum-1] = targetSite+1;
                                change = targetSite+1;
                            }
                        }
                    else if ((Lcopy[targetSite-1]==3 || Lcopy[targetSite-1]==1) && (Lcopy[targetSite+1]!=3 || Lcopy[targetSite+1]!=1))
                        {
                            //fprintf(stderr,"Second-left is different only:site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                            Lcopy[targetSite-1] = 0;
                            ++VacNum;
                            VacList[VacNum-1]=targetSite-1;
                            change = targetSite-1;
                        }
                    else if ((Lcopy[targetSite+1]==3||Lcopy[targetSite+1]==1) && (Lcopy[targetSite-1]!=3||Lcopy[targetSite-1]!=1))
                        {
                            //fprintf(stderr,"Third-right is different only :site after and before target site is %d %d\n",Lattice[targetSite+1], Lattice[targetSite-1]);
                            Lcopy[targetSite+1] = 0;
                            ++VacNum;
                            VacList[VacNum-1]=targetSite+1;
                            change = targetSite+1;
                        }
                    
                    else
                        {
                            //fprintf(stderr,"Forth-both are same or vacant:site before and after target site is %d %d\n",Lattice[targetSite-1], Lattice[targetSite+1]);
                            Lcopy[targetSite] = 2;
                            l = VacList[VacNum-1];  //last vacant site is l
                            VacList[j] = l; //update the filled site with the last vacant site
                            --VacNum;
                            change = targetSite;
                        }
                }
            }
        
        }
   
}
    
    return change;
}




// below is a NR random number generator. It generated float numbers evenly over range [0,1)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software *1(.|a. */

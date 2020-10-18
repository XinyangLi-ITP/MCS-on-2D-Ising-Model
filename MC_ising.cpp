/*
Program to calculate some physical quantities in ferromagnetic phase transition of 2D ising model.
By Xinyang Li,ITP-CAS,Oct 2020
*/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <time.h>
using namespace std;

#define Spin(i,j) spin[i+j*L]
#define PI 3.141592653589793238462643
//Initialize spin state.
void init(int L,int *spin,int &H){

 for(int i=0;i<L*L;i++){
       spin[i]=1;
   }

 int i,j;
 int ip,im,jp,jm;
 int H_single;
 for(i=0;i<L;i++){
    for(j=0;j<L;j++){
       if (i==L-1) ip = 0;
       else
          ip = i+1;
       if (j==L-1) jp = 0;
       else 
          jp = j+1;
          H_single=-Spin(i,j)*(Spin(ip,j)+Spin(i,jp));
          H=H+H_single;
    }
  }
}


//Pick a spin randomly,calculate the energy of interacting with nearest neighbor spins.Flip it according to probability. 
void flip(int L,int *spin,double &m,double T,int &H,int (*x)[2]){

 int i,j;
 int ip,im,jp,jm;
 double id,jd;	
 int H_old;
 int del_H;

 id=((double) rand())/RAND_MAX*L;
 jd=((double) rand())/RAND_MAX*L;

 i=((int) id)%L;
 j=((int) jd)%L;

 if (i==L-1)   ip = 0;
 else
    ip = i+1;
 if (j==L-1)   jp = 0;
 else
    jp =j+1;
 
 if (i==0)   im =L-1;
 else
    im = i-1;
 if (j==0)   jm =L-1;
 else
    jm = j-1;

 H_old=-Spin(i,j)*(Spin(ip,j)+Spin(im,j)+Spin(i,jp)+Spin(i,jm));
 del_H=-2*H_old;
  
 if (del_H<=0){
   Spin(i,j)=-Spin(i,j);
   H=H+del_H;
 }

 else{
   double boltz=exp(2*H_old/T);
   double rn;
   rn=((double) rand())/(double) RAND_MAX;
     if (rn<=boltz){  
        Spin(i,j)=-Spin(i,j);
        H=H+del_H;
      }
  }
 for(int i=0;i<L*L;i++){
   m+=spin[i];
  }
  m=fabs(m/(L*L));

 for (int i=0;i<L;i++){
    for (int j=0;j<L;j++){
          x[i+j*L][0]=i;
          x[i+j*L][1]=j;
    }
  }

}

//Calculate average energy and heat capacity at fixed T.
 double heat_capacity(int step,double *E,double &energy,double &Cv,int L,double T){
   int i;
   double square_H=0;
   for(i=0.5*step;i<step;i+=100){
     energy+=E[i];
     square_H+=pow(E[i],2);
    }
   energy=energy/(0.5*step/100);
   square_H=square_H/(0.5*step/100);
   Cv=(square_H-pow(energy,2))/(L*L*T*T);
}

//Calculate average magnetization,magnetic susceptibility and binder parameter at fixed T.
double magnet(int step,double *M,double &mag,double &chi,int L,double T,double &binder_parameter){
   int i;
   double square_m;
   double fourth_m;
   for(i=0.5*step;i<step;i+=100){
       mag+=M[i];
       square_m+=pow(M[i],2);
       fourth_m+=pow(M[i],4);
    }
   square_m=square_m/(0.5*step/100);
   mag=mag/(0.5*step/100);
   fourth_m=fourth_m/(0.5*step/100);
   chi=L*L*(square_m-pow(mag,2))/T;
   binder_parameter=1-(fourth_m/(3*pow(square_m,2)));
}

//Calculate spacial correlation function at fixed T.
double correlation(int L,int *spin,int (*x)[2],int nhis,int ig,double *g,int *SUM){
  double delx;
  double dely;
  double delr;
  int sum=0;
   for (int p=0;p<L*L;p++){
     for (int q=0;q<L*L;q++){
       delx=x[p][0]-x[q][0];
           while(delx>L/2)
             {delx=delx-L;}
           while(delx<0-L/2)
             {delx=delx+L;}
       dely=x[p][1]-x[q][1];
           while(dely>L/2)
             {dely=dely-L;}
           while(dely<0-L/2)
             {dely=dely+L;}

       delr=sqrt(pow(delx,2)+pow(dely,2));
       if (delr<=L/2){
         ig=round(delr);
         g[ig]=g[ig]+spin[p]*spin[q];
         SUM[ig]=SUM[ig]+1;
       }
     }
   }
}
 
//Average spacial correlation function.
double Gr(int L,int count,int nhis,int *SUM,double *g,double T,double mag){
    double r;
    char namebuffer6[256];
    sprintf(namebuffer6, "L%d_grco_T%f.dat",L,T);
    ofstream grco(namebuffer6);
    char namebuffer7[256];
    sprintf(namebuffer7, "L%d_grdisc_T%f.dat",L,T);
    ofstream grdisc(namebuffer7);
    for(int i=0;i<=nhis;i++){
        g[i]=(g[i]/count)/(L*L);
        SUM[i]=(SUM[i]/count)/(L*L);
      }
    for(int i=0;i<nhis;i++){
        r=i;
        g[i]=g[i]/(SUM[i]);
        grdisc<<r<<"  "<<g[i]<<endl;
        g[i]=g[i]-pow(mag,2);
        grco<<r<<"  "<<g[i]<<endl;
      }
   grdisc.close();
   grco.close();
 }

//Calculate autocorrelation function at fixed T.
double Autocorrelation(int L,double T,double *M,int step){
    char namebuffer8[256];
    sprintf(namebuffer8, "L%d_Autoco_T%f_tw30000.dat",L,T);
    ofstream Autoco(namebuffer8);
    double chi_0=0;
    double *chi_t;
    chi_t=new double[600]; 
    for(int j=0;j<600;j++){
        int t=0;
        int t_prime=0;
        int t_sum=0;
        double term1=0;
        double term2_1=0;
        double term2_2=0;
         t=j*L*L;
      for(int k=30000;k<step-t;k+=1){
        t_prime=k;
        t_sum=(step-t-30000)/1;
        term1+=M[t_prime]*M[t_prime+t];
        term2_1+=M[t_prime];
        term2_2+=M[t_prime+t];
      }
        term1=term1/t_sum;
        term2_1=term2_1/t_sum;
        term2_2=term2_2/t_sum;
        chi_t[j]=term1-term2_1*term2_2;
        chi_0=chi_t[0];
    }
   for(int j=0;j<600;j++){  
        chi_t[j]=chi_t[j]/chi_0;
        Autoco<<j<<" "<<chi_t[j]<<endl;
    }
   Autoco.close(); 
}

int main(int argc,char *argv[]){
    int L=(int) atof(argv[1]);
  //double T=(double) atof(argv[2]);
    char namebuffer1[256];
    char namebuffer2[256];
    char namebuffer3[256];
    char namebuffer4[256];
    char namebuffer5[256];
    
    sprintf(namebuffer1, "L%d_m_T.dat", L);
    ofstream m_T(namebuffer1);
    sprintf(namebuffer2, "L%d_chi_T.dat", L);
    ofstream chi_T(namebuffer2);
    sprintf(namebuffer3,"L%d_H_T.dat",L);
    ofstream H_T(namebuffer3);
    sprintf(namebuffer4,"L%d_Cv_T.dat",L);
    ofstream Cv_T(namebuffer4);
    sprintf(namebuffer5,"L%d_binder_T.dat",L);
    ofstream binder_T(namebuffer5);

    srand((int) time(NULL));
  
    int step=L*L*10*6000;
    int *spin;
    spin= new int[L*L];
    double m;
    double *M;
    M=new double[step];
    double *E;
    E=new double[step];
    double mag=0;
    double energy=0;
    double chi=0;
    int H=0;
    double Cv=0;
    double binder_parameter=0;

    int (*x)[2];
    x=new int[L*L][2];

    init(L,spin,H);
    //double T=2.40;
    for(double T=6.0;T>0.0;T-=0.2){
     
      double count=0;
      int nhis=0.5*L;
      int ig;
      double (*g);
      g=new double[nhis+1];
      int (*SUM);
      SUM=new int[nhis+1];
      for (int i=0;i<=nhis;i++){
        g[i]=0;
        SUM[i]=0;
      }
      for(int i=0;i<step;i++){
        m=0;   
        flip(L,spin,m,T,H,x);
        M[i]=m;
        E[i]=H;
        if(i>0.5*step&i%(200*L*L)==0){
          correlation(L,spin,x,nhis,ig,g,SUM);
          count+=1;
        }
      }
     
    magnet(step,M,mag,chi,L,T,binder_parameter);
    heat_capacity(step,E,energy,Cv,L,T);
    Gr(L,count,nhis,SUM,g,T,mag);
    Autocorrelation(L,T,M,step);
 
    H_T<<T<<"  "<<energy<<endl;
    Cv_T<<T<<"  "<<Cv<<endl;
    m_T<<T<<"  "<<mag<<endl;
    chi_T<<T<<"  "<<chi<<endl;
    binder_T<<T<<"  "<<binder_parameter<<endl;
   
    delete [] g;
   }
    delete [] spin;
    delete [] M;
    delete [] E;
    delete [] x;

    m_T.close();
    chi_T.close();
    H_T.close();
    Cv_T.close();
    binder_T.close();
    return 0;
}

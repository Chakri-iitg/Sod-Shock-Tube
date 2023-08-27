#include <iostream>
#include<fstream>
#include<cmath>
#include<bits/stdc++.h>
#define ga 1.4


using namespace std;

double max_eig(double val){    //Function for finding max eigen value

     if(val>=0){

              return val;
          }
     else{

              return 0;

 }

}


 double min_eig(double val){   //Function for finding Min eigen value

     if(val<0){

               return val;



     }
     else{

               return 0;
}

}

void eigen(int i,double **u,double *eig){     //Function for finding eigen values

      double P,U,A;

             P=((ga-1)*(u[2][i]-((pow(u[1][i],2))/(2*u[0][i]))));
             U=(u[1][i]/u[0][i]);
             A=sqrt(((ga*P)/u[0][i]));


             eig[0]=U-A;
             eig[1]=U;
             eig[2]=U+A;

}

void flux_plus(int i,double **f_1,double **u,double *eig){
                                                               //Function for finding F+
      double P,U,A,H,a;

            P=((ga-1)*(u[2][i]-((pow(u[1][i],2))/(2*u[0][i]))));
            U=(u[1][i]/u[0][i]);
            A=sqrt(((ga*P)/u[0][i]));
            H=((u[2][i]+P)/u[0][i]);
            a=((u[0][i])/(2*ga));

            f_1[0][i]=a*(max_eig(eig[0])+((2*(ga-1))*max_eig(eig[1]))+max_eig(eig[2]));
            f_1[1][i]=a*((max_eig(eig[0])*(U-A))+((2*(ga-1))*max_eig(eig[1]))+(max_eig(eig[2])*(U+A)));
            f_1[2][i]=a*((max_eig(eig[0])*(H-(U*A)))+((ga-1)*(U*U)*max_eig(eig[1]))+((H+(U*A))*max_eig(eig[2])));


}

void flux_minus(int i,double **f_1,double **u,double *eig){        //Function for finding F-

double P,U,A,H,a;

           P=((ga-1)*(u[2][i]-((pow(u[1][i],2))/(2*u[0][i]))));
           U=(u[1][i]/u[0][i]);
           A=sqrt(((ga*P)/u[0][i]));
           H=((u[2][i]+P)/u[0][i]);
           a=((u[0][i])/(2*ga));


           f_1[0][i]=a*(min_eig(eig[0])+((2*(ga-1))*min_eig(eig[1]))+min_eig(eig[2]));
           f_1[1][i]=a*((min_eig(eig[0])*(U-A))+((2*(ga-1))*min_eig(eig[1]))+(min_eig(eig[2])*(U+A)));
           f_1[2][i]=a*((min_eig(eig[0])*(H-(U*A)))+((ga-1)*(U*U)*min_eig(eig[1]))+((H+(U*A))*min_eig(eig[2])));


}

void max_lambda(double **u,double &max_val,int N){               //Function for finding max eigen value in entire domain
        double *lam;
        double P,U,A,temp;
        lam=new double[3];

         for(int i=0;i<N-1;i++){
                   P=((ga-1)*(u[2][i]-((pow(u[1][i],2))/(2*u[0][i]))));
                   U=(u[1][i]/u[0][i]);
                   A=sqrt(((ga*P)/u[0][i]));


                   lam[0]=U-A;
                   lam[1]=U;
                   lam[2]=U+A;

                     temp=0;
                         for(int j=0;j<3;j++)
                        {
                              if(fabs(lam[j]) > temp)
                              temp= lam[j];

                        }

            if(i==0)
         {
                max_val=temp;
            }

            if(temp>max_val)
         {
               max_val=temp;
            }


   }

   }





int main()
{
    ifstream in("input.txt");
    ofstream ou_1("density.txt");
    ofstream ou_3("pressure.txt");
    ofstream ou_2("velocity.txt");
    ofstream ou_4("energy.txt");
    double p_l,p_r,u_l,u_r,d_l,d_r;

   in>>d_l>>p_l>>u_l>>d_r>>p_r>>u_r;
    double ro_e_l,ro_e_r,max_lam;
    double **u,**f,*eig;
    double t=0,d_x,d_t,L=1,T=0.35;
    int N=400;
    d_x=L/N;
   eig=new double[3];
   u=new double *[3];
   f=new double *[3];

   for(int i=0;i<3;i++){

    u[i]=new double[N+1];
    f[i]=new double[N+1];


   }
    ro_e_l=((d_l*pow(u_l,2)*0.5)+(p_l/(ga-1)));
    ro_e_r=((d_r*pow(u_r,2)*0.5)+(p_r/(ga-1)));

   for(int i=0;i<N;i++){
          if(i<N/2){
               u[0][i]=d_l;
               u[1][i]=d_l*u_l;     //Calculating the u by using given initial conditions
               u[2][i]=ro_e_l;


    }
          else{
               u[0][i]=d_r;
               u[1][i]=d_r*u_r;
               u[2][i]=ro_e_r;

    }

   }


   double **f_pl,**f_mi,**f_1;
   f_pl=new double*[3];
   f_mi=new double*[3];
   f_1=new double*[3];

   for(int i=0;i<3;i++){

             f_pl[i]=new double[N+1];
             f_mi[i]=new double[N+1];
             f_1[i]=new double[N+1];

    }
    u[0][0]=u[0][1];
    u[1][0]=-u[1][1];
    u[2][0]=u[2][1];

    u[0][N]=u[0][N-1];
    u[1][N]=-u[1][N-1];
    u[2][N]=u[2][N-1];

   while(t<T){      //Time loop

    for(int i=0;i<N;i++){  //Space loop

              eigen(i,u,eig);

              flux_plus(i,f_pl,u,eig);

              eigen(i+1,u,eig);

              flux_minus(i+1,f_mi,u,eig);

                f_1[0][i]=f_pl[0][i]+f_mi[0][i+1];
                f_1[1][i]=f_pl[1][i]+f_mi[1][i+1];
                f_1[2][i]=f_pl[2][i]+f_mi[2][i+1];

                 max_lambda(u,max_lam,N);

                if(t>T)
            {
                 d_t=T-t;
            }

                else
             {


                 d_t=(0.8*d_x)/max_lam;
             }
               if(i!=0){

                       u[0][i]=u[0][i]-(((d_t/d_x))*(f_1[0][i]-f_1[0][i-1]));
                       u[1][i]=u[1][i]-(((d_t/d_x))*(f_1[1][i]-f_1[1][i-1]));
                       u[2][i]=u[2][i]-(((d_t/d_x))*(f_1[2][i]-f_1[2][i-1]));

        }
    }
     u[0][0]=u[0][1];
    u[1][0]=-u[1][1];
    u[2][0]=u[2][1];

    u[0][N]=u[0][N-1];
    u[1][N]=-u[1][N-1];
    u[2][N]=u[2][N-1];


   t+=d_t;

}


double x,p,U,e;
   for(int i=0;i<N;i++)
{           x=((2*i)+1)*(d_x*0.5);
            U=(u[1][i]/u[0][i]);
            p=((ga-1)*(u[2][i]-((pow(u[1][i],2))/(2*u[0][i]))));
            e=(p/((ga-1)*u[0][i]));

            ou_1<<x<<" "<<u[0][i]<<endl;
            ou_2<<x<<" "<<U<<endl;
            ou_3<<x<<" "<<p<<endl;
            ou_4<<x<<" "<<e<<endl;
  }


























    return 0;
}

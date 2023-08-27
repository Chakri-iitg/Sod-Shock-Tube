#include <iostream>
#include<fstream>
#include<cmath>
#include<bits/stdc++.h>
#define ga 1.4


using namespace std;

double fun_tilda(int n,int i,double **u){
                                               //Function for calculating the Tilda values of U,H,A
  double val,a,b,c;
  double h_l,h_r,p_l,p_r;



  if(n==1){                                                              //u_tilda
          a=sqrt((u[0][i+1]))*(u[1][i+1]/u[0][i+1]);
          b=sqrt((u[0][i]))*(u[1][i]/u[0][i]);
          c=sqrt(u[0][i+1])+sqrt(u[0][i]);

          val=((a+b)/c);

  }
  else if(n==2){                                                         //h_tilda
          p_l=(u[2][i]-(pow(u[1][i],2)/(2*u[0][i])))*(ga-1) ;
          p_r=(u[2][i+1]-(pow(u[1][i+1],2)/(2*u[0][i+1])))*(ga-1) ;


          h_l=(u[2][i]+p_l)/u[0][i];
          h_r=(u[2][i+1]+p_r)/u[0][i+1];

          a=sqrt(u[0][i+1])*h_r;
          b=sqrt(u[0][i])*h_l;
          c=sqrt(u[0][i+1])+sqrt(u[0][i]);
          val=((a+b)/c);
  }
  else if(n==3){
                                                                          //a_tilda

            a=fun_tilda(2,i,u);
            b=fun_tilda(1,i,u);

            val=(ga-1)*(a-(pow(b,2)*0.5));

  }

  return val;

}

void flux(int i,double **f,double **u){
                                               //Function for calculating the flux at i th cell

            f[0][i]=u[1][i] ;
            f[1][i]=(((ga-1)*u[2][i])+(((3-ga)*0.5)*(pow(u[1][i],2)/u[0][i])));
            f[2][i]=((ga*u[2][i]*u[1][i])/u[0][i])-(((ga-1)*0.5)*(pow(u[1][i],3)/pow(u[0][i],2))) ;

}

void eigen(double *eig,double **k,double a,double u,double h){

          eig[0]=u-a;
          eig[1]=u;                                                       //Function for calculating the Eigen values
          eig[2]=u+a;

          k[0][0]=k[0][1]=k[0][2]=1;

          k[1][0]=u-a;
          k[1][1]=u;
          k[1][2]=u+a;
          k[2][0]=h-(u*a);
          k[2][1]=pow(u,2)/2;
          k[2][2]=h+(u*a);

}

void Delta(int i,double *del,double **u1,double h,double u,double a){
                                                                          //Function for calculating the Delta values
 double x,y,z;
           x=u1[0][i+1]-u1[0][i];
           y=u1[1][i+1]-u1[1][i];
           z=u1[2][i+1]-u1[2][i];
           del[1]=(((ga-1)/pow(a,2))*(((h-pow(u,2))*x)+(u*y)-z)) ;

           del[0]=((1/(2*a))*((x*(u+a))-y-(a*del[1])));
           del[2]=x-del[0]-del[1];

 }

 void flux_1(int i,double *eig,double *delta,double **k,double **f_1){


   for(int j=0;j<3;j++)
        {                                                       //Function for calculating the flux at i+half cell
             if(eig[j]<0)
            {
                    f_1[0][i]=f_1[0][i]+(eig[j]*delta[j]*k[0][j]);
                    f_1[1][i]=f_1[1][i]+(eig[j]*delta[j]*k[1][j]);
                    f_1[2][i]=f_1[2][i]+(eig[j]*delta[j]*k[2][j]);
            }
        }






 }
 void max_lambda(double **u,double &max_val,int N){
double *lam;
double temp,u1,a;
lam=new double[3];

for(int i=0;i<N-1;i++){                               //This function is used to check the max eigen value in the entire domain
            u1=fun_tilda(1,i,u);
            a=fun_tilda(3,i,u);



             lam[0]=u1-a;
             lam[1]=u1;
             lam[2]=u1+a;

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




void entropy_fix(double *eig){

  double err=1e-6;                               //Entropy fix

  for(int i=0;i<3;i++){

    if(fabs(eig[i])<err){
            eig[i]=pow(eig[i],2)/err ;
            eig[i]=0.5*(eig[i]+err);

             eig[i]=fabs(eig[i]);

    }

    }






}

int main()
{
    int q;

   cout<<"Enter the number: \n 1-> Without entropy fix \n 2-> With entropy fix\n";
   cin>>q;

    ifstream in("input.txt");  //created the input text file

    ofstream ou_1("density_1.dat");
    ofstream ou_2("velocity_1.dat");
    ofstream ou_3("pressure_1.dat");
    ofstream ou_4("energy_1.dat");


    ofstream ou_5("density_2.dat");
    ofstream ou_6("velocity_2.dat");
    ofstream ou_7("pressure_2.dat");
    ofstream ou_8("energy_2.dat");

    double d_l,u_l,p_l,d_r,u_r,p_r;
    double u_tl,h_tl,a_tl,ro_e_l,ro_e_r;
    in>>d_l>>p_l>>u_l>>d_r>>p_r>>u_r;

   double t,d_x,d_t,L,max_eig,T;
    int N;
   double **u,**f,**f_1;
   double *eig,**k;
   double *delta;
   N=400;
   L=1;
   T=0.15;
   d_x=L/N;
   u=new double *[3];
   f=new double *[3];
   f_1=new double*[3];
   for(int i=0;i<3;i++){

    u[i]=new double[N+1];
    f[i]=new double[N+1];
    f_1[i]=new double[N+1];

   }
   eig=new double[3];
   k=new double*[3];
   for(int i=0;i<3;i++){

    k[i]=new double[3];

   }
   delta=new double[3];
    ro_e_l=((d_l*pow(u_l,2)*0.5)+(p_l/(ga-1)));
   ro_e_r=((d_r*pow(u_r,2)*0.5)+(p_r/(ga-1)));


   for(int i=0;i<N;i++){
    if(i<N/2){               //Calculating the intial u by using the given intial values
           u[0][i]=d_l;
           u[1][i]=d_l*u_l;
           u[2][i]=ro_e_l;


    }
    else {
           u[0][i]=d_r;
           u[1][i]=d_r*u_r;
           u[2][i]=ro_e_r;

    }

}


  t=0;



  while(t<T){       //Loop for time step

    for(int i=0;i<N-1;i++){  //Loop for space

               flux(i,f,u) ;


               u_tl=fun_tilda(1,i,u);

               h_tl=fun_tilda(2,i,u);

               a_tl=fun_tilda(3,i,u);

                eigen(eig,k,a_tl,u_tl,h_tl);

         if(q==2){
                 entropy_fix(eig);
         }

                Delta(i,delta,u,h_tl,u_tl,a_tl);

               f_1[0][i]=f[0][i];
               f_1[1][i]=f[1][i];
               f_1[2][i]=f[2][i];

               flux_1(i,eig,delta,k,f_1);


                max_lambda(u,max_eig,N);


              if(t>T)
              {
                    d_t=T-t;
              }

             else
              {


                    d_t=(0.8*d_x)/max_eig;
       }
               if(i!=0){

                         u[0][i]=u[0][i]-(((d_t/d_x))*(f_1[0][i]-f_1[0][i-1]));
                         u[1][i]=u[1][i]-(((d_t/d_x))*(f_1[1][i]-f_1[1][i-1]));
                         u[2][i]=u[2][i]-(((d_t/d_x))*(f_1[2][i]-f_1[2][i-1]));

        }






    }




                 t+=d_t;




  }


       double x,p,U,e;
    /*for(int i=0;i<N;i++)
{           x=((2*i)+1)*(d_x*0.5);
            U=(u[1][i]/u[0][i]);
            p=((ga-1)*(u[2][i]-((pow(u[1][i],2))/(2*u[0][i]))));
            e=(p/((ga-1)*u[0][i]));
       if(q==1){
            ou_1<<x<<" "<<u[0][i]<<endl;
            ou_2<<x<<" "<<U<<endl;
            ou_3<<x<<" "<<p<<endl;                   //Printing the values in Text file
            ou_4<<x<<" "<<e<<endl;
       }
      else if(q==2){
           ou_5<<x<<" "<<u[0][i]<<endl;
           ou_6<<x<<" "<<U<<endl;
           ou_7<<x<<" "<<p<<endl;
           ou_8<<x<<" "<<e<<endl;



      }






 }
 */




    return 0;
}

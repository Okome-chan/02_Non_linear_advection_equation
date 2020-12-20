
/*****************************************************************************/
//
//  burgers.cc
//  
//  ver 0.0.1   2020/12/14
//  developed by Yoneda
//
/*****************************************************************************/

#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>
#include <math.h>
#include <time.h>
#include <unistd.h>

using namespace std;

int main(){

  int typ;          //simulation type (1: Backward, 2: MacCormack)
  int t;
  int tmax;         //time step
  int collumn;      //number of collumn
  double dt;        //dt
  double dx;        //dx
  double r;         //Courant number
  double *u;        //u at time tn
  double *u_new;    //u at time tn+1
  FILE *fp;

  /**********************************/
  //Define values
  /**********************************/
  r=0.2;
  tmax=200;
  collumn=50;
  dx=4.*M_PI/double(tmax-1);
  dt=dx*r;

  
  /**********************************/
  //Setting array
  /**********************************/
  u=(double *)malloc(sizeof(double)*collumn+1);
  u_new=(double *)malloc(sizeof(double)*collumn+1);

  /**********************************/
  //Define subroutine
  /**********************************/
  void initial_array(double *u,double dx,int collumn);
  void Backward(double *u_new, double *u, double dx, double dt, int imax);
  void Maccormack(double *u_new, double *u, double dx, double dt, int imax);
  void update_value(double *u_new, double *u, int imax);
  void make_graph(double *u,double collumn,double dx,int t,FILE *fp);

  /**********************************/
  //set initial value
  /**********************************/
  initial_array(u,dx,collumn);

  /**********************************/
  //select simulation type
  /**********************************/
  cout<<"1: Backward, 2: MacCormack"<<endl;
  cout<<"Please enter simulation type"<<endl;
  cin>>typ;

  /**********************************/
  // calculate Courant number
  /**********************************/
  t=0;
  make_graph(u,collumn,dx,t,fp);
  for(t=1; t<=tmax; t++){
    if(typ==1){
      Backward(u_new,u,dx,dt,collumn);
      update_value(u_new,u,collumn);
    }
    if(typ==2){
      Maccormack(u_new,u,dx,dt,collumn);
    }
    if(t%5==0) make_graph(u,collumn,dx,t,fp);
  }

  return 0;
}


//Initiallization of array
void initial_array(double *u,double dx,int collumn){
  int i;
  for(i=1;i<=collumn;i++){
    u[i]=0.;
  }
  for(i=11;i<=20;i++){
   u[i]=sin(dx*(i-1.));
  }
  return;
}

//Backward difference scheme
void Backward(double *u_new, double *u, double dx,double dt, int imax){
  int i;
  for(i=2; i<=imax-1; i++){
    if(u[i]>=0.) u_new[i]=u[i]-dt/dx*(u[i]-u[i-1]);
    if(u[i]<0.)  u_new[i]=u[i]-dt/dx*(u[i+1]-u[i]);
  }
  //boundary condition
  u_new[imax]=u_new[2];
  u_new[1]=u_new[imax-1];

  return;
}

//Mac Cormack scheme
void Maccormack(double *u_new, double *u, double dx, double dt, int imax){
  int i;
  for(i=1; i<=imax-1; i++){
    u_new[i]=u[i]-(u[i+1]*u[i+1]-u[i]*u[i])*dt/(2.*dx);
cout<<i<<"\t"<<u[i]<<"\t"<<u_new[i]<<endl;
  }
  u_new[imax]=u_new[2];
  for(i=2; i<=imax; i++){
    u[i]=0.5*(u[i]+u_new[i])-0.5*(u_new[i]*u_new[i]-u_new[i-1]*u_new[i-1])*dt/(2.*dx);
cout<<i<<"\t"<<u[i]<<"\t"<<u_new[i]<<endl;
  }
  u[1]=u[imax-1];
//exit(1);

  return;
}



//Update simulated values
void update_value(double *u_new, double *u, int imax){
  int i;
  for(i=1;i<=imax;i++){
    u[i]=u_new[i];
  }
}

//Make graph//
void make_graph(double *u,double collumn,double dx,int t,FILE *fp)
{
  int i;
  double x;
  char data_file[256];
  FILE *data;

  fp=popen("gnuplot -persist","w");
  fprintf(fp,"set terminal png size 400,400\n");
  fprintf(fp,"set output './output/%04d.jpg'\n",t);
  fprintf(fp,"set xrange [0.0:3.2]\n");
  fprintf(fp,"set yrange [-0.2:1.2]\n");
  //sleep(1);
  sprintf(data_file,"out.dat");
  data=fopen(data_file,"w");
  for(i=1;i<=collumn;i++){
    x=(double)i*dx;
    fprintf(data,"%f\t%f\n",x,u[i]);
  }
  fprintf(fp,"plot \"%s\" using 1:2 with lines title \"sim\"\n",data_file);
  fclose(data);
  fflush(fp);
  pclose(fp);

  return;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "definicion_funciones.h"

//-----------------------------------------------------//
//             Funciones auxiliares                    //
//-----------------------------------------------------//

float dy1dx(float x,float y1,float y2,float y3){
  float ro=28.0,sigma=10.0,beta=8.0/3.0;
  float retorno=sigma*(y2-y1);
  return retorno;
}

float dy2dx(float x,float y1,float y2,float y3){
  float ro=28.0,sigma=10.0,beta=8.0/3.0;
  float retorno=y1*(ro-y3)-y2;
  return retorno;
}

float dy3dx(float x,float y1,float y2,float y3){
  float ro=28.0,sigma=10.0,beta=8.0/3.0;
  float retorno=y1*y2-beta*y3;
  return retorno;
}

float *RungeKuttaFourthOrderStep(float x_old,float y1_old,float y2_old,float y3_old,float h){
  
  float *retorno=NULL;
  float y_prime_1=0.0;
  float y_prime_2=0.0;
  float x_middle=0.0;
  float k1_y1=0.0,k2_y1=0.0,k3_y1=0.0,k4_y1=0.0;
  float k1_y2=0.0,k2_y2=0.0,k3_y2=0.0,k4_y2=0.0;
  float k1_y3=0.0,k2_y3=0.0,k3_y3=0.0,k4_y3=0.0;
  float slope=0.0,x1=0.0,x2=0.0,x3=0.0,x4=0.0;
  float y1_step1=0.0,y1_step2=0.0,y1_step3=0.0,slope_y1=0.0;
  float y2_step1=0.0,y2_step2=0.0,y2_step3=0.0,slope_y2=0.0;
  float y3_step1=0.0,y3_step2=0.0,y3_step3=0.0,slope_y3=0.0;
  float x_new=0.0,y1_new=0.0,y2_new=0.0,y3_new=0.0;
  
  retorno=generatevector_empty(retorno,4);
  
  k1_y1=dy1dx(x_old,y1_old,y2_old,y3_old);
  k1_y2=dy2dx(x_old,y1_old,y2_old,y3_old);
  k1_y3=dy3dx(x_old,y1_old,y2_old,y3_old);
  
  //Paso No.1//
  x1=x_old+(h/2.0);
  
  y1_step1=y1_old + (h/2.0) * k1_y1;
  y2_step1=y2_old + (h/2.0) * k1_y2;
  y3_step1=y3_old + (h/2.0) * k1_y3;

  k2_y1=dy1dx(x1,y1_step1,y2_step1,y3_step1);
  k2_y2=dy2dx(x1,y1_step1,y2_step1,y3_step1);
  k2_y3=dy3dx(x1,y1_step1,y2_step1,y3_step1);
  
  //Paso No.2//
  x2=x_old+(h/2.0);
  
  y1_step2=y1_old + (h/2.0) * k2_y1;
  y2_step2=y2_old + (h/2.0) * k2_y2;
  y3_step2=y3_old + (h/2.0) * k2_y3;
  
  k3_y1=dy1dx(x2,y1_step2,y2_step2,y3_step2);
  k3_y2=dy2dx(x2,y1_step2,y2_step2,y3_step2);
  k3_y3=dy3dx(x2,y1_step2,y2_step2,y3_step2);
  
  //Paso No.3//
  x3=x_old+(h/1.0);
  
  y1_step3=y1_old + (h/1.0) * k3_y1;
  y2_step3=y2_old + (h/1.0) * k3_y2;
  y3_step3=y3_old + (h/1.0) * k3_y3;
  
  k4_y1=dy1dx(x3,y1_step3,y2_step3,y3_step3);
  k4_y2=dy2dx(x3,y1_step3,y2_step3,y3_step3);
  k4_y3=dy3dx(x3,y1_step3,y2_step3,y3_step3);
  
  //Paso No.4//
  slope_y1=(1.0/6.0)*(k1_y1+2.0*k2_y1+2.0*k3_y1+k4_y1);
  slope_y2=(1.0/6.0)*(k1_y2+2.0*k2_y2+2.0*k3_y2+k4_y2);
  slope_y3=(1.0/6.0)*(k1_y3+2.0*k2_y3+2.0*k3_y3+k4_y3);
  
  x_new=x_old+h;
  y1_new=y1_old+h*slope_y1;
  y2_new=y2_old+h*slope_y2;
  y3_new=y3_old+h*slope_y3;
  
  retorno[0]=x_new;
  retorno[1]=y1_new;
  retorno[2]=y2_new; 
  retorno[3]=y3_new;
  
  return retorno;
}


float *generatevector_empty(float *vector,int n){
  int i=0;
  
  if(!(vector = malloc(sizeof(float)*n))){
    fprintf(stderr, "Problem with memory allocation");
  }
  
  for(i=0;i<n;i++){
    vector[i]=0.0;
  }
  
  return vector;
}
void intToChar(int j,char indice[]){
  if(j==0){strcat(indice,"0");}
  if(j==1){strcat(indice,"1");}
  if(j==2){strcat(indice,"2");}
  if(j==3){strcat(indice,"3");}
  if(j==4){strcat(indice,"4");}
  if(j==5){strcat(indice,"5");}
  if(j==6){strcat(indice,"6");}
  if(j==7){strcat(indice,"7");}
  if(j==8){strcat(indice,"8");}
  if(j==9){strcat(indice,"9");}
}

void escribirxy_txt(float *x,float *y1,float *y2,float *y3,int n_points,int j){
  FILE *fileout;
  int i=0;
  char filename[26]="RungeKuttaEdoDataSolution";
  char indice[2]="";
  char format[5]=".txt";
  
  intToChar(j,indice);
  
  strcat(filename,indice);
  strcat(filename,format);
  
  fileout=fopen(filename,"w");
  for(i=0;i<(n_points);i++){
    // Usando el formato t x y z
    fprintf(fileout,"%f %f %f %f\n",x[i],y1[i],y2[i],y3[i]);
  }
}

void imprimir_en_bash(float *x,float *y1,float *y2,float *y3,int n_points){
  int i=0;
  for(i=0;i<(n_points);i++){
    printf(" t[%d]=%f x[%d]=%f y[%d]=%f y[%d]=%f \n",i,x[i],i,y1[i],i,y2[i],i,y3[i]);
  }
}
void  graficarGNUPLOT(int j,int x0,int y0,int z0){

 char filename1[50]="";
 char format1[17]=".txt' using 2:3\n";
 
 char filename2[50]="";
 char format2[17]=".txt' using 2:4\n";
 
 char filename3[50]="";
 char format3[17]=".txt' using 3:4\n";
 
 char filename4[50]="";
 char format4[19]=".txt' using 2:3:4\n";
 char indice[2]="";
 

 intToChar(j,indice);

 strcat(filename1,"plot 'RungeKuttaEdoDataSolution");
 strcat(filename1,indice);
 strcat(filename1,format1);
 
 strcat(filename2,"plot 'RungeKuttaEdoDataSolution");
 strcat(filename2,indice);
 strcat(filename2,format2);
 
 strcat(filename3,"plot 'RungeKuttaEdoDataSolution");
 strcat(filename3,indice);
 strcat(filename3,format3);
 
 strcat(filename4,"splot 'RungeKuttaEdoDataSolution");
 strcat(filename4,indice);
 strcat(filename4,format4);
 
  FILE *gplot = popen("gnuplot -persist","w");
  fprintf(gplot, "set term png\n");
  fprintf(gplot, "set output 'It_%d_%d%d%d_plotTarea6.png'\n",j,x0,y0,z0);
  fprintf(gplot, "set multiplot layout 2,2 rowsfirst \n");
  
  fprintf(gplot, "set title 'Plano x-y;x0=%d, y0=%d ,z0=%d '\n",x0,y0,z0);
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "%s",filename1);
  
  fprintf(gplot, "set title 'Plano x-z;x0=%d, y0=%d ,z0=%d '\n",x0,y0,z0);
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "%s",filename2);
  
  fprintf(gplot, "set title 'Plano y-z; x0=%d, y0=%d ,z0=%d '\n",x0,y0,z0);
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "%s",filename3);
  
  fprintf(gplot, "set title 'Grafica 3D;x0=%d, y0=%d ,z0=%d '\n",x0,y0,z0);
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "%s",filename4);
  fprintf(gplot, "unset multiplot\n");

  close(gplot);
  
}

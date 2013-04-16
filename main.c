/*
Fisica Computacional EDO Runge Kutta cuarto 0rden
Nombre:Jorge Luis Mayorga-Ivan Dario Torroledo
Codigo:20111082-201116966
Parte 1/2
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "definicion_funciones.h"

int main(int argc, char **argv)
{
  /*Declaracion de variables e inicializacion nula*/
  int n_points=0,i=0;
  float*x=NULL,*y1=NULL,*y2=NULL,*y3=NULL,h=0.0,a=0.0,b=0.0;
  float y1_0=0.0,y2_0=0.0,y3_0=0.0,y1_prima_temp=0.0,y2_prima_temp=0.0,y3_prima_temp=0.0;
  float *temp=NULL,n_points_float=0.0;
  float x_old=0.0,y1_old=0.0,y2_old=0.0,y3_old=0.0;
  int N=-10,M=10,j=0,limite_superior=0,limite_inferior=0;
  /*Limites del Random*/
  limite_inferior=-10;
  limite_superior=10;
  /*Inicia el For de los 10 numeros aleatorios*/
  for(j=0;j<10;j++){
    h=0.03;
    y1_0=(rand() % limite_superior) + limite_inferior;
    y2_0=(rand() % limite_superior) + limite_inferior;
    y3_0=(rand() % limite_superior) + limite_inferior;
    a=0.0;
    b=30.0;
    n_points_float=(b-a)/h;
    n_points=(int)n_points_float;
    
    //Generar vectores vacios y mallocs
    x=generatevector_empty(x,n_points);
    temp=generatevector_empty(temp,4);
    y1=generatevector_empty(y1,n_points);
    y2=generatevector_empty(y2,n_points);
    y3=generatevector_empty(y3,n_points);
 
    //Inicio condiciones Iniciales a las variables de fase
    y1[0]=y1_0;
    y2[0]=y2_0;
    y3[0]=y3_0;
    x[0]=a;
    for(i=1;i<n_points;i++){
      
      y1_prima_temp=dy1dx(x[i-1],y1[i-1],y2[i-1],y3[i-1]);
      y2_prima_temp=dy2dx(x[i-1],y1[i-1],y2[i-1],y3[i-1]);
      y3_prima_temp=dy3dx(x[i-1],y1[i-1],y2[i-1],y3[i-1]);
      
      x_old=x[i-1];
      y1_old=y1[i-1];
      y2_old=y2[i-1];
      y3_old=y3[i-1];
      
      temp=RungeKuttaFourthOrderStep(x_old,y1_old,y2_old,y3_old,h);
      
      x[i]=temp[0];
      y1[i]=temp[1];
      y2[i]=temp[2];
      y3[i]=temp[3];
      
    }
    escribirxy_txt(x,y1,y2,y3,n_points,j);
    graficarGNUPLOT(j,y1_0,y2_0,y3_0);
  }
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//-----------------------------------------------------//
// Declaracion de funciones auxiliares                 //
//-----------------------------------------------------//

float dy1dx(float x,float y1,float y2,float y3);
float dy2dx(float x,float y1,float y2,float y3);
float dy3dx(float x,float y1,float y2,float y3);
float *RungeKuttaFourthOrderStep(float x_old,float y1_old,float y2_old,float y3_old,float h);
float *generatevector_empty(float *vector,int n);
void imprimir_en_bash(float *x,float *y1,float *y2,float *y3,int n_points);
void  escribirxy_txt(float *x,float *y1,float *y2,float *y3,int n_points,int j);
void  graficarGNUPLOT(int j,int x0,int y0,int z0);
void intToChar(int j,char indice[]);

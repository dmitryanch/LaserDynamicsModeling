//#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
double sigma, gam, r, delta, A, k, dt;
//const double dt = 0.005;

void add (double a_m[], double b_m[], double c_m[], int N)
{
 //int N = sizeof(a_m)/sizeof(a_m[0]);
 for (int i = 0; i < N; i++)
    {
     c_m[i] = a_m[i] + b_m[i];
     //c_m[2 * i + 1] = a_m[2 * i + 1] + b_m[2 * i + 1];
    }
}

void differ (double a_m[], double b_m[], double c_m[], int N)
{

 //int N = sizeof(a_m)/sizeof(a_m[0]);
 for (int i = 0; i < N; i++)
    {
     c_m[i] = a_m[i] - b_m[i];
     //c_m[2 * i + 1] = a_m[2 * i + 1] - b_m[2 * i + 1];
    }
}

void mult_c_c (double a_m[], double b_m[], double c_m[], int N)
{

 for (int i = 0; i < N; i++)
    {
     c_m[2 * i] = a_m[2 * i] * b_m[2 * i] - a_m[2 * i + 1] * b_m[2 * i + 1];

     c_m[2 * i + 1] = a_m[2 * i] * b_m[2 * i + 1] + a_m[2 * i + 1] * b_m[2 * i];
    }
}

void div_c_r (double a_m[], double b, double c_m[], int N)
{
 for (int i = 0; i < N; i++) c_m[i] = a_m[i] / b;
}

void mult_c_r (double a_m[], double b, double c_m[], int N)
{
 for (int i = 0; i < N; i++) c_m[i] = a_m[i] * b;
}

void abs_c (double a[], double abs_a[], int N)
{
 for (int i = 0; i < N; i++) abs_a[i] = sqrt( a[2 * i] * a[2 * i] + a[2 * i +1] * a[2 * i + 1]);
}

double min_m (double a[], int N)
{
 double min_a = a[0];
 for (int i = 1; i < N; i++)
    if (a[i] < min_a) min_a = a[i];
 return min_a;
}

double max_m (double a[], int N)
{
 double max_a = a[0];
 for (int i = 1; i < N; i++)
    if (a[i] > max_a) max_a = a[i];
 return max_a;
}

void matr_mul_vec (double A[][5], double b[], double c[], int N)
{
 for (int i = 0; i < N; i++)
    {
     c[i] = 0.0;
     for (int j = 0; j < N; j++) c[i]+=A[i][j] * b[j];
    }
}

// Комплексное сопряжение

void conj(double E[], double Econj[], int N)
{
for(int i = 0; i < N;i++)
	{
	Econj[2 * i]=E[ 2 * i];
	Econj[2 * i +1]=-E[2 * i + 1];
	}
}

// Норма вектора

double norma(double E[], int N)
{
double x=0;
for(int i = 0; i < N; i++) x = x + E[i] * E[i];
x=sqrt(x);
return x;
}

// diff-функция

void diff(int k[], int diffK[], int N)
{
for(int i=0;i<N-1;i++) diffK[i]=k[i+1]-k[i];
}

// Матрица линеаризации

void matlin(double E1,double E2, double P1,double P2, double D, double X[][5])
{
 extern double sigma, gam, delta, A, k;
 X[0][0]=-sigma;
 X[1][0]=-A*k*k;
 X[2][0]=D;
 X[3][0]=0;
 X[4][0]=-gam*P1;
 X[0][1]=A*k*k;
 X[1][1]=-sigma;
 X[2][1]=0;
 X[3][1]=D;
 X[4][1]=-gam*P2;
 X[0][2]=sigma;
 X[1][2]=0;
 X[2][2]=-1;
 X[3][2]=-delta;
 X[4][2]=-gam*E1;
 X[0][3]=0;
 X[1][3]=sigma;
 X[2][3]=delta;
 X[3][3]=-1;
 X[4][3]=-gam*E2;
 X[0][4]=0;
 X[1][4]=0;
 X[2][4]=E1;
 X[3][4]=E2;
 X[4][4]=-gam;
}

void next (double A[][5], double B[], double X[])
{
extern double dt;
//void mult_c_r (double a_m[], double b, double c_m[], int N);
//void matr_mul_vec (double A[][5], double b[], double c[], int N);
//void add (double a_m[], double b_m[], double c_m[], int N);
//void div_c_r (double a_m[], double b, double c_m[], int N);

double A1[5], A2[5], A3[5], A4[5], Bt[5];

matr_mul_vec (A, B, A1, 5) ;
mult_c_r (A1, dt, A1, 5) ;

div_c_r (A1,(double)2,Bt,5);
add (Bt,B,Bt, 5);
matr_mul_vec (A, Bt, A2, 5) ;
mult_c_r (A2, dt, A2, 5) ;

div_c_r (A2,(double)2,Bt,5);
add (Bt,B,Bt, 5);
matr_mul_vec (A, Bt, A3, 5);
mult_c_r (A3, dt, A3, 5) ;

add (A3,B,Bt, 5);
matr_mul_vec (A, Bt, A4, 5);
mult_c_r (A4, dt, A4, 5) ;

mult_c_r (A2, (double)2, A2, 5);
mult_c_r (A3, (double)2, A3, 5);
add (A1,A2,X, 5);
add (A3,X,X, 5);
add (A4,X,X, 5);
div_c_r (X,(double)6,X,5);
add (B,X,X, 5);
}

void f1 (double E[], double P[], double Y[])
{
 extern double sigma;
 // void differ (double a_m[], double b_m[], double c_m[], int N);
 // void mult_c_r (double a_m[], double b, double c_m[], int N);

 double PE[2];

 differ (P, E, PE, 2);
 mult_c_r (PE, sigma, Y, 2);
}

void f2 (double E[], double P[], double D, double Y[])
{
 extern double delta;
 // void mult_c_r (double a_m[], double b, double c_m[], int N);
 // void mult_c_c (double a_m[], double b_m[], double c_m[], int N);
 // void add (double a_m[], double b_m[], double c_m[], int N);

 double DE[2], deltaP[2], idelta[2];

 idelta[0] = -1.0;
 idelta[1] = -delta;

 mult_c_c (idelta, P, deltaP, 1);

 mult_c_r (E, D, DE, 2);

 add (deltaP, DE, Y, 2);
}

double f3 (double E[], double P[], double D)
{
 extern double gam, r;
 // void conj(double E[], double Econj[], int N);
 // void mult_c_c (double a_m[], double b_m[], double c_m[], int N);
 // void add (double a_m[], double b_m[], double c_m[], int N);

// Y=-gamma*(D-r+0.5*(conj(E)*P+E*conj(P)));
 double Econj[2], Pconj[2], conjEP[2], EconjP[2], sumEP[2], Y;

 conj (E, Econj, 1);
 conj (P, Pconj, 1);
 mult_c_c (Econj, P, conjEP, 1);
 mult_c_c (E, Pconj, EconjP, 1);
 add (conjEP, EconjP, sumEP, 2);

 Y = -gam * (D - r + 0.5 * sumEP[0]);
 return Y;
}

void rk (double E[], double P[], double &D)
{
 extern double dt;
 // void f1 (double E[], double P[], double Y[]);
 // void f2 (double E[], double P[], double D, double Y[]);
 // double f3 (double E[], double P[], double D);

 // void add (double a_m[], double b_m[], double c_m[], int N);
 // void div_c_r (double a_m[], double b, double c_m[], int N);
 // void mult_c_r (double a_m[], double b, double c_m[], int N);

 double a1[2], b1[2], c1, a2[2], b2[2], c2, a3[2], b3[2], c3, a4[2], b4[2], c4;
 double Et[2], Pt[2];

 f1 (E, P, a1);
 mult_c_r (a1, dt, a1, 2);
 f2 (E, P, D, b1);
 mult_c_r (b1, dt, b1, 2);
 c1 = dt * f3 (E, P, D);

 div_c_r (a1, (double)2, Et, 2);
 add (E, Et, Et, 2);
 div_c_r (b1, (double)2, Pt, 2);
 add (P, Pt, Pt, 2);

 f1 (Et, Pt, a2);
 mult_c_r (a2, dt, a2, 2);
 f2 (Et, Pt, D + c1/(double)2, b2);
 mult_c_r (b2, dt, b2, 2);
 c2 = dt * f3 (Et, Pt, D + c1/(double)2);

 div_c_r (a2, (double)2, Et, 2);
 add (E, Et, Et, 2);
 div_c_r (b2, (double)2, Pt, 2);
 add (P, Pt, Pt, 2);

 f1 (Et, Pt, a3);
 mult_c_r (a3, dt, a3, 2);
 f2 (Et, Pt, D + c2/(double)2, b3);
 mult_c_r (b3, dt, b3, 2);
 c3 = dt * f3 (Et, Pt, D + c2/(double)2);

 add (E, a3, Et, 2);
 add (P, b3, Pt, 2);

 f1 (Et, Pt, a4);
 mult_c_r (a4, dt, a4, 2);
 f2 (Et, Pt, D + c3, b4);
 mult_c_r (b4, dt, b4, 2);
 c4 = dt * f3 (Et, Pt, D + c3);

 mult_c_r (a2, (double)2, a2, 2);
 mult_c_r (a3, (double)2, a3, 2);

 add (a1, a2, Et, 2);
 add (a3, Et, Et, 2);
 add (a4, Et, Et, 2);
 div_c_r (Et, (double)6, Et, 2);
 add (E, Et, E, 2);

 mult_c_r (b2, (double)2, b2, 2);
 mult_c_r (b3, (double)2, b3, 2);

 add (b1, b2, Pt, 2);
 add (b3, Pt, Pt, 2);
 add (b4, Pt, Pt, 2);
 div_c_r (Pt, (double)6, Pt, 2);
 add (P, Pt, P, 2);

 D=D+(c1+2*c2+2*c3+c4)/6;
}



 int main ()
 {

	 ofstream f_out_y ("out.txt"); 	// Поток для записи результата
	 //ofstream f_out_e ("outE.txt"); 	// Поток для записи результата
	 //ofstream f_out_sp ("outSP.txt"); 	// Поток для записи результата
	 
	 // void rk (double a[], double b[], double &);
	 // void abs_c (double a[], double abs_a[], int N);
	 // double min_m (double a[], int N);
	 // double max_m (double a[], int N);
	 // void matlin(double E1,double E2, double P1,double P2, double D, double X[][5]);
	 // void next (double A[][5], double B[], double X[]);
	 // double norma(double E[], int N);
	 // void div_c_r (double a_m[], double b, double c_m[], int N);
	 // void diff(int k[], int diffK[], int N);

	 sigma = 0.1;
	 gam = 0.001;
	 A = 0.01;
	 dt = 0.5;
	 double r0 = 5.0;
	 delta = -0.25;

	 double t = 0.0;

	 const int Nk = 100, Nt = 3e5;
	 const double dk = 0.1;
	 double K [Nk];
	 double a[2],b[2],c,SP[Nk],Y[5],Y1[5],Yn, Q[5][5], min_E, max_E;
	 
	 int  k_positive[Nk], k_zone, k_max, nk, diffK[Nk];
	 int Ndel=10e5;
	 int i1;
	 double wRel = sqrt(2 * gam*sigma*(r0 - 1 - delta *delta) / (1 + delta * delta));
	 for (int i = 0; i < Nk; i++) K[i] = i * dk;

	 const int NM = 10, NW = 10;
	 double M[NM], W[NW];
	 /*double XM[NM*NW], YW[NM*NW], max_SPs[NM*NW];
	 bool K0[NM*NW];*/
	 /*int nextIndex = 0;
	 for (int i = 0; i < NM*NW; i++)
	 {
		 XM[i] = -1;
		 YW[i] = -1;
	 }*/

	 for (int i = 0; i < NM; i++) M[i] = 1.0 / NM*(i + 1);
	 for (int i = 0; i < NW; i++) W[i] = 0.5 + 2.0 / (NW - 1)*i;
	
	 double *E = new double [2 * Nt];
	 double *P = new double [2 * Nt];
	 double *D = new double [Nt];
	 double *abs_E = new double [Nt];

	 E[0] = 0.1;
	 E[1] = 0;
	 P[0] = 0;
	 P[1] = 0;
	 D[0] = 0;

	 a[0] = E[0];
	 a[1] = E[1];
	 b[0] = P[0];
	 b[1] = P[1];
	 c = D[0];

	 // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 for (int i = 0; i < NM; i++)
	 {
		 double m = M[i];
		 for (int j = 0; j < NW; j++)
		 {
			 double w = W[j] * wRel;

			 a[0] = 0.1; a[1] = 0;
			 b[0] = 0; b[1] = 0;
			 c = 0;
			 for (i1 = 0; i1 < (Nt + Ndel); i1++)
			 {
				 r = r0*(1. + m*sin(w*(double)i1*dt));
				 rk(a, b, c);
				 if (i1 >= Ndel)
				 {
					 E[2 * (i1 - Ndel)] = a[0];
					 E[2 * (i1 - Ndel) + 1] = a[1];
					 P[2 * (i1 - Ndel)] = b[0];
					 P[2 * (i1 - Ndel) + 1] = b[1];
					 D[i1 - Ndel] = c;
				 }
			 }
			 abs_c(E, abs_E, Nt);
			 min_E = min_m(abs_E, Nt);
			 max_E = max_m(abs_E, Nt);
			 /*if (max_E-min_E>0.1)
				y[j+NR*i]=2;*/
			 for (int ik = 0; ik < Nk; ik++) SP[ik] = 0.0;
			 for (int ik = 0; ik < Nk; ik++)
			 {
				 for (int iy = 1; iy < 5; iy++) Y[iy] = 0.0;
				 Y[0] = 1.0;
				 k = K[ik];
				 for (int iE = 0; iE < Nt; iE++)
				 {
					 matlin(E[2 * iE], E[2 * iE + 1], P[2 * iE], P[2 * iE + 1], D[iE], Q);
					 next(Q, Y, Y1);
					 for (int iy = 0; iy < 5; iy++) Y[iy] = Y1[iy];
					 Yn = norma(Y, 5);
					 SP[ik] = SP[ik] + log(Yn);
					 div_c_r(Y, Yn, Y, 5);
				 }
				 SP[ik] = SP[ik] / (Nt*dt);
			 }
			 // Обработка спектра старшего ляпуновского показателя
			 double max_SP = max_m(SP, Nk);
			 //if (max_SP > 1e-10)
			 {
				/* XM[nextIndex] = m;
				 YW[nextIndex] = w;
				 max_SPs[nextIndex] = max_SP;
				 K0[nextIndex] = max_SP == SP[0];*/

				 f_out_y << m << ' ' << w / wRel << ' ' << max_SP << ' ' << SP[0] << endl; 	//вывод в файл
				 
				 //nextIndex++;
			 }
		 }
	 }
	 
	 //for (int i = 0; i < nextIndex; i++)
	 //{
		// f_out_y << XM[i] << ' ' << YW[i] / wRel << ' ' << max_SPs[i] << K0[i] << endl; 	//вывод в файл
	 //}
	
	/*for (int i = 0; i < Nt; i++)
		f_out_e << sqrt(E[2*i]*E[2*i] + E[2*i+1]*E[2*i+1]) << endl;*/

	/*for (int i = 0; i < Nk; i++)
		f_out_sp << SP[i] << endl;*/

	return 0;


}
//---------------------------------------------------------------------------

#include "nrutil.c"
#define NRANSI
//#include "C:\libreria\nrutil.h"
#include "nrutil.h"
//void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
//	void (*derivs)(float, float [], float [])) ;

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double []))
{
	int i;
	double xh,hh,h6,*dym,*dyt,*yt;

	dym=dvector(0,n);
	dyt=dvector(0,n);
	yt=dvector(0,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	
	(*derivs)(x,y,dydx);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_dvector(yt,0,n);
	free_dvector(dyt,0,n);
	free_dvector(dym,0,n);
}
#undef NRANSI

	
void rk4d0(double y[], double dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double []))
{
	int i;
	double xh,hh,h6,*dym,*dyt,*yt;

	dym=dvector(0,n);
	dyt=dvector(0,n);
	yt=dvector(0,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_dvector(yt,0,n);
	free_dvector(dyt,0,n);
	free_dvector(dym,0,n);
}
#undef NRANSI


void rk4d1(double y[], double dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double []))
{
	int i;
	double xh,hh,h6,*dym,*dyt,*yt;

	dym=dvector(1,n);
	dyt=dvector(1,n);
	yt=dvector(1,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=1;i<=n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_dvector(yt,1,n);
	free_dvector(dyt,1,n);
	free_dvector(dym,1,n);
}
#undef NRANSI






void rk4f0(float y[], float dydx[], int n, float x, float h, float yout[],
	void (*derivs)(float, float [], float []))
{
	int i;
	float xh,hh,h6,*dym,*dyt,*yt;

	dym=vector(0,n);
	dyt=vector(0,n);
	yt=vector(0,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_vector(yt,0,n);
	free_vector(dyt,0,n);
	free_vector(dym,0,n);
}
#undef NRANSI

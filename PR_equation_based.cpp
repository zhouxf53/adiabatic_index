// Copyright [2010] <Xunfei Zhou>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
using namespace std;
double pc = 4.24766e6;
double Tc = 3 69.825;
double vc = 0.00454545455;
double M = 44.097;  // g/mol
double Rm = 8.3142;     // j/mol/K
double omega = 0.1524;

double vT_p(double v, double T);
double pT_v(double p, double T);
double pv_T(double p, double v);
double kappa(double p, double v, double T);


void main(){
	double p,v,T,n,s;
	T = 350;	
	p = 2e6;
	v = pT_v(p,T);
	n = kappa(p,v,T);	
	
	cout<<"P="<<p<<endl;
	cout<<"v="<<v<<endl;
	cout<<"T="<<T<<endl;
	
	cout<<"n="<<n<<endl;	
	cin>>s;
}


double vT_p(double v,double T){
	double a,b,alpha,k,p;
	double R = Rm/M*1000;
	k = 0.37464+1.54226*omega-0.26992*omega*omega;
	alpha = (1.0+k*(1.0-pow(T/Tc,0.5)))*(1.0+k*(1.0-pow(T/Tc,0.5)));
	a = 0.45724*alpha*R*R*Tc*Tc/pc;
	b = 0.07780*R*Tc/pc;
	p = R*T/(v-b)-a/(v*v+2*b*v-b*b);
	return p;
}

double pT_v(double p,double T){
	double v = 0;
	double v0 = 0.00000000000000001,v1=1;
	double vmid;
	double fv0, fvmid;
	while(fabs(v1 - v0) > 1.0e-8){
		vmid = (v0 + v1) / 2.0;
		fv0 = vT_p(v0, T)-p;
		fvmid = vT_p(vmid,T)-p;
		if(fv0 * fvmid <0)
			v1 = vmid;
		else
			v0 = vmid;
	}
	 v = v0;
	return v;
}

double pv_T(double p,double v){
	double T = 0;
	double T0 = 250,T1 = 360;
	double Tmid;
	double fT0, fTmid;
	while(fabs(T1-T0)>1.0e-8){
		Tmid = (T0+T1) / 2.0;
		fT0 = vT_p(v,T0)  - p;
		fTmid = vT_p(v,Tmid) - p;
		if(fT0 * fTmid < 0)
			T1 = Tmid;
		else
			T0 = Tmid;
	}
	 T = T0;
	
	return T;
}

double kappa(double p,double v,double T){
	//double R=8.31451;
	double Z;
    double R=Rm/M*1000;
	double Tr=T/Tc,pr=p/pc;
	double cp_R=3.20946+5.40615*Tr+1.89922*Tr*Tr;	
	double a,b,A,B,alpha,k,E;
	k=0.37464+1.54226*omega-0.26992*omega*omega;	
	alpha=pow((1.0+k*(1.0-pow(Tr,0.5))),2);
	a=0.45724*alpha*R*R*Tc*Tc/pc;	
	b=0.07780*R*Tc/pc;
	A=a*p/R/R/T/T;
	B=b*p/R/T;
	Z=p*v/R/T;
	E=k*sqrt(Tr/alpha);

	double n,m;
	double zt,zp;
	double pdzdp,tdzdt;
	pdzdp=(2*A*B-2*B*B-3*B*B*B-Z*Z*B-Z*(A-6*B*B-2*B))/(3*Z*Z-2*(1-B)*Z+A-3*B*B-2*B);
	tdzdt=(B*Z*Z+(-6*B*B-2*B+(E+2)*A)*Z-(E+3)*A*B+2*B*B+3*B*B*B)/(3*Z*Z-2*(1-B)*Z+A-3*B*B-2*B);
	zp=Z-pdzdp;
	zt=Z+tdzdt;

	double modify,modify1,modify2,modify3;
    modify1=-A*b/2/B/alpha;
	modify2=k*(k+1)*sqrt(Tr);
    modify3=log((v+b-sqrt(2.0)*b)/(v+b+sqrt(2.0)*b))/2/b/sqrt(2.0);
	modify=modify1*modify2*modify3;
   cp_R=cp_R-(1-zt*zt/zp)-modify;

	m=zt/cp_R;
	n=Z/(zp-zt*zt/cp_R);                            

	return n;
}

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void function(double* k,double* f,double a,double b,double c)
{
	k[0]=a*(f[1]-f[0]);
        k[1]=f[0]*(b-f[2])-f[1];
        k[2]=f[0]*f[1]-c*f[2];
}


void rk4(double*,int,int,double,const double,const double,const double);

int main()
{
	const double a=10.0;
	const double b=28.0;
	const double c=8.0/3.0;
	const double tmax=100.0;
	double dt;

	cout << "Enter a stepsize dt<1(should also be a fraction of 100): " << endl;
	cin >> dt;
	
	int n=tmax/dt;
	int d=3;
	double f[3];
	f[0]=1;f[1]=1;f[2]=1;

	rk4(f,d,n,dt,a,b,c);
	
}

void rk4(double* f,int d, int n,double dt, const double a, const double b, const double c)
{
	double k1[3];double k2[3];double k3[3];double k4[3];
	double ftemp[3];
	
	ofstream output("data.txt");

	for (int i=1;i<=n;i++)
	{
	function(k1,f,a,b,c);
        
        ftemp[0]=f[0]+(dt/2.)*k1[0];
        ftemp[1]=f[1]+(dt/2.)*k1[1];
        ftemp[2]=f[2]+(dt/2.)*k1[2];
        function(k2,ftemp,a,b,c);
        
        ftemp[0]=f[0]+(dt/2.)*k2[0];
        ftemp[1]=f[1]+(dt/2.)*k2[1];
        ftemp[2]=f[2]+(dt/2.)*k2[2];
        function(k3,ftemp,a,b,c);
        
        ftemp[0]=f[0]+dt*k3[0];
        ftemp[1]=f[1]+dt*k3[1];
        ftemp[2]=f[2]+dt*k3[2];
        function(k4,ftemp,a,b,c);
	
	
	f[0]=f[0]+(dt/6.)*(k1[0]+2.*k2[0]+2.*k3[0]+k4[0]);
	f[1]=f[1]+(dt/6.)*(k1[1]+2.*k2[1]+2.*k3[1]+k4[1]);
	f[2]=f[2]+(dt/6.)*(k1[2]+2.*k2[2]+2.*k3[2]+k4[2]);
	
	
	
	output << i*dt << "\t" << f[0] << "\t" << f[1] << "\t" << f[2] << endl;
	
	}
	output.close();

}

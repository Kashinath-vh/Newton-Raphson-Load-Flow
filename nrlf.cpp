#include<iostream>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;

complex<double> sum(complex<double> * x,int n)
{
int i;
complex<double> y;
for(i=1;i<=n;i++)
y=y + x[i];
return y;
}

double max(double* x,int n)
{int i;
double y;
y=abs(x[1]);
for(i=2;i<=n;i++)
if(abs(x[i])>y)
y=abs(x[i]);
return y;
}

int min(complex<double>* x,int n)
{int i;
double y;
y=abs(x[1]);
int index=1;
for(i=2;i<=n;i++)
if(abs(x[i])<y)
{
y=abs(x[i]);
index=i;
}
return index;
}

//	Invert the Jacobian by Shipley's method
int invshipley(double a[240][240],int n )
{
int i,j,k;
for(k=1;k<=n;k++)
{
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
{
if(i==k || j==k)
continue;
a[i][j]=a[i][j] - (a[i][k]*a[k][j])/a[k][k];
}
for(i=1;i<=n;i++)
{if (i==k)
continue;
a[i][k]=(-1*a[i][k])/a[k][k];
}
for(j=1;j<=n;j++)
{if (j==k)
continue;
a[k][j]=(-1.0*a[k][j])/a[k][k];
}
a[k][k]=(-1.0/a[k][k]);
}
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
a[i][j]=-1.0*a[i][j];
return 0;
}

int main()
{
int i,j,l,from_bus,to_bus,k,no_iterations,no_buses,no_lines,temp,precision_no=4,no_PQ,no_PV,linedata[120][3];
double tolerance,base_kV,base_MVA,r,x,J[240][240],delP[240],delV[240],Vang,Vmag,Ploss;
complex<double> z,Z_base,SL[120],Y[120][120],Scal[120],V[120],I[120],SLtotal,Il;
string comment;
fstream fp;

//	Open input file and read data
fp.open("69bus.txt",ios::in);
if(fp.is_open())
{
getline(fp,comment);
fp>>no_buses;
fp>>no_lines;
fp>>tolerance;
fp>>no_iterations;
fp>>base_kV;
fp>>base_MVA;
getline(fp,comment);
getline(fp,comment);
Z_base=(base_kV*base_kV)/base_MVA;

for(i=1;i<=no_lines;i++)
{
fp>>temp;
fp>>from_bus;
fp>>to_bus;
fp>>r;
fp>>x;
z=(complex<double>(r,x));
fp>>r;
fp>>x;

Y[from_bus][to_bus]=Y[from_bus][to_bus] - Z_base/z;
Y[from_bus][from_bus]=Y[from_bus][from_bus] + Z_base/z;
Y[to_bus][to_bus]=Y[to_bus][to_bus] + Z_base/z;
Y[to_bus][from_bus]=Y[from_bus][to_bus];
SL[to_bus]=(complex<double>(r,x))/(base_MVA*1000.0);
linedata[i][1]=from_bus;
linedata[i][2]=to_bus;
}
}

else
{
cout<<"Unable to open file\n";
return 0;
}
fp.close();

//	Open output file
fp.open("nrlfop.txt",ios::out);

V[1]=1.0 + 0i;
for(i=2;i<=no_buses;i++)
V[i]= 1. + 0i;

for(k=1;k<=no_iterations;k++)
{
for(i=1;i<=no_buses;i++)
{complex<double> temp=0;
for(j=1;j<=no_buses;j++)
temp+=Y[i][j]*V[j];
I[i]=temp;
}

for(i=1;i<=no_buses;i++)
Scal[i]=(V[i]*conj(I[i]));

//    J11
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
{
if(i==j)
J[i-1][j-1]=-(imag(Scal[i]) + abs(V[i])*abs(V[i])*imag(Y[i][i]));
else
J[i-1][j-1]=-(abs(V[i])*abs(V[j])*abs(Y[i][j])*sin( arg(Y[i][j]) + arg(V[j]) - arg(V[i])));
}
//J21
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
J[i-2+no_buses][j-1]= real(Scal[i]) - abs(V[i])*abs(V[i])*real(Y[i][i]);
else
J[i-2+no_buses][j-1]= -abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J12
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
J[i-1][j-2+no_buses]= real(Scal[i]) + abs(V[i])*abs(V[i])*real(Y[i][i]);
else
J[i-1][j-2+no_buses]=abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J22
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
J[i-2+no_buses][j-2+no_buses]= imag(Scal[i]) - abs(V[i])*abs(V[i])*imag(Y[i][i]);
else
J[i-2+no_buses][j-2+no_buses]=-abs(V[i])*abs(V[j])*abs(Y[i][j])*sin(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

invshipley(J,2*(no_buses - 1));
 
for(i=2;i<=no_buses;i++)
delP[i-1]=-(real(Scal[i]) + real(SL[i]));
for(i=no_buses;i<=2*(no_buses - 1);i++)
delP[i]=-(imag(Scal[i-no_buses+2]) + imag(SL[i-no_buses+2]));

for(i=1;i<=2*(no_buses - 1);i++)
{
double temp=0;
for(j=1;j<=2*(no_buses - 1);j++)
temp+=J[i][j]*delP[j];
delV[i]=temp;
}

for(i=1;i<=(no_buses - 1);i++)
{
Vang=arg(V[i+1])+delV[i];
Vmag=abs(V[i+1])*(1.0 + delV[i+ no_buses - 1]);
V[i+1]=polar(Vmag,Vang);
}

if(tolerance>max(delP,(2*no_buses-2)))  
break;
}
fp<<std::fixed<<std::setprecision(precision_no);
if(k<=no_iterations)
fp<<"\n\nSolution converged after "<<k<<" iterations";
else
fp<<"\n\nSolution did not converge";
fp<<"\nBus_no\t\tVoltage(magnitude)\t\tAngle(radian)\n";
fp<<"____________________________________________________________\n";
for(i=1;i<=no_buses;i++)
fp<<i<<"\t\t"<<abs(V[i])<<"\t\t\t\t"<<arg(V[i])<<"\n";
fp<<"____________________________________________________________\n";

fp<<"\nLine_no\t\tCurrent(magnitude)\t\tAngle(radian)\n";
fp<<"____________________________________________________________\n";
for(i=1;i<=no_lines;i++)
{
from_bus=linedata[i][1];
to_bus=linedata[i][2];
Il=(-V[from_bus]+V[to_bus])*Y[from_bus][to_bus];
fp<<i<<"\t\t"<<abs(Il)*base_MVA*1000/base_kV<<"\t\t\t\t"<<arg(Il)<<"\n";
}
fp<<"____________________________________________________________\n";

SLtotal=sum(Scal,no_buses);
Ploss=real(SLtotal)*100000;
fp<<"\n\nTotal P loss ="<<real(SLtotal)*100000<<" kW\t\tTotal Q loss ="<<imag(SLtotal)*100000<<" kVAR";
fp<<"\n\nVoltage is minimum at bus "<<min(V,no_buses);
fp.close();
return 0;
}

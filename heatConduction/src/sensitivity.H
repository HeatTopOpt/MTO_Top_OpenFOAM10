
Info<<"sensitivity analysis"<<endl;
fsensh=-3*(k-kmin)*xh*xh*(fvc::grad(T)&fvc::grad(T));
if (objFunction==2)
{
   fsensh=-3*(k-kmin)*xh*xh*(fvc::grad(T)&fvc::grad(Ta));
}
fsensh=fsensh/gMax(mag(fsensh.primitiveField()));
gsensh.primitiveFieldRef()=mesh.V()/gMax(mesh.V());
#include "filter_sens.H"

gx[0]=V;//constraint function

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=fsens[i]/N;//sensitivities of objective function
 dgdx[0][i]=gsens[i]/N;//sensitivities of constraint function
}


Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}


#include "filter_x.H"
//***********************************//update mu and lambda
opt++;

if(solid_area)
{
   forAll(cells, i)
   {
      xh[cells[i]]=1.0;
   }
}
xh.correctBoundaryConditions();
if(runTime.writeTime())
{
   xh.write();      
}
DT=kmin+(k-kmin)*xh*xh*xh;

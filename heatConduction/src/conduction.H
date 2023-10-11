Info<< "\nslove conduction\n" << endl;

solve(fvm::laplacian(DT, T)+Q);

if (objFunction==2)
{
   sumT=0;
   maxT=0;
   for(i=0;i<n;i++)
   {
      sumT=sumT+std::pow(T[i],Pnorm);
      maxT=std::max(T[i],maxT);
   }      
   reduce(sumT, sumOp<scalar>());
   reduce(maxT, maxOp<scalar>());
   sumT=std::pow(sumT,1.0/Pnorm-1);
   solve(fvm::laplacian(DT, Ta)+b2*sumT*Foam::pow(T,Pnorm-1));
}
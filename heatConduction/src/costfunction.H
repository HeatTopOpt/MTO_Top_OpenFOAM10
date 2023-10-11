Info<<"cost function"<<endl;

meanT=fvc::domainIntegrate(T).value()/area;
V=fvc::domainIntegrate(xh).value()/area-voluse; 

if (Pstream::master())
{ 
   ofstream outfile1("Voluse.txt",std::ios::app);
   outfile1<<V<<"\n";
   outfile1.close();

   ofstream outfile2("meanT.txt",std::ios::app);
   outfile2<<meanT<<"\n";
   outfile2.close(); 

   ofstream outfile3("Time.txt",std::ios::app);
   outfile3<<runTime.elapsedClockTime()<<"\n";
   outfile3.close();    

   ofstream outfile4("maxT.txt",std::ios::app);
   outfile4<<maxT<<"\n";
   outfile4.close();   
}
if(runTime.writeTime())
{
   T.write();     
}
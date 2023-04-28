clear all

%get photon number ampltiudes of kerr coherent states

i=sqrt(-1);


alp=0.7;

phi=0.5;

N=5;

for j=1:N
    
    %photon number
     p=j-1;
     
    
    c(j)= exp(i*phi*p^2) *alp^p/sqrt(factorial(p)) ;
    
    
    c(j);
end


c=c*exp(-abs(alp)^2/2);

c*c'

cp=c;
    


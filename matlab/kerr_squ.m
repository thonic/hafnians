clear all

%get photon number AMPLITUDEs of kerr-squeezed states

i=sqrt(-1);


%squeezing para
r=asinh(0.5);


%strength of kerr
phi=0.1;

N=10;

c=zeros(1,N);

for j=1:2:N
    
    
    %photon number
     p=j-1;
     
    
    c(j)= exp(i*phi*p^2)*sqrt(factorial(p)*tanh(r)^p);
    
    c(j)=c(j)/factorial(p/2)/2^(p/2);
    
    
    c(j);
end


c=c*sqrt(sech(r));
sum(c*c')



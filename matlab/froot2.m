clear all


%this program will find the roots and the 

i=sqrt(-1);

cp=[1; 1];

% cp=[1; 0.5; .25];
% cp=[1;1;1];

%kerr coherent state 
% cp=[0.7827;   0.4808 + 0.2627i;  -0.1129 + 0.2466i;  -0.0231 - 0.1071i;  -0.0056 + 0.0380i]

%kerr squeezed state 0.9973
% cp=[0.9457 + 0.0000i;   0.0000 + 0.0000i;   0.2755 + 0.1165i;   0.0000 + 0.0000i;  -0.0034 + 0.1158i];

% cp=rand(3,1); %+i*rand(3,1);

% cp=[1; 1; 1; 1];





cp=cp/sqrt(cp'*cp);
t=0.99999999;
phi=acos(t);
sq=asinh(1);

%size of state vector
     n=max(size(cp));
     
     
%max photon number
nm=n-1;
     
N=2*n;
maxn=N;

am=zeros(maxn+2);

ad=zeros(maxn+2);

cvac=zeros(n,1);
cvac(1)=1;

for j=1:maxn+1;
    am(j,j+1)=sqrt(j);
    ad(j+1,j)=sqrt(j);
end
    

ba=(cosh(sq)*am+ad);

bam=zeros(n,n,n+1);

for j=0:n
   
    ba2=ba^j;
    
    bam(:,:,j+1)=ba2(1:n,1:n);
    
end


% bam=ba(1:n,1:n);
bam2=ba^2;
bam2=bam2(1:n,1:n);

hvec=zeros(n,1);
hv2=hvec;

% cp2=cp;
cp3=cp;
% 
% hvec(1) = cp(n)/sqrt(factorial(nm));
% 
% cp2=cp2-hvec(1)*bam(:,:,n)*cvac
% 
% hvec(2) = cp(n-1)/sqrt(factorial(nm-1));
% 
% cp2=cp2-hvec(2)*bam(:,:,n-1)*cvac



for j=1:n
    
    hv2(j)=cp3(n-j+1)/sqrt(factorial(nm-j+1));
    
    hv2(j);
    
    cp3=cp3-hv2(j)*bam(:,:,n-j+1)*cvac;
    cp3;
    
    
end
    



% hvec(3)=cp2(n-2)/sqrt(factorial(nm-2));

% cp2=cp-hvec(1)*bam2*cvac-hvec(2)*bam(:,:,1)*cvac-hvec(3)*cvac;


cp4=cp;
for j=1:n
    
    cp4=cp4-hv2(j)*bam(:,:,n-j+1)*cvac;
    
%     rv(1,j)=hv2(n-j+1);
    
end

% rv=flip(hv2');

% r=roots(hv2);

hv2f=flip(hv2);


% r(1,1)=hv2(3);
% r(1,2)=hv2(2);
% r(1,3)=hv2(1);


% r=roots(hv2f);
r=roots(hv2);


m=zeros(n-1);

% beta = m * alp
% alp vec goes from alp(N+1) to alp (2)
% one more number than beta

for j=n-1:-1:1
    
    for k=n-1:-1:j
        
      
    m(j,k)=t^((n-1)-k);
    
    end
end

im=inv(m);

% alp vec = im  * r (beta )


alp=[0;im*r];



s1=0;
s2=0;
for j=2:n
    j
    s1=s1+alp(j)*t^(n+1-j);
    
    s2=s2+conj(alp(j))*t^(j-n-1);
    
end

s1=s1*cosh(sq);

sdiff = (s2-s1)/cosh(sq);

x = real(sdiff)/(t^n-t^-n/cosh(sq) );

y = imag(sdiff)/(t^n+t^-n/cosh(sq) );

alp(1)=x+i*y;

% alp(1)=(s2-s1)/(sqrt(2)*t^2-1/t^2);




s3=0;
s4=0;
for j=1:n
    s3=s3+alp(j)*t^(n+1-j);
    
    s4=s4+conj(alp(j))*t^(j-n-1);
    
end

s3=s3*cosh(sq);

norm(s3-s4)


%%%%%%%%
%%%%%%%
%%%%%





break



alp=zeros(n,1);

alp(n)=r(n-1);

% alp(n) = r(n-2)-alp(



break

alp=zeros(n,1);

alp(3)=r(2);

alp(2)=(r(1)-alp(3))/t;




s1=0;
s2=0;
for j=2:3
    s1=s1+alp(j)*t^(n-j);
    
    s2=s2+conj(alp(j))*t^(j-n);
    
end

s1=s1*cosh(sq);

alp(1)=(s2-s1)/(sqrt(2)*t^2-1/t^2);


a1=tanh(sq)*(alp(3)+alp(2)/t) - (alp(3)+t*alp(2));
a1=a1/(t^2-tanh(sq)/t^2);


s3=0;
s4=0;
for j=1:3
    s3=s3+alp(j)*t^(n+1-j);
    
    s4=s4+conj(alp(j))*t^(j-n-1);
    
end

s3=s3*cosh(sq);

norm(s3-s4)





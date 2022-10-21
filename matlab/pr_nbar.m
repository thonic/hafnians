%this program will calculate the probabilty 

% of photon pattern for all modes given a CV matrix, Displacment vector
%and photon pattern, pp


%takes the cv matrix (2Mx2M) and Big DV (2M) and photon pattern (M)


function x = pr_nbar(cv,bdv,pp)


N=max(size(cv))/2;
qcv = cv + eye(2*N)/2;
iq=inv(qcv);
f=eye(2*N)-iq;
A=[zeros(N) eye(N);eye(N) zeros(N)]*f; 

bigf = iq*bdv;

%norm fac
nf= exp(-1/2*bdv.'*iq*bdv)/sqrt(det(qcv));

%norm fac 2 - factorial factor bar(n)!
nf2=1;
for j=1:N
    nf2=nf2*factorial(pp(j));
end
nf2;
%photon pattern is repeated for alpha, alpha*
pp2= [pp pp];

%vector of variables 
av  = sym('a', [1 2*N]);
aexp4=av*A*av.'+2*av*bigf;
%multimode gaussian fucntion 
f1b(av)=exp(1/2*aexp4);

%copy of f 
g(av)=f1b;

for j=1:2*N
    
    p_ind=pp2(j);
%     
    for k=1:p_ind
%         
    j;
    av(j);
    g(av)=diff(g,av(j));

    end
end

c=cell(1,2*N);
for j=1:2*N
    c{j} = av(j);
end

v=zeros(1,2*N);

x = eval(subs(g,c,{v}))*nf/nf2;





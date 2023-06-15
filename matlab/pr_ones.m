%this program will calculate the probabilty 

% of 1 1 1 1 for all modes given a CV matrix 


%takes the cv matrix (2Mx2M) and Big DV (2M)

function x = pr_ones(cv,bdv)


%dimnesion of system
N=max(size(cv))/2;


qcv = cv + eye(2*N)/2;

iq=inv(qcv);

f=eye(2*N)-iq;
A=[zeros(N) eye(N);eye(N) zeros(N)]*f; 



bigf = iq*bdv;


%norm fac
nf= exp(-1/2*bdv.'*iq*bdv)/sqrt(det(qcv));

av  = sym('a', [1 2*N]);
 
aexp4=av*A*av.'+2*av*bigf;


f1(av)=exp(1/2*aexp4);


g(av)=f1;

for j=1:2*N
g(av)=diff(g,av(j));
end

c=cell(1,2*N);
for j=1:2*N
    c{j} = av(j);
end

v=zeros(1,2*N);

x = eval(subs(g,c,{v}))*nf;

% x=pr1;


%  av  = sym('a', [1; 2])
% 
% m=[2 1;1 3];
% 
% f(av) = exp( av*m*av.');
% g(av) = f;
% for j=1:2
% g(av)=diff(g,av(j));
% end

% eval(g(0,0))

% c=cell(1,2*N);
% for j=1:2*N
%     c{j} = av(j);
% end
% 
% v=zeros(1,2*N);
% 
% subs(g,c,{v});











% this program will take the alp's and squeezing
% and compare the output probabilities from the defn of the state
% with derivatives of the Q function
% and  haf/loop hafnian form 
% for a single mode

% USES froot2.m as an input to give alp parameter and cp initial state

i=sqrt(-1);

N=max(size(cp)); %total number of modes


%permutation matrix 
p=[zeros(N) eye(N);eye(N) zeros(N)];



%squeezing parameter 
sq=asinh(1);
%squeezing matrices
cm=eye(N);
cm(1,1)=cosh(sq);
sm=zeros(N);
sm(1,1)=sinh(sq);


%squeezing transformatriom symplectic 
s=[cm sm;sm cm];
s2=[cm -sm;-sm cm];

%initial covariance matrix 
cv=s*s'/2;


%transmission coeficitent 
% t=0.99999999;
phi=acos(t);



totbs=eye(N);
%beamsplitter matrix 
U=[cos(phi) sin(phi); sin(phi) -cos(phi)];

cv1=cv; %initial state 


%loop for transforming the covariance matrix ONLY 
for j=2:N
    
    bs=create_bs(U,1,j,N); %create beamsplitter between modes 1 and j
    
    bt=blkdiag(bs,conj(bs)); %symplectic beamsplitter trans
    
    totbs=bs*totbs;
    
    cv=bt*cv*bt'; %new covariance matrix 
    
end


tbs2=blkdiag(totbs,conj(totbs));



cv=s2*cv*s2';
cv1=cv;
qcv=cv+eye(2*N)/2;
iq=inv(qcv);
f=eye(2*N)-iq;
A=p*f; 

%displacement vector 
dv=zeros(N,1);
dv(1)=dv(1)+alp(1);


% loop for evolving displacement vector 
for j=2:max(size(alp));


bs=create_bs(U,1,j,N);

% beamsplitter transformation of dv
dv=bs*dv;

% add in new displacement 
dv(1)=dv(1)+alp(j);


end

% total displacment vector for a and a^*
bdv = [dv;conj(dv)];


% final parameters for Q function-  exp 1/2[ xAx + 2Fx ]
bigA = A; 
bigf = iq*bdv;

%norm const 
nf=exp(-1/2*bdv.'*iq*bdv)/sqrt(det(qcv));



%%%%%% stops here %%%%%%%%%%%

break





%prob of 1 photon in herald mode 
% is approx just prob 1 photon from coherent state

%reduced cv matrix 
red_bdv = delete_vec(bdv,[1,N+1]);
red_cv = delete_cv(cv,[1]);
red_qcv = delete_cv(qcv,[1]);
 red_iq = inv(red_qcv);
red_f=eye(2*N-2)-red_iq;
p_red=[zeros(N-1) eye(N-1);eye(N-1) zeros(N-1)];
 red_A=p_red*red_f;
 
 redf= red_iq*red_bdv;


rnf= exp(-1/2*red_bdv.'*red_iq*red_bdv)/sqrt(det(red_qcv));

ph12b = pr_ones(red_cv,red_bdv);

ph1=real(ph12b);


break


%1st way to calculate the prob of 1 photon in mode 2 (detection mode)
% ph1 = abs(bdv(2))^2*exp(-abs(bdv(2))^2);

% ph1 = ph1*abs(bdv(3))^2*exp(-abs(bdv(3))^2);




% pr_nbar(red_cv,red_bdv,[1 1 1 1])
% 
% p0= pr_nbar(cv,bdv,[0 1 1 1 1]) ;
% p1= pr_nbar(cv,bdv,[1 1 1 1 1]);
% p2= pr_nbar(cv,bdv,[2 1 1 1 1]);
% p3=pr_nbar(cv,bdv,[3 1 1 1 1]);
% p4=pr_nbar(cv,bdv,[4 1 1 1 1]);
% p5=pr_nbar(cv,bdv,[5 1 1 1 1]);
% 
% 
% e0=norm(p0/ph1-abs(cp(1)).^2)
% e1=norm(p1/ph1-abs(cp(2)).^2)
% e2=norm(p2/ph1-abs(cp(3)).^2)
% e3=norm(p3/ph1-abs(cp(4)).^2)
% e4=norm(p4/ph1-abs(cp(5)).^2)
e5=norm(p5/ph1-0)

break


syms x1 x2 x3 x4
xv=[x1 x2 x3 x4];


rxv=[x2 x4];

rxvd=rxv*[0 1;1 0];
rexp = rxv*red_A*rxv.'+ 2*rxv*redf;

%2nd way to calculate the prob of 1 photon in mode 2 (detection mode)
rf  =exp(1/2*rexp);
rf2(x2,x4) = diff(rf,x2,x4);
ph12 = eval(rf2(0,0))*rnf;

% pr_ones(cv,bdv)

xvd=xv*p;

aexp4=xv*bigA*xv.'+2*xv*bigf;
% aexp4=xv*bigA*xv.'-2*xv*bf2;

f1(x1,x2,x3,x4)=exp(1/2*aexp4);

f2(x1,x2,x3,x4)=diff(f1,x2,x4);



% 'Pr(1,0)'
% 
% pr10(x1,x2,x3,x4)=diff(f1,x1,x3);
% p10=eval(pr10(0,0,0,0))*nf



% 'Pr(0,1)'
pr01(x1,x2,x3,x4)=diff(f1,x2,x4);
% eval(pr01(0,0,0,0))*nf
% c(1,2)

p01=eval(pr01(0,0,0,0))*nf;

% p01/ph1

% pr_nbar(cv,bdv,[0 1])


% 'Pr(1,1)'
pr11(x1,x2,x3,x4)=diff(f1,x1,x2,x3,x4);


p11=eval(pr11(0,0,0,0))*nf;

% p11/ph1


% 'Pr(2,1)'
pr21(x1,x2,x3,x4)=diff(pr11,x1,x3);
%remember factorial 2!
p21=eval(pr21(0,0,0,0))*nf/2;


%3rd way to calculate the prob of 1 photon in mode 2 (detection mode)
%using full state Pr(0,1) + Pr(1,1) + Pr(1,2)...
ph13=p01+p11+p21;

 abs(cp).^2

'p01'
p01/ph13

pr_nbar(cv,bdv,[0 1])/ph13

norm(cp(1)^2-p01/ph13)

'p11'
p11/ph13

pr_nbar(cv,bdv,[1 1])/ph13


norm(abs(cp(2))^2-p01/ph13)

'p21'
p21/ph13

pr_nbar(cv,bdv,[2 1])/ph13


norm(0-p21/ph13)

break


%this will test to find the Q function of the states
%for a state of the form c_0|0> + c_1|1>


x2=0;
x4=0;


f3(x1,x3)=subs(f2);

f1b(x1,x3)=subs(f1);

f4(x1,x3)=f3/f1b;



vpa(f3,2)
vpa(f1b,2)
vpa(f3/f1b,2)
vpa(expand(f4),4)


eval(f3(0,0))\ph13


f4=(f3*inv(f1b))

vpa(f3,3)
%pretty(f3)

break

syms t r 

bs1=[t r 0 0;r -t 0 0; 0 0 1 0;0 0 0 1];

bs2=[t 0 r  0;0 1 0 0; r 0 -t 0;0 0 0 1];

bs3=[t 0 0 r;0 1 0 0; 0 0 1 0;r 0 0 -t];

bs4=[]









% this program will create the Non-gaussian output state

%it will create the CV matrix and Dv vector 



% and sample from it
% uses brute force diff of symbolic variables
%
% takes the alp values from froot2.m

i=sqrt(-1);


M=max(size(cp)); %number of modes per NG state, 
        %should be equal to number of photon basis states used i.e.
        %size(cp)

K=2; %number of NG states  AND number of interferometer modes
% NOTE this should be changed to general interferometers 


N=M*K; %total number of modes

%how to number modes 
% interferometer at M*(K-1) to M*K 

%squeezing parameter
sq=asinh(1);


%CV mat for single mode
cm=eye(M);
cm(1,1)=cosh(sq);
sm=zeros(M);
sm(1,1)=sinh(sq);

%squeezing matrix
s=[cm sm;sm cm];
%inverse squeezing matrix
s2=[cm -sm;-sm cm];


%initial squeezing transform
cv=s*s'/2;
phi=acos(t);
% phi=1e-01;
totbs=eye(M);
%beamsplitter for state generation, phi should be small
U=[cos(phi) sin(phi); sin(phi) -cos(phi)];

cv1=cv;
for j=2:M    
    bs=create_bs(U,1,j,M);    
    bt=blkdiag(bs,conj(bs));    
    totbs=bs*totbs;    
    cv=bt*cv*bt';    
end
%final squeezing transform
cv1=s2*cv*s2';




%displacement vector 
dv=zeros(M,1);
dv(1)=dv(1)+alp(1);
for j=2:max(size(alp));
bs=create_bs(U,1,j,M);
dv=bs*dv;
dv(1)=dv(1)+alp(j);
end
dv1=dv;
%big DV vector
bdv1 = [dv1;conj(dv1)];


%%%%%%%%%%%
%Full CV matrix for N=KxM modes
%%%%%%%%%%

cm=eye(N);
sm=zeros(N);
c=1;
for j=1:K
    cm(c,c)=cosh(sq);
    
    sm(c,c)=sinh(sq);
    
    c=c+M;
    
end

s=[cm sm;sm cm];
s2=[cm -sm;-sm cm];
%initial cv
cv=s*s'/2;
%transmission coeficitent 
% t=0.99999999;
phi=acos(t)
% phi=1e-05;
U=[cos(phi) sin(phi); sin(phi) -cos(phi)];
%beamsplitter for state generation is block-diagonal

%displacement vector 
dv=zeros(M,1);
dv(1)=dv(1)+alp(1);

bs1=eye(M);
for j=2:M
    
    bs=create_bs(U,1,j,M);
    
    bs1=bs*bs1;
    
    %for displacement vector
    dv=bs*dv;

    dv(1)=dv(1)+alp(j);
    
end

totbs=bs1;

bigdv = dv;
for j=1:K-1
    bigdv=[ bigdv; dv];
end


totdv = [bigdv;conj(bigdv)];

for j=1:K-1
    
    totbs=blkdiag(totbs,bs1);
    
end

totbs=blkdiag(totbs,conj(totbs));


cv=totbs*cv*totbs';
cv= s2*cv*s2';
cv_in=cv;


%interferometer uni 
%random interferometer
% this can be chosen to be anything 
% uint = RandOrthMat(K);

%50/50 beamsplitter
uint = [1 1;-1 1]/sqrt(2);

%reshape changes the 2x2 matrix to a larger matrix, KMxKM, between the modes 
% specified in the argument

bigu = reshape_mat(uint,M,K,N);
% bigu = create_bs(uint,1,K+1,N);

% bigt = bigu*totbs;

bigt=blkdiag(bigu,conj(bigu));

cv= bigt*cv*bigt';

return

%reduced CV (delete row & colummn) and DV (delete row) for herald modes

rcv1=   delete_cv(cv1,[1]);
rdv1=  delete_vec(dv1,[1]);
rbv1 = [rdv1; conj(rdv1)];

% pr of ones in one NG state 
pr1a = pr_ones(rcv1,rbv1);




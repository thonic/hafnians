

%this program will reshape matrix nxn to NxN

%small mat is m

%no of NG states :  K = N/M  

%modes per NG state is M

%total modes is N (M*K)

function x = reshape_mat(m,M,K,N)

   if M*K ~= N
    
    'error, M, K and N not sized properly'
    return
   end
   
   if max(size(m)) ~= K
    
    'error, matrix and M not sized properly'
    return
   end
   

   new_mat=eye(2*N);
   
   c1=1;

   j1=1;
   for j1=1:max(size(m))
        c2=1;
       for j2=1:max(size(m))
            
           new_mat(c1,c2)=m(j1,j2);
       
            new_mat(c1+M,c2+M) = conj(m(j1,j2));
           
           c2=c2+2*M;
       end 
       c1=c1+2*M;
   end
   

   
x=new_mat;
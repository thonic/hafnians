%this program will remove the rows/colums of the matrix 
%from the vector v

%v is the n_pattern of photons from the nmulticombs function 
%i.e. v(i,:)=[1 1 0 2.. ] is n_1=n_2=1, n_3=0, n_4=2...

%we now delete RC that n_j=0


function x = delete_rc2(m,v)

%     if max(size(v)) ~= max(size(m))
%         'error, m and v not sized properly'
%         return
%     end


%     for j=max(size(v)):-1:1
%     m(:,j)=[];
%     m(j,:)=[];
%     end
%    
% x=m;
    
for j=max(size(v)):-1:1
    
        if v(j)==0
        m(:,j)=[]; 
        m(j,:)=[];
        end
end

x=m;
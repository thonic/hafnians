

%this program will remove the entries of the displacement vector 
%given the vector v

function x = delete_vec(m,v)

%     if max(size(v)) ~= max(size(m))
%         'error, m and v not sized properly'
%         return
%     end
    
   for j=max(size(v)):-1:1
    
        if v(j)==0
        m(j)=[]; 
        
        end
end

x=m;
  
function x = find_partition(n)

%n is size of matrix

if rem(n,2) ~= 0
    'error-n not even'
    return;
end

s2=factd(n-1);


c=SetPartition(n,n/2);%total partitions of 2*n numbers into n groups

s1=max(size(c));%bell or stirlings number = total number of partitions
 
b=vertcat(c{:}); 
  
a=zeros(s2,n); %array that will hold the new partitions 


    t2=0;%counter for a perfect matching permutation

for j=1:s1
    
   
    
    t1=0;%counter for partition size
      
    for k=1:n/2  
      if max(size(b{j,k})) ~= 2
          t1=1;
          break
      end        
    end
 
           
    if t1 == 0
        t2=t2+1; 
       for k=1:n/2         
           a(t2,2*k-1)=b{j,k}(1,1);
            a(t2,2*k)=b{j,k}(1,2);           
       end
       
       
    end
    
    
end

if s2~=t2
'error, s2~=t2'
return;
end

for j=1:s2
	
	for k=1:2:n-1
        
	if a(j,k) > a(j,k+1)
        'error- a not order properly'
        return;
    end
    
	end
	

end

% dlmwrite('f2.txt',a,' ');

x=a;

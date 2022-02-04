clear all


%this function will find the perfect matching partitions

n=3;

s2=factd(2*n-1);


c=SetPartition(2*n,n);%total partitions of 2*n numbers into n groups

s1=max(size(c));%bell or stirlings number = total number of partitions
 
b=vertcat(c{:}); 
  
a=zeros(s2,2*n); %array that will hold the new partitions 


    t2=0;%counter for a perfect matching permutation

for j=1:s1
    
   
    
    t1=0;%counter for partition size
      
    for k=1:n  
      if max(size(b{j,k})) ~= 2
          t1=1;
          break
      end        
    end
 
           
    if t1 == 0
        t2=t2+1; 
       for k=1:n         
           a(t2,2*k-1)=b{j,k}(1,1);
            a(t2,2*k)=b{j,k}(1,2);           
       end
       
       
    end
    
    
end
    


t2

s2
if s2~=t2
'error, s2~=t2'
end

break

%check the order of the partitions 
for j=1:s2
	
	for k=1:2:2*n-1
        
	if a(j,k) > a(j,k+1)
        'error'
    end
    
	end
	

end


%write 'a' matrix to a file

A = magic(6);

fileID = fopen('myfile.txt','w');
nbytes = fprintf(fileID,'%5d %5d %5d %5d\n',a)
fclose(fileID);

type('myfile.txt')



break

The following should do the trick:

fid = fopen('coeffs.txt','wt');  % Note the 'wt' for writing in text mode
fprintf(fid,'%f\n',a);  % The format string is applied to each element of a
fclose(fid);

break

dlmwrite('f2.txt',a,' ');



fid=fopen('myfile.txt','w')

for j=1:15
    for k=1:6
    fprintf(fid,'%d',a(j,k));
    end
end

fclose(fid)






function x = hafnian(m)

s1=max(size(m));
n=max(size(m));


s2=factd(n-1);

a=find_partition(n);


tot=0;
for j=1:s2
	j;
	pr2=1;
	for k=1:2:s1-1
        k;
	pr2=pr2*m(a(j,k),a(j,k+1));
	end
	tot=tot+pr2;

end

x=tot;
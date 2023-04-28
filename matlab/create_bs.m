%this will create a beamsplitter matrix u
%between modes n and m given a system size N

%assume n<m

function y = create_bs(U,n,m,N)


if n == m
    'error, n = m'
    return
end

y=eye(N);

y(n,n)=U(1,1);
y(m,m)=U(2,2);

y(n,m)=U(1,2);
y(m,n)=U(2,1);


end
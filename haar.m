%Function to provide haar random unitary matrix.
function u=haar(n)
z=(normrnd(0,1,[n,n])+1i*normrnd(0,1,[n,n]))/sqrt(2);
[q,r]=qr(z);
lambda=zeros(n,n);
for i=1:n
    lambda(i,i)=r(i,i)/abs(r(i,i));
end
u=q*lambda; 
end 

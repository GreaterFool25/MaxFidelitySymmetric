%Function takes a vector and arranges it into a diagonal matrix.
function q= Diagm(Q)
q=0*Q*Q'; 
for i=1:size(Q,1)
    q(i,i)=Q(i,1); 
end 

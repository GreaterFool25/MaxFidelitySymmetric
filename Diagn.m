%Function called in CoherenceTest.m
function q= Diagn(Q)
q=0*Q; 
for i=1:size(Q,1)
    q(i,i)=Q(i,i); 
end 

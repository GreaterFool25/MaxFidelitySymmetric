%function to produce random density matrix. haar.m provides a haar random 
%unitary matrix.
function rho=rho_Rnd(d)
u=haar(d);dg=diag(rand(d,1));dg=dg/trace(dg);
rho=u*dg*u';
end


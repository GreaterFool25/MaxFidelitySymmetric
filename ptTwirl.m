%Function called in ptTwirlTest code. Rho is a bipartite state, with d-dim
%subsystems.
function rhoT=ptTwirl(rho)
d=sqrt(size(rho,1));
rhoT=kron(eye(d)/d,tracea(rho,d));
end 
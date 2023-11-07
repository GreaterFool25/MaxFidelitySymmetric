%Function called in ptTwirlTest code. Performs the twirl channel
%corresponding to the H_max problem.
% Rho is a bipartite state, with d-dim
%subsystems. Uses tracea.m that provides the partial trace (over
%subsystem A) for a given bipartite state on systems AB.

function rhoT=ptTwirl(rho)
d=sqrt(size(rho,1));
rhoT=kron(eye(d)/d,tracea(rho,d)); 
end 
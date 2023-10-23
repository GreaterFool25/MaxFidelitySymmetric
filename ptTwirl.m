%Function called in ptTwirlTest code.
function rhoT=ptTwirl(rho)
d=sqrt(size(rho,1));
rhoT=kron(eye(d)/d,tracea(rho,d));
end 
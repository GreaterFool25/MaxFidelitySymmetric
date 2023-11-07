%Function called in CoherenceTest code. Performs the twirl channel
%corresponding to the fidelity of asymmetry problem on a given state of any
%dimension.

function rhoT=Cohtwirl(rho)
    d=size(rho,1);
    rhoT=zeros(d,d);
    for z=0:(d-1)
        Z=zeros(d,d);
        for j=0:(d-1)
            Z(j+1,j+1)=exp(2*pi*j*z*1i/d);  %constructing generalised Z gate.
        end
        rhoT=rhoT+(1/d)*Z*rho*Z'; %constructing the operation of twirling gate 
                                  %by applying generalised Z gate 
    end 
end 


%Code to compare SDP's and iterative map algorithm for calculating fidelity of coherence. 
% Save function Cohtwirl.m (separate files) and then execute
% rest of the code. 
%The rho_Rnd.m (and Haar.m which it calls) function also needs to be 
% saved- it provides random density matrices 
zz=0;
TimeDataSDP=[];TimeDataIter=[];  %array to store runtime data.
x=4:4:36; %dimensions to be tested.

for d=x
    
T=100;     %no. of instances the runtimes of the iterative map and SDP are averaged over.
xSDP=[];xIter=[];
for t=1:T  % no. of cases to be averaged over.
rho=rho_Rnd(d); %random density matrix whose fidelity of coherence we intend to calculate.
tic
cvx_begin sdp    % defining the SDP to calculate the fidelity of coherence
                 %using CVX in matlab. The solver is set to SDPT3.
variable Z(d,d) complex
variable Q(d,1) nonnegative
maximize(0.5*trace(Z + Z'))
[rho Z; Z' Diagm(Q)] >= 0   % Other manners of expressing the same constraints 
                            %exist but this one is most efficient.
trace(Diagm(Q))==1          %diagonalises vector to a diagonal matrix
cvx_end
xSDP=[xSDP toc];

sig=(Cohtwirl(rho^.5))^2;sig=sig/trace(sig);%setting initial point to the
                                            %soln given for petz-renyi
                                            %entropies

vec=[];  % sig is random state that the 
                                           % iterative map will work on.
tic
while abs( abs( Fidelity(rho,sig/trace(sig)) )-cvx_optval )>10^(-5)  %code to run iterative map till it beats the SDP solution.
    sig=sig^(-1/2)*( Cohtwirl( (sig^.5*rho*sig^.5)^.5 ) )^2*sig^(-1/2);
end 
xIter=[xIter toc];
abs( Fidelity(rho,sig/trace(sig)) )-cvx_optval

end
TimeDataSDP=[TimeDataSDP; [mean(log(xSDP)) std(log(xSDP))]];
TimeDataIter=[TimeDataIter; [mean(log(xIter)) std(log(xIter))]];
zz=zz+1;
end

%Plotting code.
errorbar(x',TimeDataSDP(:,1),TimeDataSDP(:,2))
hold on 
errorbar(x',TimeDataIter(:,1),TimeDataIter(:,2))
title('Runtime comparison for SDP and Iterative Algorithm ')
xlabel('Dimension') 
xlim([1 x(end)+6]) 
ylim([-9,4])
ylabel('Log of runtime') 
legend({"SDP solver",'Iterative Algorithm'},'Location','southeast')

%Code to compare SDP's and iterative map algorithm for H_max. Save
%ptTwirl.m as separate function and execute code in command line.
%The rho_Rnd.m (and Haar.m which it calls) function also needs to be 
% saved- it provides random density matrices 
zz=0;
TimeDataSDP=[];TimeDataIter=[];
x=2:7; % dimensions to be tested.
X_Hmax_AllData_SDP=[];X_Hmax_AllData_Iter=[];

for d=x 
    
T=100; % no. of cases to be averaged over.
xSDP=[];xIter=[];
for t=1:T
rho=rho_Rnd(d^2);  %H_max calculated for this random state, can be modified 
% to state of choice.
tic
cvx_begin sdp   % defining the SDP to calculate the H_max
                 %using CVX in matlab. The solver is set to SDPT3.
variable Z(d^2,d^2) complex
variable Q(d,d) complex semidefinite
maximize(0.5*trace(Z + Z'))
[rho Z; Z' kron(eye(d)/d,Q)] >= 0 
trace(Q)==1
%Diagn(Q)-Q==0
cvx_end
xSDP=[xSDP toc];

sig=(ptTwirl(rho^.5))^2;sig=sig/trace(sig); %setting initial point to the
                                            %soln given for petz-renyi
                                            %entropies

vec=[];  %sig is the random seed for iterative map.
tic
%the while condition is to stop the iterative map once it is sufficiently
%close to the soln given by the SDP.
while abs( abs( Fidelity(rho,sig/trace(sig)) )-Fidelity(rho,kron(eye(d)/d,Q)/trace(kron(eye(d)/d,Q))) )>10^(-5)
    sig=sig^(-1/2)*( ptTwirl( (sig^.5*rho*sig^.5)^.5 ) )^2*sig^(-1/2);
    vec=[vec Fidelity(rho,sig)];   
end 
xIter=[xIter toc];
abs( Fidelity(rho,sig/trace(sig)) )-cvx_optval

end

X_Hmax_AllData_SDP=[X_Hmax_AllData_SDP;xSDP];
X_Hmax_AllData_Iter=[X_Hmax_AllData_Iter;xIter];
TimeDataSDP=[TimeDataSDP; [mean(xSDP) std((xSDP))]];
TimeDataIter=[TimeDataIter; [(mean(xIter)) std((xIter))]];
zz=zz+1;
end

%Plotting code.
errorbar(x',TimeDataSDP(:,1),TimeDataSDP(:,2))
hold on 
errorbar(x',TimeDataIter(:,1),TimeDataIter(:,2))
title('Runtime comparison between SDP and Iterative Algorithm ')
xlabel('Dimension') 
xlim([1 x(end)+1]) 
ylim([-9,5])
ylabel('Log of runtime') 
legend({"SDP solver",'Iterative Algorithm'},'Location','southeast')
%Code to compare SDP's and iterative map algorithm for Cohtwirl. Save
%function Cohtwirl.m (separate file) and then execute rest of the code. 
zz=0;
TimeDataSDP=[];TimeDataIter=[];
x=4:4:36; 
%+U5*rho*U5'+U6*rho*U6'+U7*rho*U7'+U8*rho*U8'+U9*rho*U9'+U10*rho*U10'
for d=x
    
T=100;
xSDP=[];xIter=[];
for t=1:T
rho=rho_Rnd(d);
tic
cvx_begin sdp
variable Z(d,d) complex
variable Q(d,d) complex semidefinite
maximize(0.5*trace(Z + Z'))
[rho Z; Z' Diagn(Q)] >= 0
trace(Q)==1
%Diagn(Q)-Q==0
cvx_end
xSDP=[xSDP toc];

sig=rho_Rnd(d);vec=[];sig=Cohtwirl(sig);
tic
while abs( abs( Fidelity(rho,sig/trace(sig)) )-cvx_optval )>10^(-5)
    sig=sig^(-1/2)*( Cohtwirl( (sig^.5*rho*sig^.5)^.5 ) )^2*sig^(-1/2);
    vec=[vec Fidelity(rho,sig)];
end 
xIter=[xIter toc];
abs( Fidelity(rho,sig/trace(sig)) )-cvx_optval

end
TimeDataSDP=[TimeDataSDP; [mean(log(xSDP)) std(log(xSDP))]];
TimeDataIter=[TimeDataIter; [mean(log(xIter)) std(log(xIter))]];
zz=zz+1;
end


errorbar(x',TimeDataSDP(:,1),TimeDataSDP(:,2))
hold on 
errorbar(x',TimeDataIter(:,1),TimeDataIter(:,2))
title('Runtime comparison for SDP and Iterative Algorithm ')
xlabel('Dimension') 
xlim([1 x(end)+6]) 
ylim([-7,4])
ylabel('Log of runtime') 
legend({"SDP solver",'Iterative Algorithm'},'Location','southeast')
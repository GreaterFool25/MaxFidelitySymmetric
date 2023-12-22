%This programs compares the fixed-point iterative algorithm with the SDP
%solver SDPT3 for the fidelity of coherence. We test 6 values of the dimension ranging
%form 4 to 44. For each dimension we record 100 samples.

clear all
format long


x=[4 12 20 28 36 44]; % dimensions to be tested.

% We collect all data in a matrix. 
%The rows label the dimensions, the columns label the samples for fixed dimension
X_coh_AllData_SDP=[]; %data SDP
X_coh_AllData_Iter=[]; %data of the fixed-point algorithm 

%Auxiliaries arrays
xSDP=[];
xIter=[];

errSDP=0;%counter for the Inaccurate solutions of the SDP

for d=x %loop for each dimension
    xSDP=[];xIter=[];
    T=100;  %numer of sample for each value of the dimension
    t=0;
    while t<T %loop for the samples
        t=t+1;
        R=randRho(d); %random matrix generation
        
        lambdamin=min(eig(R)); %minimum eigenvalue
        lambdamax=max(eig(R)); %maximum eigenvalue
        mu=lambdamin^(1/2)/(4*lambdamax^(3/2)); %strong convexity parameter
        
        %SDP solver
        tic
        cvx_begin sdp quiet
        variable Z(d,d) complex
        variable Q(d,1) nonnegative
        maximize(0.5*trace(Z + Z'))
        [R Z; Z' diag(Q)] >= 0   
        trace(diag(Q))==1          
        cvx_end
        
        %If the solution of the SDP is accurate, solve the fixed-point algorithm
        if strcmp(cvx_status, 'Solved')
            xSDP=[xSDP toc]; %record the SDP runtime sample
            
            
            S=(diag(diag((R^.5))))^2; %initial state for the fixed-point algorithm
            Grad=eye(d)-S^(-1/2)*diag(diag(((S^.5*R*S^.5)^.5)))*S^(-1/2); %gradient
            
            %We use the PL inequality ocndition to have a guarantee on the closeness of
            %each iteration to the optimal value
            
            %Fixed-point algorithm
            tic
            while (1/(2*mu))*trace(Grad*Grad')>10^(-9)
                S=S^(-1/2)*(diag(diag((S^.5*R*S^.5)^.5)))^2*S^(-1/2); 
                Grad=eye(d)-S^(-1/2)*diag(diag((S^.5*R*S^.5)^.5))*S^(-1/2);
            end 
            xIter=[xIter toc]; %record the fixed-point algorithm runtime sample
       else 
            errSDP=errSDP+1; %count the number of times the SDP retuns an inaccurate solution
            t=t-1;  %repeat the iteration if the solution of the SDP was inaccuarate
        end
    end
    %record all data
    X_coh_AllData_SDP=[X_coh_AllData_SDP; xSDP]
    X_coh_AllData_Iter=[X_coh_AllData_Iter; xIter]
end







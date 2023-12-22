%This programs compares the fixed-point iterative algorithm with the SDP
%solver SDPT3 for the max-conditional entropy. We test marginal dimensions ranging
%form 2 to 6. For each dimension we record 100 samples.

clear all
format long

x=[2 3 4 5 6 7]; % dimensions of the marginal to be tested.

% We collect all data in a matrix. 
%The rows label the dimensions, the columns label the samples for fixed dimension
X_Hmax_AllData_SDP=[]; %data SDP
X_Hmax_AllData_Iter=[]; %data of the fixed-point algorithm 

%Auxiliaries arrays
xSDP=[]; 
xIter=[];

errSDP=0; %counter for the Inaccurate solutions of the SDP
for d=x %loop for each dimension
    xSDP=[];xIter=[];
    T=100; %numer of sample for each value of the dimension
    t=0;
    while t<T  %loop for the samples 
        t=t+1;
        R=randRho(d^2); %random matrix generation
        
        lambdamin=min(eig(R)); %minimum eigenvalue
        lambdamax=max(eig(R)); %maximum eigenvalue
        mu=lambdamin^(1/2)/(4*lambdamax^(3/2)); %strong convexity parameter
        
        %SDP solver
        tic
        cvx_begin sdp quiet
        variable Z(d^2,d^2) complex
        variable Q(d,d) complex semidefinite
        maximize(0.5*trace(Z + Z'))
        [R Z; Z' kron(eye(d)/d,Q)] >= 0 
        trace(Q)==1
        cvx_end
      
        %If the solution of the SDP is accurate, solve the fixed-point algorithm
         
        if strcmp(cvx_status, 'Solved')
            xSDP=[xSDP toc]; %record the SDP runtime sample
            
            S=(ptTwirl(R^.5))^2; %initial state for the fixed-point algorithm
            Grad=eye(d^2)-S^(-1/2)*ptTwirl((S^.5*R*S^.5)^.5)*S^(-1/2); %gradient
            
            %We use the PL inequality ocndition to have a guarantee on the closeness of
            %each iteration to the optimal value
            
            %Fixed-point algorithm
            tic
            while (1/(2*mu))*trace(Grad*Grad')>10^(-9) 
                S=S^(-1/2)*(ptTwirl((S^.5*R*S^.5)^.5))^2*S^(-1/2); 
                Grad=eye(d^2)-S^(-1/2)*ptTwirl((S^.5*R*S^.5)^.5)*S^(-1/2);
            end 
            xIter=[xIter toc]; %record the fixed-point algorithm runtime sample
       else 
            errSDP=errSDP+1; %count the number of times the SDP retuns an inaccurate solution
            t=t-1; %repeat the iteration if the solution of the SDP was inaccuarate
        end
    end
    %record all data
    X_Hmax_AllData_SDP=[X_Hmax_AllData_SDP; xSDP]
    X_Hmax_AllData_Iter=[X_Hmax_AllData_Iter; xIter]
end







function flag=LMI_Aut16_th3(A,B,C,K,L,h,r0,etaM,alpha)
% This MATLAB program checks the feasibility of LMIs from Theorem 3 of the paper 
% A. Selivanov and E. Fridman, "Observer-based input-to-state stabilization of networked control systems with large uncertain delays," Automatica, 2016. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B,C - the parameters of the system (20); 
% K     - the controller gain from (22); 
% L     - the observer gain from (7); 
% h     - the sampling period; 
% r0    - the sensors-to-controller known constant delay; 
% etaM  - a bound for the sensors-to-controller time-varying delay uncertainty; 
% alpha - a desired decay rate; 

% Output: 
% flag==1 if LMIs are feasible, flag==0 otherwise. 

n=size(A,1); 

%% Decision variables 
P1=sdpvar(n); 
P2=sdpvar(n); 
S=sdpvar(n); 
R=sdpvar(n); 
W=sdpvar(n); 
G=sdpvar(n,n,'f'); 

%% Notations 
D=expm(A*r0)*L*C*expm(-A*r0); 
rhoM=exp(-2*alpha*etaM); 

%% The main LMI (N<=0)
N=blkvar; 
N(1,1)=P1*A+A'*P1+2*alpha*P1+S-rhoM*R; 
N(1,2)=-P1*D; 
N(1,3)=P1*B*K+rhoM*(R-G); 
N(1,4)=rhoM*G; 
N(1,5)=-P1*D; 
N(1,6)=etaM*A'*R;
N(2,2)=P2*(A+D)+(A+D)'*P2+2*alpha*P2; 
N(2,5)=P2*D; 
N(2,6)=-etaM*D'*R; 
N(2,7)=h*exp(alpha*h)*(A+D)'*W; 
N(3,3)=-rhoM*(R-G)-rhoM*(R-G)'; 
N(3,4)=rhoM*(R-G); 
N(3,6)=etaM*(B*K)'*R;
N(4,4)=-rhoM*(S+R); 
N(5,5)=-pi^2/4*W; 
N(5,6)=-etaM*D'*R; 
N(5,7)=h*exp(alpha*h)*D'*W; 
N(6,6)=-R; 
N(7,7)=-W; 
N=sdpvar(N); 

%% Park's condition 
Park=[R G; G' R]; 

%% Solution of LMIs
LMIs=[N<=0,P1>=0,P2>=0,S>=0,Park>=0,R>=0,W>=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    flag=(min(primal)>=0 && all(primal(1:3)>0)); 
else
    yalmiperror(sol.problem); 
end
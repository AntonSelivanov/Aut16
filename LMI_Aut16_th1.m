function OmegaVal=LMI_Aut16_th1(A,B,C,K,L,h,r0,etaM,r1,muM,alpha,sigma)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov and E. Fridman, "Observer-based input-to-state stabilization of networked control systems with large uncertain delays," Automatica, 2016. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B,C - the parameters of the system (6); 
% K     - the controller gain from (10); 
% L     - the observer gain from (7); 
% h     - the maximum sampling period: s_{k+1}-s_k<=h; 
% r0    - the sensors-to-controller known constant delay; 
% etaM  - a bound for the sensors-to-controller time-varying delay uncertainty; 
% r1    - the controller-to-actuators known constant delay; 
% muM   - a bound for the controller-to-actuators time-varying delay uncertainty; 
% alpha - a desired decay rate; 
% sigma - the event-triggering parameter from (5); 

% Output: 
% OmegaVal - the value of Omega from (5). If OmegaVal is empty, the LMIs are not feasible

n=size(A,1); 
m=size(B,2); 

%% Decision variables 
P1=sdpvar(n); 
P2=sdpvar(n); 
S=sdpvar(n); 
S0=sdpvar(n); 
S1=sdpvar(n); 
R0=sdpvar(n); 
R1=sdpvar(n); 
W=sdpvar(n); 
Omega=sdpvar(m); 
G0=sdpvar(n,n,'f');
G1=sdpvar(n,n,'f');
G2=sdpvar(n,n,'f');
G3=sdpvar(n,n,'f'); 

%% Notations
tauBar=h+etaM; 
tauM=r0+r1+etaM+muM+h; 
D=expm(A*(r0+r1))*L*C*expm(-A*(r0+r1)); 
H=tauBar^2*R0+(tauM-r0-r1)^2*R1; 
rhoA=expm(A*(r0+r1)); 
rhoBar=exp(-2*alpha*tauBar); 
rhoM=exp(-2*alpha*tauM); 

%% The main LMI (Phi<=0)
Phi=blkvar; 
Phi(1,1)=P1*A+A'*P1+2*alpha*P1+S0-rhoBar*R0; 
Phi(1,2)=P1*B*K+rhoBar*(R0-G0); 
Phi(1,3)=rhoBar*G0; 
Phi(1,8)=-P1*D; 
Phi(1,9)=-P1*D; 
Phi(1,11)=A'*H;
Phi(2,2)=-rhoBar*(R0-G0)-rhoBar*(R0-G0)'; 
Phi(2,3)=rhoBar*(R0-G0); 
Phi(2,11)=(B*K)'*H;
Phi(3,3)=rhoBar*(S-S0-R0); 
Phi(4,4)=exp(-2*alpha*(r0+r1))*(S1-S)-rhoM*R1; 
Phi(4,5)=rhoM*(R1-G1); 
Phi(4,6)=rhoM*(G1-G2); 
Phi(4,7)=rhoM*G2; 
Phi(5,5)=-rhoM*(R1-G1)-rhoM*(R1-G1)'; 
Phi(5,6)=rhoM*(R1-G1)-rhoM*(G3-G2); 
Phi(5,7)=rhoM*(G3-G2); 
Phi(5,8)=-(rhoA*B*K)'*P2; 
Phi(5,12)=-h*exp(alpha*h)*(rhoA*B*K)'*W; 
Phi(6,6)=-rhoM*(R1-G3)-rhoM*(R1-G3)'+sigma*K'*Omega*K;
Phi(6,7)=rhoM*(R1-G3); 
Phi(6,8)=(rhoA*B*K)'*P2; 
Phi(6,12)=h*exp(alpha*h)*(rhoA*B*K)'*W; 
Phi(7,7)=-rhoM*(S1+R1); 
Phi(8,8)=P2*(A+D)+(A+D)'*P2+2*alpha*P2; 
Phi(8,9)=P2*D; 
Phi(8,10)=P2*rhoA*B; 
Phi(8,11)=-D'*H; 
Phi(8,12)=h*exp(alpha*h)*(A+D)'*W; 
Phi(9,9)=-pi^2/4*W; 
Phi(9,11)=-D'*H; 
Phi(9,12)=h*exp(alpha*h)*D'*W; 
Phi(10,10)=-Omega; 
Phi(10,12)=h*exp(alpha*h)*(rhoA*B)'*W; 
Phi(11,11)=-H; 
Phi(12,12)=-W; 
Phi=sdpvar(Phi); 

%% Park's conditions 
Park0=[R0 G0; G0' R0]; 
Park1=[R1 G1; G1' R1]; 
Park2=[R1 G2; G2' R1]; 
Park3=[R1 G3; G3' R1]; 

%% Solution of LMIs
LMIs=[Phi<=0,P1>=0,P2>=0,S>=0,S0>=0,S1>=0,R0>=0,R1>=0,W>=0,Omega>=0,Park0>=0,Park1>=0,Park2>=0,Park3>=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    [primal,~]=check(LMIs); % Checking that the solver returned a proper solution
    if (min(primal)>=0 && all(primal(1:3)>0))
        OmegaVal=double(Omega); 
    end
else
    yalmiperror(sol.problem); 
end
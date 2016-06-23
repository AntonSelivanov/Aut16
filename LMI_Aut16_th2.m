function OmegaVal=LMI_Aut16_th2(A,B,C,K,L,h,r1,muM,alpha,sigma)
% This MATLAB program checks the feasibility of LMIs from Theorem 2 of the paper 
% A. Selivanov and E. Fridman, "Observer-based input-to-state stabilization of networked control systems with large uncertain delays," Automatica, 2016. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B,C - the parameters of the system (16); 
% K     - the controller gain from (17); 
% L     - the observer gain from (7); 
% h     - the waiting time from (15); 
% r1    - the controller-to-actuators known constant delay; 
% muM   - a bound for the controller-to-actuators time-varying delay uncertainty; 
% alpha - a desired decay rate; 
% sigma - the event-triggering parameter from (15); 

% Output: 
% OmegaVal - the value of Omega from (15). If OmegaVal is empty, the LMIs are not feasible

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
Omega=sdpvar(m); 
G0=sdpvar(n,n,'f');
G1=sdpvar(n,n,'f');

%% Notations
D=expm(A*r1)*L*C*expm(-A*r1); 
H=muM^2*R0+h^2*R1; 
rhoA=expm(A*r1); 
rhoM=exp(-2*alpha*(r1+muM)); 
rhoBar=exp(-2*alpha*(r1+h+muM)); 

%% The main LMIs (Xi<=0, Psi<=0)
Xi=blkvar; % t\in[t_k,t_k+h)
Xi(1,1)=P2*(A+D)+(A+D)'*P2+2*alpha*P2; 
Xi(1,2)=-D'*P1; 
Xi(1,3)=-P2*rhoA*B*K; 
Xi(1,5)=P2*rhoA*B*K; 
Xi(1,7)=-D'*H; 
Xi(2,2)=P1*(A+B*K)+(A+B*K)'*P1+2*alpha*P1+S; 
Xi(2,7)=(A+B*K)'*H; 
Xi(3,3)=exp(-2*alpha*r1)*(S0-S)-rhoM*R0; 
Xi(3,4)=rhoM*R0; 
Xi(4,4)=rhoM*(S1-S0-R0)-rhoBar*R1; 
Xi(4,5)=rhoBar*(R1-G1); 
Xi(4,6)=rhoBar*G1; 
Xi(5,5)=-rhoBar*(R1-G1)-rhoBar*(R1-G1)'; 
Xi(5,6)=rhoBar*(R1-G1); 
Xi(6,6)=-rhoBar*(S1+R1); 
Xi(7,7)=-H; 
Xi=sdpvar(Xi); 

Psi=blkvar; % t\in[t_k+h,t_{k+1})
Psi(1,1)=P2*(A+D)+(A+D)'*P2+2*alpha*P2; 
Psi(1,2)=-D'*P1; 
Psi(1,3)=-P2*rhoA*B*K; 
Psi(1,4)=P2*rhoA*B*K; 
Psi(1,7)=P2*rhoA*B; 
Psi(1,8)=-D'*H; 
Psi(2,2)=P1*(A+B*K)+(A+B*K)'*P1+2*alpha*P1+S; 
Psi(2,8)=(A+B*K)'*H; 
Psi(3,3)=exp(-2*alpha*r1)*(S0-S)-rhoM*R0; 
Psi(3,4)=rhoM*(R0-G0); 
Psi(3,5)=rhoM*G0; 
Psi(4,4)=-rhoM*(R0-G0)-rhoM*(R0-G0)'+sigma*K'*Omega*K; 
Psi(4,5)=rhoM*(R0-G0); 
Psi(5,5)=rhoM*(S1-S0-R0)-rhoBar*R1; 
Psi(5,6)=rhoBar*R1; 
Psi(6,6)=-rhoBar*(S1+R1); 
Psi(7,7)=-Omega; 
Psi(8,8)=-H; 
Psi=sdpvar(Psi); 

%% Park's conditions 
Park0=[R0 G0; G0' R0]; 
Park1=[R1 G1; G1' R1]; 

%% Solution of LMIs
LMIs=[P1>=0,P2>=0,S>=0,S0>=0,S1>=0,R0>=0,R1>=0,Omega>=0,Park0>=0,Park1>=0,Xi<=0,Psi<=0]; 
options=sdpsettings('solver','sedumi','verbose',0); 
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    if (min(primal)>=0 && all(primal(1:2)>0)) 
        OmegaVal=double(Omega); 
    end
else
    yalmiperror(sol.problem); 
end
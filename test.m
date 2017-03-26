% This MATLAB program checks the feasibility of LMIs from Theorems 1-3 of the paper 
% A. Selivanov and E. Fridman, "Observer-based input-to-state stabilization of networked control systems with large uncertain delays," Automatica, vol. 74, pp. 63–70, 2016
% for the inverted pendulum on a cart given in Section 5. 

%% System parameters
M=10;   % the cart mass
m=1;    % the pendulum mass
l=3;    % the length of the pendulum arm
g=10;   % the gravitational acceleration

A=[0 1 0 0; 0 0 -m*g/M 0; 0 0 0 1; 0 0 g/l 0]; 
B=[0; 1/M; 0; -1/(M*l)]; 
C=[1 0 0 0; 0 0 1 0]; 
K=[2 12 378 210]; 
L=-(place(A',C',[-4.2 -6 -7.1 -8]))'; 

%% LMIs of Theorem 1
h=.039; r0=.1; etaM=.005; r1=.1; muM=.005; alpha=.001; sigma=.01; 
display(['Omega=' num2str(LMI_Aut16_th1(A,B,C,K,L,h,r0,etaM,r1,muM,alpha,sigma))]); 

%% LMIs of Theorem 2
h=.088; r1=.1; muM=.005; alpha=.001; sigma=.01; 
display(['Omega=' num2str(LMI_Aut16_th2(A,B,C,K,L,h,r1,muM,alpha,sigma))]); 

%% LMIs of Theorem 3
h=.112; r0=.1; etaM=.005; alpha=.001;
if LMI_Aut16_th3(A,B,C,K,L,h,r0,etaM,alpha)
    display('Feasible'); 
else
    display('Not Feasible'); 
end 

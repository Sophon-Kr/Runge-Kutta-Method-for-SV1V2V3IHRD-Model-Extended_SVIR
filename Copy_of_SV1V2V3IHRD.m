clc;                                   % Clears the screen
clear all;

%%%% #specifying Parameters used# %%%%
rho = 66186727; 
zeta = 0.015;  % Death rate of COVID-19 (March, 2020)
eta=0.09;
omega1=0.04;
omega2=0.001;
epsilon1=0.641;
epsilon2=0.704;
myu=3.6529e-5;
alpha=0.2;
lambda=0.1;
beta = 0.02;
h = 0.0001;  
tfinal = 10;
N = ceil(tfinal/h); % Calculates upto y(10000) ceil(tfinal/h)

%%% #initial conditions# %%%
t(1) = 0;
S(1) = 10000;
V_1(1)= 10000;
V_2(1)= 10000;
I(1)= 10000;
H(1)= 10000;
R(1)=10000;
D=10000;
S0 = 10000;
V_10= 10000;
V_20= 10000;
I0= 10000;
H0= 10000;
R0=10000;
D0=10000;
%step size


% #initial condition equations# %
f1 = @ (t, S,V_1,V_2,I,H,R,D) eta*R-beta*S*I-omega1*S-myu*S; %s node 
f2 = @ (t, S,V_1,V_2,I,H,R,D) omega1*S-omega2*V_1-beta*(1-epsilon1)*V_1*I-myu*V_1; %v1 node 
f3 = @ (t, S,V_1,V_2,I,H,R,D) omega2*V_1-beta*(1-epsilon2)*V_2*I-myu*V_2; %v2 node
f4 = @ (t, S,V_1,V_2,I,H,R,D) beta*S*I+beta*(1-epsilon1)*V_1*I+beta*(1-epsilon2)*V_2*I-alpha*I-lambda*I-zeta*I ; %i node 
f5 = @ (t, S,V_1,V_2,I,H,R,D) alpha*I-lambda*I-zeta*H; %h node 
f6 = @ (t, S,V_1,V_2,I,H,R,D) lambda*I+lambda*I-myu*R-eta*R ; %r node 
f7 = @ (t, S,V_1,V_2,I,H,R,D) zeta*I+zeta*H ; %d node


for i=1:N  % update loop
%update time
t(i+1) = t(i) + h; 

  %%%%% Stage One %%%%%
 K1S = f1( t(i), S(i),V_1(i),V_2(i),I(i),H(i),R(i),D(i));
 K1V_1 = f2( t(i), S(i),V_1(i),V_2(i),I(i),H(i),R(i),D(i));
 K1V_2 = f3( t(i), S(i),V_1(i),V_2(i),I(i),H(i),R(i),D(i));
 K1I = f4( t(i), S(i),V_1(i),V_2(i),I(i),H(i),R(i),D(i));
 K1H = f5( t(i), S(i),V_1(i),V_2(i),I(i),H(i),R(i),D(i));
 K1R = f6( t(i), S(i),V_1(i),V_2(i),I(i),H(i),R(i),D(i));
 K1D = f7( t(i), S(i),V_1(i),V_2(i),I(i),H(i),R(i),D(i));
 
 %%%%% Stage Two %%%%%%
 K2S = f1( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2, + I(i)+0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R + D(i)+ 0.5*h*K1D);
 K2V_1 = f2 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2, + I(i)+0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R + D(i)+ 0.5*h*K1D);
 K2V_2 = f3 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2, + I(i)+0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R + D(i)+ 0.5*h*K1D);
 K2I = f4 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2, + I(i)+0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R + D(i)+ 0.5*h*K1D);
 K2H = f5 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2, + I(i)+0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R + D(i)+ 0.5*h*K1D);
 K2R = f6 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2, + I(i)+0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R + D(i)+ 0.5*h*K1D);
 K2D = f7 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2, + I(i)+0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R + D(i)+ 0.5*h*K1D);

 %%%%% Stage Three %%%%%
 K3S = f1( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2, + I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R + D(i)+ 0.5*h*K2D);
 K3V_1 = f2 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2, + I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R +D(i)+ 0.5*h*K2D);
 K3V_2 = f3 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2, + I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R +D(i)+ 0.5*h*K2D);
 K3I = f4 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2, + I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R +D(i)+ 0.5*h*K2D);
 K3H = f5 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2, + I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R +D(i)+ 0.5*h*K2D);
 K3R = f6 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2, + I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R +D(i)+ 0.5*h*K2D);
 K3D = f7 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2, + I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R +D(i)+ 0.5*h*K2D);

 %%%%%  Stage Four %%%%%
 K4S = f1( t(i)+ h, S(i)+ h*K3S,V_1(i)+ h*K3V_1,V_2(i)+ h*K3V_2, + I(i)+ h*K3I,H(i)+ h*K3H,R(i)+ h*K3R+D(i)+ h*K3D);
 K4V_1 = f2 ( t(i)+ h, S(i)+ h*K3S,V_1(i)+ h*K3V_1,V_2(i)+ h*K3V_2, + I(i)+ h*K3I,H(i)+ h*K3H,R(i)+ h*K3R+D(i)+ h*K3D);
 K4V_2 = f3 ( t(i)+ h, S(i)+ h*K3S,V_1(i)+ h*K3V_1,V_2(i)+ h*K3V_2, + I(i)+ h*K3I,H(i)+ h*K3H,R(i)+ h*K3R+D(i)+ h*K3D);
 K4I = f4 ( t(i)+ h, S(i)+ h*K3S,V_1(i)+ h*K3V_1,V_2(i)+ h*K3V_2, + I(i)+ h*K3I,H(i)+ h*K3H,R(i)+ h*K3R+D(i)+ h*K3D);
 K4H = f5 ( t(i)+ h, S(i)+ h*K3S,V_1(i)+ h*K3V_1,V_2(i)+ h*K3V_2, + I(i)+ h*K3I,H(i)+ h*K3H,R(i)+ h*K3R+D(i)+ h*K3D);
 K4R = f6 ( t(i)+ h, S(i)+ h*K3S,V_1(i)+ h*K3V_1,V_2(i)+ h*K3V_2, + I(i)+ h*K3I,H(i)+ h*K3H,R(i)+ h*K3R+D(i)+ h*K3D);
 K4D = f7 ( t(i)+ h, S(i)+ h*K3S,V_1(i)+ h*K3V_1,V_2(i)+ h*K3V_2, + I(i)+ h*K3I,H(i)+ h*K3H,R(i)+ h*K3R+D(i)+ h*K3D);
    
    
    %%%%% Now, the main equations %%%%%
    S(i+1) = S(i) +(1/6)*(K1S+ 2*K2S+2*K3S+K4S )*h;
    V_1(i+1)= V_1(i) +(1/6)*(K1V_1+ 2*K2V_1+2*K3V_1+K4V_1 )*h;
    V_2(i+1)= V_2(i) +(1/6)*(K1V_2+ 2*K2V_2+2*K3V_2+K4V_2 )*h;
    I(i+1)= I(i)+(1/6)*(K1I+ 2*K2I+2*K3I+K4I )*h;
    H(i+1)= H(i) +(1/6)*(K1H+ 2*K2H+2*K3H+K4H )*h;
    R(i+1)= R(i) +(1/6)*(K1R+ 2*K2R+2*K3R+K4R )*h;
    D(i+1)= D(i) +(1/6)*(K1D+ 2*K2D+2*K3D+K4D )*h;
  
    
   
 
    
end

%plot the solutions
%figure(1); clf(1)
%plot(t,S)

%hold on
%plot(t,E)

%plot(t,I)
%plot(t,R)


%legend('S(t)')
%legend('E(t)')
%legend('I(t)')
%legend('R(t)')

%  The function to plots all together
% plot(t,S,t,E,t,I,t,R)
% legend('S(t)', 'E(t)', 'I(t)', 'R(t)')
% xlabel('Time(years)')
% ylabel('Populations')
% set(gca, 'Fontsize', 12)

%  The function to plots all together
% plot(t,V_1,t,V_2,t,V_3,t,I)
% legend('V_1(t)','V_2(t)','V_3(t)','I(t)')
 plot(t, S,t,V_1,t,V_2,t,I,t,H,t,R,t,D)
 legend('S(t)','V_1(t)','V_2(t)','I(t)','H(t)','R(t)','D(t)')

xlabel('Time(years)')
ylabel('Populations')
set(gca, 'Fontsize', 12)

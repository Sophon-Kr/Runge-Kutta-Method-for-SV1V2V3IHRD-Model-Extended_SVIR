clc;                                   % Clears the screen
clear all;

%%%% #specifying Parameters used# %%%%
rho = 66186727; 
eta = 0.09;
beta = 2.27;
omega1 = 0.004;
omega2 = 0.001;
omega3 = 0.0001;
epsilon1 = 0.641; 
epsilon2 = 0.704;
epsilon3 = 0.767;
alpha = 0.2;
lamda = 0.1;
zeta = 0.015;
mu = 0.000036529;
h = 0.0001;  
tfinal = 10;
N = ceil(tfinal/h); % Calculates upto y(10000) ceil(tfinal/h)

%%% #initial conditions# %%%
t(1) = 0;
S(1) = 800000;
V_1(1)= 150000;
V_2(1)= 150000;
V_3(1)= 150000;
I(1)= 150000;
H(1)= 150000;
R(1)=150000;
S0 = 800000;
V_10= 150000;
V_20= 150000;
V_30= 150000;
I0= 150000;
H0= 150000;
R0=150000;
%step size


% #initial condition equations# %
f1 = @ (t, S,V_1,V_2,V_3,I,H,R) rho+(eta*R)-(beta*S*I)-(omega1*S)-(mu*S); %s node 
f2 = @ (t, S,V_1,V_2,V_3,I,H,R) (omega1*S)-(omega2*V_1)-(beta*(1-epsilon1)*V_1*I)-(mu*V_1); %v1 node 
f3 = @ (t, S,V_1,V_2,V_3,I,H,R) (omega2*V_1)-(omega3*V_2)-(beta*(1-epsilon2)*V_2*I)-(mu*V_2); %v2 node 
f4 = @ (t, S,V_1,V_2,V_3,I,H,R) (omega3*V_2)-(beta*(1-epsilon3)*V_3*I)-(mu*V_3); %v3 node 
f5 = @ (t, S,V_1,V_2,V_3,I,H,R) (beta*S*I)+(beta*(1-epsilon1)*V_1*I)+(beta*(1-epsilon2)*V_2*I)+(beta*(1-epsilon3)*V_3*I)-(alpha*I)-(lamda*I)-(zeta*I) ; %i node 
f6 = @ (t, S,V_1,V_2,V_3,I,H,R) (alpha*I)-(lamda*I)-(zeta*H) ; %h node 
f7 = @ (t, S,V_1,V_2,V_3,I,H,R) (lamda*I)+(lamda*I)-(mu*R)-(eta*R) ; %r node 



for i=1:N  % update loop
%update time
t(i+1) = t(i) + h; 

  %%%%% Stage One %%%%%
 K1S = f1( t(i), S(i),V_1(i),V_2(i),V_3(i),I(i),H(i),R(i));
 K1V_1 = f2( t(i), S(i),V_1(i),V_2(i),V_3(i),I(i),H(i),R(i));
 K1V_2 = f3( t(i), S(i),V_1(i),V_2(i),V_3(i),I(i),H(i),R(i));
 K1V_3 = f4( t(i), S(i),V_1(i),V_2(i),V_3(i),I(i),H(i),R(i));
 K1I = f5( t(i), S(i),V_1(i),V_2(i),V_3(i),I(i),H(i),R(i));
 K1H = f6( t(i), S(i),V_1(i),V_2(i),V_3(i),I(i),H(i),R(i));
 K1R = f7( t(i), S(i),V_1(i),V_2(i),V_3(i),I(i),H(i),R(i));
 
 %%%%% Stage Two %%%%%%
 K2S = f1( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2,V_3(i)+ 0.5*h*K1V_3,I(i)+ 0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R);
 K2V_1 = f2 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2,V_3(i)+ 0.5*h*K1V_3,I(i)+ 0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R);
 K2V_2 = f3 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2,V_3(i)+ 0.5*h*K1V_3,I(i)+ 0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R);
 K2V_3 = f4 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2,V_3(i)+ 0.5*h*K1V_3,I(i)+ 0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R);
 K2I = f5 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2,V_3(i)+ 0.5*h*K1V_3,I(i)+ 0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R );
 K2H = f6 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2,V_3(i)+ 0.5*h*K1V_3,I(i)+ 0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R);
 K2R = f7 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K1S,V_1(i)+ 0.5*h*K1V_1,V_2(i)+ 0.5*h*K1V_2,V_3(i)+ 0.5*h*K1V_3,I(i)+ 0.5*h*K1I,H(i)+ 0.5*h*K1H,R(i)+ 0.5*h*K1R);

 %%%%% Stage Three %%%%%
 K3S = f1( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2,V_3(i)+ 0.5*h*K2V_3,I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R);
 K3V_1 = f2 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2,V_3(i)+ 0.5*h*K2V_3,I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R);
 K3V_2 = f3 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2,V_3(i)+ 0.5*h*K2V_3,I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R);
 K3V_3 = f4 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2,V_3(i)+ 0.5*h*K2V_3,I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R);
 K3I = f5 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2,V_3(i)+ 0.5*h*K2V_3,I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R);
 K3H = f6 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2,V_3(i)+ 0.5*h*K2V_3,I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R);
 K3R = f7 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K2S,V_1(i)+ 0.5*h*K2V_1,V_2(i)+ 0.5*h*K2V_2,V_3(i)+ 0.5*h*K2V_3,I(i)+ 0.5*h*K2I,H(i)+ 0.5*h*K2H,R(i)+ 0.5*h*K2R);

 %%%%%  Stage Four %%%%%
 K4S = f1( t(i)+ 0.5*h, S(i)+ 0.5*h*K3S,V_1(i)+ 0.5*h*K3V_1,V_2(i)+ 0.5*h*K3V_2,V_3(i)+ 0.5*h*K3V_3,I(i)+ 0.5*h*K3I,H(i)+ 0.5*h*K3H,R(i)+ 0.5*h*K3R);
 K4V_1 = f2 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K3S,V_1(i)+ 0.5*h*K3V_1,V_2(i)+ 0.5*h*K3V_2,V_3(i)+ 0.5*h*K3V_3,I(i)+ 0.5*h*K3I,H(i)+ 0.5*h*K3H,R(i)+ 0.5*h*K3R);
 K4V_2 = f3 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K3S,V_1(i)+ 0.5*h*K3V_1,V_2(i)+ 0.5*h*K3V_2,V_3(i)+ 0.5*h*K3V_3,I(i)+ 0.5*h*K3I,H(i)+ 0.5*h*K3H,R(i)+ 0.5*h*K3R);
 K4V_3 = f4 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K3S,V_1(i)+ 0.5*h*K3V_1,V_2(i)+ 0.5*h*K3V_2,V_3(i)+ 0.5*h*K3V_3,I(i)+ 0.5*h*K3I,H(i)+ 0.5*h*K3H,R(i)+ 0.5*h*K3R);
 K4I = f5 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K3S,V_1(i)+ 0.5*h*K3V_1,V_2(i)+ 0.5*h*K3V_2,V_3(i)+ 0.5*h*K3V_3,I(i)+ 0.5*h*K3I,H(i)+ 0.5*h*K3H,R(i)+ 0.5*h*K3R);
 K4H = f6 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K3S,V_1(i)+ 0.5*h*K3V_1,V_2(i)+ 0.5*h*K3V_2,V_3(i)+ 0.5*h*K3V_3,I(i)+ 0.5*h*K3I,H(i)+ 0.5*h*K3H,R(i)+ 0.5*h*K3R);
 K4R = f7 ( t(i)+ 0.5*h, S(i)+ 0.5*h*K3S,V_1(i)+ 0.5*h*K3V_1,V_2(i)+ 0.5*h*K3V_2,V_3(i)+ 0.5*h*K3V_3,I(i)+ 0.5*h*K3I,H(i)+ 0.5*h*K3H,R(i)+ 0.5*h*K3R);
    
    
    %%%%% Now, the main equations %%%%%
    S(i+1) = S(i) +(1/6)*(K1S+ 2*K2S+2*K3S+K4S )*h;
    V_1(i+1)= V_1(i) +(1/6)*(K1V_1+ 2*K2V_1+2*K3V_1+K4V_1 )*h;
    V_2(i+1)= V_2(i) +(1/6)*(K1V_2+ 2*K2V_2+2*K3V_2+K4V_2 )*h;
    V_3(i+1)= V_3(i) +(1/6)*(K1V_3+ 2*K2V_3+2*K3V_3+K4V_3 )*h;
    I(i+1)= I(i) +(1/6)*(K1I+ 2*K2I+2*K3I+K4I )*h;
    H(i+1)= H(i) +(1/6)*(K1H+ 2*K2H+2*K3H+K4H )*h;
    R(i+1)= R(i) +(1/6)*(K1R+ 2*K2R+2*K3R+K4R )*h;
  
    
   
 
    
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
 plot(t, S,t,V_1,t,V_2,t,V_3,t,I,t,H,t,R)
 legend('S(t)','V_1(t)','V_2(t)','V_3(t)','I(t)','H(t)','R(t)')

xlabel('Time(years)')
ylabel('Populations')
set(gca, 'Fontsize', 12)

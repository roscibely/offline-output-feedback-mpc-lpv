clear; clc; close all
%% Model parameters
B = [0.0935; 0.00478];
C = [0.333 -1];
D = [0];
N = 5;
alf=[1 5];
bet=[0.1 1];
alpha=alf(2)*rand(1,N) +alf(1);
beta=bet(2)*rand(1,N)+bet(1);
A1 = [0.872 -0.0623*alf(1); 0.0935 0.997];
A2 = [0.872 -0.0623*alf(2); 0.0935 0.997];
B1= beta(1)*B;
B2=beta(N)*B;
n=size(B1,1); m=size(B1,2);
%Weighting matrix
Le = 1*eye(2);
R = 1;
%Constrain
umax = 1;
%Initial states
x = [-1.5; -0.2]; %initial
xo = [-0.5; 1]; %observer
u=1;
%% x_set
xset=x;
for k=1:N
    A = [0.872 -0.0623*alpha(k); 0.0935 0.997];
    xset(:,k+1)=A*xset(:,k);
end
%% LPV matrices
% LPV variables
alpha_=sdpvar(1);
beta_=sdpvar(1);
A_ = [0.872 -0.0623*alpha_; 0.0935 0.997];
B_ = beta_*B;
%% Off-line robust observer design
p = sqrt(0.6);
Ge = sdpvar(2,2, 'full');
Pe = sdpvar(2,2, 'symmetric');
Ye = sdpvar(2,1);
Lmi= [Pe>=0, [p^2*(Ge+Ge'-Pe)-Le (Ge*A_-Ye*C)'; Ge*A_-Ye*C Pe]>=0];
Lmi = [Lmi, alf(1)<=alpha_<=alf(2), uncertain(alpha_)];
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
optimize(Lmi,-trace(Ge),ops);
Lp = inv(value(Ge))*value(Ye);
e = x-xo;% estimation error
%% AW With relaxation
Xa=sdpvar(n,n, 'full');
Qa=sdpvar(n,n,'symmetric');
La =  sdpvar(m,n, 'full');
Ua = sdpvar(m,m,  'symmetric');
xp = sdpvar(n,1);
mua=sdpvar(1);
LMI1=[[-(Xa+Xa'-Qa) -La' zeros(n,m) (C*Xa+D*La)' (A1*Xa+B1*La)';
    -La -2*Ua eye(m) (D*Ua)' (B1*Ua)';
    zeros(m,n) eye(m) -mua*eye(m) zeros(m,n) zeros(m,m);
    (C*Xa+D*La) (D*Ua) zeros(m,m) -eye(m) zeros(m,n);
    (A1*Xa+B1*La) (B1*Ua) zeros(n,m) zeros(n,m) -Qa]<=0];
LMI1=[LMI1, [-(Xa+Xa'-Qa) -La' zeros(n,m) (C*Xa+D*La)' (A2*Xa+B2*La)';
    -La -2*Ua eye(m) (D*Ua)' (B2*Ua)';
    zeros(m,n) eye(m) -mua*eye(m) zeros(m,n) zeros(m,m);
    (C*Xa+D*La) (D*Ua) zeros(m,m) -eye(m) zeros(m,n);
    (A2*Xa+B2*La) (B2*Ua) zeros(n,m) zeros(n,m) -Qa]<=0];
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
AW = optimizer(LMI1,mua,ops,xp,{Xa,La,mua});
sol = AW{[0;0]};
Fa1 = sol{2}*inv(sol{1});
disp('AW gain'); fprintf('%f ', Fa1);
%%  OFF-line output feedback  MPC
%LMI variables
Q = sdpvar(2,2, 'symmetric');
gamma = sdpvar(1,1);
Y0 = sdpvar(1,2, 'full');
Y1 = sdpvar(1,2, 'full');
Y2 = sdpvar(1,2, 'full');
G = sdpvar(2,2, 'full');
X = sdpvar(1,1, 'full');
xp = sdpvar(2,1);
Y=Y0+alpha_*Y1+beta_*Y2;
% constraints and objective
objective = gamma;
%optimization object
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
gammav(1)=45;
% Optimization with LPV variables
LA_ = [G+G'-Q G*A_'+Y'*B_' G*sqrtm(Le) Y'*sqrtm(R);
    A_*G+B_*Y Q zeros(2,2) zeros(2,1);
    sqrtm(Le)*G zeros(2,2) gamma*eye(2) zeros(2,1);
    sqrtm(R)*Y zeros(1,2) zeros(1,2) gamma*eye(1)];
L3 = [[X Y; Y' G+G'-Q]>=0, X<=umax.^2];
Gn_=1*eye(2);
Qn_=1*eye(2);
u=0
for i = 1:N
    LMIs = [LA_ >= 0, L3, gamma*(G+G'-Q)-0.001*(G'*G)>=0];
    LMIs = [LMIs, alf(1)<=alpha_<=alf(2), uncertain(alpha_)];
    LMIs = [LMIs, bet(1)<=beta_<=bet(2), uncertain(beta_)];
    LMIs = [LMIs, G-Q>=0, G>=Gn_,Q>=Qn_, gamma<=gammav(i)];
    y(:,i)= C*xset(:,i);
    % minimization
    T(:,i)=A_*xset(:, i) + Lp*(y(:,i) -C*xset(:,i));
    L4 = [1  T(:,i)'         (sqrtm(Le)*xset(:,i))' (sqrtm(R)*u);
        T(:,i)       Q             zeros(2,2)    zeros(2,1);
        sqrtm(Le)*xset(:,i) zeros(2,2) gamma*eye(2)  zeros(2,1);
        sqrtm(R)*u zeros(1,2)    zeros(1,2)  gamma*eye(1)];
    LMIs = [LMIs, L4>=0]
    controller = optimizer(LMIs,objective,ops,xp,{G,Y0,Y1,Y2,gamma,Q});
    sol= controller{{xset(:,i)}};
    Gn(:,:,i) = sol{1};
    Yn0(:,:,i) = sol{2};
    Yn1(:,:,i) = sol{3};
    Yn2(:,:,i) = sol{4};
    gammav(i+1) = sol{5};
    Qn(:,:,i) = sol{6};
    Gn_=Gn(:,:,i);
    Qn_=Qn(:,:,i);
end
N_=5;
F0 = Yn0(:,:,N_)*inv(Gn(:,:,N_));
F1 = Yn1(:,:,N_)*inv(Gn(:,:,N_));
F2 = Yn2(:,:,N_)*inv(Gn(:,:,N_));
%% Simulation
xo = [-0.5; 1]; %observer
x = [-1.5; -0.2]; %initial
x1=x;
alpha=alf(2)*rand(1,100) +alf(1);
beta=bet(2)*rand(1,100)+bet(1);
u(1)=1; u_til=1;
time=100
r(1:20)=0; r(21:35)=1; r(36:time)=0; 
for i=1:time
    A = [0.872 -0.0623*alpha(i); 0.0935 0.997];
    % anti-windup
    x1(:,i+1) = A*x1(:,i)+beta(i)*B*u_til(i); % Tentar utilizar o xo
    y1(:,i)= C*x1(:,i);
    u1(i+1) = Fa1*x1(:,i+1);
    u_ant=u(i);
    u(i+1) =  -(F0+alpha(i)*F1+beta(i)*F2)*xo(:,i) + (r(i)-C*x(:,i)) -u1(i+1);
    if u(i+1)>1
        u(i+1) =1; end
    if u(i+1)<0
        u(i+1) =0; end
    u_dep=u(i+1);
    % Anti windup application
    u_til(i+1)=u_ant-u_dep;
    E(:,i) = eig(A+(F0+alpha(i)*F1+beta(i)*F2)*beta(i)*B);
    x(:,i+1) = A*x(:,i)+beta(i)*B*u(i+1);
    y(:,i)= C*x(:,i);
    %estimate state
    xo(:,i+1) = A*xo(:, i) + beta(i)*B*u(i+1) + Lp*(y(:,i) -C*xo(:,i));
    %dynamics error
    e(:,i+1) = (A-Lp*C)*e(:,i);
end

figure(1)
subplot(211);  stairs(x(1,:),'r-','LineWidth',1.5); hold on;
stairs(x(2,:),'k--','LineWidth',1.5); grid on; legend('x_1', 'x_2'); 
title('States');
subplot(212);  stairs(xo(1,:),'r-','LineWidth',1.5); hold on; 
stairs(xo(2,:),'k--','LineWidth',1.5); grid on; legend('x_1', 'x_2'); 
title('Observer states');
figure(3);
re=real(E);
im=imag(E);
plot(re',im','kx'); zgrid;

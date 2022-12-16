
%INPUT DATA

% % COMMON INPUT

%include global variables;
include_global;

%preprocessor phase
[d,v,mel,ke,fd,LM]=preprocessor2;

%initiating the figure for the later plot
% figure('Color',[1 1 1])

% % Calculation and ASSEMBLY of element matrices
 
%Stiffness Matrix
[K]=assembly_stiffness_matrix(ke);

% %Mass Matrix
[M]=assembly_mass_matrix(LM,mel);


%% STABILITY: here I call the stability subroutine for ALL the discretizations
% 
% % NB: the BCs are applied! That's why I reduce "red"
% 
% % Newmark
% 
% M_newk_red=M(2:3,2:3);
% K_newk_red=K(2:3,2:3);
% 
% [Eigen_A_newk,wh_vec_newk] = Stability(M_newk_red,K_newk_red);
% 
% % Jacobi
% 
% M_Jac_red=M(2:3,2:3);
% K_Jac_plus_red=diag(diag(K(2:3,2:3)));
% K_Jac_minus_red=K_Jac_plus_red-K(2:3,2:3);
% 
% [Eigen_A_Jac,wh_vec_Jac] = Stability(M_Jac_red,K_Jac_plus_red);
% 
% % Gauss_Seidel
% 
% M_GSeid_red=M(2:3,2:3);
% K_Gseid_plus_red=tril(K(2:3,2:3));
% K_Gseid_minus_red=K_Gseid_plus_red-K(2:3,2:3);
% 
% [Eigen_A_GSeid,wh_vec_GSeid] = Stability(M_GSeid_red,K_Gseid_plus_red);
% 
% 
% % Compute rho(M^-1*K) for time step
% 
% rho=max(abs(eig(inv(M_newk_red)*K_newk_red)))
% rhoJ=max(abs(eig(inv(M_Jac_red)*K_Jac_plus_red)))
% rhoGS=max(abs(eig(inv(M_GSeid_red)*K_Gseid_plus_red)))



%%

% %SOLUTION PHASE

% %SOLUTION of the STATIC PROBLEM to find the initial displacement d

% %initialize a0
a0=M\(fd(:,1)-K*d);
a=a0;
%a(1:nnp,1)=0;

%Extract K_E, K_F and K_EF matrices from K:
K_E = K(1:nd,1:nd);
K_F = K(nd+1:nnp,nd+1:nnp);
K_EF= K(1:nd,nd+1:nnp);

%Extract f_F vector (prescribed forces)
f_F = fd(1+nd:nnp,1);

%Extract d_E vector (prescribed displacements)
d_E = d(1:nd,1);

%FIND THE UNKNOWN (FREE) DISPLACEMENTS VECTOR d_F
d_F = K_F\(f_F-K_EF'*d_E);

%now let's reconstruct the global displacement vector by putting together the
%prescribed displacements d_E and the newly found d_F vector
d = [d_E
     d_F];
d=[0;6;12];
%COMPUTE THE UNKNOWN REACTION r
f_E = K_E*d_E+K_EF*d_F;

%so we take the initial displacement and 
%the associated initial istantaneous acceleration as I.C. (not for the
%analytical solution
a0=M\(-K*d);
% a0=(M+(dt^2)*K)\(-K*d);
a0(1)=0;

% % NEWMARK PREDICTOR-CORRECTOR SOLUTION - then used as a reference solution

[d1N,v1N,a1N,U_dN,U_vN,U_aN]=predictor_corrector(d,v,a0,M,K,fd);

% %plot of the NEWMARK P.C. solution
 x=linspace(0,dt*nt,nt+1);
 y2=U_dN(nnp,:);
% plot(x,y2,'g');
% hold on
% %

% % JACOBI WAVEFORM RELAXATION SOLUTION

%%start counting cpu time for JACOBI WR
cpu1=cputime;

% Partitioning

Kplus=diag(diag(K));
Kminus=Kplus-K;


% Kplus=tril(K);
% Kminus=Kplus-K;
% Starting solving

d_0=d;
v_0=v;
a_0=a0;


U_d=zeros(nnp,nt+1);
U_v=zeros(nnp,nt+1);
U_a=zeros(nnp,nt+1);
 
U_d(:,1)=d;
U_v(:,1)=v;
U_a(:,1)=a0;

%Iniziatise WR time-discrete matrix (initial WRs=initial conditions)
WR=zeros(nnp,nt+1);
WR2=zeros(nnp,nt+1);
for qq=1:nnp
WR(qq,:)=d(qq,1);
end

WRa=zeros(nnp,nt+1);
for qqa=1:(nnp)
    WRa(qqa,:)=a_0(qqa,1);
end

WRv=zeros(nnp,nt+1);
for qqa=1:(nnp)
    WRv(qqa,:)=v_0(qqa,1);
end

%Initialize first row of WR_STOR=ITERATION 0, THEN COMPARED WITH ITERATION
%1!

WR_STOR(1,:)=WR(nnp,:);
WR_STOR2(1,:)=WR2(nnp,:);
STOR_line=1;

%let's initialize e (error) and the counting i
e=1;
e2=1;
i=0;

%ni=2;
while ((e>10^-14||e2>10^-14))
    
%for iter_ations=1:1  
d=d_0;
v=zeros(nnp,1);
a=a0;

%SOLUTION OF THE MATRIX SYSTEM
%vector force
fd1=zeros(nnp,nt+1);
for jj=1:nt+1
fd1(:,jj)=Kminus*WR(:,jj); %per ogni colonna di WR (spostamento a un t moltiplico per Kminus che � come moltiplicare per l'1/6 di prima (coeff riduttivo/moltiplicativo))
end

 
%SOLUTION

A=M+beta_b*(dt^2)*Kplus;
[L,U]=lu(A);

for n=1:nt
%PREDICTOR PHASE
d1p=d(:,1)+dt*v(:,1)+((dt^2)/2)*(1-2*beta_b)*a(:,1);
v1p=v(:,1)+(1-gamma_b)*dt*a(:,1);

% ***mytag*** if break here, same results . 
disp('d1p');
disp(d1p);
 
%SOLUTION PHASE

%SOLUTION
% disp('b')
b=fd1(:,n+1)-Kplus*d1p;
disp('b')
disp(b)
break

%LU_decomposition
z=L\b;
a1=U\z;
a1(1,1)=0;



%CORRECTOR PHASE
d1=d1p+beta_b*(dt^2)*a1;
v1=v1p+(1-gamma_b)*dt*a1;

%SUBSTITUTING
U_d(:,n+1)=d1;
U_v(:,n+1)=v1;
U_a(:,n+1)=a1;
 
d(:,1)=d1;
v(:,1)=v1;
a(:,1)=a1;
 
 
WR(:,n+1)=d1;
WRv(:,n+1)=v1;
WRa(:,n+1)=a1;


end


i=i+1;

WR;
%disp('WR');
%disp(WR);



%If the counter is even e=1 so you do another cycle; if the counter is odd
%you store the WR for the considered point in the matrix WR_STOR as a row
%and you compare it with the previous row (previous iteration) to find the
%error and check if it's bigger than the desired limit.

% if mod(i,2) == 0
%   %iteration number is even
%   e=1;
% else
  %iteration number is odd
  STOR_line=STOR_line+1;
  WR_STOR(STOR_line,:)=WR(nnp,:);
  cc=WR_STOR(STOR_line,:)-WR_STOR(STOR_line-1,:);
  %disp('cc');
  %disp(cc);


  dd=abs(cc);
  e=max(dd);

  %disp('dd');
  %disp(dd);
  disp('e');
  disp(e);
  
% end
   
  WR_STOR2(STOR_line,:)=WR(nnp-1,:);
  cc2=WR_STOR2(STOR_line,:)-WR_STOR2(STOR_line-1,:);
  dd2=abs(cc2);
  e2=max(dd2);
  disp('e2');
  disp(e2);
  

% % %plot
% % C = {'c','r','g','y','k',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.
% x=linspace(0,dt*nt,nt+1);
% y=WR(nnp,:);
% % plot(x,y,'color',C{i})
% plot(x,y)
% hold on

end

disp('the CPU TIME of the Waveform Relaxation (Jacobi) is:');
cpu2=cputime-cpu1;
disp(cpu2);

disp('the number of iterations necessary for convergence by the Waveform Relaxation (Jacobi) is:');
disp(' ');         %insterted just to have a space in the plot of the answer
disp(i);

% %plot
 x=linspace(0,dt*nt,nt+1);
 yG2=WR(nnp,:);
% plot(x,yG2,'b')
% hold on
 
% PLOT EXACT SOLUTION: displacement over time of the last point of the rod
 
tf=dt*nt;
ti=0;
 
%displacement
ue=@(t) 0; 
ue=@(t) ue(t)+(1/(sqrt(6)))*(6*(sqrt(3))+6*(sqrt(6)))*cos((sqrt((1-((sqrt(2))/2))*(1/18)))*t)...
             -(1/(sqrt(6)))*(6*(sqrt(3))-6*(sqrt(6)))*cos((sqrt((1+(sqrt(2))/2)*(1/18)))*t);
t=ti:dt:tf;

%velocity
ve=@(t) 0; 
ve=@(t) ve(t)-((1/(sqrt(6)))*(sqrt((1-(sqrt(2))/2)*(1/18))))*(6*(sqrt(3))+6*(sqrt(6)))*sin((sqrt((1-((sqrt(2))/2))*(1/18)))*t)...
             +((1/(sqrt(6)))*(sqrt((1+(sqrt(2))/2)*(1/18))))*(6*(sqrt(3))-6*(sqrt(6)))*sin((sqrt((1+(sqrt(2))/2)*(1/18)))*t);
t=ti:dt:tf;
%acceleration
ae=@(t) 0; 
ae=@(t) ae(t)-((1/(sqrt(6)))*((1-(sqrt(2))/2)*(1/18)))*(6*(sqrt(3))+6*(sqrt(6)))*cos((sqrt((1-((sqrt(2))/2))*(1/18)))*t)...
             +((1/(sqrt(6)))*((1+(sqrt(2))/2)*(1/18)))*(6*(sqrt(3))-6*(sqrt(6)))*cos((sqrt((1+(sqrt(2))/2)*(1/18)))*t);
t=ti:dt:tf;
      
% %PLOT    
% plot(t,ue(t),'r')
% 
% xlabel('time');
% ylabel('horizontal displacement');
% title('horizontal displacement of last point (tip) of the rod');
% %
% legend('Newmark','Gauss-Seidel WR','Exact analytical solution with 2 elements/masses');
clc
clear all
close all
%% Using this algorithm , you can get beamforming and nullforming coding configurations!
%% 
global theta1 phi1 lambda NA NE Nx Ny  eps phi theta beamwidth_lower_theta beamwidth_lower_phi theta_desired1 phi_desired1 d f c k R D0 DD D thetam thetan phim phin
%%  Desired beam parameter set
% Note: The metasurface usually has Field of View (FoV) about [-60 60]!
theta_desired1=30;    %desired_theta 
phi_desired1=40;             %desired_phi

beamwidth_lower_theta = 10;   %
beamwidth_lower_phi = 10;
theta1 = deg2rad(theta_desired1);
phi1 = deg2rad(phi_desired1);
f=5.25*1e9;           %center_frequency
c=3*1e8;           %
lambda=c/f;  
k=2*pi/lambda;

thetam = find([-90:90] == (theta_desired1 - beamwidth_lower_theta/2));
thetan = find([-90:90] == (theta_desired1 + beamwidth_lower_theta/2));
phim = find([-90:90] == (phi_desired1 - beamwidth_lower_phi/2));
phin = find([-90:90] == (phi_desired1 + beamwidth_lower_phi/2));
disp([thetam thetan; phim phin]);
%%  Metasurface parameter set
d=0.024;
Nx = 16; % unit-cell number
Ny=16; 
NA = 181;
NE = 181;
Nx=16;
eps=0.0001;
R=1.5;
phi=linspace(-pi/2,pi/2,181);
theta=linspace(-pi/2,pi/2,181);
[D0,D,DD] = arriveMTSphase;
[P,P1,u,udeg] = idealandcomphase(theta1,phi1,DD,D);
%% 
[actual_BF,code_BF,const_BF,~] = SupersurfaceCode_jchfit(u,D);
%% 
[patternall,patterntheta] = Actualpatternjch(actual_BF);
figure;
% pattern_norm = pattern_dbw/max(max(pattern_dbw));
mesh(theta*180/pi,phi*180/pi,patternall);
% patternCustom(pattern_dbw,theta*180/pi,phi*180/pi);
xlabel('Elevation angles (deg)');
ylabel('Azimuth angles (deg)');
title('The beam pattern of Beamforming');
view([0,90]);
% plot(theta*180/pi,patterntheta);
title('Actual Beamforming pattern');
disp(sum(sum(patternall)))
%% 
[mainindex] = find([-90:1:90] == theta_desired1);
patternthetaMain = patterntheta(mainindex - beamwidth_lower_theta/2 : mainindex + beamwidth_lower_theta/2);
uuu  =reshape(D+u,1,256);
%% PSO Process
%% 
W1=1;%wight1
W2=1;%wight2
M=500;  %pop size 
NUM=256;      %
iter_max=100;  %Max Iteration
limt = [0,100];
% vlimt=[-180, 180]*pi/180;
vlimt=[-100, 100];  %
threshold=1e-100;   %

w=0.9;    
c1=2;      
c2=2;       
%% Population initialization


for i=1                                     %
pop(i,:)=uuu; 
      V(i,:)=vlimt(2)*rand(1,NUM);         
    [fitness(i),actualmix{i},paroutput(i,:)]=FitObjjch_dir_detaBF_threeOBJ(pop(i,:),DD,patternthetaMain);
end
%% 
[patternall,patterntheta] = Actualpatternjch(actualmix{1});
%% 
for i=2:M                                
pop(i,:)=uuu; 
    V(i,:)=vlimt(2)*rand(1,NUM);       
    [fitness(i),actualmix{i},paroutput(i,:)] = FitObjjch_dir_detaBF_threeOBJ(pop(i,:),DD,patternthetaMain);
end

%% Calculate Initial Value


[bestfitness, bestindex]=min(fitness);   
pBest=pop;                                          
gBest=pop(bestindex,:);  
gBactual =actualmix{i}; 
fitness_pbest=fitness;                          
fitness_gbest=bestfitness;                     
BestPopPar = paroutput(bestindex,:);
disp(BestPopPar);
%% 
x = gBest;
xactual = gBactual;
% Popt = [x;x;x;x;x;x;x;x;x;x;x;x;x;x;x;x];
Popt = reshape(x,Nx,Ny);
% [pattern,temp1] = Actualpatternjch(Popt);
figure
[ttttt,pattern_optimized0_discrete] = Actualpatternjch(xactual);
mesh(theta*180/pi,phi*180/pi,ttttt);
% patternCustom(pattern_dbw,theta*180/pi,phi*180/pi);
xlabel('Elevation angles (deg)');
ylabel('Azimuth angles (deg)');
view([0,90]);
title('The beam pattern of Null_forming');
%% PSO Main Loop
for i=1: iter_max           
    for j=1:M                  
        

        V(j,:)=w*V(j,:)+c1*rand*(pBest(j,:)-pop(j,:))+c2*rand*(gBest-pop(j,:));   %rand直接生成1个（0，1）的数据
        V(j,find(V(j,:)>vlimt(2)))=vlimt(2);
        V(j,find(V(j,:)<vlimt(1)))=vlimt(1);
        

        pop(j,:)=pop(j,:)+V(j,:);

%         pop(j, find(pop(j,:)>limt(2)))=limt(2);
%         pop(j,find(pop(j,:)<limt(1)))=limt(1);
%         pop(j,:) = vlimt(2)*rand(1,NUM);%%%%%

        [fitness(j),actualmix{j},paroutput(j,:)]=FitObjjch_dir_detaBF_threeOBJ(pop(j,:),DD,patternthetaMain);
        

        if  fitness(j) < fitness_pbest(j)
            pBest(j,:)=pop(j,:);
            fitness_pbest(j)=fitness(j);
        end
        
 
        if fitness(j) < fitness_gbest
            gBest=pop(j,:);
            gBactual = actualmix{j};

            fitness_gbest=fitness(j);
            BestPopPar = paroutput(j,:);
        end
        
    end
    allpop{i} = gBest;
    result(i) =fitness_gbest;  %储存历代全局历史最优适应度值
    ParResult(i,:) = BestPopPar;
    
   if result(i) < threshold
        phase=gBest;
        break
   end
  w=0.5*(iter_max-i)/iter_max+0.4;
   picture(i)=result(i);
   disp(i) 
   disp(fitness_gbest)
   actualMixtest{i} = gBactual;
end
%% 

x = gBest;
xactual = gBactual;
Popt = reshape(x,Nx,Ny);
%% The Nulling beam pattern in disired direction
[Patt,~] = Actualpatternjch(xactual);
figure
% pattern_norm = pattern_dbw/max(max(pattern_dbw));
mesh(theta*180/pi,phi*180/pi,Patt);
% patternCustom(pattern_dbw,theta*180/pi,phi*180/pi);
xlabel('Elevation angles (deg)');
ylabel('Azimuth angles (deg)');
view([0,90]);
title('The beam pattern of Null_forming');
sum(sum(Patt))
%% Evaluation and Observe Nulling beam coding configurations.
phasecom = xactual/lambda*2*pi - D ;
phasecom = rem(phasecom,2*pi);

phasedeg = rad2deg(phasecom);
codevali = zeros(16,16);
actualvali = zeros(16,16);
for i = 1:1:16
    for j = 1:1:16
        if phasedeg(i,j) < 0
            phasedeg(i,j) = phasedeg(i,j)+360;
        end
        if round(phasedeg(i,j)) == 0
            codevali(i,j) = 1;
            actualvali(i,j)=D(i,j)+deg2rad(0);
        elseif round(phasedeg(i,j)) == 90
            codevali(i,j) = 0;
            actualvali(i,j)=D(i,j)+deg2rad(90);
        elseif round(phasedeg(i,j)) == 180
            codevali(i,j) = 3;
            actualvali(i,j)=D(i,j)+deg2rad(180);
        elseif round(phasedeg(i,j)) == 270
            codevali(i,j) = 2;
            actualvali(i,j)=D(i,j)+deg2rad(270);
        end
    end
end
%% 
actualvali = mod(actualvali,2*pi);
actualvali=actualvali*lambda/(2*pi);
dlmwrite('Beamforming_Configuration.txt',code_BF);
dlmwrite('Nullforming_Configuration.txt',codevali);
%% 
for i = 1:5:100
    tt = actualMixtest{i};
    [aa,~] = Actualpatternjch(tt);
    figure
% pattern_norm = pattern_dbw/max(max(pattern_dbw));
mesh(theta*180/pi,phi*180/pi,aa);
% patternCustom(pattern_dbw,theta*180/pi,phi*180/pi);
xlabel('Elevation angles (deg)');
ylabel('Azimuth angles (deg)');
view([0,90]);
title('The beam pattern of Null_forming');
sum(sum(aa))
end
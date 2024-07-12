%% 
function [loss,actual_new,paroutput]= FitObjjch_dir_detaBF(xx,DD,pattern)
global D theta_desired1  phi_desired1 lambda Nx Ny beamwidth_lower_theta beamwidth_lower_phi thetam thetan phim phin
% actual=zeros(Ny,Nx);
% AF=zeros(length(phi),length(theta));
% x=reshape(x,Nx,Ny);
% for i=1:Ny
%     for j=1:Nx
%         actual(i,j)=phase1(i,j)+x(i,j);
%     end 
% end
% actual=actual*lambda/(2*pi);
xx = reshape(xx,Nx,Ny);
% u = reshape(u,Nx,Ny);
thetarem = [theta_desired1 - beamwidth_lower_theta/2:theta_desired1 + beamwidth_lower_theta/2];
phirem = [phi_desired1 - beamwidth_lower_phi/2:phi_desired1 + beamwidth_lower_phi/2];
theta = -90:1:90;
phi = -90:1:90;
% for i =1:1:length(thetarem)
%     for j =1:1:length(phirem)
%          AF0=exp(sqrt(-1) *2*pi/lambda*(sin(deg2rad(thetarem(i)))*sin(deg2rad(phirem(j)))*real(DD)+sin(deg2rad(thetarem(i)))*cos(deg2rad(phirem(j)))*imag(DD)...
%               -xx));
%          AF(j,i)=sum(sum(AF0));
%     end
% end
comphi = xx-D;
u = mod(comphi,2*pi);
udeg = rad2deg(u);
phase_dis = zeros(Nx,Ny);
for i = 1:1:16
    for j = 1:1:16
        if udeg(i,j) < 45 || udeg(i,j) > 315
            phase_dis(i,j) = 0;
        elseif udeg(i,j) < 135 && udeg(i,j) > 45
            phase_dis(i,j) = 90;
        elseif udeg(i,j) >135 && udeg(i,j) < 225
            phase_dis(i,j) = 180;
        elseif udeg(i,j) > 225 && udeg(i,j) < 315
            phase_dis(i,j) = 270;
        end
    end
end
actual_new = D + deg2rad(phase_dis);
actual_new = mod(actual_new,2*pi);
actual_new = actual_new*lambda/(2*pi);
otherrem = [-90:(thetarem(1)-1)  (thetarem(end)+1):90];
% for i =1:1:length(theta)
%     for j =1:1:length(phi)
%          AF0all=exp(sqrt(-1) *2*pi/lambda*(sin(deg2rad(theta(i)))*sin(deg2rad(phi(j)))*real(DD)+sin(deg2rad(theta(i)))*cos(deg2rad(phi(j)))*imag(DD)...
%               -actual_new));
%          AFall(j,i)=sum(sum(AF0all));
%     end
% end
% AFallSum = sum(sum(abs(20*log10(abs(AFall)))));
% AFallStabelPower = 20*log10(abs(AFall))/AFallSum*15e4;
% for i =1:1:length(otherrem)
%     for j =1:1:length(phirem)
%          AF0other=exp(sqrt(-1) *2*pi/lambda*(sin(deg2rad(otherrem(i)))*sin(deg2rad(phirem(j)))*real(DD)+sin(deg2rad(otherrem(i)))*cos(deg2rad(phirem(j)))*imag(DD)...
%               -actual_new));
%          AFother(j,i)=sum(sum(AF0other));
%     end
% end
for i =1:1:length(thetarem)
    for j =1:1:length(phirem)
         AF0=exp(sqrt(-1) *2*pi/lambda*(sin(deg2rad(thetarem(i)))*sin(deg2rad(phirem(j)))*real(DD)+sin(deg2rad(thetarem(i)))*cos(deg2rad(phirem(j)))*imag(DD)...
              -actual_new));
         AF(j,i)=sum(sum(AF0));
    end
end
% thetam = find(theta == (theta_desired1 - beamwidth_lower_theta/2));
% thetan = find(theta == (theta_desired1 + beamwidth_lower_theta/2));
% phim = find(phi == (phi_desired1 - beamwidth_lower_phi/2));
% phin = find(phi == (phi_desired1 + beamwidth_lower_phi/2));
% disp([thetam thetan; phim phin]);
% AF = AFallStabelPower(phim:phin , thetam:thetan);
% AFother = AFallStabelPower([1:phim-1 phin+1:end] , [1:thetam-1 thetan+1:end]);

% max_p=max(max(abs(AF) ) ) ;
AF=abs(20*log10(abs(AF)));
AFBF = max(max(pattern));% BF波束在期望none点平均值
% AFother = 20*log10(abs(AFother)); 
% afomin = min(min(AFother));
% afomax = max(max(AFother));
% afodeta = afomax - afomin;
% BFwidth = var(pattern);% BF波束在期望none点平整度（方差）
loss1=mean(mean(AF));% 期望none点平均值
loss2 = var(var(AF));% 期望none点平整度（方差）
% loss2 = max(max(AF)) - min(min(AF))
% loss2 = var(reshape(AF,1,(length(thetarem)*length(phirem))));
% % loss3 = var(var(AFother));% 期望其他方向平整度（方差）
% [e1,e2] = size(AFother);
% loss3 = var(reshape(AFother,1,e1*e2));
% corloss = [loss2,loss1,loss3];
% loss = 1/afovar;
% loss = sqrt((1/(abs(loss2-BFwidth))^2)+abs(loss1 - AFBF)^2+(1/abs((loss3-0)))^2);
% distance1 = abs(abs(loss1 - AFBF) - 60); % none点方向增益差与阈值 = 70的距离
% distance1  = abs(loss1 - AFBF);
% distance1  = abs(loss1);

distance2 = abs(loss2 - 0);% none点方向波束宽度与BF波束距离
% 
% distance3 = abs(loss3 - 0); % 其他方向与阈值 = 0方差的距离

paroutput = [loss1 loss2];

% loss = sqrt((loss1)^2 + distance2^2);
loss = loss1;
% if loss2 > 0.5
%     loss = 1e9;
% end

% loss = distance1;

% loss = sqrt(distance3^2);
% if loss2 > BFwidth
%     loss = 1e-5;
% end
% 
% if loss1 > AFBF
%     loss = 10e5;
% end
% 


% if afomax > mean(mean(AF))
%     loss = 1e5;
% end
end
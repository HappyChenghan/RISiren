function[pattern_dbw,temp1] = Actualpatternjch(B)
% This function can calculate the beam pattern of each metasurface's phase
global lambda d Nx Ny theta_desired1 phi_desired1 DD
eps = 0.0001;
NA = 181;
NE = 181;
theta_indx = find([-90:90] == theta_desired1);
phi_indx = find([-90:90] == phi_desired1);
%%
%求个角度相位
phi = linspace(-pi/2,pi/2,NA);
theta = linspace(-pi/2,pi/2,NE);
aa = d/2:d:(Ny-1)*d+d/2;
DD1 = repmat(aa',1,Nx);
bb = d/2:d:(Nx-1)*d+d/2;
DD2 = repmat(bb,Ny,1);
% DD = DD1+sqrt(-1).*DD2;
 for jj = 1:length(phi) 
    for ii = 1:length(theta)    
        pattern0 = exp(sqrt(-1) *2*pi/lambda*(sin(theta(ii))*sin(phi(jj))*real(DD)+sin(theta(ii))*cos(phi(jj))*imag(DD)...
              -B));
        pattern(jj,ii) = sum(sum(pattern0));
    end
end
%%
%% 画图
AFallSum = sum(sum(abs(20*log10(abs(pattern)))));
AFallStabelPower = 20*log10(abs(pattern))/AFallSum*15e4;
pattern_dbw = AFallStabelPower;
% pattern_dbw=abs(pattern);
%% 三维方向图

temp1 = pattern_dbw(phi_indx,:);
% figure
for n=1:length(theta)
   for m=1:length(phi)
       if n==91
       temp1(m)=pattern_dbw(n,m);
       end
end
end
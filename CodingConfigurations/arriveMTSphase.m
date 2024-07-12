function [D0,D,DD] = arriveMTSphase
global lambda d Nx Ny  R 
position_x=Nx*d/2+0;
position_y=Ny*d/2;%馈源位置
D0=zeros(Ny,Nx); 
for i=1:Ny
    for j=1:Nx
        dy=d/2+(i-1)*d-position_y;
        dx=d/2+(j-1)*d-position_x;
        D0(i,j)=sqrt(dy^2+dx^2+R^2);
       
    end
end
D=2*pi/lambda*(D0);%发送波束的实际相位
%% 
% [actual,code] =SupersurfaceCode1(theta1,phi1);
% [u1,D1,P1] = CompAndActual(theta1,phi1);
aa=d/2:d:(Ny-1)*d+d/2;
DD1=repmat(aa',1,Nx);
bb=d/2:d:(Nx-1)*d+d/2;
DD2=repmat(bb,Ny,1);
DD=DD1+sqrt(-1).*DD2;
end
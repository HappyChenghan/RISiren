function [P,P1,u,udeg] =idealandcomphase(theta1,phi1,DD,D)
global lambda
P1=sin(phi1)*sin(theta1)*real(DD)+sin(theta1)*cos(phi1)*imag(DD);
  
P=2*pi/lambda*P1;%目标波束相位
comphi=P-D;
% u1=mod(comphi,2*pi);
u=mod(comphi,2*pi);
udeg =rad2deg(u);
end
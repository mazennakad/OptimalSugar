 function [nu] = Dynamic_Viscosity_T_Sucrose_CW(T,S)
%--------------------------------------------------------------------------
%P.A. Bouchard and C. Grandjean (1995) A neural network correlation of the 
%variation of viscosity of sucrose aqueous solutions with temperature and 
%concentration. Lebemsm-Wis.u.-Tech.28, 157-159.
% Input - c in mol/m3; T = K
% Equations below apply for c = g/kg; T=C; nu=mPa s;
% 1 mol sucrose = 342.3 gm; %H2O density = 1g/ml
%--------------------------------------------------------------------------
T1=T-273.15; C=S*342.3/1000;
N=length(S);

U(1,1:N)=C/900; U(2,1:N)=T1/90; U(3,1:N)=1;

w(1,1)=-0.8232639;
w(1,2)=-6.2633560;
w(1,3)=0.9544210;
w(1,4)=-9.930766;    
w(1,5)=-46.32714;

w(2,1)=0.1552180; 
w(2,2)=3.655903; 
w(2,3)=0.8726018; 
w(2,4)=1.631972;         
w(2,5)=14.32810;

w(3,1)= 2.121665; 
w(3,2)=9.890403; 
w(3,3)=3.996007; 
w(3,4)=13.76070; 
w(3,5)=44.07883; 

w(4,1)=-26.77806; 
w(4,2)=-85.06680;  
w(4,3)=-62.71757; 
w(4,4)=-225.7701; 
w(4,5)=-345.4294; 
w(4,6)=740.5184; 


WW=w';
X=(1+exp(-WW(1:5,1:3)*U)).^(-1);
X(6,1:N)=1;
Y=(1+exp(-w(4,:)*X)).^(-1);
nu=(10.^(5.6*Y-1))*(1e-3);        % now in Pa s


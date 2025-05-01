function Phi = poly_matrix(p,r)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
Phi_temp=zeros(p^2,p^(r+1));  
coeffs=zeros(r+1,p^(r+1));
str=dec2base(0,p,r+1);
for j_1= 1:p^(r+1)
    coeffs_temp=dec2base(j_1-1,p,r+1)-str;
    for k_1=1:r+1
        if coeffs_temp(1,k_1)>p
            coeffs_temp(1,k_1)=coeffs_temp(1,k_1)-7;
        end
    end
    coeffs(:,j_1)=coeffs_temp';
end
for j_1= 1:p^(r+1)
   construction_M=zeros(p,p);
   for k_1=1:p
       construction_M(mod(polyval(coeffs(:,j_1),k_1-1),p)+1,k_1)=1;
   end
   construction_v=zeros(p^2,1);
   for i_1=1:p^2
       construction_v(i_1)=construction_M(i_1);
   end
   Phi_temp(:,j_1)=construction_v;
end
Phi=1/sqrt(p).*Phi_temp; 


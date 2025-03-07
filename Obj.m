function f = Obj(x,l1R,l2R,l1D,l2D)

% s1 = [x(1);x(2)]; s2 = [x(3);x(4)];
% 
phi_1 = x(5);phi_2 = x(6);phi_3 = x(7);phi_4 = x(8);
s11 = x(1);s12 = x(2);s21 = x(3);s22 = x(4);

[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4] = Para(s11,s12,s21,s22,phi_1,phi_2,phi_3,phi_4,l1R,l2R,l1D,l2D);


% f =  -(phi_1 - phi_3);
% f = x(7) - x(5);

% f =    (EnergyBarrier(x,l1R,l2R,l1D,l2D) - 0.01).^2;
f =  - EnergyBarrier(x,l1R,l2R,l1D,l2D);

% f =  norm(s1 + s2).^2;


end
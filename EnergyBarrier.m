function [Eb] = EnergyBarrier(x,l1R,l2R,l1D,l2D)
phi_1 = x(5);phi_2 = x(6);
eta_1 = phi_1;
phi_minus = min(phi_1,phi_2);
phi_plus = max(phi_1,phi_2);

interval_of_eta_2 = 100;

E_eta2 = zeros(interval_of_eta_2,2);
[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4] = Para(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),l1R,l2R,l1D,l2D);


for k = 1:interval_of_eta_2
    eta_2 = -pi + 2*pi/interval_of_eta_2*k;
    myenergyfun = @(t)Energy_spr(eta_1,eta_2,t(1),t(2),s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4);
    [xx,pe] = fminsearch(@(t) myenergyfun(t),[0,0]);
    E_eta2(k,1) = eta_2;
    E_eta2(k,2) = pe;


end


xx = E_eta2(:,1); yy = E_eta2(:,2);
kkk = find(xx < phi_plus & xx > phi_minus);
Eb = max(yy(kkk));



end
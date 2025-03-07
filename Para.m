function[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4,x_r,y_r,x_d,y_d] = Para(s11,s12,s21,s22,phi_1,phi_2,phi_3,phi_4,l1R,l2R,l1D,l2D)

s1 = [s11;s12];
s2 = [s21;s22];
Ident = Rot(0);

R1 = Rot(phi_1);R2 = Rot(phi_2); R3 = Rot(phi_3); R4 = Rot(phi_4);
delta_12 = R1 - R2;
delta_23 = R2 - R3;
delta_34 = R3 - R4;
delta_41 = R4 - R1;

% Matrix that solves all the equality constraints
s3 = inv(delta_34)*delta_41*s1 + (-Ident - inv(delta_34)*delta_23)*s2;
s4 = (-Ident - inv(delta_34)*delta_41)*s1 + inv(delta_34)*delta_23*s2;
t1 = Ident*s1 + inv(delta_12)*R2*l1R - inv(delta_12)*l1D;
t2 = Ident*s2 - inv(delta_12)*R1*l1R + inv(delta_12)*l1D;
t3 = inv(delta_34)*delta_41*s1 - (Ident + inv(delta_34)*delta_23)*s2 - inv(delta_34)*R4*l1R + inv(delta_34)*l1D;
t4 = (-Ident - inv(delta_34)*delta_41)*s1 + inv(delta_34)*delta_23*s2 + inv(delta_34)*R3*l1R - inv(delta_34)*l1D;
u1 = Ident*s1 + inv(delta_41)*R4*l2R - inv(delta_41)*l2D;
u2 = Ident*s2 + inv(delta_23)*R3*l2R - inv(delta_23)*l2D;
u3 = inv(delta_34)*delta_41*s1 + (-Ident - inv(delta_34)*delta_23)*s2 - inv(delta_23)*R2*l2R + inv(delta_23)*l2D;
u4 = -(Ident + inv(delta_34)*delta_41)*s1 + inv(delta_34)*delta_23*s2 - inv(delta_41)*R1*l2R + inv(delta_41)*l2D;
v1 = Ident*s1 + inv(delta_41)*R4*l2R + inv(delta_12)*R2*l1R - inv(delta_41)*l2D - inv(delta_12)*l1D;
v2 = Ident*s2 + inv(delta_23)*R3*l2R - inv(delta_12)*R1*l1R - inv(delta_23)*l2D + inv(delta_12)*l1D;
v3 = inv(delta_34)*delta_41*s1 - (Ident + inv(delta_34)*delta_23)*s2 - inv(delta_23)*R2*l2R - inv(delta_34)*R4*l1R + inv(delta_23)*l2D + inv(delta_34)*l1D;
v4 = -(Ident + inv(delta_34)*delta_41)*s1 + inv(delta_34)*delta_23*s2 - inv(delta_41)*R1*l2R + inv(delta_34)*R3*l1R + inv(delta_41)*l2D - inv(delta_34)*l1D;

x_11 = [0;0]; x_21 = s1; x_31 = s1 - u1; x_41 = t1;
x_12 = s1; x_22 = s1 + s2; x_32 = s1 + s2 - t2; x_42 = s1 + u2;
x_13 = s1 + s2; x_23 = s1 + s2 + t3; x_33 = s1 + s2 + t3 - v3; x_43 = -s4;
x_14 = [0;0]; x_24 = -s4; x_34 = -s4 + u4; x_44 = -t4;

x_r = [x_11(1) x_12(1) x_13(1) x_14(1);x_21(1) x_22(1) x_23(1) x_24(1);x_31(1) x_32(1) x_33(1) x_34(1);x_41(1) x_42(1) x_43(1) x_44(1)];
y_r = [x_11(2) x_12(2) x_13(2) x_14(2);x_21(2) x_22(2) x_23(2) x_24(2);x_31(2) x_32(2) x_33(2) x_34(2);x_41(2) x_42(2) x_43(2) x_44(2)];

y_11 = [0;0]; y_21 = R1*s1; y_31 = R1*(s1 - u1); y_41 = R1*t1;
y_12 = R1*s1; y_22 = R1*s1 + R2*s2; y_32 = R1*s1 + R2*(s2 - t2); y_42 = R1*s1 + R2*u2;
y_13 = R1*s1 + R2*s2; y_23 = R1*s1 + R2*s2 + R3*t3; y_33 = R1*s1 + R2*s2 + R3*(t3 - v3); y_43 = -R4*s4;
y_14 = [0;0]; y_24 = -R4*s4; y_34 = R4*(-s4 + u4); y_44 = -R4*t4;

x_d = [y_11(1) y_12(1) y_13(1) y_14(1);y_21(1) y_22(1) y_23(1) y_24(1);y_31(1) y_32(1) y_33(1) y_34(1);y_41(1) y_42(1) y_43(1) y_44(1)];
y_d = [y_11(2) y_12(2) y_13(2) y_14(2);y_21(2) y_22(2) y_23(2) y_24(2);y_31(2) y_32(2) y_33(2) y_34(2);y_41(2) y_42(2) y_43(2) y_44(2)];


end
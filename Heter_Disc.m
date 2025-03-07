%% initialization
ccc; tic
m = 40;
n = 40;
x_ref = zeros(2,1,m + 1,n + 1);
x_def = zeros(2,1,m + 1,n + 1);
Mat_of_coor_ref = cell(m,n);
Mat_of_coor_def = cell(m,n);
for j = 1:n
    for i = 1:m
        Mat_of_coor_def{i,j} = zeros(12,2);
        Mat_of_coor_ref{i,j} = zeros(12,2);
    end
end
for j = 1:(n + 1)
    for i = 1:(m + 1)
        x_ref(:,:,i,j) = [i - 1;j - 1]/m*2 - [1;1];
    end
end

%% sqaure to disc - v2 elliptical grid mapping
for jj = 1:(n + 1)
    for ii = 1:(m + 1)
        xtemp = x_ref(:,:,ii,jj); xxx = xtemp(1); yyy = xtemp(2);
        x_def(:,:,ii,jj) = [xxx*sqrt(1 - yyy.^2/2);yyy*sqrt(1 - xxx.^2/2)];
    end
end

fac = zeros(m*n,4);

vert_1 = zeros((m+1)*(n+1),2);
vert_1(:,1) = reshape(x_ref(1,:,:,:),[1,(m + 1)*(n + 1)]);
vert_1(:,2) = reshape(x_ref(2,:,:,:),[1,(m + 1)*(n + 1)]);
vert_2 = zeros((m+1)*(n+1),2);
vert_2(:,1) = reshape(x_def(1,:,:,:),[1,(m + 1)*(n + 1)]);
vert_2(:,2) = reshape(x_def(2,:,:,:),[1,(m + 1)*(n + 1)]);

for jjj = 1:n
    for iii = 1:m
        fac((jjj-1)*m+iii,1) = (m+1)*(jjj-1)+iii;
        fac((jjj-1)*m+iii,2) = (m+1)*(jjj-1)+iii+1;
        fac((jjj-1)*m+iii,3) = (m+1)*(jjj-1)+iii+1+m+1;
        fac((jjj-1)*m+iii,4) = (m+1)*(jjj-1)+iii+m+1;
    end
end

tcolor = repmat([1,0.89,0.88],[m*n,1]);
figure(1);axis equal;axis off;
title('Quad Mesh of Square');
patch('Faces',fac,'Vertices',vert_1,'FaceVertexCData',tcolor,'FaceColor','flat','LineWidth',1.5);
drawnow;


figure(2);axis equal;axis off;
title('Quad Mesh of Disc');
patch('Faces',fac,'Vertices',vert_2,'FaceVertexCData',tcolor,'FaceColor','flat','LineWidth',1.5);
drawnow;




%% start optimization
e1 = [1;0]; e2 = [0;1];
l1R = e1./m*2; l2R = e2./n*2;
aic = zeros(m*n,8);
ic0 = [0.45/m;-0.05/m;0.45/n;0.05/n;-0.2*pi;0.15*pi;-0.2*pi;0.2*pi];
aic(1,:) = ic0;

indexmatrix = zeros(m,n);
for j = 1:n
    for i = 1:m
        if j == 1
            indexmatrix(i,j) = i - 1;
        elseif i == 1
            indexmatrix(i,j) = 1 + (j - 2)*m;
        else
            indexmatrix(i,j) = (j - 1)*m + i - 1;

        end

    end
end


for jjjj = 1:n
    for iiii = 1:m

        l1D = (x_def(1,1,iiii+1,jjjj) - x_def(1,1,iiii,jjjj))*e1 + (x_def(2,1,iiii+1,jjjj) - x_def(2,1,iiii,jjjj))*e2; l2D = (x_def(1,1,iiii,jjjj+1) - x_def(1,1,iiii,jjjj))*e1 + (x_def(2,1,iiii,jjjj+1) - x_def(2,1,iiii,jjjj))*e2;
        indexx = indexmatrix(iiii,jjjj);
        ic = aic(indexx+1,:);
        options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',10000,'MaxIterations',1000,'Display','iter');
        if iiii == 1 && jjjj == 1

            [x,y] = fmincon(@(x)Obj(x,l1R,l2R,l1D,l2D),ic,[],[],[],[],[],[],@(x)Constraint(x,l1R,l2R,l1D,l2D),options);

        else
            x_0 = ic;
            [x,y] = fmincon(@(x)Obj_2(x,x_0),ic,[],[],[],[],[],[],@(x)Constraint(x,l1R,l2R,l1D,l2D),options);

        end

        aic((jjjj-1)*m+iiii+1,:) = x;

        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,x_r,y_r,x_d,y_d] = Para(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),l1R,l2R,l1D,l2D);

        tt = x_ref(:,:,iiii,jjjj);
        x_r = x_r + tt(1); y_r = y_r + tt(2);

        x_node1 = x_r(1,1); x_node2 = x_r(4,1); x_node3 = x_r(3,1); x_node4 = x_r(2,1);
        x_node5 = x_r(4,2); x_node6 = x_r(3,2); x_node7 = x_r(2,2); x_node8 = x_r(2,3);
        x_node9 = x_r(3,3); x_node10 = x_r(2,4); x_node11 = x_r(3,4); x_node12 = x_r(4,4);
        y_node1 = y_r(1,1); y_node2 = y_r(4,1); y_node3 = y_r(3,1); y_node4 = y_r(2,1);
        y_node5 = y_r(4,2); y_node6 = y_r(3,2); y_node7 = y_r(2,2); y_node8 = y_r(2,3);
        y_node9 = y_r(3,3); y_node10 = y_r(2,4); y_node11 = y_r(3,4); y_node12 = y_r(4,4);

        xref_temp = [x_node1,y_node1;x_node2,y_node2;x_node3,y_node3;x_node4,y_node4;x_node5,y_node5;x_node6,y_node6;x_node7,y_node7;x_node8,y_node8;x_node9,y_node9;x_node10,y_node10;x_node11,y_node11;x_node12,y_node12];
        Mat_of_coor_ref{iiii,jjjj} = xref_temp;


        ttt = x_def(:,:,iiii,jjjj);
        x_d = x_d + ttt(1); y_d = y_d + ttt(2);

        x_node1_d = x_d(1,1); x_node2_d = x_d(4,1); x_node3_d = x_d(3,1); x_node4_d = x_d(2,1);
        x_node5_d = x_d(4,2); x_node6_d = x_d(3,2); x_node7_d = x_d(2,2); x_node8_d = x_d(2,3);
        x_node9_d = x_d(3,3); x_node10_d = x_d(2,4); x_node11_d = x_d(3,4); x_node12_d = x_d(4,4);
        y_node1_d = y_d(1,1); y_node2_d = y_d(4,1); y_node3_d = y_d(3,1); y_node4_d = y_d(2,1);
        y_node5_d = y_d(4,2); y_node6_d = y_d(3,2); y_node7_d = y_d(2,2); y_node8_d = y_d(2,3);
        y_node9_d = y_d(3,3); y_node10_d = y_d(2,4); y_node11_d = y_d(3,4); y_node12_d = y_d(4,4);
        xdef_temp = [x_node1_d,y_node1_d;x_node2_d,y_node2_d;x_node3_d,y_node3_d;x_node4_d,y_node4_d;x_node5_d,y_node5_d;x_node6_d,y_node6_d;x_node7_d,y_node7_d;x_node8_d,y_node8_d;x_node9_d,y_node9_d;x_node10_d,y_node10_d;x_node11_d,y_node11_d;x_node12_d,y_node12_d];
        Mat_of_coor_def{iiii,jjjj} = xdef_temp;
        
    end
end

for j = 1:n
    for i = 1:m
        if i == 1 && j == 1
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempi1j = Mat_of_coor_ref{i + 1,j}; ref_tempij1 = Mat_of_coor_ref{i,j + 1};
            ref_tempij(6,:) = 1/2*(ref_tempij(6,:) + ref_tempi1j(2,:));
            ref_tempij(8,:) = 1/2*(ref_tempij(8,:) + ref_tempi1j(12,:));
            ref_tempij(9,:) = 1/2*(ref_tempij(9,:) + ref_tempij1(5,:));
            ref_tempij(11,:) = 1/2*(ref_tempij1(3,:) + ref_tempij(11,:));
            
            def_tempij = Mat_of_coor_def{i,j}; def_tempi1j = Mat_of_coor_def{i + 1,j}; def_tempij1 = Mat_of_coor_def{i,j + 1};
            def_tempij(6,:) = 1/2*(def_tempij(6,:) + def_tempi1j(2,:));
            def_tempij(8,:) = 1/2*(def_tempij(8,:) + def_tempi1j(12,:));
            def_tempij(9,:) = 1/2*(def_tempij(9,:) + def_tempij1(5,:));
            def_tempij(11,:) = 1/2*(def_tempij1(3,:) + def_tempij(11,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        elseif j == 1 && i < m
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempim1j = Mat_of_coor_ref{i - 1,j}; ref_tempi1j = Mat_of_coor_ref{i + 1,j}; ref_tempij1 = Mat_of_coor_ref{i,j + 1};
            ref_tempij(2,:) = ref_tempim1j(6,:); ref_tempij(12,:) = ref_tempim1j(8,:);
            ref_tempij(6,:) = 1/2*(ref_tempij(6,:) + ref_tempi1j(2,:));
            ref_tempij(8,:) = 1/2*(ref_tempij(8,:) + ref_tempi1j(12,:));
            ref_tempij(9,:) = 1/2*(ref_tempij(9,:) + ref_tempij1(5,:));
            ref_tempij(11,:) = 1/2*(ref_tempij1(3,:) + ref_tempij(11,:));

            def_tempij = Mat_of_coor_def{i,j}; def_tempim1j = Mat_of_coor_def{i - 1,j}; def_tempi1j = Mat_of_coor_def{i + 1,j}; def_tempij1 = Mat_of_coor_def{i,j + 1};
            def_tempij(2,:) = def_tempim1j(6,:); def_tempij(12,:) = def_tempim1j(8,:);
            def_tempij(6,:) = 1/2*(def_tempij(6,:) + def_tempi1j(2,:));
            def_tempij(8,:) = 1/2*(def_tempij(8,:) + def_tempi1j(12,:));
            def_tempij(9,:) = 1/2*(def_tempij(9,:) + def_tempij1(5,:));
            def_tempij(11,:) = 1/2*(def_tempij1(3,:) + def_tempij(11,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        elseif i == m && j == 1
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempim1j = Mat_of_coor_ref{i - 1,j}; ref_tempij1 = Mat_of_coor_ref{i,j + 1};
            ref_tempij(2,:) = ref_tempim1j(6,:); ref_tempij(12,:) = ref_tempim1j(8,:);
            ref_tempij(9,:) = 1/2*(ref_tempij(9,:) + ref_tempij1(5,:));
            ref_tempij(11,:) = 1/2*(ref_tempij1(3,:) + ref_tempij(11,:));

            def_tempij = Mat_of_coor_def{i,j}; def_tempim1j = Mat_of_coor_def{i - 1,j}; def_tempij1 = Mat_of_coor_def{i,j + 1};
            def_tempij(2,:) = def_tempim1j(6,:); def_tempij(12,:) = def_tempim1j(8,:);
            def_tempij(9,:) = 1/2*(def_tempij(9,:) + def_tempij1(5,:));
            def_tempij(11,:) = 1/2*(def_tempij1(3,:) + def_tempij(11,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        elseif i == 1 && j < n
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempijm1 = Mat_of_coor_ref{i,j - 1}; ref_tempi1j = Mat_of_coor_ref{i + 1,j}; ref_tempij1 = Mat_of_coor_ref{i,j + 1};
            ref_tempij(5,:) = ref_tempijm1(9,:); ref_tempij(3,:) = ref_tempijm1(11,:);
            ref_tempij(6,:) = 1/2*(ref_tempij(6,:) + ref_tempi1j(2,:));
            ref_tempij(8,:) = 1/2*(ref_tempij(8,:) + ref_tempi1j(12,:));
            ref_tempij(9,:) = 1/2*(ref_tempij(9,:) + ref_tempij1(5,:));
            ref_tempij(11,:) = 1/2*(ref_tempij1(3,:) + ref_tempij(11,:));

            def_tempij = Mat_of_coor_def{i,j}; def_tempijm1 = Mat_of_coor_def{i,j - 1}; def_tempi1j = Mat_of_coor_def{i + 1,j}; def_tempij1 = Mat_of_coor_def{i,j + 1};
            def_tempij(5,:) = def_tempijm1(9,:); def_tempij(3,:) = def_tempijm1(11,:);
            def_tempij(6,:) = 1/2*(def_tempij(6,:) + def_tempi1j(2,:));
            def_tempij(8,:) = 1/2*(def_tempij(8,:) + def_tempi1j(12,:));
            def_tempij(9,:) = 1/2*(def_tempij(9,:) + def_tempij1(5,:));
            def_tempij(11,:) = 1/2*(def_tempij1(3,:) + def_tempij(11,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        elseif i < m && j < n
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempijm1 = Mat_of_coor_ref{i,j - 1}; ref_tempim1j = Mat_of_coor_ref{i - 1,j}; ref_tempi1j = Mat_of_coor_ref{i + 1,j}; ref_tempij1 = Mat_of_coor_ref{i,j + 1};
            ref_tempij(5,:) = ref_tempijm1(9,:); ref_tempij(3,:) = ref_tempijm1(11,:);ref_tempij(2,:) = ref_tempim1j(6,:); ref_tempij(12,:) = ref_tempim1j(8,:);
            ref_tempij(6,:) = 1/2*(ref_tempij(6,:) + ref_tempi1j(2,:));
            ref_tempij(8,:) = 1/2*(ref_tempij(8,:) + ref_tempi1j(12,:));
            ref_tempij(9,:) = 1/2*(ref_tempij(9,:) + ref_tempij1(5,:));
            ref_tempij(11,:) = 1/2*(ref_tempij1(3,:) + ref_tempij(11,:));

            def_tempij = Mat_of_coor_def{i,j}; def_tempijm1 = Mat_of_coor_def{i,j - 1}; def_tempim1j = Mat_of_coor_def{i - 1,j}; def_tempi1j = Mat_of_coor_def{i + 1,j}; def_tempij1 = Mat_of_coor_def{i,j + 1};
            def_tempij(5,:) = def_tempijm1(9,:); def_tempij(3,:) = def_tempijm1(11,:);def_tempij(2,:) = def_tempim1j(6,:); def_tempij(12,:) = def_tempim1j(8,:);
            def_tempij(6,:) = 1/2*(def_tempij(6,:) + def_tempi1j(2,:));
            def_tempij(8,:) = 1/2*(def_tempij(8,:) + def_tempi1j(12,:));
            def_tempij(9,:) = 1/2*(def_tempij(9,:) + def_tempij1(5,:));
            def_tempij(11,:) = 1/2*(def_tempij1(3,:) + def_tempij(11,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        elseif i == m && j < n
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempijm1 = Mat_of_coor_ref{i,j - 1}; ref_tempim1j = Mat_of_coor_ref{i - 1,j}; ref_tempij1 = Mat_of_coor_ref{i,j + 1};
            ref_tempij(5,:) = ref_tempijm1(9,:); ref_tempij(3,:) = ref_tempijm1(11,:);ref_tempij(2,:) = ref_tempim1j(6,:); ref_tempij(12,:) = ref_tempim1j(8,:);
            ref_tempij(9,:) = 1/2*(ref_tempij(9,:) + ref_tempij1(5,:));
            ref_tempij(11,:) = 1/2*(ref_tempij1(3,:) + ref_tempij(11,:));

            def_tempij = Mat_of_coor_def{i,j}; def_tempijm1 = Mat_of_coor_def{i,j - 1}; def_tempim1j = Mat_of_coor_def{i - 1,j}; def_tempij1 = Mat_of_coor_def{i,j + 1};
            def_tempij(5,:) = def_tempijm1(9,:); def_tempij(3,:) = def_tempijm1(11,:);def_tempij(2,:) = def_tempim1j(6,:); def_tempij(12,:) = def_tempim1j(8,:);
            def_tempij(9,:) = 1/2*(def_tempij(9,:) + def_tempij1(5,:));
            def_tempij(11,:) = 1/2*(def_tempij1(3,:) + def_tempij(11,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        elseif i == 1 && j == n
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempijm1 = Mat_of_coor_ref{i,j - 1}; ref_tempi1j = Mat_of_coor_ref{i + 1,j};
            ref_tempij(5,:) = ref_tempijm1(9,:); ref_tempij(3,:) = ref_tempijm1(11,:);
            ref_tempij(6,:) = 1/2*(ref_tempij(6,:) + ref_tempi1j(2,:));
            ref_tempij(8,:) = 1/2*(ref_tempij(8,:) + ref_tempi1j(12,:));

            def_tempij = Mat_of_coor_def{i,j}; def_tempijm1 = Mat_of_coor_def{i,j - 1}; def_tempi1j = Mat_of_coor_def{i + 1,j};
            def_tempij(5,:) = def_tempijm1(9,:); def_tempij(3,:) = def_tempijm1(11,:);
            def_tempij(6,:) = 1/2*(def_tempij(6,:) + def_tempi1j(2,:));
            def_tempij(8,:) = 1/2*(def_tempij(8,:) + def_tempi1j(12,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        elseif j == n && i ~= m
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempijm1 = Mat_of_coor_ref{i,j - 1}; ref_tempim1j = Mat_of_coor_ref{i - 1,j}; ref_tempi1j = Mat_of_coor_ref{i + 1,j};
            ref_tempij(5,:) = ref_tempijm1(9,:); ref_tempij(3,:) = ref_tempijm1(11,:);ref_tempij(2,:) = ref_tempim1j(6,:); ref_tempij(12,:) = ref_tempim1j(8,:);
            ref_tempij(6,:) = 1/2*(ref_tempij(6,:) + ref_tempi1j(2,:));
            ref_tempij(8,:) = 1/2*(ref_tempij(8,:) + ref_tempi1j(12,:));

            def_tempij = Mat_of_coor_def{i,j}; def_tempijm1 = Mat_of_coor_def{i,j - 1}; def_tempim1j = Mat_of_coor_def{i - 1,j}; def_tempi1j = Mat_of_coor_def{i + 1,j};
            def_tempij(5,:) = def_tempijm1(9,:); def_tempij(3,:) = def_tempijm1(11,:);def_tempij(2,:) = def_tempim1j(6,:); def_tempij(12,:) = def_tempim1j(8,:);
            def_tempij(6,:) = 1/2*(def_tempij(6,:) + def_tempi1j(2,:));
            def_tempij(8,:) = 1/2*(def_tempij(8,:) + def_tempi1j(12,:));

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        else
            ref_tempij = Mat_of_coor_ref{i,j}; ref_tempijm1 = Mat_of_coor_ref{i,j - 1}; ref_tempim1j = Mat_of_coor_ref{i - 1,j};
            ref_tempij(5,:) = ref_tempijm1(9,:); ref_tempij(3,:) = ref_tempijm1(11,:);ref_tempij(2,:) = ref_tempim1j(6,:); ref_tempij(12,:) = ref_tempim1j(8,:);

            def_tempij = Mat_of_coor_def{i,j}; def_tempijm1 = Mat_of_coor_def{i,j - 1}; def_tempim1j = Mat_of_coor_def{i - 1,j};
            def_tempij(5,:) = def_tempijm1(9,:); def_tempij(3,:) = def_tempijm1(11,:);def_tempij(2,:) = def_tempim1j(6,:); def_tempij(12,:) = def_tempim1j(8,:);

            Mat_of_coor_ref{i,j} = ref_tempij;
            Mat_of_coor_def{i,j} = def_tempij;
        end
    end
end

matrix_plot_xr = zeros(4,4*m*n);
matrix_plot_yr = zeros(4,4*m*n);
matrix_plot_xd = zeros(4,4*m*n);
matrix_plot_yd = zeros(4,4*m*n);

for j = 1:n
    for i = 1:m
        for k = 1:4
            inindex = [1,2,3,4;4,5,6,7;7,8,9,10;10,11,12,1];
            Mrtemp = Mat_of_coor_ref{i,j}; Mdtemp = Mat_of_coor_def{i,j};
            intemp = inindex(k,:); index_1 = intemp(1);index_2 = intemp(2);index_3 = intemp(3);index_4 = intemp(4);
            matrix_plot_xr(:,((j - 1)*m + i - 1)*4 + k) = [Mrtemp(index_1,1),Mrtemp(index_2,1),Mrtemp(index_3,1),Mrtemp(index_4,1)];
            matrix_plot_yr(:,((j - 1)*m + i - 1)*4 + k) = [Mrtemp(index_1,2),Mrtemp(index_2,2),Mrtemp(index_3,2),Mrtemp(index_4,2)];
            matrix_plot_xd(:,((j - 1)*m + i - 1)*4 + k) = [Mdtemp(index_1,1),Mdtemp(index_2,1),Mdtemp(index_3,1),Mdtemp(index_4,1)];
            matrix_plot_yd(:,((j - 1)*m + i - 1)*4 + k) = [Mdtemp(index_1,2),Mdtemp(index_2,2),Mdtemp(index_3,2),Mdtemp(index_4,2)];


        end
    end
end

figure(3);hold on;
patch(matrix_plot_xr,matrix_plot_yr,[0 0 0],'LineStyle','none');axis equal;axis off;
title('First Stable State as a Square');

figure(4);hold on;
patch(matrix_plot_xd,matrix_plot_yd,[0 0 0],'LineStyle','none');axis equal;axis off;
title('Second Stable State as a Disc');








toc
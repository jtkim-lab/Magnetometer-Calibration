% magnetometer calibration
% 
% Jungtaek Kim
% jungtaek.kim@jt-inc.net
% 

function [A_i, B]=magCali(data)

    [n,m]=size(data);

    if m ~= 3
        exit;
    end
    
    D = zeros(n, 10);
    
    for i = 1:n
        D(i, 1) = data(i,1)*data(i,1);
        D(i, 2) = data(i,2)*data(i,2);
        D(i, 3) = data(i,3)*data(i,3);
        D(i, 4) = 2.0*data(i,2)*data(i,3);
        D(i, 5) = 2.0*data(i,1)*data(i,3);
        D(i, 6) = 2.0*data(i,1)*data(i,2);
        D(i, 7) = 2.0*data(i,1);
        D(i, 8) = 2.0*data(i,2);
        D(i, 9) = 2.0*data(i,3);
        D(i, 10) = 1.0;
    end
    
    S = D' * D;
    
    C = [0.0,  0.5,  0.5,  0.0,  0.0,  0.0;
         0.5,  0.0,  0.5,  0.0,  0.0,  0.0;
         0.5,  0.5,  0.0,  0.0,  0.0,  0.0;
         0.0,  0.0,  0.0,-0.25,  0.0,  0.0;
         0.0,  0.0,  0.0,  0.0,-0.25,  0.0;
         0.0,  0.0,  0.0,  0.0,  0.0,-0.25];
         
     
    CC = [-1.0,  1.0,  1.0,  0.0,  0.0,  0.0;
           1.0, -1.0,  1.0,  0.0,  0.0,  0.0;
           1.0,  1.0, -1.0,  0.0,  0.0,  0.0;
           0.0,  0.0,  0.0, -4.0,  0.0,  0.0;
           0.0,  0.0,  0.0,  0.0, -4.0,  0.0;
           0.0,  0.0,  0.0,  0.0,  0.0, -4.0];
    
    S11 = S(1:6, 1:6);
    S12 = S(1:6, 7:10);
    S21 = S(7:10, 1:6);
    S22 = S(7:10, 7:10);
    I4 = eye(4);
    
    L_S22 = chol(S22, 'lower');
    
    L_S22_i = inv(L_S22);
    L_S22_t_i = inv(L_S22');
    
    S22_i = L_S22_t_i * L_S22_i;
    
%     S22_i = pinv(S22);%I4\S22;
    
    SS = S11 - S12 * S22_i * S21;
    
    E = C * SS;
    
    [V,L] = eig(E);
    
    max = L(1, 1);
    ind = 1;
    
    for i = 1:6
        if max < L(i,i)
            max = L(i,i);
            ind = i;
        end
    end
    
    v1 = V(:, ind);
    
    if v1(1) < 0.0
        v1 = - v1;
    end
    
    norm_v1 = norm(v1);
    v1 = v1/norm_v1;
    
    v2 = S22_i * S21 * v1;
    
    v = zeros(10);
    v = complex(v);
    
    v(1) = v1(1);
    v(2) = v1(2);
    v(3) = v1(3);
    v(4) = v1(4);
    v(5) = v1(5);
    v(6) = v1(6);
    v(7) = -v2(1);
    v(8) = -v2(2);
    v(9) = -v2(3);
    v(10) = -v2(4);
    
    Q = [v(1), v(6), v(5);
         v(6), v(2), v(4);
         v(5), v(4), v(3)];
    U = [v(7); v(8); v(9)];
    
    I3 = eye(3);
    
    L_Q = chol(Q, 'lower');
    
    L_Q_i = inv(L_Q);
    L_Q_t_i = inv(L_Q');
    
    Q_i = L_Q_t_i * L_Q_i;
%     Q_i = pinv(Q);%I3\Q;
    
    B = Q_i * U;
    
    B = - B;
    QB = Q * B;
    btqb = B' * QB;
    
    J = v(10);
    hmb = sqrt(btqb - J);
    
    [VV, LL] = eig(Q);
    
    VV1 = VV(:, 1);
    VV2 = VV(:, 2);
    VV3 = VV(:, 3);
    
    norm_VV1 = norm(VV1);
    norm_VV2 = norm(VV2);
    norm_VV3 = norm(VV3);
    
    VV1 = VV1 / norm_VV1;
    VV2 = VV2 / norm_VV2;
    VV3 = VV3 / norm_VV3;
    
    VV = [VV1, VV2, VV3];
    
    vdz = VV * LL;
    
    SQ = vdz * VV';
    
    hm = 0.569;
    
    A_i = SQ * hm / hmb;
    
    C = [data(:, 1) + B(1), data(:, 2) + B(2), data(:, 3) + B(3)]';
    
    
% estimate k by calculating SE
%%%%%%%% chaning theta %%%%%%%%%%%
theta = 75/180*pi;
x_pos = [29 32 102 99];    %75
z_pos = [26 34 18 11];     %75
% x_pos = [36 28 95 104];    %60
% z_pos = [47 38 8 18];     %60
% x_pos = [35 43 85 94];    %45
% z_pos = [50 58 10 17];  %45
% x_pos = [78 40 54 89];      %30
% z_pos = [7 59 63 13];    %30
% x_pos = [89 60 44 71];     %20
% z_pos = [4 67 63 5];   %20
%%%%%%%%%%%%%
k_start = 1.5;
k_end = 2.5;
step = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = lambda;
Ns = round((k_end-k_start)/step+1);
SE_sum = zeros(1,Ns);
[aa bb] = size(v_est_z);


rotate_theta = atan(tan(theta)*1.2889); % (depth of z / pixel of z) / (width of x / pixel of x)
x_v = round(sort(cos(rotate_theta)*x_pos+sin(rotate_theta)*z_pos));
z_v = round(sort(-sin(rotate_theta)*x_pos+cos(rotate_theta)*z_pos));
x_test = [x_v(2):x_v(3)];
z_test = [z_v(2):z_v(3)];

tmp = zeros(max(size(z_test)),max(size(x_test)));
for kk = 1:Ns
    kkk = k_start+(kk-1)*step;
    for row = 1:aa
        for col = 1:bb
          [angle(row,col), bw(row,col), se(row,col)] = Doppler_angle(squeeze(store_dopplersig(row,col,:))',v_est_z_vfilter(row,col),theta,lambda,kkk);
        end
    end
    kk
    A = isnan(se);
    for zz = z_test
        for xx = x_test
            x_test_new = round(cos(rotate_theta)*xx - sin(rotate_theta)*zz);
            z_test_new = round(sin(rotate_theta)*xx + cos(rotate_theta)*zz);
            if z_test_new==0
                z_test_new = z_test_new+1;
            end
            tmp(zz-min(z_test)+1,xx-min(x_test)+1) = v_est_z(z_test_new,x_test_new);
            SE_sum(kk)
            if ~A(z_test_new,x_test_new)
                SE_sum(kk) =  SE_sum(kk)+ se(z_test_new,x_test_new);
            end
        end
    end
end

cc = find(SE_sum == min(SE_sum));
MSE = min(SE_sum)/max(size(z_test))/max(size(x_test));
est_k = (cc-1).*step+k_start;
sqrt(MSE)/pi
est_k

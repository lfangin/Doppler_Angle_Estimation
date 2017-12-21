
clear all
%close all
create_setting
%load bf_data_0_deg
%%%%%%%%% change this when changing theta %%%%%%%%%%%%%%%
load bf_data_0_deg_firing_128_75_0.05.mat % velocity = 0.1 % firing = 128
theta 	= 75/180*pi;  			% The angle between the flow direction and the z axis.
%bf_data_0_deg_firing_128 = bf_data_0_deg_firing_128;
vel_thres_min = 0.001;
energy_bound = 1e-49;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate a suitable Tprf
v0 					= 0.05; 				% Largest velocity of scatterers [m/s]
% velocity 			= v0.*within_vessel';
c = 1540;
fc = 5e6;
f_max 	= 2*v0*cos(theta)/c*fc;
PRF 	= 2*10^3;				%4*f_max;  (?????)
PRI 	= 1/PRF; 				% Time between pulse emissions [s]
firing 	= 128;					% Number of shoots   
Np = 128;

dz = 1.54/100/1000;
lambda = c/fc;
dz1 = dz*5;
%%%%%% Blood flow estimation 
%%%% Parameters for velocity flow
Vmax=lambda/4/PRI; % in m/s
window_len = round(5*lambda/(dz1)); % axial window length for summation
if mod(window_len,2) == 0 
    window_len = window_len +1;
end
overlap_window_len = (window_len-1)/2;   % overlap window length axially


%%%%%
cut_wall = 0.2;
min_cutoff_v = cut_wall*(1/PRI/2)*lambda/2;
b = fir1(20,cut_wall,'high');

a = 1;
%store_dopplersig = zeros(5,5,8);

for ii = 1 : 128    
    
    bbflow = squeeze(bf_data_0_deg_firing_128(:,ii,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% Color Doppler estimation 
    [m n] = size(bbflow);  % m: fast-time , n:slow-time, q: scan line  
    sig_len = Np;  % the number of samples used for calculating the velocity in the slow-time axis
    size_r = floor((m-window_len)/overlap_window_len)+1;   % number of steps in the fast-time direction
 
   
    for kk = 1:size_r 
    
       sv_sum = sum(bbflow((kk-1)*overlap_window_len+1 : (kk-1)*overlap_window_len + window_len,:),1)/window_len; % after wall filtering  
       %sv_sum_fil = filter(b,a,sv_sum);
       sv_sum_fil = conv2(1,b,sv_sum,'same');
       R0_wf = sum(sv_sum_fil.*conj(sv_sum_fil) )/sig_len; % power of x(t)
       wf_energy(kk,ii) = R0_wf;
       if R0_wf < energy_bound
           mv0(kk,ii) = 0;
           vel0(kk,ii) = 0;
           sv_sum_fil = zeros(1,128);
       else
            [mv0(kk,ii) vel0(kk,ii)] = Doppler_1Dauto(sv_sum_fil,c,PRI,fc);   % using 1-D autocorrelation method  (in units of m/s)
       end
       store_dopplersig(kk,ii,:) = sv_sum_fil;

       
%        %%% Pre-clutter-filtering
%        doppler_sig_pre = sv_sum;
%        R0_nwf = sum(doppler_sig_pre.*conj(doppler_sig_pre) )/sig_len; % power of x(t)
%        nwf_energy(kk,ii) = R0_nwf; 
%       
%        %%% Post-clutter-filtering
%        doppler_sig = sv_sum_fil; % after wall filtering
%        R0_wf = sum(doppler_sig.*conj(doppler_sig) )/sig_len; % power of x(t)
%        wf_energy(kk,ii) = R0_wf; 
       
       %%% Autocorrelation function
       
    end

end

v_est_z = medfilt2(vel0(:,:,1),[3,3]);
[aa bb] = size(v_est_z);
v_larger = abs(v_est_z) >= vel_thres_min;
v_est_z_vfilter = v_larger.*v_est_z;
%find(abs(v_est_z)< vel_thres_min)
store_dopplersig = repmat(v_larger,[1,1,128]).*store_dopplersig;

%
kkk = 2.0;
for row = 1:aa
    for col = 1:bb
      [angle(row,col), bw(row,col), se(row,col)] = Doppler_angle(squeeze(store_dopplersig(row,col,:))',v_est_z_vfilter(row,col),theta,lambda,kkk);
    end
end

%%%%%%% figures %%%%%%%%%%%%%%%%%%%
aaaa    =  max(size(bbflow));
%Rmin    = toff(1)*c/2; % minimum range to be imaged
Rmin    = 1.0000e-05*c/2;
dz2     = 1/(fc*4)*c/2;
Rmax    = Rmin + aaaa*dz2;
range   = Rmin:dz2:Rmax-dz2; 
x_axis_label = [-(128-1)/2*lambda:lambda : (128-1)/2*lambda]*1000;
%imagesc(x_axis_label,range*1000,abs(mv_interp1))

% 
% figure;
% imagesc(x_axis_label,range*1000,abs(bw))
% colormap(gray)
% colorbar
% xlabel('Azimuth (mm)')
% ylabel('Range (mm)')
% title('Bandwidth')

% figure;
% imagesc(x_axis_label,range*1000,abs(k))
% colormap(gray)
% colorbar
% xlabel('Azimuth (mm)')
% ylabel('Range (mm)')
% title('Estimated k')

figure;
imagesc(abs(v_est_z_vfilter))
%colormap(gray)
colorbar
xlabel('Azimuth (mm)')
ylabel('Range (mm)')
title('Velocity')

figure;
imagesc(abs(bw))
%colormap(gray)
colorbar
xlabel('Azimuth (mm)')
ylabel('Range (mm)')
title('bw')



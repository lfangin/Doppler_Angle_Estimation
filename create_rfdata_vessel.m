

clear all


eval('load /Users/jacquelinelin/NTUEE_106_1/US_project/dopper/data/setting/setting.mat')


    
%%%%%%%%% Flow data setting %%%%%%%%%%%%%%%%%%%%

% Set the seed of the random number generator
randn('seed',sum(100*clock))

% Initialize the ranges for the scatterers
% Notice that the coordinates are in meters

x_range = 0.025; % x range for the scatterers [m]
y_range = 0.005; % y range for the scatterers [m]
%y_range = 0;   %for test
z_range = 0.010; % z range for the scatterers [m]
z_offset = 0.020; % Offset of the mid-point of the scatterers [m]

R = 0.003; % Radius of blood vessel [m]

% Set the number of scatterers. It should be roughly
% 10 scatterers per resolution cell
N=round(10*x_range/(5*lambda)*y_range/(5*lambda)*z_range/(lambda*2)); %8556
%N = 40000;
%N = 20; %for test
disp([num2str(N), ' Scatterers'])

% Generate the coordinates and amplitude
% Coordinates are rectangular within the range.
% The amplitude has a Gaussian distribution.
x = x_range*(rand(1,N)-0.5);
y = y_range*(rand(1,N)-0.5);
z = z_range*(rand(1,N)-0.5);

% Find which scatterers that lie within the blood vessel
r 				= (y.^2+z.^2).^0.5;
within_vessel 	= (r < R)';

% Assign an amplitude and a velocity for each scatterer
v0 					= 0.05; 				% Largest velocity of scatterers [m/s]
velocity 			= v0.*within_vessel';
blood_to_stationary = 0.1; 				% Ratio between amplitude of blood to stationary tissue   (-20dB)
amp 				= randn(N,1).*((1-within_vessel) + within_vessel*blood_to_stationary);

aaaaa = [75];

for ttttt = aaaaa
    % Calculate a suitable Tprf
    theta 	= ttttt/180*pi;  			% The angle between the flow direction and the z axis.
    f_max 	= 2*v0*cos(theta)/c*fc;
    PRF 	= 2*10^3;				%   4*f_max;  (?????)
    PRI 	= 1/PRF; 				% Time between pulse emissions [s]
    firing 	= 128;					% Number of shoots 
    tilted_angle = 0;
    thetas_tx  = tilted_angle*pi/180;
    thetas_rx = thetas_tx;

    % for pp = 1 : length(tilted_angel)

    % Find the response by calling field
    % for ii = 1:SL


    % xdc_center_focus (rx, [0 0 0]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Tx Beamforming
    %     sin_th = sin(-1*max_th/180*pi) + (ii-1)*d_th;
    %     th = asin(sin_th);

        focus_tx = [100000000*sin(thetas_tx) 0 100000000*cos(thetas_tx)];
    % 
        xdc_focus(tx, 0, focus_tx);
        %%% Tx apodization (No apodization)
        %xdc_apodization (xmit_aperture, 0, apo_vector);

        %recordXZ = [];  %for test
        SL_data = [];
        for jj = 1:firing
          %%%%%% Generate the rotated and offset block of sample
          jj
          xnew = x*sin(theta) + z*cos(theta);
          znew = z*sin(theta) - x*cos(theta) + z_offset;
          scatterers = [xnew; y; znew;]' ;

          %recordXZ(1:2,jj) = [x(1),z(1)] ; %for test

    %       scatter(x,z);
    %       if jj~= firing
    %         hold on;
    %       end
          % Calculate the received response
          [rf_data, toff(jj)] = calc_scat_multi(tx,rx, scatterers, amp);
          SL_data(1:max(size(rf_data)),:,jj) = rf_data;


          % Propagate the scatterers and alias them to lie within the correct range
           x = x + velocity*PRI;
           outside_range = (x > x_range/2);
           x = x - x_range*outside_range;


        end


         cmd=['save /Users/jacquelinelin/NTUEE_106_1/US_project/dopper/data/phantom/vessel/CFM_uBF',num2str(ttttt),'_degree_firing_',num2str(firing),'_',num2str(v0),'.mat SL_data toff'];
         eval(cmd);
end

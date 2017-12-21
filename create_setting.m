
path(path, '/Users/jacquelinelin/Documents/MATLAB/Field_II_ver_3_24_mac')
%path('D:\�ի�\Vector_flow\fieldii\fieldiiwindows')

% set initial parameters
% unit: Hz, m, rad, s, rad
c                   = 1540;             % sound velocity
fc                  = 5e6;              % center frequency
fs                  = 100e6;            % sampling frequency
nelex               = 128;              % number of elements in x direction, min = 2
neley               = 1;                % number of elements in y direction, min = 1
nsample             = 4096;             % number of samples
offset              = 0e-3;             % offset in z direction

lambda              = c/fc;
ele_size_x          = lambda*0.95;      % element size in x direction
if neley == 1                           % element size in y direction
    ele_size_y      = 5e-3;
else
    ele_size_y      = lambda*0.95;  
end
kerfx               = ele_size_x/20;    % kerf in x direction
kerfy               = ele_size_y/20;    % kerf in y direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



% parameters don't need to be changed
npadding            = nsample;                              % number of zero padding
nbeamx              = nelex;                                % number of beams in x direction
nbeamy              = neley;                                % number of beams in y direction

dft                 = fs/(nsample + npadding);
dt                  = 1/fs;                                 % delta t in fast time direction
dz                  = c*dt/2;                               % delta z in fast time direction

offset_index        = floor( offset / (dz) );

pitch_x             = ele_size_x + kerfx;                   % pitch in x direction
pitch_y             = ele_size_y + kerfy;                   % pitch in y direction
transducer_size_x   = nelex*ele_size_x + (nelex-1)*kerfx;   % transducer siez in x direction
transducer_size_y   = neley*ele_size_y + (neley-1)*kerfy;   % transducer size in y direction

depth               = nsample*dz; 

ex                  = -(nelex-1)*pitch_x/2 + pitch_x*(0:nelex-1); % x value of each element
ey                  = -(neley-1)*pitch_y/2 + pitch_y*(0:neley-1); % y value of each element

zmid                = offset + 20/1000; 
fnumber             = zmid/transducer_size_x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



% set steering angle: thetas & phis
switch 0
    case 0
        thetas_tx   = [0]*pi/180;        
        phis_tx     = [0]*pi/180;
        thetas_rx   = [0]*pi/180;
        phis_rx     = [0]*pi/180;    
    case 1
        % same focusing quality as the optimal multifocus
        delta_angle     = lambda/transducer_size_x;
        n_compounding   = ceil(transducer_size_x/lambda/fnumber);
        steering_range  = delta_angle*n_compounding*180/pi;
        for i = ceil(-n_compounding/2) : ceil(n_compounding/2)-1
            try
                thetas_tx = [thetas_tx delta_angle*i];
            catch
                thetas_tx = delta_angle*i;
            end
        end
        thetas_rx=[0]*pi/180;
        phis_tx = [0]*pi/180;
        phis_rx = [0]*pi/180;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% field II function
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path(path,'fieldii/fieldiiwindows');
    field_init(0);
% set the sampling frequency
set_sampling(fs);

% set the transducer
if neley==1
    tx      = xdc_linear_array(nelex, ele_size_x, ele_size_y, kerfx, 1, 10, [0 0 1000000000000000000000000]);    
    rx      = xdc_linear_array(nelex, ele_size_x, ele_size_y, kerfx, 1, 10, [0 0 1000000000000000000000000]);
else
    enabled = ones(nelex, neley);
    tx      = xdc_2d_array(nelex, neley, ele_size_x, ele_size_y, kerfx, kerfy, enabled, 1, 1, [0 0 0]);
    rx      = xdc_2d_array(nelex, neley, ele_size_x, ele_size_y, kerfx, kerfy, enabled, 1, 1, [0 0 0]);
end

% set the impulse response
impulse_response = sin(2*pi*fc*(0:1/fs:2/fc));
impulse_response = impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse(tx, impulse_response);
xdc_impulse(rx, impulse_response);

% set the excitation
excitation = sin(2*pi*fc*(0:1/fs:2/fc));
xdc_excitation(tx, excitation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Store the setting
cmd = 'save /Users/jacquelinelin/NTUEE_106_1/US_project/dopper/data/setting/setting.mat';
cmd = [cmd, ' c fc fs nbeamx nbeamy nelex neley nsample npadding'];
cmd = [cmd, ' dft dt dz ele_size_x ele_size_y kerfx kerfy lambda offset offset_index pitch_x pitch_y transducer_size_y transducer_size_x'];
cmd = [cmd, ' ex ey depth zmid'];
cmd = [cmd, ' tx rx impulse_response excitation'];
cmd = [cmd, ' thetas_tx thetas_rx phis_tx phis_rx'];
eval(cmd);


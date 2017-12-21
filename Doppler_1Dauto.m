function [freq_shift vel] = Doppler_1Dauto(bfmot_sum,c,PRI,fo);

%%% in unit of m/s
%%% c: sound velocity
%%% PRI: pulse repetition interval (seconds)
%%% fo: center frequency

        I_comp = real(bfmot_sum);
        Q_comp = imag(bfmot_sum);
        I1Q = sum([I_comp 0].*[0 Q_comp]);   % correlation between I(n+1)*Q(n)
        IQ1 = sum([Q_comp 0].*[0 I_comp]);   %                     I(n)  *Q(n+1)
        I1I = sum([I_comp 0].*[0 I_comp]);   %                     I(n+1)*I(n)
        Q1Q = sum([Q_comp 0].*[0 Q_comp]);   %                     Q(n+1)*Q(n)
        
        corr_theta1 = atan2((I1Q-IQ1),(I1I+Q1Q));
        freq_shift  = 2*corr_theta1/(2*pi*PRI);
        vel= c/4/pi/fo/PRI*corr_theta1;  % in unit of m/s
%         freq_shift = vel/c*fo;


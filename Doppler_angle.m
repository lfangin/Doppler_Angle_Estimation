function [angle, bw, se] = Doppler_angle(bfmot_sum,vel,theta,lambda,k)
% vel = v_est_z(kk,ii)
% k = constant
% angle = the estimated flow angle // 
% R(T) = 1/(N-1)sum(conj(S(+1))*S())
% R(0) = 1/2(N-1)sum(|S()|^2+|S(+1)|^2)
        %store_sig = bfmot_sum;
        bfmot_sum = bfmot_sum(5:120);
        
        Ns = max(size(bfmot_sum)); 
        w = lambda;
        PRI =    5.0000e-04;
        I_comp = real(bfmot_sum);
        Q_comp = imag(bfmot_sum);
        I1Q = sum([I_comp 0].*[0 Q_comp]);   % correlation between I(n+1)*Q(n)
        IQ1 = sum([Q_comp 0].*[0 I_comp]);   %                     I(n)  *Q(n+1)
        I1I = sum([I_comp 0].*[0 I_comp]);   %                     I(n+1)*I(n)
        Q1Q = sum([Q_comp 0].*[0 Q_comp]);   %                     Q(n+1)*Q(n)
        
        S_square = abs(bfmot_sum).^2;
        
        R_0 = (2*sum(S_square)-S_square(1)-S_square(Ns))/(2*(Ns-1));        
        R_T = complex(I1I+Q1Q,IQ1-I1Q)/(Ns-1);
        
        sigma_s = 2/(PRI^2)*(1-(abs(R_T)/R_0));     % sigma^2 = 2/(PRI)^2*(1-|(R(T)|/R(0)))
        bw = sqrt(sigma_s);
        angle = atan(w*bw/vel/k);
        se = (abs(theta)-abs(angle))^2;
        
        
        

Fs = PRF;
%y = drawsig;
zstart = 19;
xstart = 64;

figure;
for i = 1:5
    y = squeeze(store_dopplersig(zstart+i,xstart,5:120))';
    NFFT = length(y);
    Y = fftshift(fft(y,NFFT));
    F = ((-0.5:1/NFFT:0.5-1/NFFT)*Fs).';
    magnitudeY = 20*log10(abs(Y));
    %phaseY = unwrap(imag(log(Y)));
    
    subplot(1,5,i)
    plot(F,repelem([max(magnitudeY)-20],NFFT),F,magnitudeY);
    xlabel('Frequency (in hertz)');
    ylabel('Amplitude (in dB)');
    title(['Magnitude Response of (',num2str(zstart+i),' ', num2str(xstart) ')']);
end


slow_time_x = 5:120;
slow_time_magn = squeeze((real(store_dopplersig(20,xstart,slow_time_x))));
figure;plot(slow_time_x,slow_time_magn);
function H = filterV(V)
x=1:numel(V); clc; close all

d = fdesign.lowpass('Fp,Fst',0.001,.15);
H = design(d, 'equiripple');  % Design an FIR equiripple filter
H.Arithmetic = class(V);
%fvtool(H); info(H)
%V(isnan(V))=0; fig; plot(x,V,'Display','data'); plot(x,filter(H,V),'Display','FIR'); xyzlabel('t (sample)','amplitude (mV)'); fcntight()


%d = fdesign.lowpass('N,F3dB',20,50/1024*2); SAME AS fdesign.lowpass('N,F3dB',20,50,1024);
%d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.10,0.15,3,60);
%designmethods(d)
%H = design(d,'cheby2');
%fvtool(H); info(H)
%V(isnan(V))=0;  plot(x,filter(H,V),'Display','IIR'); legend show; fcnlinewidth(2)


%fcnfft(double(V'),x,1)

%F=lowPassFilter(.10,.15,.5); plot(x,filter(F,V))


%Create a Butterworth bandstop filter for data sampled at 10 kHz. The stopband is [1,1.5] kHz. The order of the filter is 20.
%d = fdesign.bandstop('N,F3dB1,F3dB2',20,1e3,1.5e3,1e4);

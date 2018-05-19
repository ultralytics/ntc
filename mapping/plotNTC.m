function [] = plotNTC(input0,output)

i = MTCSRCCH2pixelID(NTCIRSDmap(1:128,'SRCCH'));
Cdata = output(1).A(i);
AlphaData = ones(size(Cdata))*.1;  

i=Cdata>1; AlphaData(i)=Cdata(i)/max3(fcnsigmarejection(Cdata(i),1.5,3));
%AlphaData(i) = 0.9;

input=[];
load NTCIRSDinputcube.mat
fig(1,1,3,2.4); fcnPlotDetector(input,Cdata,AlphaData);
xyzlabel('','','',sprintf('%s event %.0f',str_(input0.MTC.A.filename),input0.MTC.A.events(input0.eventNumber)));
h=colorbar('East'); h.Label.String='V (samples)';

 
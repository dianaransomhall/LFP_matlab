






% load the needed variables
load 'F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\11_19_2013_all.mat' ;


EPA=0 ;

if EPA
    % EPA Computer
    addpath(genpath('\\AA.AD.EPA.GOV\ORD\RTP\USERS\A-D\dhall05\Net MyDocuments\MATLAB\NeuroShare')) ;
else
    %personal computer
    addpath(genpath('C:\Users\Diana\Documents\MATLAB\Cina_Herr_Diana'));
end

%set library
%personal computer
if EPA
    ns_SetLibrary('\\AA.AD.EPA.GOV\ORD\RTP\USERS\A-D\dhall05\Net MyDocuments\MATLAB\find-1.1\nsMCDLibrary.dll')
else
    ns_SetLibrary('C:\Users\Diana\Documents\MATLAB\SpikeSorting\find-1.1\nsMCDLibrary64.dll')
end




if EPA
    % EPA computer
file1='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\September 4, 2013  BIC\090413 20607 BIC FPSPK0001.mcd';
file2='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\September5, 2013 BIC  CAR\September 5, 2013 BIC\090513 8395 BIC\090513 8395 BIC FPSPK0001.mcd'
file3='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\September5, 2013 BIC  CAR\September 5, 2013 BIC\090513 15152 BIC\090513 15152 BIC FPSPK0001.mcd'
file4='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\September5, 2013 BIC  CAR\September 5, 2013 BIC\090513 18405 BIC\090513 18405 BIC FPSPK0002.mcd'
else
% personal computer
 file1='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September 4, 2013 #20607 BIC\090413 20607 BIC FPSPK0001.mcd';
 file2='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 8395 BIC\090513 8395 BIC FPSPK0001.mcd'
 file3='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 15152 BIC\090513 15152 BIC FPSPK0001.mcd'
 file4='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 18405 BIC\090513 18405 BIC FPSPK0002.mcd'
end

files.bic = {file1, file2, file3, file4 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

curFile = 1;  ent = 1;

% get file
file = files.bic{curFile} ;
    
[nsresult, hfile]=ns_OpenFile(file)
[nsresult,info]=ns_GetFileInfo(hfile)


 analogEnt = ent +60 ;
        
[nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
%[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
[nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);


size(data) % 7632500

Fs = 25000; %sampling frequency /s      

% create spectrum object, segment length is 10,000, overlap = 50%
h = spectrum.welch;
hopts = psdopts(h, data )
set(hopts, 'Fs', Fs, 'SpectrumType', 'Onesided', 'ConfLevel', 0.95 ) ;

Hpsd = psd(h, data ) ;

plot(Hpsd)





















% generate_derived_LFP_data.m
% diana hall
% november 5, 2013
%%%Use NeuroShare
%you need the FIND toolbox and the NeuroShare

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


ns_GetLibraryInfo


% Ctrl+Z to undo
%file='F:\Hall\09-2013 BIC & CAR\September 4, 2013 #20607 BIC\090413 20607 BIC FPSPK.mcd';

%[nsresult, hfile]=ns_OpenFile(file)
%[nsresult,info]=ns_GetFileInfo(hfile)

% Find out about the Entity types
% Then read specific entity info and data

%mcs_Info(hfile)

%stream       type      infotext        n       TimeSpan 
%     elec0001      2     Analog Entity      60        301.600000  
%     spks0001      3     Analog Entity      60        301.600000
% 'spks0001 0059 0059   57', is the 59th channel for (60 total)
% 57 is name of channel in matrix notation
% "spks0001 0000 0000" - "spks0001 0059 0059"
% spks is a "segment" type in file 1-60, while "elec" are "analog entities"
% and are 61-120





%%%%%%%%%%%%%%%+++++  BIC   ++++++++++++%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%BIC%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files

file1='F:\EXTRAP\FP SPIKE\BIC\090513 8395 BIC\090513 8395 BIC FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\BIC\090513 15152 BIC\090513 15152 BIC FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\BIC\090513 18405 BIC\090513 18405 BIC FPSPK0002.mcd';
file4='F:\EXTRAP\FP SPIKE\BIC\December 18, 2013 18395 BIC\121813 18395 BIC FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\BIC\November 19, 2013 BIC\111913 9858 BIC FPSPK0001.mcd';

bic.filepath = {file1, file2, file3, file4, file5 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
 % get params for FFT
 fs=25000; %sampling rate
 L = 255*fs ;% take 255 s of data = 255*25000
 NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
 w_kaiser = kaiser(L, 0.5); 
 f = fs/2*linspace(0,1,NFFT/2+1); 
 index_50Hz = find(f>=50,1);

for curFile=1:length(bic.filepath)    
    file = bic.filepath{curFile} ;   
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    temp = mcs_Info(hfile);   

    % loop through the ent's with sufficient high frequency activity
    for ent= 1:60
        %analog ent is 60 more than current
        analogEnt = ent +60 ;
        [nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);

            curData_temp=data(1:L);

            %apply kaiser window to attenuate signal
            curData = w_kaiser.*curData_temp ; %element-wise operations

            X = fft(curData, NFFT)/L;

            % find what index of f>=50 (Hz)
            X_mag = abs(X);
            % X_mag = real(X);

            bic.file(curFile).channel(ent).psd = X_mag(1:index_50Hz);
        
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.bic)
    




% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_bic.mat';
% save(filename, 'bic' ) ;















%%%%%%%%%%%%%%%+++++  PERM   ++++++++++++%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files
file1='F:\EXTRAP\FP SPIKE\PER\11-19-2013 PERM\112113 18413 PERM FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM  19854\121813 19854 PERM FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 19398\121813 19398 PERM FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 20608\121813 20608 PERM FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\PER\13079\111913 13079 PERM FPSPK0001.mcd';

per.filepath = {file1, file2, file3, file4, file5 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
 % get params for FFT
 fs=25000; %sampling rate
 L = 255*fs ;% take 255 s of data = 255*25000
 NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
 w_kaiser = kaiser(L, 0.5); 
 f = fs/2*linspace(0,1,NFFT/2+1); 
 index_50Hz = find(f>=50,1);

for curFile=1:length(per.filepath)    
    file = per.filepath{curFile} ;   
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    temp = mcs_Info(hfile);   

    % loop through the ent's with sufficient high frequency activity
    for ent= 1:60
        %analog ent is 60 more than current
        analogEnt = ent +60 ;
        [nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);

            curData_temp=data(1:L);

            %apply kaiser window to attenuate signal
            curData = w_kaiser.*curData_temp ; %element-wise operations

            X = fft(curData, NFFT)/L;

            % find what index of f>=50 (Hz)
            X_mag = abs(X);
            % X_mag = real(X);

            per.file(curFile).channel(ent).psd = X_mag(1:index_50Hz);
        
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.bic)
    


% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_29_2014_per.mat';
% save(filename, 'per' ) ;

















%%%%%%%++++   CAR  30uM   +++++++++++++
 
%%%%%%%%%%%%%%%%%%%CAR%%%%%%%%%%%%%%%%%%%%%%%%%%%
%personal computer

file1='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19855\112113 19855 CAR FPSPK0002.mcd';
file2='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\20600\112113 20600 CAR FPSPK0002.mcd';
file3='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 18412\121713 18412 CAR FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 19858\121713 19858 CAR FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 20599\121713 20599 CAR FPSPK0001.mcd';

% 30uM CAR
car.filepath={file1, file2, file3, file4, file5} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
 % get params for FFT
 fs=25000; %sampling rate
 L = 255*fs ;% take 255 s of data = 255*25000
 NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
 w_kaiser = kaiser(L, 0.5); 
 f = fs/2*linspace(0,1,NFFT/2+1); 
 index_50Hz = find(f>=50,1);
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_car.mat');
 
for curFile=1:length(car.filepath)    
    file = car.filepath{curFile} ;   
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    temp = mcs_Info(hfile);   

    % loop through the ent's with sufficient high frequency activity
    for ent= 1:60
        %analog ent is 60 more than current
        analogEnt = ent +60 ;
        [nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);

            curData_temp=data(1:L);

            %apply kaiser window to attenuate signal
            curData = w_kaiser.*curData_temp ; %element-wise operations

            X = fft(curData, NFFT)/L;

            % find what index of f>=50 (Hz)
            X_mag = abs(X);
            % X_mag = real(X);

            car.file(curFile).channel(ent).psd = X_mag(1:index_50Hz);
        
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.bic)

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_29_2014_car.mat';
% save(filename, 'car' ) ;













%%%%%%%%%%%%%%%%%%%   CONTROL            %%%%%%%%%%%%%%%%%%%%%%%%%%%
 file1='F:\EXTRAP\FP SPIKE\BIC\090513 8395 BIC\090513 8395 BIC FPSPK.mcd';
 file2='F:\EXTRAP\FP SPIKE\BIC\090513 15152 BIC\090513 15152 BIC FPSPK.mcd';
 file3='F:\EXTRAP\FP SPIKE\BIC\090513 18405 BIC\090513 18405 BIC FPSPK.mcd';
 file4='F:\EXTRAP\FP SPIKE\BIC\December 18, 2013 18395 BIC\121813 18395 BIC FPSPK.mcd';
 file5='F:\EXTRAP\FP SPIKE\BIC\November 19, 2013 BIC\111913 9858 BIC FPSPK.mcd';
 % perm files
 file6='F:\EXTRAP\FP SPIKE\PER\11-19-2013 PERM\112113 18413 PERM FPSPK.mcd';
 file7='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM  19854\121813 19854 PERM FPSPK.mcd';
 file8='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 19398\121813 19398 PERM FPSPK.mcd';
 file9='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 20608\121813 20608 PERM FPSPK.mcd';
 file10='F:\EXTRAP\FP SPIKE\PER\13079\111913 13079 PERM FPSPK.mcd';
 % car files
 file11='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19855\112113 19855 CAR FPSPK.mcd';
 file12='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\20600\112113 20600 CAR FPSPK.mcd';
 file13='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 18412\121713 18412 CAR FPSPK.mcd';
 file14='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 19858\121713 19858 CAR FPSPK.mcd';
 file15='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 20599\121713 20599 CAR FPSPK.mcd';
 
% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4, file5, file6, file7, file8,...
    file9, file10, file11, file12, file13, file14, file15 } ;
carcon.filepath={file11,file12,file13,file14,file15} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
 % get params for FFT
 fs=25000; %sampling rate
 L = 255*fs ;% take 255 s of data = 255*25000
 NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
 w_kaiser = kaiser(L, 0.5); 
 f = fs/2*linspace(0,1,NFFT/2+1); 
 index_50Hz = find(f>=50,1);
 % load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_28_2014_con.mat');

for curFile=1:length(carcon.filepath)    
    file = con.filepath{curFile} ;   
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    temp = mcs_Info(hfile);   

    % loop through the ent's with sufficient high frequency activity
    for ent= 1:60
        %analog ent is 60 more than current
        analogEnt = ent +60 ;
        [nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);

            curData_temp=data(1:L);

            %apply kaiser window to attenuate signal
            curData = w_kaiser.*curData_temp ; %element-wise operations

            X = fft(curData, NFFT)/L;

            % find what index of f>=50 (Hz)
            X_mag = abs(X);
            % X_mag = real(X);

            carcon.file(curFile).channel(ent).psd = X_mag(1:index_50Hz);
        
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.bic)

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_29_2014_carcon.mat';
% save(filename, 'carcon' ) ;


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

for curFile=1:length(files.bic)
    
    file = files.bic{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    bic.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot

    spikes.n=0 ;
    for ent = 1:60

        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;

        bic.file(curFile).channel(ent).nspikes = entity.ItemCount ;

        % indicator of ae or not
        bic.file(curFile).channel(ent).ae =...
            (bic.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount

            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            bic.file(curFile).channel(ent).train(nspike) = timestamp ;
            
        end; %end of loop through timestamp in given entity

    end; %end of loop through all entities

    L = 25000; %length of 1 second signal
    %make L-point Kaiser window
    %beta=0.5, default value accoding to matlab documentation
    w_kaiser = kaiser(L, 0.5); 
    

    % loop through the ent's with sufficient high frequency activity
    for ent= 1:60
        %analog ent is 60 more than current
        analogEnt = ent +60 ;
        
        [nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);

        for sec=1:nseconds

            %take a chunk of length 1 second
            fs = 25000;
            curData_temp=data((sec-1)*fs+1:sec*fs);

            %apply kaiser window to attenuate signal
            curData = w_kaiser.*curData_temp ; %element-wise operations
            
            
            L = length(curData) ; %7,540,000
            NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
            

            X = fft(curData, NFFT)/L;

            f = fs/2*linspace(0,1,NFFT/2+1); 

            X_mag = abs(X);
            % X_mag = real(X);
            
            bic.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);


        end; %end of loop through sec in an entity
        

        %for band=1:20
        for sec = 1:nseconds
        
            tempA = cell(1,nbands ) ; %used for bands

             for band=1:nbands
               %get indices of psd corresponding to current freq band
             index_t=( bands_needed(band) <= f(1:200) & f(1:200) <= bands_needed(band+1)) ;

                % get the particular band power at each band
                % this is a column of power for freq. from band 
                % defined by the column
                
                tempA{band} = bic.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = bic(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                bic.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( bic.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %bic(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( bic(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            bic.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %bic(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     

    end;%for loop through entitites
    
end;  % for curFile=1:length(files.bic)
    




for curFile=1:length(files.bic)
    
    % start a matrix
    bandM = zeros(bic.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:bic.file(curFile).nseconds
                bandM(sec, ent, band) = bic.file(curFile).channel(ent).sec(sec).avpw(band) ;
                %bandM(sec,entCol,band) = bic(sec).file(curFile).ent(ent).avpw(band);
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    bic.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;

    
end; %loop through all files


% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(bic.file(curFile).nseconds/60)*60 ;

bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
fs = 25000; %sampling rate        
L=fs ; %length of 1 sec worth of data
NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
f = fs/2*linspace(0,1,NFFT/2+1); %hertz
%index corresponding to freq between 1 and 50 Hz 
index_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
% f(index_freq) becomes actual 


 % get power for specific bands
f_temp = cell(1,5) ;
for band=1:length(bands_needed)-1
  % for each band get the corresponding frequencies
  index_freq = ( f >=bands_needed(band)& f <=bands_needed(band+1) ) ;

  f_temp{band} = f( index_freq ) ;
end;
        
bic.fbands = f_temp ;

%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 

for curFile=1:length(files.bic)

    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(bic.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd

           temp(curSec, : ) = transpose( bic.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;

        end;
        bic.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        bic.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(bic.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd

                temp(curSec, : ) = transpose( bic.file(curFile).channel(ent).sec(curSec).pw{band} ) ;

            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            bic.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files





% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\11_19_2013_bic.mat';
% save(filename, 'bic', 'files' ) ;

















%%%%%%%++++   CAR     +++++++++++++
 
%%%%%%%%%%%%%%%%%%%BIC%%%%%%%%%%%%%%%%%%%%%%%%%%%
%personal computer

if ~EPA
 % lower concentration of CAR deemed too low to show effect
 % file1='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK0001.mcd';
 % file2='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK0001.mcd';
 % file3='F:\EXTRAP\FP SPIKE\September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK0001.mcd' ;
 % file4='F:\EXTRAP\FP SPIKE\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK0001.mcd';
 % file5='F:\EXTRAP\FP SPIKE\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK0002.mcd' ;
 file1='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\19855\112113 19855 CAR FPSPK0001.mcd' ;%15 conc
 file2='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\19857\112113 19857 CAR FPSPK0001.mcd' ;%15 conc
 file3='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\20600\112113 20600 CAR FPSPK0001.mcd' ;%15 conc
 file4='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\19855\112113 19855 CAR FPSPK0002.mcd' ;%30 conc
 file5='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\20600\112113 20600 CAR FPSPK0002.mcd' ;%30 conc
else 
 %file1='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK0001.mcd';
 %file2='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK0001.mcd';
 %file3='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK0001.mcd';
 %file4='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK0001.mcd';
 %file5='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK0002.mcd';
end;


% file5 is not working for some reason
% files.car={file1, file2, file3, file4} ;
files.car={file1, file2, file3, file4, file5} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(files.car)
    
    file = files.car{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    car.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot

    spikes.n=0 ;
    for ent = 1:60

        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;

        car.file(curFile).channel(ent).nspikes = entity.ItemCount ;

        % indicator of ae or not
        car.file(curFile).channel(ent).ae =...
            (car.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount

            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            car.file(curFile).channel(ent).train(nspike) = timestamp ;
            
        end; %end of loop through timestamp in given entity

    end; %end of loop through all entities
    
    
    L = 25000; %length of 1 second signal
    %make L-point Kaiser window
    %beta=0.5, default value accoding to matlab documentation
    w_kaiser = kaiser(L, 0.5); 


    % loop through the ent's with sufficient high frequency activity
    for ent= 1:60
        %analog ent is 60 more than current
        analogEnt = ent +60 ;
        
        [nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);

        for sec=1:nseconds

            
            %take a chunk of length 1 second
            fs = 25000;
            curData_temp=data((sec-1)*fs+1:sec*fs);
      
            %apply kaiser window to attenuate signal
            curData = w_kaiser.*curData_temp ; %element-wise operations
            
            
            L = length(curData) ; %7,540,000
            NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
      

            X = fft(curData, NFFT)/L;

            f = fs/2*linspace(0,1,NFFT/2+1); 

            X_mag = abs(X);
            % X_mag = real(X); %just get real part

            car.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);


        end; %end of loop through sec in an entity
        

        %for band=1:20
        for sec = 1:nseconds
        
            tempA = cell(1,nbands ) ; %used for bands

             for band=1:nbands
               %get indices of psd corresponding to current freq band
             index_t=( bands_needed(band) <= f(1:200) & f(1:200) <= bands_needed(band+1)) ;

                % get the particular band power at each band
                % this is a column of power for freq. from band 
                % defined by the column
                
                tempA{band} = car.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = bic(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                car.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( car.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %bic(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( bic(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            car.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %bic(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     

    end;%for loop through entitites
    
end;  % for curFile=1:length(files.bic)
% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\Matlab_workspace\11_05_2013_vars';
% save(filename, 'spikes', 'bic', 'car', 'files', 'apwBIC', 'apwCAR' )



for curFile=1:length(files.car)
    
    % start a matrix
    bandM = zeros( car.file(curFile).nseconds , 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:car.file(curFile).nseconds
                bandM(sec, ent, band) = car.file(curFile).channel(ent).sec(sec).avpw(band) ;
                %bandM(sec,entCol,band) = bic(sec).file(curFile).ent(ent).avpw(band);
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    car.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;

    
end; %loop through all files



% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
nbands = length( bands_needed) - 1 ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(car.file(curFile).nseconds/60)*60 ;

bands_needed=[1 4 8 14 30 50];

fs = 25000; %sampling rate        
L=fs ; %length of 1 sec worth of data
NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
f = fs/2*linspace(0,1,NFFT/2+1); %hertz
%index corresponding to freq between 1 and 50 Hz 
index_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
% f(index_freq) becomes actual 


 % get power for specific bands
f_temp = cell(1,5) ;
for band=1:length(bands_needed)-1
  % for each band get the corresponding frequencies
  index_freq = ( f >=bands_needed(band)& f <=bands_needed(band+1) ) ;

  f_temp{band} = f( index_freq ) ;
end;
        
car.fbands = f_temp ;

%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 

for curFile=1:length(files.car)
    for ent = 1:60
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(car.file(curFile).channel(ent).sec(1).psd(index_all_freq ) ) ) ;
        
        for curSec=1:nsecondsRd

           temp(curSec, : ) = transpose( car.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;

        end;
        car.file(curFile).channel(ent).f = f(index_all_freq ) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        car.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(car.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd

                temp(curSec, : ) = transpose( car.file(curFile).channel(ent).sec(curSec).pw{band} ) ;

            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            car.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_car.mat';
% save(filename, 'car', 'files' ) ;













%%%%%%%%%%%%%%%%%%%   CONTROL            %%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~EPA
 file1='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September 4, 2013 #20607 BIC\090413 20607 BIC FPSPK.mcd';
 file2='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 8395 BIC\090513 8395 BIC FPSPK.mcd';
 file3='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 15152 BIC\090513 15152 BIC FPSPK.mcd' ;
 file4='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 18405 BIC\090513 18405 BIC FPSPK.mcd';
 %file5='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK.mcd';
 %file6='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK.mcd';
 %file7='F:\EXTRAP\FP SPIKE\September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK.mcd' ;
 %file8='F:\EXTRAP\FP SPIKE\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK.mcd' ;
 %file9='F:\EXTRAP\FP SPIKE\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK.mcd' ;
 file5='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\19855\112113 19855 CAR FPSPK.mcd' ;%15 conc
 file6='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\19857\112113 19857 CAR FPSPK.mcd' ;%15 conc
 file7='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\20600\112113 20600 CAR FPSPK.mcd' ;%15 conc
 file8='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\19855\112113 19855 CAR FPSPK.mcd' ;%30 conc
 file9='F:\EXTRAP\FP SPIKE\November 21, 2013\CAR\20600\112113 20600 CAR FPSPK.mcd' ;%30 conc
 
else 
%EPA computer
   % file1='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September 4, 2013 #20607 BIC\090413 20607 BIC FPSPK.mcd';
   % file2='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 8395 BIC\090513 8395 BIC FPSPK.mcd';
   % file3='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 15152 BIC\090513 15152 BIC FPSPK.mcd' ;
   % file4='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 18405 BIC\090513 18405 BIC FPSPK.mcd';
   % file5='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK.mcd';
   % file6='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK.mcd';          
   % file7=strcat( 'L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\',...
   %             'September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK.mcd')
   % file8 = 'L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK.mcd' ;
   % file9='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK.mcd';

end;

%new files
files.con={file1, file2, file3, file4, file5, file6, file7, file8, file9} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;


for curFile=1:length(files.con)
    
    file = files.con{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    con.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot

    spikes.n=0 ;
    for ent = 1:60

        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;

        con.file(curFile).channel(ent).nspikes = entity.ItemCount ;

        % indicator of ae or not
        con.file(curFile).channel(ent).ae =...
            (con.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount

            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            con.file(curFile).channel(ent).train(nspike) = timestamp ;
            
        end; %end of loop through timestamp in given entity

    end; %end of loop through all entities
    
    
    
    
    L = 25000; %length of 1 second signal
    %make L-point Kaiser window
    %beta=0.5, default value accoding to matlab documentation
    w_kaiser = kaiser(L, 0.5); 


    % loop through the ent's with sufficient high frequency activity
    for ent= 1:60
        %analog ent is 60 more than current
        analogEnt = ent +60 ;
        
        [nsresult,entity] = ns_GetEntityInfo(hfile,analogEnt)
        
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,analogEnt,1,entity.ItemCount);

        for sec=1:nseconds


           %take a chunk of length 1 second
            fs = 25000;
            curData_temp=data((sec-1)*fs+1:sec*fs);
    
            %apply kaiser window to attenuate signal
            curData = w_kaiser.*curData_temp ; %element-wise operations
            
            
            L = length(curData) ; 
            NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
        

            X = fft(curData, NFFT)/L;

            f = fs/2*linspace(0,1,NFFT/2+1); 

             X_mag = abs(X);
            % X_mag = real(X); % just get the real part to avoid averaging phase
            
            con.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);


        end; %end of loop through sec in an entity
        

        %for band=1:20
        for sec = 1:nseconds
        
            tempA = cell(1,nbands ) ; %used for bands

             for band=1:nbands
               %get indices of psd corresponding to current freq band
             index_t=( bands_needed(band) <= f(1:200) & f(1:200) <= bands_needed(band+1)) ;

                % get the particular band power at each band
                % this is a column of power for freq. from band 
                % defined by the column
                
                tempA{band} = con.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = bic(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                con.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( con.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;

            end; %end through bands
            
            con.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %bic(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     

    end;%for loop through entitites
    
end;  % for curFile=1:length(files.con)
    




for curFile=1:length(files.con)
    
    % start a matrix
    bandM = zeros(con.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:con.file(curFile).nseconds
                bandM(sec, ent, band) = con.file(curFile).channel(ent).sec(sec).avpw(band) ;
                %bandM(sec,entCol,band) = bic(sec).file(curFile).ent(ent).avpw(band);
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    con.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;

    
end; %loop through all files


% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
nbands = length( bands_needed) - 1 ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(con.file(curFile).nseconds/60)*60 ;

bands_needed=[1 4 8 14 30 50];

fs = 25000; %sampling rate        
L=fs ; %length of 1 sec worth of data
NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
f = fs/2*linspace(0,1,NFFT/2+1); %hertz
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 


 % get power for specific bands
f_temp = cell(1,5) ;
for band=1:length(bands_needed)-1
  % for each band get the corresponding frequencies
  index_freq = ( f >=bands_needed(band)& f <=bands_needed(band+1) ) ;

  f_temp{band} = f( index_freq ) ;
end;
        
con.fbands = f_temp ;


for curFile=1:length(files.con)
    for ent = 1:60
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(con.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd

           temp(curSec, : ) = transpose( con.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;

        end;
        con.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        con.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(con.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd

                temp(curSec, : ) = transpose( con.file(curFile).channel(ent).sec(curSec).pw{band} ) ;

            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            con.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files




% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_con.mat';
% save(filename, 'con', 'files' ) ;


% now combine all three files
% load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\11_19_2013_bic.mat')
% load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_car.mat')
% load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_con.mat')


% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_all.mat';
% save(filename, 'con', 'car', 'bic', 'files' ) ;




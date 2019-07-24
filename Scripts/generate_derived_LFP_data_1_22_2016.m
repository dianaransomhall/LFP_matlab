% generate_derived_LFP_data.m
% diana hall
% november 5, 2013
%%%Use NeuroShare
%you need the FIND toolbox and the NeuroShare
EPA=0 ;
mac=1;
% now all the find tools are avialable
if 0
    if EPA
        % EPA Computer
        addpath(genpath('\\AA.AD.EPA.GOV\ORD\RTP\USERS\A-D\dhall05\Net MyDocuments\MATLAB\NeuroShare')) ;
    else
        if mac
            addpath(genpath('/Users/dh2744/Dropbox/epa/LFP/find-1.1_64'));

        else    
            %personal computer
            addpath(genpath('C:\Users\diana_000\Documents\MATLAB\Cina_Herr_Diana'));
            endl
        end
    end
end
%set library
%personal computer
if EPA
    ns_SetLibrary('\\AA.AD.EPA.GOV\ORD\RTP\USERS\A-D\dhall05\Net MyDocuments\MATLAB\find-1.1\nsMCDLibrary.dll')
else
    if mac
        %ns_SetLibrary('/Users/dh2744/Dropbox/epa/LFP/NeuroExplorerNeuroShareLibrary.dll' );
        addpath(genpath('/Users/dh2744/Dropbox/epa/LFP/nsMCDLibrary_MacOSX_3.7b') );
        addpath('/Users/dh2744/Dropbox/epa/LFP/nsMCDLibrary_MacOSX_3.7b/Matlab/Matlab-Import-Filter/Matlab_Interface');
        [ns_RESULT] = ns_SetLibrary('/Users/dh2744/Dropbox/epa/LFP/nsMCDLibrary_MacOSX_3.7b/nsMCDLibrary/nsMCDLibrary.dylib' );
    else
        ns_SetLibrary('C:\Users\diana_000\\Documents\MATLAB\SpikeSorting\find-1.1\nsMCDLibrary64.dll')
    end;
end;

ns_GetLibraryInfo

% Ctrl+Z to undo
%file='F:\EXTRAP\FP SPIKE\DOM\Februaru 6, 2014 DA 15155\020614 15155 DA FPSPK0001.mcd' ;
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

file1='/Volumes/BACKUP/EXTRAP/FP SPIKE/BIC/090513 8395 BIC/090513 8395 BIC FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\BIC\090513 15152 BIC\090513 15152 BIC FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\BIC\090513 18405 BIC\090513 18405 BIC FPSPK0002.mcd';
file4='F:\EXTRAP\FP SPIKE\BIC\December 18, 2013 18395 BIC\121813 18395 BIC FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\BIC\November 19, 2013 BIC\111913 9858 BIC FPSPK0001.mcd';


bic.filepath = {file1, file2, file3, file4, file5 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(bic.filepath)
    
    file = bic.filepath{curFile} ;
    
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
            
            % curData2=data((fs*12):(fs*13));
            % L=length(curData2);NFFT=2^nextpow2(L);
            % X2=fft(curData2, NFFT)/L;
            L = length(curData) ; %25k
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
    



for curFile=1:length(bic.filepath)
    
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
for curFile=1:length(bic.filepath)
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




% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_bic.mat';
% save(filename, 'bic' ) ;




%%%%%%%%%%%%%%%%%%%   CONTROL   BIC         %%%%%%%%%%%%%%%%%%%%%%%%%%%
% bic files
file1='F:\EXTRAP\FP SPIKE\BIC\090513 8395 BIC\090513 8395 BIC FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\BIC\090513 15152 BIC\090513 15152 BIC FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\BIC\090513 18405 BIC\090513 18405 BIC FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\BIC\December 18, 2013 18395 BIC\121813 18395 BIC FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\BIC\November 19, 2013 BIC\111913 9858 BIC FPSPK.mcd';

% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4, file5 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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



% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_biccon.mat';
% save(filename, 'con' ) ;













%%%%%%%%%%%%%%%+++++  PERM25 0uM  ++++++++++++%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files
file1='F:\EXTRAP\FP SPIKE\PER\11-19-2013 PERM\112113 18413 PERM FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM  19854\121813 19854 PERM FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 19398\121813 19398 PERM FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 20608\121813 20608 PERM FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\PER\13079\111913 13079 PERM FPSPK0001.mcd';
per.filepath = {file1, file2, file3, file4, file5 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(per.filepath )
    
    file = per.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    per.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        per.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        per.file(curFile).channel(ent).ae =...
            (per.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            per.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            per.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %per(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

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
                
                tempA{band} = per.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = per(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                per.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( per.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %per(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( per(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            per.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %per(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.per)
    



for curFile=1:length(per.filepath)
    
    % start a matrix
    bandM = zeros(per.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:per.file(curFile).nseconds
                bandM(sec, ent, band) = per.file(curFile).channel(ent).sec(sec).avpw(band) ;
               
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    per.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(per.file(curFile).nseconds/60)*60 ;
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
        
per.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(per.filepath)
    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(per.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( per.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        per.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        per.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(per.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( per.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            per.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_per25.mat';
% save(filename, 'per' ) ;








%%%%%%%%%%%%%%%%%%%   CONTROL   PERM25        %%%%%%%%%%%%%%%%%%%%%%%%%%%
% per files
file1='F:\EXTRAP\FP SPIKE\PER\11-19-2013 PERM\112113 18413 PERM FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM  19854\121813 19854 PERM FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 19398\121813 19398 PERM FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\PER\12-18-13 PERM 20608\121813 20608 PERM FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\PER\13079\111913 13079 PERM FPSPK.mcd';

% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4, file5 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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



% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_per15con.mat';
% save(filename, 'con' ) ;















%%%%%%%%%%%%%%%+++++  PERM 50uM  ++++++++++++%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files
file1='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\February 11, 2014  PERM  19857\021114 19857 PERM FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\February 11, 2014 PERM 15158\021114 15158 PERM FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\February 11, 2014 PERM 20607\021114 20607 PERM FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 9858\041014 9858 PERM FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 18400\041014 18400 PERM FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 19398\041014 19398 PERM FPSPK0001.mcd';
file7='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 18397  10-30-14\103014 18397 PERM FPSPK0001.mcd';
file8='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\September 9, 2014 20607\090914 20607 PERM FPSPK0001.mcd';
per.filepath = {file1, file2, file3, file4, file5, file6,file7,file8 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(per.filepath )
 
    file = per.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    per.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        per.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        per.file(curFile).channel(ent).ae =...
            (per.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            per.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            per.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %per(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

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
                
                tempA{band} = per.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = per(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                per.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( per.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %per(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( per(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            per.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %per(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.per)
    



for curFile=1:length(per.filepath)
    
    % start a matrix
    bandM = zeros(per.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:per.file(curFile).nseconds
                bandM(sec, ent, band) = per.file(curFile).channel(ent).sec(sec).avpw(band) ;
               
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    per.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(per.file(curFile).nseconds/60)*60 ;
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
        
per.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(per.filepath)
    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(per.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( per.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        per.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        per.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(per.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( per.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            per.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_per50.mat';
% save(filename, 'per' ) ;






%%%%%%%%%%%%%%%%%%%   CONTROL  PERM50   %%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\February 11, 2014  PERM  19857\021114 19857 PERM FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\February 11, 2014 PERM 15158\021114 15158 PERM FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\February 11, 2014 PERM 20607\021114 20607 PERM FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 9858\041014 9858 PERM FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 18400\041014 18400 PERM FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 19398\041014 19398 PERM FPSPK.mcd';
file7='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\PERM 18397  10-30-14\103014 18397 PERM FPSPK.mcd';
file8='F:\EXTRAP\FP SPIKE\PER\50uM_Permethrin\September 9, 2014 20607\090914 20607 PERM FPSPK.mcd';
% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4,file5,file6,file7,file8 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=7:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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



% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conPER50.mat';
% save(filename, 'con' ) ;
























%%%%%%%++++   CAR5 uM   +++++++++++++
 
%%%%%%%%%%%%%%%%%%%CAR5 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%personal computer
 file1='F:\EXTRAP\FP SPIKE\CAR\5uM\09-19-13 CAR  20613\091913 20613 CAR FPSPK0001.mcd';
 file2='F:\EXTRAP\FP SPIKE\CAR\5uM\October 1 2013 CAR 15157\100113 15157 CAR FPSPK0001.mcd';
 file3='F:\EXTRAP\FP SPIKE\CAR\5uM\October 1 2013 CAR 20596\100113 20596 CAR FPSPK0002.mcd';
 file4='F:\EXTRAP\FP SPIKE\CAR\5uM\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK0001.mcd';
 file5='F:\EXTRAP\FP SPIKE\CAR\5uM\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK0001.mcd';

% 5uM CAR
car.filepath={file1, file2, file3, file4, file5 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(car.filepath)
    
    file = car.filepath{curFile} ;
    
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


for curFile=1:length(car.filepath)
    
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
for curFile=1:length(car.filepath)
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
% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car5.mat';
% save(filename, 'car' ) ;





%%%%%%%%%%%%%%%%%%%   CONTROL   CAR5         %%%%%%%%%%%%%%%%%%%%%%%%%%%
% car files
 file1='F:\EXTRAP\FP SPIKE\CAR\5uM\09-19-13 CAR  20613\091913 20613 CAR FPSPK.mcd';
 file2='F:\EXTRAP\FP SPIKE\CAR\5uM\October 1 2013 CAR 15157\100113 15157 CAR FPSPK.mcd';
 file3='F:\EXTRAP\FP SPIKE\CAR\5uM\October 1 2013 CAR 20596\100113 20596 CAR FPSPK.mcd';
 file4='F:\EXTRAP\FP SPIKE\CAR\5uM\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK.mcd';
 file5='F:\EXTRAP\FP SPIKE\CAR\5uM\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK.mcd';
% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4, file5 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file);
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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



% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2014_car5con.mat';
% save(filename, 'con' ) ;



























%%%%%%%++++   CAR  15uM   +++++++++++++
 
%%%%%%%%%%%%%%%%%%%CAR%%%%%%%%%%%%%%%%%%%%%%%%%%%
%personal computer
file1='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19855\112113 19855 CAR FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\20600\112113 20600 CAR FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19857_15uMOnly\112113 19857 CAR FPSPK0001.mcd';
% 15uM CAR
car.filepath={file1, file2, file3 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(car.filepath)
    
    file = car.filepath{curFile} ;
    
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


for curFile=1:length(car.filepath)
    
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
for curFile=1:length(car.filepath)
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
% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car15.mat';
% save(filename, 'car' ) ;



%%%%%%%%%%%%%%%%%%%   CONTROL   CAR15         %%%%%%%%%%%%%%%%%%%%%%%%%%%
% car files
 file1='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19855\112113 19855 CAR FPSPK.mcd';
 file2='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\20600\112113 20600 CAR FPSPK.mcd';
 file3='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19857_15uMOnly\112113 19857 CAR FPSPK.mcd';
% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file);
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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



% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2014_car15con.mat';
% save(filename, 'con' ) ;

















%%%%%%%++++   CAR  30uM   +++++++++++++
 
%%%%%%%%%%%%%%%%%%%CAR%%%%%%%%%%%%%%%%%%%%%%%%%%%
%personal computer
file1='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19855\112113 19855 CAR FPSPK0002.mcd';
file2='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\20600\112113 20600 CAR FPSPK0002.mcd';
file3='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 18412\121713 18412 CAR FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 19858\121713 19858 CAR FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 20599\121713 20599 CAR FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\October 15, 2014 18416\101514 18416 CAR FPSPK0001.mcd';
file7='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\October 16, 2014 CAR 26018\101614 26018 CAR FPSPK0001.mcd';
file8='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\September 11, 2014 20613\091114 20613 CAR FPSPK0001.mcd';
file9='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\18408 Sept 9 2014\090914 18408 CAR FPSPK0001.mcd';
% 30uM CAR
car.filepath={file1, file2, file3, file4, file5, file6,file7,file8,file9} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=9:length(car.filepath)
    
    file = car.filepath{curFile} ;
    
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


for curFile=1:length(car.filepath)
    
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
for curFile=1:length(car.filepath)
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
% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30.mat';
% save(filename, 'car' ) ;



%%%%%%%%%%%%%%%%%%%   CONTROL   CAR30         %%%%%%%%%%%%%%%%%%%%%%%%%%%
% car files
file1='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\19855\112113 19855 CAR FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\20600\112113 20600 CAR FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 18412\121713 18412 CAR FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 19858\121713 19858 CAR FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\December 17, 2013 CAR 20599\121713 20599 CAR FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\October 15, 2014 18416\101514 18416 CAR FPSPK.mcd';
file7='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\October 16, 2014 CAR 26018\101614 26018 CAR FPSPK.mcd';
file8='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\September 11, 2014 20613\091114 20613 CAR FPSPK.mcd';
file9='F:\EXTRAP\FP SPIKE\CAR\15uM_30uM\18408 Sept 9 2014\090914 18408 CAR FPSPK.mcd';

% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4, file5, file6, file7, file8,file9 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file);
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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



% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30con.mat';
% save(filename, 'con' ) ;







%%%%%%%++++   DA:=demoic acid 30uM   +++++++++++++
 
%%%%%%%%%%%%%%%%%%%CAR%%%%%%%%%%%%%%%%%%%%%%%%%%%
%personal computer
file1='F:\EXTRAP\FP SPIKE\DOM\Februaru 6, 2014 DA 15155\020614 15155 DA FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\DOM\DA 18412\040914 18412 DA FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\DOM\February 11, 2014 DA 18413\021114 18413 DA FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\DOM\February 11, 2014 DA 18416\021114 18416 DA FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\DOM\February 11, 2014 DA 19860\021114 19860 DA FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\DOM\DA 20608\040914 20608 DA FPSPK0001.mcd';

% 300uM DA
da.filepath={file1, file2, file3, file4, file5, file6} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(da.filepath)
    
    file = da.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    da.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        da.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        da.file(curFile).channel(ent).ae =...
            (da.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            da.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            da.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
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
                
                tempA{band} = da.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = bic(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                da.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( da.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %bic(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( bic(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            da.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %bic(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.bic)


for curFile=1:length(da.filepath)
    
    % start a matrix
    bandM = zeros( da.file(curFile).nseconds , 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:da.file(curFile).nseconds
                bandM(sec, ent, band) = da.file(curFile).channel(ent).sec(sec).avpw(band) ;
                %bandM(sec,entCol,band) = bic(sec).file(curFile).ent(ent).avpw(band);
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    da.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files


% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
nbands = length( bands_needed) - 1 ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(da.file(curFile).nseconds/60)*60 ;
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
        
da.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(da.filepath)
    for ent = 1:60
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(da.file(curFile).channel(ent).sec(1).psd(index_all_freq ) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( da.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        da.file(curFile).channel(ent).f = f(index_all_freq ) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        da.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(da.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( da.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            da.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files
% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_da.mat';
% save(filename, 'da' ) ;








%%%%%%%%%%%%%%%%%%%   CONTROL   DA 300nM   %%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\DOM\Februaru 6, 2014 DA 15155\020614 15155 DA FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\DOM\DA 18412\040914 18412 DA FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\DOM\February 11, 2014 DA 18413\021114 18413 DA FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\DOM\February 11, 2014 DA 18416\021114 18416 DA FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\DOM\February 11, 2014 DA 19860\021114 19860 DA FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\DOM\DA 20608\040914 20608 DA FPSPK.mcd';

% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4, file5, file6 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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
% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_conDA.mat';
% save(filename, 'con' ) ;



















%%%%%%%%%%   H20
%%%%%%%%%%%%%%%+++++  H20 (not control) ++++++++++++%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files
file1='F:\EXTRAP\FP SPIKE\H2O\15158 H2O\052914 15158 CON FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\H2O\18401 H2O\050114 18401 CON FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\H2O\19855 H2O\052214 19855 CON FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\H2O\19858 H2O 05-21-14\052114 19858 CON FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\H2O\April 16, 2014 H2O 18408\041614 18408 CON FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\H2O\April 16, 2014 H2O 19857\041614 19857 CON FPSPK0001.mcd';
file7='F:\EXTRAP\FP SPIKE\H2O\May 20, 2014 151223\052014 151223 CON FPSPK0001.mcd';
file8='F:\EXTRAP\FP SPIKE\H2O\May 22, 2014 18395\052214 18395 CON FPSPK0001.mcd';
file9='F:\EXTRAP\FP SPIKE\H2O\May 22,2014 20613\052214 20613 CON FPSPK0001.mcd';
file10='F:\EXTRAP\FP SPIKE\H2O\May 29, 2014 19861\052914 19861 CON FPSPK0001.mcd';
file11='F:\EXTRAP\FP SPIKE\H2O\October 15, 2014 19860\101514 19860 CON FPSPK0001.mcd';
file12='F:\EXTRAP\FP SPIKE\H2O\October 16, 2014 H2O 26014\101614 26014 CON FPSPK0001.mcd';
file13='F:\EXTRAP\FP SPIKE\H2O\October 23, 2014 H2O\102314 19861 CON FPSPK0001.mcd';
file14='F:\EXTRAP\FP SPIKE\H2O\July 1, 2014 20613\070114 20613 CON FPSPK0001.mcd';
file15='F:\EXTRAP\FP SPIKE\H2O\20608 H2O 121014\121014 20608 CON FPSPK0001.mcd';
file16='F:\EXTRAP\FP SPIKE\H2O\19861 H2O 19861\121014 19861 CON FPSPK0001.mcd';
h20.filepath = {file1, file2, file3, file4,file5,file6,file7,file8,file9,...
    file10, file11, file12, file13, file14, file15, file16 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=14:length(h20.filepath )
    
    file = h20.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    h20.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        h20.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        h20.file(curFile).channel(ent).ae =...
            (h20.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            h20.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            h20.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %per(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

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
                
                tempA{band} = h20.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = per(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                h20.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( h20.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %per(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( per(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            h20.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %per(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.per)
    



for curFile=14:length(h20.filepath)
    
    % start a matrix
    bandM = zeros(h20.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:h20.file(curFile).nseconds
                bandM(sec, ent, band) =h20.file(curFile).channel(ent).sec(sec).avpw(band) ;
               
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    h20.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(h20.file(curFile).nseconds/60)*60 ;
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
        
h20.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=11:length(h20.filepath)
    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(h20.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( h20.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        h20.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        h20.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(h20.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( h20.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            h20.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_h20.mat';
% save(filename, 'h20' ) ;


%%%%%%%%%%%%%%%%%%%   CONTROL  H20   %%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\H2O\15158 H2O\052914 15158 CON FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\H2O\18401 H2O\050114 18401 CON FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\H2O\19855 H2O\052214 19855 CON FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\H2O\19858 H2O 05-21-14\052114 19858 CON FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\H2O\April 16, 2014 H2O 18408\041614 18408 CON FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\H2O\April 16, 2014 H2O 19857\041614 19857 CON FPSPK.mcd';
file7='F:\EXTRAP\FP SPIKE\H2O\May 20, 2014 151223\052014 151223 CON FPSPK.mcd';
file8='F:\EXTRAP\FP SPIKE\H2O\May 22, 2014 18395\052214 18395 CON FPSPK.mcd';
file9='F:\EXTRAP\FP SPIKE\H2O\May 22,2014 20613\052214 20613 CON FPSPK.mcd';
file10='F:\EXTRAP\FP SPIKE\H2O\May 29, 2014 19861\052914 19861 CON FPSPK.mcd';
file11='F:\EXTRAP\FP SPIKE\H2O\October 15, 2014 19860\101514 19860 CON FPSPK.mcd';
file12='F:\EXTRAP\FP SPIKE\H2O\October 16, 2014 H2O 26014\101614 26014 CON FPSPK.mcd';
file13='F:\EXTRAP\FP SPIKE\H2O\October 23, 2014 H2O\102314 19861 CON FPSPK.mcd';
file14='F:\EXTRAP\FP SPIKE\H2O\July 1, 2014 20613\070114 20613 CON FPSPK.mcd';
file15='F:\EXTRAP\FP SPIKE\H2O\20608 H2O 121014\121014 20608 CON FPSPK.mcd';
file16='F:\EXTRAP\FP SPIKE\H2O\19861 H2O 19861\121014 19861 CON FPSPK.mcd';
% last control file doesn't work well
% we remove file9 because it's corrupted as of now
con.filepath={file1, file2, file3, file4, file5,file6,file7,file8,file9,...
    file10, file11, file12, file13, file14, file15, file16 } ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=13:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
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
    



for curFile=11:length(con.filepath)
    
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

for curFile=11:length(con.filepath)
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

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conh20.mat';
% save(filename, 'con' ) ;




















%%%%%%%%%%  ACE
%%%%%%%%%%%%%%%+++++  ACE (not control) ++++++++++++%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files
file1='F:\EXTRAP\FP SPIKE\ACE\#20613 June 4, 2014\060414 20613 ACE FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\ACE\June 3, 2014 13102\060314 13102 ACE FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\ACE\June 3, 2014 19861\060314 19861 ACE FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\ACE\May 20, 2014 20596\052014 20596 ACE FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\ACE\May 22, 2014 151223\052214 151223 CON FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\ACE\October 23, 2014 ACE\102314 26014 ACE FPSPK0001.mcd';
file7='F:\EXTRAP\FP SPIKE\ACE\November 12, 2014 15158 ACE\111214 15158 ACE FPSPK0001.mcd';
file8='F:\EXTRAP\FP SPIKE\ACE\ACE 19861 10-30-14\103014 19861 ACE FPSPK0001.mcd';
file9='F:\EXTRAP\FP SPIKE\ACE\October 23, 2014 ACE\102314 26014 ACE FPSPK0001.mcd';
ace.filepath = {file1, file2, file3, file4,file5,file6,file7,file8,file9 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=7:length(ace.filepath )
    
    file = ace.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    ace.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        ace.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        ace.file(curFile).channel(ent).ae =...
            (ace.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            ace.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            ace.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %per(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

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
                
                tempA{band} = ace.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = per(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                ace.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( ace.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %per(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( per(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            ace.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %per(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.per)
    



for curFile=1:length(ace.filepath)
    
    % start a matrix
    bandM = zeros(ace.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:ace.file(curFile).nseconds
                bandM(sec, ent, band) =ace.file(curFile).channel(ent).sec(sec).avpw(band) ;
               
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    ace.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(ace.file(curFile).nseconds/60)*60 ;
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
        
h20.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(ace.filepath)
    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(ace.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( ace.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        ace.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        ace.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(ace.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( ace.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            ace.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_ace.mat';
% save(filename, 'ace' ) ;


%%%%%%%%%%%%%%%%%%%   CONTROL  ACE   %%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\ACE\#20613 June 4, 2014\060414 20613 ACE FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\ACE\June 3, 2014 13102\060314 13102 ACE FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\ACE\June 3, 2014 19861\060314 19861 ACE FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\ACE\May 20, 2014 20596\052014 20596 ACE FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\ACE\May 22, 2014 151223\052214 151223 CON FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\ACE\October 23, 2014 ACE\102314 26014 ACE FPSPK.mcd';
file7='F:\EXTRAP\FP SPIKE\ACE\November 12, 2014 15158 ACE\111214 15158 ACE FPSPK.mcd';
file8='F:\EXTRAP\FP SPIKE\ACE\ACE 19861 10-30-14\103014 19861 ACE FPSPK.mcd';
file9='F:\EXTRAP\FP SPIKE\ACE\October 23, 2014 ACE\102314 26014 ACE FPSPK.mcd';

con.filepath={file1, file2, file3 , file4, file5, file6,file7,file8, file9} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=7:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conace.mat';
% save(filename, 'con' ) ;


















%%%%%%%%%%  TRI
%%%%%%%%%%%%%%%+++++  TRI (not control) ++++++++++++%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files
file1='F:\EXTRAP\FP SPIKE\TRI\December 18, 2013 TRI\121813 18405 TRI FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\TRI\Nov 21, 2013 18416\112113 18416 TRI FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\TRI\November 13, 2013 19861 TRI\111313 19861 TRI FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\TRI\November 14, 2013 TRI\15158\111413 15158 TRI FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\TRI\October 1 2013 TRI 18406\100113 18406 TRI FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\TRI\October 1 2013 TRI 20608\100113 20608 TRI FPSPK0001.mcd';
file7='F:\EXTRAP\FP SPIKE\TRI\Sept 26, 2013 TRI 18412\092613 18412 TRI FPSPK0001.mcd';
%file8='F:\EXTRAP\FP SPIKE\TRI\November 14, 2013 TRI\111413 20598 TRI FPSPK0001.mcd';
tri.filepath = {file1, file2, file3, file4,file5,file6,file7 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(tri.filepath )
    
    file = tri.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    tri.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        tri.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        tri.file(curFile).channel(ent).ae =...
            (tri.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            tri.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            tri.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %per(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

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
                
                tempA{band} = tri.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = per(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                tri.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( tri.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %per(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( per(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            tri.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %per(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.per)
    



for curFile=1:length(tri.filepath)
    
    % start a matrix
    bandM = zeros(tri.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:tri.file(curFile).nseconds
                bandM(sec, ent, band) =tri.file(curFile).channel(ent).sec(sec).avpw(band) ;
               
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    tri.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(tri.file(curFile).nseconds/60)*60 ;
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
        
h20.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(tri.filepath)
    for ent = 1:60
        
        % nsecondsRd is nseconds rounded down to nearest 60 second increment
          nsecondsRd = floor(tri.file(curFile).nseconds/60)*60 ;
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(tri.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( tri.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        tri.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        tri.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(tri.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( tri.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            tri.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_tri.mat';
% save(filename, 'tri' ) ;


%%%%%%%%%%%%%%%%%%%   CONTROL  TRI   %%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\TRI\December 18, 2013 TRI\121813 18405 TRI FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\TRI\Nov 21, 2013 18416\112113 18416 TRI FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\TRI\November 13, 2013 19861 TRI\111313 19861 TRI FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\TRI\November 14, 2013 TRI\15158\111413 15158 TRI FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\TRI\October 1 2013 TRI 18406\100113 18406 TRI FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\TRI\October 1 2013 TRI 20608\100113 20608 TRI FPSPK.mcd';
file7='F:\EXTRAP\FP SPIKE\TRI\Sept 26, 2013 TRI 18412\092613 18412 TRI FPSPK.mcd';
con.filepath = {file1, file2, file3, file4,file5,file6,file7 }  ;


for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile);
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
    for ent = 1:60
        % nsecondsRd is nseconds rounded down to nearest 60 second increment
        nsecondsRd = min(300, floor(tri.file(curFile).nseconds/60)*60) ;
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

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_contri.mat';
% save(filename, 'con' ) ;








%%%%%%%%%%%%%%%+++++  GLY   ++++++++++++%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%GLY%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\GLY\August 6, 2014 GLY\080614 26016 GLY FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\GLY\August 21, 2014 GLY 18416\082114 18416 GLY FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\GLY\August 21, GLY 151223\082114 151223 GLY FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\GLY\June 18, 2014 19860\061814 19860 GLY FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\GLY\June 18, 2014 20608\061814 20608 GLY FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\GLY\October 16, 2014 GLY 26022\101614 26022 GLY FPSPK0001.mcd';
file7='F:\EXTRAP\FP SPIKE\GLY\September 11, 2014 18412\091114 18412 GLY FPSPK0001.mcd';
file8='F:\EXTRAP\FP SPIKE\GLY\September 25, 2014\092514 20608 GLY FPSPK0001.mcd';
file9='F:\EXTRAP\FP SPIKE\GLY\October 23, 2014 19860 GLY\102314 19860 GLY FPSPK0001.mcd';
gly.filepath = {file1, file2, file3, file4, file5,file6,file7,file8, file9 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=9:length(gly.filepath)
    
    file = gly.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    gly.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        gly.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        gly.file(curFile).channel(ent).ae =...
            (gly.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            gly.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            % curData2=data((fs*12):(fs*13));
            % L=length(curData2);NFFT=2^nextpow2(L);
            % X2=fft(curData2, NFFT)/L;
            L = length(curData) ; %25k
            NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L
            
            X = fft(curData, NFFT)/L;
            f = fs/2*linspace(0,1,NFFT/2+1); 
            X_mag = abs(X);
            % X_mag = real(X);
            
            gly.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
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
                
                tempA{band} = gly.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = bic(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                gly.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( gly.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %bic(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( bic(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            gly.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %bic(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.gly)
    

for curFile=1:length(gly.filepath)
    
    % start a matrix
    bandM = zeros(gly.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:gly.file(curFile).nseconds
                bandM(sec, ent, band) = gly.file(curFile).channel(ent).sec(sec).avpw(band) ;
                %bandM(sec,entCol,band) = bic(sec).file(curFile).ent(ent).avpw(band);
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    gly.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(gly.file(curFile).nseconds/60)*60 ;
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
        
gly.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(gly.filepath)
    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(gly.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( gly.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        gly.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        gly.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(gly.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( gly.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            gly.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files




% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_gly.mat';
% save(filename, 'gly' ) ;




%%%%%%%%%%%%%%%+++++  GLY   ++++++++++++%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%GLY%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\GLY\August 6, 2014 GLY\080614 26016 GLY FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\GLY\August 21, 2014 GLY 18416\082114 18416 GLY FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\GLY\August 21, GLY 151223\082114 151223 GLY FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\GLY\June 18, 2014 19860\061814 19860 GLY FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\GLY\June 18, 2014 20608\061814 20608 GLY FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\GLY\October 16, 2014 GLY 26022\101614 26022 GLY FPSPK.mcd';
file7='F:\EXTRAP\FP SPIKE\GLY\September 11, 2014 18412\091114 18412 GLY FPSPK.mcd';
file8='F:\EXTRAP\FP SPIKE\GLY\September 25, 2014\092514 20608 GLY FPSPK.mcd';
file9='F:\EXTRAP\FP SPIKE\GLY\October 23, 2014 19860 GLY\102314 19860 GLY FPSPK.mcd';
con.filepath = {file1, file2, file3, file4,file5,file6,file7, file8, file9 }  ;


for curFile=9:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile);
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
    for ent = 1:60
        % nsecondsRd is nseconds rounded down to nearest 60 second increment
        nsecondsRd = min(300, floor(gly.file(curFile).nseconds/60)*60) ;
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

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_congly.mat';
% save(filename, 'con' ) ;













%%%%%%%%%%  DMSO
%%%%%%%%%%%%%%%+++++  DMSO (not control) ++++++++++++%%%%%%%%%%%%%%%%
%the FPSSPK0001 are the treatment files
file1='F:\EXTRAP\FP SPIKE\DMSO\December 3, 2014\120314 20600 DMET FPSPK0001.mcd';
file2='F:\EXTRAP\FP SPIKE\DMSO\December 11, 2014 18406\121114 18406 DMET FPSPK0001.mcd';
file3='F:\EXTRAP\FP SPIKE\DMSO\December 11, 2014 18416\121114 18416 DMET FPSPK0001.mcd';
file4='F:\EXTRAP\FP SPIKE\DMSO\December 11, 2014 19861\121114 19861 DMET FPSPK0001.mcd';
file5='F:\EXTRAP\FP SPIKE\DMSO\December 17, 2014 26014\121714 26014 DMET FPSPK0001.mcd';
file6='F:\EXTRAP\FP SPIKE\DMSO\December 18, 2014  26016\121814 26016 DMET FPSPK0001.mcd';
file7='F:\EXTRAP\FP SPIKE\DMSO\DMET 20608 Dec 18-14\121814 20608 DMET FPSPK0001.mcd';
file8='F:\EXTRAP\FP SPIKE\DMSO\November 6, 2014 18415\110614 18415 DMET FPSPK0001.mcd';
file9='F:\EXTRAP\FP SPIKE\DMSO\November 12, 2014 26024 DMET\111214 26024 DMET FPSPK0001.mcd';
file10='F:\EXTRAP\FP SPIKE\DMSO\November 13, 2014 13102 DMET\111314 13102 DMET FPSPK0002.mcd';
dmso.filepath = {file1,file2,file3,file4,file5,file6,file7,file8,file9,file10 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=7:length(dmso.filepath )
    
    file = dmso.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    dmso.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        dmso.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        dmso.file(curFile).channel(ent).ae =...
            (dmso.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            dmso.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            dmso.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %per(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

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
                
                tempA{band} = dmso.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = per(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                dmso.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( dmso.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %per(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( per(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            dmso.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %per(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.per)
    



for curFile=1:length(dmso.filepath)
    
    % start a matrix
    bandM = zeros(dmso.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:dmso.file(curFile).nseconds
                bandM(sec, ent, band) =dmso.file(curFile).channel(ent).sec(sec).avpw(band) ;
               
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    dmso.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(dmso.file(curFile).nseconds/60)*60 ;
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
        
dmso.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(dmso.filepath)
    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(dmso.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( dmso.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        dmso.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        dmso.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(dmso.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( dmso.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            dmso.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_dmso.mat';
% save(filename, 'dmso' ) ;







%%%%%%%%%%%%%%%%%%%   CONTROL  DMSO   %%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='F:\EXTRAP\FP SPIKE\DMSO\December 3, 2014\120314 20600 DMET FPSPK.mcd';
file2='F:\EXTRAP\FP SPIKE\DMSO\December 11, 2014 18406\121114 18406 DMET FPSPK.mcd';
file3='F:\EXTRAP\FP SPIKE\DMSO\December 11, 2014 18416\121114 18416 DMET FPSPK.mcd';
file4='F:\EXTRAP\FP SPIKE\DMSO\December 11, 2014 19861\121114 19861 DMET FPSPK.mcd';
file5='F:\EXTRAP\FP SPIKE\DMSO\December 17, 2014 26014\121714 26014 DMET FPSPK.mcd';
file6='F:\EXTRAP\FP SPIKE\DMSO\December 18, 2014  26016\121814 26016 DMET FPSPK.mcd';
file7='F:\EXTRAP\FP SPIKE\DMSO\DMET 20608 Dec 18-14\121814 20608 DMET FPSPK.mcd';
file8='F:\EXTRAP\FP SPIKE\DMSO\November 6, 2014 18415\110614 18415 DMET FPSPK.mcd';
file9='F:\EXTRAP\FP SPIKE\DMSO\November 12, 2014 26024 DMET\111214 26024 DMET FPSPK.mcd';
file10='F:\EXTRAP\FP SPIKE\DMSO\November 13, 2014 13102 DMET\111314 13102 DMET FPSPK.mcd';


con.filepath={file1, file2, file3 ,file4, file5, file6,file7,file8,file9,file10} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile);
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_condmso.mat';
% save(filename, 'con' ) ;











%+++++++++++++++++++++++++++++LINDANE+++++++++++++++++++++++++++++++++++++++++
%%%%%%%%%%  LINDANE
%%%%%%%%%%%%%%%+++++  LINDANE ++++++++++++%%%%%%%%%%%%%%%%
%
file1='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/011315 26011 LIN FPSPK0002.mcd';
file2='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/011515 19857 LIN FPSPK0002.mcd';
file3='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/012115 20599 LIN FPSPK0002.mcd';
file4='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/012215 26022 LIN FPSPK0002.mcd';
file5='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/020515 13079 LIN FPSPK0002.mcd';
file6='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/021215 26016 LIN FPSPK0002.mcd';
file7='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031215 19858 LIN FPSPK0002.mcd';
file8='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031815 20594 LIN FPSPK0002.mcd';
file9='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031915 26024 LIN FPSPK0002.mcd';
file10='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/043015 20594 LIN FPSPK0002.mcd';
file11='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/050615 19858 LIN FPSPK0002.mcd';
file12='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/062515 18419 LIN FPSPK0002.mcd';
file13='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031915 26024 LIN FPSPK0002.mcd';
file14='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031915 26024 LIN FPSPK0002.mcd';
lin.filepath = {file1,file2,file3,file4,file5,file6,file7,...
    file8,file9,file10 }  ;
        
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;
for curFile=1:length(lin.filepath )
    
    file = lin.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)
    % Find out about the Entity types
    % Then read specific entity info and data
    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    
    lin.file(curFile).nseconds = nseconds ;
    
    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot
    spikes.n=0 ;
    for ent = 1:60
        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;
        lin.file(curFile).channel(ent).nspikes = entity.ItemCount ;
        % indicator of ae or not
        lin.file(curFile).channel(ent).ae =...
            (lin.file(curFile).channel(ent).nspikes >floor( (nseconds*(5/60)) ) ) ;   
        
        for nspike = 1:entity.ItemCount
            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            lin.file(curFile).channel(ent).train(nspike) = timestamp ;
            
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
            
            lin.file(curFile).channel(ent).sec(sec).psd = X_mag(1:200);
            %per(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

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
                
                tempA{band} = lin.file(curFile).channel(ent).sec(sec).psd(index_t) ;
                %tempA{band} = per(sec).file(curFile).ent(ent).psd(index_t) ;
                
                %now get mean power
                lin.file(curFile).channel(ent).sec(sec).avpw(band) = ...
                    mean( lin.file(curFile).channel(ent).sec(sec).psd(index_t) ) ;
                %per(sec).file(curFile).ent(ent).avpw(band) = ...
                %     mean( per(sec).file(curFile).ent(ent).psd(index_t) ) ;
        
            end; %end through bands
            
            lin.file(curFile).channel(ent).sec(sec).pw = tempA ;
           %per(sec).file(curFile).ent(ent).pw = tempA ;
        
        end;  %end of loop through secs     
    end;%for loop through entitites
    
end;  % for curFile=1:length(files.per)
    



for curFile=1:length(lin.filepath)
    
    % start a matrix
    bandM = zeros(lin.file(curFile).nseconds, 60 ,nbands) ;
    
    % make a multi-array, one matrix for each band
    for band=1:nbands
        for ent = 1:60
            for sec=1:lin.file(curFile).nseconds
                bandM(sec, ent, band) =lin.file(curFile).channel(ent).sec(sec).avpw(band) ;
               
            end; %sec loop
        end; %ent loop
    end; %band loop
    
    lin.file(curFile).bandM = bandM ;
    %apwBIC(curFile).bandM = bandM ;
    
end; %loop through all files

% matrix at all the psd for each band in each channel
allPw = cell(1,5) ;
% nsecondsRd is nseconds rounded down to nearest 60 second increment
nsecondsRd = floor(lin.file(curFile).nseconds/60)*60 ;
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
        
lin.fbands = f_temp ;
%index corresponding to freq between 1 and 50 Hz 
index_all_freq = ( 1 <= f(1:200) & f(1:200) <= 50 ) ; 
for curFile=1:length(lin.filepath)
    for ent = 1:60
        
        
         %get matrix of power from all bands 1-50Hz
        clear temp ;
        temp = -ones(nsecondsRd,...
              length(lin.file(curFile).channel(ent).sec(1).psd(index_all_freq) ) ) ;
        
        for curSec=1:nsecondsRd
           temp(curSec, : ) = transpose( lin.file(curFile).channel(ent).sec(curSec).psd(index_all_freq) ) ;
        end;
        lin.file(curFile).channel(ent).f = f(index_all_freq) ;
        % nsecondsRd by nobs in all bands (length( f ) )
        lin.file(curFile).channel(ent).allPwallBands = temp ;
       
        
        for band = 1:nbands
            clear temp ;
            temp = -ones(nsecondsRd,...
              length(lin.file(curFile).channel(ent).sec(1).pw{band}) ) ;
            for curSec=1:nsecondsRd
                temp(curSec, : ) = transpose( lin.file(curFile).channel(ent).sec(curSec).pw{band} ) ;
            end;
            tempAllBands = temp ;
            % allPw cell of matrices, 1 matrix for each band
            %  matrix is nSecondsRd by number of obs. in each freq band
            lin.file(curFile).channel(ent).allPw{band} = temp ;
            
        end;
        
       
    end; % for ent
    
end; % loop through files

 filename='/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/02_15_2016_lindane.mat';
 save(filename, 'lin' ) ;







%%%%%%%%%%%%%%%%%%%   CONTROL  LINDANE   %%%%%%%%%%%%%%%%%%%%%%%%%%%
file1='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/011315 26011 LIN FPSPK.mcd';
file2='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/011515 19857 LIN FPSPK.mcd';
file3='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/012115 20599 LIN FPSPK.mcd';
file4='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/012215 26022 LIN FPSPK.mcd';
file5='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/020515 13079 LIN FPSPK.mcd';
file6='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/021215 26016 LIN FPSPK.mcd';
file7='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031215 19858 LIN FPSPK.mcd';
file8='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031815 20594 LIN FPSPK.mcd';
file9='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031915 26024 LIN FPSPK.mcd';
file10='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/043015 20594 LIN FPSPK.mcd';
file11='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/050615 19858 LIN FPSPK.mcd';
file12='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/062515 18419 LIN FPSPK.mcd';
file13='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031915 26024 LIN FPSPK.mcd';
file14='/Volumes/BACKUP/EXTRAP/FP SPIKE/LIN/031915 26024 LIN FPSPK.mcd';
con.filepath = {file1,file2,file3,file4,file5,file6,file7,...
    file8,file9,file10 }  ;


con.filepath={file1, file2, file3 ,file4, file5, file6,file7,file8,file9,file10} ;
bands_needed=[1 4 8 14 30 50];
nbands = length( bands_needed) - 1 ;

for curFile=1:length(con.filepath)
    
    file = con.filepath{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile);
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
    



for curFile=1:length(con.filepath)
    
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

for curFile=1:length(con.filepath)
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

 filename='/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/02_15_2016_conlin.mat';
 save(filename, 'con' ) ;
%LFP_Data_analysis.m
%diana hall
%Oct 23, 2013
% purpose: to read in data that was collected on MCD by cina,
%  that contain low frequency and high frequency signals
%  we do a PSD and then apply an ANOVA to average power
%  over a given frequency (control-baseline



%%%Use NeuroShare
%you need the FIND toolbox and the NeuroShare

EPA=0 ;

if EPA
    % EPA Computer
    addpath(genpath('\\AA.AD.EPA.GOV\ORD\RTP\USERS\A-D\dhall05\Net MyDocuments\MATLAB\NeuroShare')) ;
else
    %personal computer
%
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

for curFile=1:length(files.bic)
    
    file = files.bic{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;
    

    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot

    spikes.n=0 ;
    for ent = 1:60

        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;

        spikes.bic(ent).file(curFile).n = entity.ItemCount;
        index_ent(ent) =(spikes.bic(ent).file(curFile).n >25 ) ;   
        
        for nspike = 1:entity.ItemCount

            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            spikes.bic(ent).file(curFile).train(nspike) = timestamp ;

        end; %end of loop through timestamp in given entity

    end; %end of loop through all entities


    %get wells that are at least 5spikes/minute = 25spikes in 5mins

    
    %print
    graphTitle='Graph_BIC.ps' ;
    
    printFig=false;
    a = 61:120 ;
    for ent=a(index_ent)
        %ent=64;
        %ent=64 is best, since spike train 4 is best
        [nsresult,entity] = ns_GetEntityInfo(hfile,ent)
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,ent,1,entity.ItemCount);

        for sec=1:nseconds


            %take a chunk of length 1 second
            fs = 25000;
            curData=data((sec-1)*fs+1:sec*fs);

            L=length(curData) ; %7,540,000
            NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L

            X = fft(curData, NFFT)/L;

            f = fs/2*linspace(0,1,NFFT/2+1); 

            X_mag = abs(X);


            bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

            
            
           if printFig; 
                %open new figure
                figure; %this opens a blank figure
                hold off;

                subplot(3,1,1)

                index=(spikes.bic(ent-60).train<=sec)&(spikes.bic(ent-60).train>sec-1);
                in=zeros(1,1)
                in(1)=sec-1;
                in(2+ sum(index) )=sec;

                if sum(index)~=0
                    in(2:sum(index)+1 )=spikes.bic(ent-60).train(index);
                end;


                axis([sec-1 sec -1 2])
                plot([in;in],[ones(size(in));zeros(size(in))],'k-')



                subplot(3,1,2)
                plot(25000/2*curData);
                ylim([-0.25 0.3])
                hold on;

                %we're plotting the low frequency 1-40Hz over top
                %
                t = 0 : 1/fs : 1 - 1/fs ;
                co=8; %make the cut off 8Hz
                sig=2*co*sinc(2*co*t);
                %taking first 40 frequency is equivalent to
                %convolving original signal with sig above
                v=conv(sig, curData);

                plot(v(1:25000), 'r')
                ylim([-0.25 0.3])
                %V = fft(v, NFFT)/L;  V_mag = abs(V);

                %plot(f(2:120), 2*V_mag(2:120 ))
                %ylim([0 0.1])

                % plot(f, 2*X_mag(1:NFFT/2 + 1)) %frequency (Hz)
                subplot(3,1,3)



                plot(f(2:120), 2*X_mag(2:120 )) %frequency (Hz)
                ylim([0 0.000002])
                hold on;
                xlim([1,100])
                title(strcat('PSD of 090413 20607 BIC Second ',...
                    num2str(sec),'Entity',num2str(ent) ) )
                xlabel('Frequency (Hz)')
                ylabel('Magnitude = amplitude ^2')

                if sec==1
                    %print('-dpdf',  graphTitle );
                    print('-dpsc2',  graphTitle );
                else
                    %print('-dpdf', '-append', graphTitle)
                    print('-dpsc2', '-append', graphTitle)
                end;
                close;
            end; %end of if printFig

        end; %end of loop through sec in an entity
        
        
        
        bands_needed=[1 4 8 14 30 50];
        %for band=1:20
        for band=1:(length( bands_needed)-1)

             %get indices of psd corresponding to current freq band
             index_t=( bands_needed(band) <= f(1:200) & f(1:200) <= bands_needed(band+1)) ;
            for sec = 1:nseconds

                %32768, vector means add 0 elements along first dim, X alond 2nd
                %pad index_t with zeros, since psd is 2-sided
                %index = padarray(index_t,[0 length(bic(sec).psd)-length(index_t) ], 'post' );
                %only one bic channel, compute average pw at each freq. band
                %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);
                if ~isempty( bic(sec).file(curFile).ent(ent).psd(index_t) )
                    bic(sec).file(curFile).ent(ent).avpw(band) =...
                        mean( bic(sec).file(curFile).ent(ent).psd(index_t)  );
                else
                    bic(sec).file(curFile).ent(ent).avpw(band) = -1 ;
                end;
                
            end;

        end;  %end of loop through bands

        
        

    end;%for loop through entitites
    
    
    
    
    %now we want to loop through all channels 1 sec at a time
    bands_needed=[1 4 8 14 30 50];
    numBands = length( bands_needed)-1 ;
    bandM=cell(1,5);     

        %for band=1:20
        for band=1:numBands
            bandM{band}=zeros(nseconds, 60) ;
             %get indices of psd corresponding to current freq band
             index_t=( bands_needed(band) <= f(1:200) & f(1:200) <= bands_needed(band+1)) ;
             
          for ent= a(1:60) ;
            for sec = 1:nseconds

                %32768, vector means add 0 elements along first dim, X alond 2nd
                %pad index_t with zeros, since psd is 2-sided
                %index = padarray(index_t,[0 length(bic(sec).psd)-length(index_t) ], 'post' );
                %only one bic channel, compute average pw at each freq. band
                %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);
                bandM{band}(sec,(ent-60)) =...
                    bic(sec).file(curFile).ent(ent).avpw(band);
            end;
          end; %loop thru ents
        end;  %end of loop through bands


    
end; %loop through all files

%filename='F:\EXTRAP\Matlab_workspace\9_18_2013';
%save(filename)

















%%%%%%%++++   CAR     +++++++++++++
 
%%%%%%%%%%%%%%%%%%%BIC%%%%%%%%%%%%%%%%%%%%%%%%%%%
%personal computer
% file1='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK0001.mcd';
% file2='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK0001.mcd';
% file3=strcat( 'L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\',...
%                 'September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK0001.mcd')
% file4='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK0001.mcd';
% file5='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK0002.mcd';

% EPA computer
 file1='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK0001.mcd';
 file2='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK0001.mcd';
 file3=strcat( 'L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\',...
                 'September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK0001.mcd')
 file4='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK0001.mcd';
 file5='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK0002.mcd';




files.car={file1, file2, file3, file4, file5} ;

for curFile=1:length(files.car)
    
    file = files.car{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    mcs_Info(hfile)

    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot

    spikes.n=0 ;
    for ent = 1:60

        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;

        spikes.car(ent).file(curFile).n = entity.ItemCount;
        index_ent(ent) =(spikes.car(ent).file(curFile).n >25 ) ;      
        
        for nspike = 1:entity.ItemCount

            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            spikes.car(ent).file(curFile).train(nspike) = timestamp ;

        end; %end of loop through timestamp in given entity

    end; %end of loop through all entities


    %get wells that are at least 5spikes/minute = 25spikes in 5mins

    
    %print
    graphTitle='Graph_BIC.ps' ;
    
    printFig=false;
    a = 61:120 ;
    for ent=a(index_ent)
        %ent=64;
        %ent=64 is best, since spike train 4 is best
        [nsresult,entity] = ns_GetEntityInfo(hfile,ent)
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,ent,1,entity.ItemCount);

        for sec=1:30


            %take a chunk of length 1 second
            fs = 25000;
            curData=data((sec-1)*fs+1:sec*fs);

            L=length(curData) ; %7,540,000
            NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L

            X = fft(curData, NFFT)/L;

            f = fs/2*linspace(0,1,NFFT/2+1); 

            X_mag = abs(X);


            car(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

            
            
           if printFig; 
                %open new figure
                figure; %this opens a blank figure
                hold off;

                subplot(3,1,1)

                index=(spikes.bic(ent-60).train<=sec)&(spikes.bic(ent-60).train>sec-1);
                in=zeros(1,1)
                in(1)=sec-1;
                in(2+ sum(index) )=sec;

                if sum(index)~=0
                    in(2:sum(index)+1 )=spikes.bic(ent-60).train(index);
                end;


                axis([sec-1 sec -1 2])
                plot([in;in],[ones(size(in));zeros(size(in))],'k-')



                subplot(3,1,2)
                plot(25000/2*curData);
                ylim([-0.25 0.3])
                hold on;

                %we're plotting the low frequency 1-40Hz over top
                %
                t = 0 : 1/fs : 1 - 1/fs ;
                co=8; %make the cut off 8Hz
                sig=2*co*sinc(2*co*t);
                %taking first 40 frequency is equivalent to
                %convolving original signal with sig above
                v=conv(sig, curData);

                plot(v(1:25000), 'r')
                ylim([-0.25 0.3])
                %V = fft(v, NFFT)/L;  V_mag = abs(V);

                %plot(f(2:120), 2*V_mag(2:120 ))
                %ylim([0 0.1])

                % plot(f, 2*X_mag(1:NFFT/2 + 1)) %frequency (Hz)
                subplot(3,1,3)



                plot(f(2:120), 2*X_mag(2:120 )) %frequency (Hz)
                ylim([0 0.000002])
                hold on;
                xlim([1,100])
                title(strcat('PSD of 090413 20607 BIC Second ',...
                    num2str(sec),'Entity',num2str(ent) ) )
                xlabel('Frequency (Hz)')
                ylabel('Magnitude = amplitude ^2')

                if sec==1
                    %print('-dpdf',  graphTitle );
                    print('-dpsc2',  graphTitle );
                else
                    %print('-dpdf', '-append', graphTitle)
                    print('-dpsc2', '-append', graphTitle)
                end;
                close;
            end; %end of if printFig

        end; %end of loop through sec in an entity
        
        
        bands_needed=[1 4 8 14 30 50];
        %for band=1:20
        for band=1:(length( bands_needed)-1)

             %get indices of psd corresponding to current freq band
             index_t=( bands_needed(band) <= f(1:200) & f(1:200) <= bands_needed(band+1)) ;
            for sec = 1:30

                %32768, vector means add 0 elements along first dim, X alond 2nd
                %pad index_t with zeros, since psd is 2-sided
                %index = padarray(index_t,[0 length(bic(sec).psd)-length(index_t) ], 'post' );
                %only one bic channel, compute average pw at each freq. band
                %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);
                car(sec).file(curFile).ent(ent).avpw(band) =...
                    mean( car(sec).file(curFile).ent(ent).psd(index_t)  );
            end;

        end;  %end of loop through bands


    end;%for loop through entitites


    
end; %loop through all files

%filename='F:\EXTRAP\Matlab_workspace\9_18_2013';
%save(filename)













%%%%%%%%%%%%%%%%%%%   CONTROL            %%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONTROL
% file1='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September 4, 2013 #20607 BIC\090413 20607 BIC FPSPK.mcd';
% file2='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 8395 BIC\090513 8395 BIC FPSPK.mcd';
% file3='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 15152 BIC\090513 15152 BIC FPSPK.mcd' ;
% file4='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 18405 BIC\090513 18405 BIC FPSPK.mcd';
% file5='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK.mcd';
% file6='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK.mcd';
% 

%EPA computer
file1='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September 4, 2013 #20607 BIC\090413 20607 BIC FPSPK.mcd';
file2='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 8395 BIC\090513 8395 BIC FPSPK.mcd';
file3='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 15152 BIC\090513 15152 BIC FPSPK.mcd' ;
file4='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5, 2013 BIC\090513 18405 BIC\090513 18405 BIC FPSPK.mcd';
file5='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK.mcd';
file6='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK.mcd';


            
file7=strcat( 'L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\',...
                'September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK.mcd')
% 
file8 = 'L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK.mcd' ;
file9='L:\Lab\NHEERL_MEA\Extrap\EXTRAP\FP SPIKE\09-10  2013\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK.mcd';


files.con={file1, file2, file3, file4, file5, file6, file7, file8 } ;

for curFile=1:length(files.con)
    
    file = files.con{curFile} ;
    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    mcs_Info(hfile)

    %go through Segment data and get the 1. number of spikes a.k.a
    %entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
    %to get 'timestamp', compose these to get one raster plot

    spikes.n=0 ;
    for ent = 1:60

        % "segment", entitly type is 3 we use the below tools
        [nsresult,entity] = ns_GetEntityInfo(hfile, ent);
        spikes.n=spikes.n+1 ;

        spikes.con(ent).file(curFile).n = entity.ItemCount;
        index_ent(ent) =(spikes.con(ent).file(curFile).n >25 ) ;      
        
        for nspike = 1:entity.ItemCount

            %sample count is # of samples in waveform
            [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
            spikes.con(ent).file(curFile).train(nspike) = timestamp ;

        end; %end of loop through timestamp in given entity

    end; %end of loop through all entities


    %get wells that are at least 5spikes/minute = 25spikes in 5mins

    
    %print
    graphTitle='Graph_BIC.ps' ;
    
    printFig=false;
    a = 61:120 ;
    for ent=a(index_ent)
        %ent=64;
        %ent=64 is best, since spike train 4 is best
        [nsresult,entity] = ns_GetEntityInfo(hfile,ent)
        %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
        [nsresult,count,data]=ns_GetAnalogData(hfile,ent,1,entity.ItemCount);

        for sec=1:30


            %take a chunk of length 1 second
            fs = 25000;
            curData=data((sec-1)*fs+1:sec*fs);

            L=length(curData) ; %7,540,000
            NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L

            X = fft(curData, NFFT)/L;

            f = fs/2*linspace(0,1,NFFT/2+1); 

            X_mag = abs(X);


            con(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);

            
            
           if printFig; 
                %open new figure
                figure; %this opens a blank figure
                hold off;

                subplot(3,1,1)

                index=(spikes.con(ent-60).train<=sec)&(spikes.con(ent-60).train>sec-1);
                in=zeros(1,1)
                in(1)=sec-1;
                in(2+ sum(index) )=sec;

                if sum(index)~=0
                    in(2:sum(index)+1 )=spikes.con(ent-60).train(index);
                end;


                axis([sec-1 sec -1 2])
                plot([in;in],[ones(size(in));zeros(size(in))],'k-')



                subplot(3,1,2)
                plot(25000/2*curData);
                ylim([-0.25 0.3])
                hold on;

                %we're plotting the low frequency 1-40Hz over top
                %
                t = 0 : 1/fs : 1 - 1/fs ;
                co=8; %make the cut off 8Hz
                sig=2*co*sinc(2*co*t);
                %taking first 40 frequency is equivalent to
                %convolving original signal with sig above
                v=conv(sig, curData);

                plot(v(1:25000), 'r')
                ylim([-0.25 0.3])
                %V = fft(v, NFFT)/L;  V_mag = abs(V);

                %plot(f(2:120), 2*V_mag(2:120 ))
                %ylim([0 0.1])

                % plot(f, 2*X_mag(1:NFFT/2 + 1)) %frequency (Hz)
                subplot(3,1,3)



                plot(f(2:120), 2*X_mag(2:120 )) %frequency (Hz)
                ylim([0 0.000002])
                hold on;
                xlim([1,100])
                title(strcat('PSD of 090413 20607 BIC Second ',...
                    num2str(sec),'Entity',num2str(ent) ) )
                xlabel('Frequency (Hz)')
                ylabel('Magnitude = amplitude ^2')

                if sec==1
                    %print('-dpdf',  graphTitle );
                    print('-dpsc2',  graphTitle );
                else
                    %print('-dpdf', '-append', graphTitle)
                    print('-dpsc2', '-append', graphTitle)
                end;
                close;
            end; %end of if printFig

        end; %end of loop through sec in an entity
        
        
        
        bands_needed=[1 4 8 14 30 50];
        %for band=1:20
        for band=1:(length( bands_needed)-1)

             %get indices of psd corresponding to current freq band
             index_t=( bands_needed(band) <= f(1:200) & f(1:200) <= bands_needed(band+1)) ;
            for sec = 1:30

                %32768, vector means add 0 elements along first dim, X alond 2nd
                %pad index_t with zeros, since psd is 2-sided
                %index = padarray(index_t,[0 length(bic(sec).psd)-length(index_t) ], 'post' );
                %only one bic channel, compute average pw at each freq. band
                %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);
                con(sec).file(curFile).ent(ent).avpw(band) =...
                    mean( con(sec).file(curFile).ent(ent).psd(index_t)  );
            end;

        end;  %end of loop through bands
        
        

    end;%for loop through entitites


    
end; %loop through all files



%filename='F:\EXTRAP\Matlab_workspace\9_18_2013';
%filename='C:\Users\Diana\Dropbox\Extrap\lfp_9_18_2013' ;
%save(filename)
























%%%%%%%%%%%%%%++++++++ ANALYSIS  ++++++++++++++

%difference in avpw
%question becomes is average power different between Con, Car and BIC
% or is it just different depending on #spikes
% i could try regression with #spikes, and factor variables for CON,
% CAR and BIC
% or ANCOVA

%files.con{5:6} %are for CAR, 1:4 are for BIC
for band=1:5
    
    Avpw.band(band).car.data=zeros(30,2);
    row.band(band)=0;
    %car files
    for curFile=1:2 %note that car has 2 files: files.car
        
        
        for ent=61:120;
                    %determine which ents to use
        %index_bic(ent) =(spikes.bic(ent-60).file(curFile).n >25 ) ;
        index_car(ent) =(spikes.car(ent-60).file(curFile).n >25 ) ;
        %note: you add 4 to curFile since curFile=5 and 6 are the car
        %control
        index_con(ent) =(spikes.con(ent-60).file(curFile+4).n >25 ) ;
        if (index_car(ent) & index_con(ent) )
                  
            for sec=1:30
                
                %get current row
                row.band(band)=row.band(band)+1;
                
                %con(sec).file(curFile).ent(ent).avpw(band)
                %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);
                %con(sec).file(curFile).ent(ent).avpw(band)
                
                temp = [con(sec).file(curFile+4).ent(ent).avpw(band)...
                    car(sec).file(curFile).ent(ent).avpw(band) ] ;
                Avpw.band(band).car.data(row.band(band),:) = temp ;
            end;%loop through seconds
            
        end; %if all have valid ents
        end; %loop through entities
    end; %loop through files
    
    
    
    Avpw.band(band).bic.data=zeros(30,2);
    row.band(band)=0;
    
    %car files
    for curFile=1:4 %note that car has 2 files: files.car
        
        
        for ent=61:120;
                    %determine which ents to use
        index_bic(ent) =(spikes.bic(ent-60).file(curFile).n >25 ) ;
        %index_car(ent) =(spikes.car(ent-60).file(curFile).n >25 ) ;
        %note: you add 4 to curFile since curFile=5 and 6 are the car
        %control
        index_con(ent) =(spikes.con(ent-60).file(curFile).n >25 ) ;
        if (index_bic(ent) & index_con(ent) )
                  
            for sec=1:30
                
                %get current row
                row.band(band)=row.band(band)+1;
                
                %con(sec).file(curFile).ent(ent).avpw(band)
                %bic(sec).file(curFile).ent(ent).psd=2*X_mag(1:200);
                %con(sec).file(curFile).ent(ent).avpw(band)
                
                temp = [con(sec).file(curFile).ent(ent).avpw(band)...
                    bic(sec).file(curFile).ent(ent).avpw(band) ] ;
                Avpw.band(band).bic.data(row.band(band),:) = temp ;
            end;%loop through seconds
            
        end; %if all have valid ents
        end; %loop through entities
    end; %loop through files
    
    
end; %end of loop through bands


%filename='F:\EXTRAP\Matlab_workspace\9_18_2013';
%filename='C:\Users\Diana\Dropbox\Extrap\lfp_9_18_2013' ;
%save(filename)



%%Anova average difference between CON CAR and BIC

%bands_needed= 1-4Hz; 4-8Hz; 8- 14Hz; 14-30Hz; 30-50Hz
% BIC
for band=1:5
    %take the log 
    [p, anovatab, stats] = anova1( log( Avpw.band(band).bic.data(:,1:2) ) );
    stats  
    pause ;
    close all;
    
end;


%%Anova average difference between CAR and BIC
for band=1:20
    [p, anovatab, stats] = anova1( log( Avpw.band(band).car.data(:,1:2) ) );
    stats  
    pause ;
    close all;
    
end;








%now, do ancova taking into account spikes
%aoctool(x,y,group,alpha,xname,yname,gname)
% x=ind variable is = # spikes in each sec
% y=dependent variable = avg. power
% group= categorical variable = CON, CAR or BIC

%we need to gather #spikes in each sec interval
%which channel we chose for each treatment
conEnt=27; carEnt=53; bicEnt=4;
for sec=1:30
    
    %get # spikes at each second
    %CAR
    index=(spikes.car(carEnt).train<=sec)&(spikes.car(carEnt).train>sec-1);    
    nspikes.car(sec)=sum(index) ;
    
    %CON
    index=(spikes.con(conEnt).train<=sec)&(spikes.con(conEnt).train>sec-1);    
    nspikes.con(sec)=sum(index) ;
    
    %BIC
    index=(spikes.bic(bicEnt).train<=sec)&(spikes.bic(bicEnt).train>sec-1);    
    nspikes.bic(sec)=sum(index) ;
    
end;


%join #spikes at a given second with other data
for band=1:20
        
    %x=dependent variable=number of spikes
    temp1 = [nspikes.con  nspikes.car  nspikes.bic ]' ;
    
    %y=temp2=dep. variable= average power at each second 1-30
    temp2 =  vertcat(Avpw.band(band).data(:,1),...
        Avpw.band(band).data(:,2),Avpw.band(band).data(:,3) ) ;

    % CON=1, CAR=2, BIC=3
    temp3 =[repmat(1,1,30) repmat(2,1,30) repmat(3,1,30 )]' ;
    
    %joing data in column x, y , group
    Avpw.band(band).ancovaData = horzcat(temp1, temp2, temp3);
    
      
    
end;


%ancova is significant at the 0.05 level at band=20;
band=1;
aoctool(Avpw.band(band).ancovaData(:,1) ,...
    Avpw.band(band).ancovaData(:,2),...
    Avpw.band(band).ancovaData(:,3),0.05,'spikes',...
    'average power',horzcat('Control','Carboryl', 'BIC') );

%try permutation test?
%also, your results may be different just because average voltage
%was somehow a bit different, but not actual chemicals
%try standardizing each variable to get homogeneous firing rates
























%%%%%%%%%%%%%+++++++ Old Code for Graphing ++++++%%%%%%%%%%%%%%%%%%%%%%%%

[nsresult, hfile]=ns_OpenFile(file)
[nsresult,info]=ns_GetFileInfo(hfile)

% Find out about the Entity types
% Then read specific entity info and data

mcs_Info(hfile)

%go through Segment data and get the 1. number of spikes a.k.a
%entity.ItemCount 2. use this to loop through ns_GetSegmentData(hfile,3,i)
%to get 'timestamp', compose these to get one raster plot

spikes.n=0 ;
for ent = 1:60
    
% "segment", entitly type is 3 we use the below tools
[nsresult,entity] = ns_GetEntityInfo(hfile, ent);
spikes.n=spikes.n+1 ;
spikes.con(ent).n = 0 ;

spikes.con(ent).n = entity.ItemCount;

for nspike = 1:entity.ItemCount
    
    %sample count is # of samples in waveform
    [nsresult,timestamp,data,samplecount,unitid] =ns_GetSegmentData(hfile,ent,nspike);
    spikes.con(ent).train(nspike) = timestamp ;

end; %end of loop through timestamp in given entity

end; %end of loop through all entities




%print
graphTitle='Graph_Control.ps' ;

%for ent=61:80
    ent=44+60;
    %many spike trains are all zeros,
    %  ent=27+60; % is best
    [nsresult,entity] = ns_GetEntityInfo(hfile,ent)
    %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
    [nsresult,count,data]=ns_GetAnalogData(hfile,ent,1,entity.ItemCount);
     
    
    for sec=1:30
        %open new figure
        figure; %this opens a blank figure
        hold off;

        %take a chunk of length 1 second
        fs = 25000;
        curData=data((sec-1)*fs+1:sec*fs);

        L=length(curData) ; %25,000
        NFFT = 2^nextpow2(L) ; % gives you the power of 2 >=L

        X = fft(curData, NFFT)/L;

        f = fs/2*linspace(0,1,NFFT/2+1); 

        X_mag = abs(X);
        
        
        con(sec).psd = 2*X_mag ; 
        
        
        
        
        subplot(3,1,1)
        
        index=(spikes.con(ent-60).train<=sec)&(spikes.con(ent-60).train>sec-1);
        in=zeros(1,1)
        in(1)=sec-1;
        in(2+ sum(index) )=sec;
        
        if sum(index)~=0
            in(2:sum(index)+1 )=spikes.con(ent-60).train(index);
        end;
        
        
        axis([sec-1 sec -1 2])
        plot([in;in],[ones(size(in));zeros(size(in))],'k-')
     
        
        
        subplot(3,1,2)
        plot(25000/2*curData);
        ylim([-0.25 0.3])
        hold on;
        
        %we're plotting the low frequency 1-40Hz over top
        %
        t = 0 : 1/fs : 1 - 1/fs ;
        co=8; %make the cut off 8Hz
        sig=2*co*sinc(2*co*t);
        %taking first 40 frequency is equivalent to
        %convolving original signal with sig above
        v=conv(sig, curData);
       
        plot(v(1:25000), 'r')
        ylim([-0.25 0.3])
        %V = fft(v, NFFT)/L;  V_mag = abs(V);

        %plot(f(2:120), 2*V_mag(2:120 ))
        %ylim([0 0.1])

        % plot(f, 2*X_mag(1:NFFT/2 + 1)) %frequency (Hz)
        subplot(3,1,3)
        
        
        
        plot(f(2:120), 2*X_mag(2:120 )) %frequency (Hz)
        ylim([0 0.000002])
        hold on;
        xlim([1,100])
        title(strcat('PSD of Control Second ',...
            num2str(sec),'Entity',num2str(ent) ) )
        xlabel('Frequency (Hz)')
        ylabel('Magnitude = amplitude ^2')

        if sec==1
            print('-dpsc2',  graphTitle );
        else
            print('-dpsc2', '-append', graphTitle)
        end;
        close;


    end; %for loop
    
    
%end;






%compute average power over 1-100Hz, in 5Hz increments,
%total of 20 bands per second
for band=1:20
    
     %get indices of psd corresponding to current freq band
     index_t=((band-1)*5+1<= f & f <= band*5) ;
    for sec = 1:30
       
        %32768, vector means add 0 elements along first dim, X alond 2nd
        %pad index_t with zeros, since psd is 2-sided
        index = padarray(index_t,[0 length(con(sec).psd)-length(index_t) ], 'post' );
        %only one bic channel, compute average pw at each freq. band
        con(sec).avpw(band) = mean( con(sec).psd(index) );
    end;
    
end;











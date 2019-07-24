% make_raw_spike_psd_figure.m
% Diana Hall
% 11/05/2013
% purpose: to make raw, spike, psd figure


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


ns_GetLibraryInfo


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

% choose parameters for plot
curFile=1; ent=64; % this has a lot of activity
file = files.bic{curFile} ;

% choose file name for graph
[pathstr,name,ext]=fileparts(files.bic{1});
graphTitle= strcat('3PartFig_BIC_Channel', num2str(ent-60),'_Date',num2str(name(1:6) ),'.ps') ;
graphFile = strcat('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\',graphTitle) ;  
% spikes.bic(ent).file(curFile).train(nspike) gives you info on activity

    
    [nsresult, hfile]=ns_OpenFile(file)
    [nsresult,info]=ns_GetFileInfo(hfile)

    % Find out about the Entity types
    % Then read specific entity info and data

    
    temp = mcs_Info(hfile);
    nseconds = floor( str2num( temp{1,5} ) ) ;


  % load data,       
 [nsresult,entity] = ns_GetEntityInfo(hfile,ent)
 %[nsresult,analog] = ns_GetAnalogInfo(hfile,ent)
 [nsresult,count,data]=ns_GetAnalogData(hfile,ent,1,entity.ItemCount);

%make a plot for each second
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

 end; %end of loop through sec in an entity
        
  

















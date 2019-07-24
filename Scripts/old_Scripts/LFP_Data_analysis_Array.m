%LFP_Data_analysis.m
%diana hall
%Oct 23, 2013
% purpose: to read in data that was collected on MCD by cina,
%  that contain low frequency and high frequency signals
%  we do a PSD and then apply an ANOVA to average power
%  over a given frequency (control-baseline



% summary of code
% read in all files, 2. get timestamp data, determine which channels
% have greater than 5spikes/minute of data 3. get LFP data for those 
% channels with higher than 5spikes/minute including FFT
%  4. do anova



%%%%%++++++++++  Electrode level analysis +++++++++++++++++
load 'F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\11_19_2013_all.mat' ;


 for band=1:5
     %take the log 
     [p, anovatab, stats] = anova1( log( Avpw.band(band).bic.data(:,1:2) ) );
     stats  
     pause ;
     close all;
     
 end;
% 
% 
 %%Anova average difference between CAR and CON
 for band=1:5
     [p, anovatab, stats] = anova1( log( Avpw.band(band).car.data(:,1:2) ) );
     stats  
     pause ;
     close all;
     
 end;





 
 
 
 
 
%%++++++++++++++++++band M analysis SECOND LEVEL ++++++++++++++++ 
 
 
%data prep
% File 1-4 are BIC,
% elimate channels below 5spikes/min
% bicANOVA is a num Files X bands*2(1 first cols for Control,2nd BIC)
% matrix
format short eng %use to see power with many 0s
% make BICanova_sec 

 BICanova_sec_temp = cell(1,5) ;
 BICanova_sec = cell(1,5) ;
 for band=1:(length(bands_needed)-1 ) ;
      temp= zeros(1,2) ;
      for nfile=1:length(files.bic) 
         for ch=1:60
             % all control (1st column), BIC 2nd column
             temp = vertcat(temp,...
                 horzcat( apwCON(nfile).bandM(:,ch,band),...
                  apwBIC(nfile).bandM(:,ch,band) ) ) ;
         end; % go through channels
      end; % go through nfile
     BICanova_sec_temp{band} = temp ;
     % get non-zero entries
     %for control
     [rowCON colCON vCON] = find( BICanova_sec_temp{band}(:,1)) ;
     %for BIC
     [rowCAR colCAR vCAR] = find( BICanova_sec_temp{band}(:,2)) ;
     %indexNZ is what was non-zero in both control and BIC
     indexNZ = intersect(rowCON, rowCAR) ;
     BICanova_sec{band} = BICanova_sec_temp{band}(indexNZ,:) ;
     
 end; %end of loop through band


% make CARanova_sec 
 CARanova_sec_temp = cell(1,5) ;
 CARanova_sec = cell(1,5) ;
 nCARfiles = length(files.car)-1 ;
 nBICfiles = length(files.bic) ;
 for band=1:(length(bands_needed)-1 ) ;
      temp= zeros(1,2) ;
      
      for nfile=1:nCARfiles 
         for ch=1:60
             % all control (1st column), CAR 2nd column
            % we have to offset control files since BIC control 
             % preceed CAR control
             nseconds = min ( size(apwCON(nfile+nBICfiles).bandM(:,ch,band),1 ),... 
                 size( apwCAR(nfile).bandM(:,ch,band), 1 ) ) ;
             temp = vertcat(temp,...
                 horzcat( apwCON(nfile+nBICfiles).bandM(1:nseconds,ch,band),...
                 apwCAR(nfile).bandM(1:nseconds,ch,band) ) ) ;
         end; % go through channels
      end; % go through nfile
     CARanova_sec_temp{band} = temp ;
     % get non-zero entries
     %for control
     [rowCON colCON vCON] = find( CARanova_sec_temp{band}(:,1)) ;
     %for BIC
     [rowCAR colCAR vCAR] = find( CARanova_sec_temp{band}(:,2)) ;
     %indexNZ is what was non-zero in both control and BIC
     indexNZ = intersect(rowCON, rowCAR) ;
     CARanova_sec{band} = CARanova_sec_temp{band}(indexNZ,:) ;
     
 end; %end of loop through band



%%Anova average difference between CAR and CON

% %bands_needed= 1-4Hz; 4-8Hz; 8- 14Hz; 14-30Hz; 30-50Hz
% % BIC
 for band=1:5
     %take the log 
     [p, anovatab, stats] = anova1( log( BICanova_sec{band} ) );
     stats  
     pause ;
     close all;
     
 end;
% 
% 
 %%Anova average difference between CAR and CON
 for band=1:5
     [p, anovatab, stats] = anova1( log( CARanova_sec{band} ) );
     stats  
     pause ;
     close all;
     
 end;
% 
% 
%%++++++++++ CAR Log plot +++++++++++
figure() %opens a figure
for band=1:5    

    subplot(1,5,band)
    boxplot( log( CARanova_sec{band}) ,...
        'label', {'Control', 'CAR'} ) 
    
    hold on ;  
    if band==3
       title({'Log(Average Power) in Carboryl vs. Control by Sec',...
                'Electrodes >= 5 spikes/min'}  ) ;
    elseif band==1
        ylabel ('log power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxplotLogCAR_AllBands_seconds.pdf');  
print( '-dpdf',  plotName) ;

plotName=strcat('F:\EXTRAP\figures\matlab\boxplotLogCAR_AllBands_seconds.tiff');  
print( '-dtiffn',  plotName) ;


%%++++++++++ BIC Log plot +++++++++++
figure() %opens a figure
for band=1:5    

    subplot(1,5,band) ;
    
    boxplot( log( BICanova_sec{band} ) ,...
        'label', {'Control', 'BIC'} ) 
    hold on ;  
    if band==3
       title({'Log(Average Power) in BIC vs. Control by Sec',...
                'Electrodes >= 5 spikes/min'}  ) ;
    elseif band==1
        ylabel ('log power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxplotLogBIC_AllBands_seconds.pdf');  
print( '-dpdf',  plotName) ;

plotName=strcat('F:\EXTRAP\figures\matlab\boxplotLogBIC_AllBands_seconds.tiff');  
print( '-dtiffn',  plotName) ;











%%%%% ++++++++++++++++chip level analysis +++++++++++++++++++++++++
% now do average across all channels anova

%data prep
% File 1-4 are BIC,
% elimate channels below 5spikes/min
% bicANOVA is a num Files X bands*2(1 first cols for Control,2nd BIC)
% matrix
format short eng %use to see power with many 0s

 BICanova = zeros(length(files.bic), (length(bands_needed)-1)*2  ) ;
for nfile=1:length(files.bic)
     for band=1:(length(bands_needed)-1 ) ;

        %returns all non-zero entries, since entries at the end of some
        % non-zero columns are zeros
        % r is vectors of row lengths, c corresponding column indices
        %  v is all non-zero entries
        [r,c, v ]=find( apwBIC(nfile).bandM(:,:,band) ) ;
        wellBIC(band) = mean( v ) ;
    
    	[rCon,cCon, vCon ]=find( apwCON(nfile).bandM(:,:,band) ) ;
        wellCON(band) = mean( vCon ) ;
     end; 
     % 2=BIC, 1=control
     BICanova(nfile,:) = horzcat( wellCON, wellBIC)  ;
end;

nameAnova = [strcat(bandChar,' CON'), strcat(bandChar, ' BIC') ] ;


%%%%  anova data set for CAR %%%%
 CARanova = zeros(length(files.car)-1, (length(bands_needed)-1)*2  ) ;
 % last CAR file is distorted
for nfile=1:length(files.car)-1
     for band=1:(length(bands_needed)-1 ) ;

        %returns all non-zero entries, since entries at the end of some
        % non-zero columns are zeros
        % r is vectors of row lengths, c corresponding column indices
        %  v is all non-zero entries
        [r,c, v ]=find( apwCAR(nfile).bandM(:,:,band) ) ;
        wellCAR(band) = mean( v ) ;
    
    	[rCon,cCon, vCon ]=find( apwCON(nfile).bandM(:,:,band) ) ;
        wellCON(band) = mean( vCon ) ;
     end; 
     % 2=BIC, 1=control
     CARanova(nfile,:) = horzcat( wellCON, wellCAR)  ;
end;

nameAnova = [strcat(bandChar,' CON'), strcat(bandChar, ' CAR') ] ;



%%Anova average difference between CAR and CON
bandChar={'1-4Hz','4-8Hz', '8-14Hz', '14-30Hz','30-50Hz'} ;
for band=1:5

    %control is first 5 columns, trt is the 2nd 5 columns
    [p, anovatab, stats] = anova1( CARanova(:, [band, band+5]) ) ;
 
    
    
    pause ;
    close all;
    
end;

%%++++++++++ CAR Log plot +++++++++++
figure() %opens a figure
for band=1:5    

    subplot(1,5,band) ;
    names=repmat({'Control', 'CAR'},1,4) ;
    %scatter( [1 1 1 1 2 2 2 2],...
   %     vertcat(CARanova(:,1), log( CARanova(:,[1+band] ) )) )  ,...
     %   'label', {'Control', 'CAR'} ) ;
    
    
    boxplot(  [CARanova(:,1), log( CARanova(:,[1+band] ) )],...
        'label', {'Control', 'CAR'} ) 
    hold on ;
    
    if band==3
       title({'Log(Average Power) in Carboryl vs. Control n=4 chips',...
                'Electrodes >= 5 spikes/min'}  ) ;
    elseif band==1
        ylabel ('log power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxplotLogCAR_AllBands.pdf');  
print( '-dpdf',  plotName) ;


    
%+++++++++++++++++++++ CAR no log
figure() %opens a figure
for band=1:5    

    subplot(1,5,band) ;
    names=repmat({'Control', 'CAR'},1,4) ;
    boxplot(  BICanova(:,[1, 1+band]  ),...
        'label', {'Control', 'CAR'} ) 
    hold on ;
    
    if band==3
       title({'Average Power in Carboryl vs. Control n=4 chips',...
                'Electrodes >= 5 spikes/min'}  ) ;
    elseif band==1
        ylabel ('power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxCAR_AllBands.pdf');  
print( '-dpdf',  plotName) ;




    
    
    
 %%++++++++++++++Anova average difference between BIC and CON
bandChar={'1-4Hz','4-8Hz', '8-14Hz', '14-30Hz','30-50Hz'} ;
for band=1:5

    group = regexprep( cellstr( int2str(CARanova(:,1) )), '1', 'Control' )
    group = regexprep( group, '2', 'BIC' )
    [p, anovatab, stats] = anova1( [BICanova(:,1) log( BICanova(:,[ 1+band]) )] ) ;
    stats  

    
    pause ;
    close all;
    
end;



% +++++++++ log BIC ++++++++++++++
figure() %opens a figure
for band=1:5    

    subplot(5,1,band) ;
    boxplot(  [BICanova(:,1) log( BICanova(:,[ 1+band]) )] ) )
    
    hold on ;
    if band==1
    title({'Average Power in Bicuculline vs. Control n=4 chips',...
                'n= 4 wells'}  ) ;
    elseif band==5
           xlabel(cellstr( ['Control ';'Bicuculline'] ) ) ; 
    end;
  
    ylabel ('log power' ) ;

end;
    
plotName=strcat('F:\EXTRAP\figures\matlab\boxLogBIC_AllBands.pdf');  
print( '-dpdf',  plotName) ;
    
    
    
 %+++++++++++++++++++++ BIC no log
figure() %opens a figure
for band=1:5    

    subplot(1,5,band) ;
    names=repmat({'Control', 'BIC'},1,4) ;
    boxplot(  BICanova(:,[1, 1+band]  ),...
        'label', {'Control', 'BIC'} ) 
    hold on ;
    
    if band==3
       title({'Average Power in Bicuculline vs. Control n=4 chips',...
                'Electrodes >= 5 spikes/min'}  ) ;
    elseif band==1
        ylabel ('power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxBIC_AllBands.pdf');  
print( '-dpdf',  plotName) ;














    
    
    
    
    

%%%%% ++++++++++++++Electrode level analysis +++++++++++++++++++++++++
% channel level analysis

%data prep
% File 1-4 are BIC,
% elimate channels below 5spikes/min
% allBIC is a 3 dim array (nrows=chips, ncols=first 60:CON, 
% 61-120:BIC, 3rd dim = bands 1-5
format short eng %use to see power with many 0s

for nfile=1:length(files.bic)
     for band=1:(length(bands_needed)-1 ) ;

        %returns all non-zero entries, since entries at the end of some
        % non-zero columns are zeros
        % r is vectors of row lengths, c corresponding column indices
        %  v is all non-zero entries
        temp = size(apwBIC(nfile).bandM) ;
        numElect= temp(2) ;

        for elect=1:numElect
            %get entries without trailing 0s
            [r,c, v ] = find( apwBIC(nfile).bandM(:, elect, band) ) ;
            tempBIC(elect) = mean(v) ;
            
            [rCon,cCon, vCon ]=find( apwCON(nfile).bandM(:,elect,band) ) ;
            tempCON(elect) = mean( vCon ) ;
            
            
        end; %end loop through electrodes
        tempAllBIC( (1+(nfile-1)*numElect) : ( nfile*numElect) ,:,band) =...
            horzcat(transpose( tempCON),transpose( tempBIC ) ) ;
    
      end;  %end loop through bands

end;

% we have to make it into a list becuase of the NaN cause variable 
% data lengths
ylimBIC = 0;
for band=1:5
    index{band} = ~isnan(tempAllBIC(:,1,band ) ) & ~isnan(tempAllBIC(:,2,band)) ;
    allBIC{band}= tempAllBIC(index{band} ,:,band ) ;
    temp = max(allBIC{band} ) ;
    ylimBIC = max( max( ylimBIC, temp) ) ;

end;



%%%%  anova data set for CAR %%%%

 % last CAR file is distorted

for nfile=1:length(files.car)-1
     for band=1:(length(bands_needed)-1 ) ;

        %returns all non-zero entries, since entries at the end of some
        % non-zero columns are zeros
        % r is vectors of row lengths, c corresponding column indices
        %  v is all non-zero entries
        temp = size(apwCAR(nfile).bandM) ;
        numElect= temp(2) ;

        for elect=1:numElect
            %get entries without trailing 0s
            [r,c, v ] = find( apwCAR(nfile).bandM(:, elect, band) ) ;
            tempCAR(elect) = mean(v) ;
            
            [rCon,cCon, vCon ]=find( apwCON(nfile).bandM(:,elect,band) ) ;
            tempCON(elect) = mean( vCon ) ;
            
            
        end; %end loop through electrodes
        tempAllCAR( (1+(nfile-1)*numElect) : ( nfile*numElect) ,:,band) =...
            horzcat(transpose( tempCON),transpose( tempCAR ) ) ;
    
      end;  %end loop through bands

end;

%filter out rows (Control and BIC or CAR) where NaN appears

% we have to make it into a list becuase of the NaN cause variable 
% data lengths
ylimCAR = 0;
for band=1:5
    index{band} = ~isnan(tempAllCAR(:,1,band ) ) & ~isnan(tempAllCAR(:,2,band)) ;
    allCAR{band}= tempAllCAR(index{band} ,:,band ) ;
    temp = max(allCAR{band} ) ;
    ylimCAR = max( max( ylimCAR, temp) ) ;

end;


%%% The 148th is fucked up since it's a 1 in CAR
% I manually changed it to 0

%%Anova average difference between CAR and CON
bandChar={'1-4Hz','4-8Hz', '8-14Hz', '14-30Hz','30-50Hz'} ;
for band=1:5

    %control is first 5 columns, trt is the 2nd 5 columns
    [p, anovatab, stats] = anova1( allCAR{band} ) ;
 
    pause ;
    close all;
    
end;

%%++++++++++ CAR Log plot +++++++++++
figure() %opens a figure
for band=1:5    

    subplot(1,5,band) ;
    names=repmat({'Control', 'CAR'},1,4) ;
    boxplot(  log( allCAR{band} ) ,...
        'label', {'Control', 'CAR'} ) 
    hold on ;
    
    if band==3
       title( strcat( 'Log(Average Power) by Electrode in Carboryl vs. Control',...
           'n=', num2str(nElect) , ' electrodes > 5 spikes/min')  ) ;
    elseif band==1
        ylabel ('log power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxplotLogCAR_AllBands_Electrodes.pdf');  
print( '-dpdf',  plotName) ;


    
%+++++++++++++++++++++ CAR no log
figure() %opens a figure
for band=1:5    

    
    subplot(1,5,band) ;
    names=repmat({'Control', 'CAR'},1,4) ;
    boxplot(  allCAR{band} ,...
        'label', {'Control', 'CAR'} ) 
    axis([0 3 0 ylimCAR] ) ;
    hold on ;
    
    if band==3
       nElect = max( [ size(allBIC{1}),size(allBIC{1}),...
           size(allBIC{3}),size(allBIC{4}),size(allBIC{5}) ] ) ;
       title( strcat( 'Average Power by Electrode in Carboryl vs. Control',...
           'n=', num2str(nElect) , ' electrodes > 5 spikes/min')  ) ;
    elseif band==1
        ylabel ('power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxCAR_AllBands_Electrodes.pdf');  
print( '-dpdf',  plotName) ;




    
    
    
 %%++++++++++++++Anova average difference between BIC and CON
bandChar={'1-4Hz','4-8Hz', '8-14Hz', '14-30Hz','30-50Hz'} ;
for band=1:5

    [p, anovatab, stats] = anova1( log(allBIC{band} ) ) ;
    stats    
    pause ;
    close all;
    
end;



%%++++++++++ CAR Log plot +++++++++++
figure() %opens a figure
for band=1:5    

    subplot(1,5,band) ;
    names=repmat({'Control', 'BIC'},1,4) ;
    boxplot(  log( allBIC{band} ) ,...
        'label', {'Control', 'BIC'} ) 
    hold on ;
    
    if band==3
       title( strcat( 'Log(Average Power) by Electrode in BIC vs. Control',...
           'n=', num2str(nElect) , ' electrodes > 5 spikes/min')  ) ;
    elseif band==1
        ylabel ('log power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxLogBIC_AllBands_Electrodes.pdf');  
print( '-dpdf',  plotName) ;
    
    
    
 %+++++++++++++++++++++ BIC no log
figure() %opens a figure
for band=1:5    

    subplot(1,5,band) ;
    names=repmat({'Control', 'BIC'},1,4) ;
    boxplot(  BICanova(:,[1, 1+band]  ),...
        'label', {'Control', 'BIC'} ) 
    hold on ;
    
    if band==3
       title( strcat( 'Average Power by Electrode in BIC vs. Control',...
           'n=', num2str(nElect) , ' electrodes > 5 spikes/min')  ) ;
    elseif band==1
        ylabel ('power' ) ;
    end;
    xlabel(bandChar{band} ) ; 

end;

plotName=strcat('F:\EXTRAP\figures\matlab\boxBIC_AllBands_Electrodes.pdf');  
print( '-dpdf',  plotName) ;














    
    

% ++++++++++++ ANCOVA +++++++++++++++++++++++++++
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











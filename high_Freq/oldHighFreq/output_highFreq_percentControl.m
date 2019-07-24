%% output_highFreq_percentControl.m
% diana hall
% november 5, 2013
%%%Use NeuroShare
%you need the FIND toolbox and the NeuroShare



%%
mac=1;
if ~mac
    xlsfilestr='F:\EXTRAP\matlab_Cina_Herr_Diana\high_Freq\high_freq_03_01_2016_percentBL.xls';
    title='High Frequency Information for All Files';
    xlswrite(xlsfilestr,title, 1,'A1' ) ;


    file1='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_bic.mat';
    file2='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_biccon.mat';

    file3='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car5.mat';
    file4='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car5con.mat';

    file5='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car15.mat' ;
    file6='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car15con.mat' ;

    file7='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30.mat' ;
    file8='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30con.mat' ;

    file9='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\02_11_2014_per25.mat';
    file10='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_per25con.mat';

    file11='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_per50.mat' ;
    file12='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conPER50.mat' ;

    file13='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_da.mat';
    file14='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_conDA.mat';

    file15='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_tri.mat';
    file16='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_contri.mat';

    file17='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_ace.mat';
    file18='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conace.mat';

    file19='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_h20.mat' ;
    file20='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conh20.mat' ;

    file21='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_gly.mat' ;
    file22='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_congly.mat' ;

    file23='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_dmso.mat';
    file24='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_condmso.mat';

    all.filepath = {file1, file2, file3, file4, file5,...
        file6,file7,file8,file9,file10,file11,file12,file13,file14,...
        file15,file16, file17, file18, file19, file20, file21,file22,file23,file24}  ;

    startCell=5;

    % 
    for curFile=[2 4 6 8 10 12 14 16 18 20 22 24]
        %load control file
        s=load( all.filepath{curFile}  );
        name=fieldnames(s);

        chem = s.(name{1}) ;
        oChunk=[]; fileTitle=[];

        trtemp=load( all.filepath{curFile-1});
        tname=fieldnames(trtemp) ;
        trt=trtemp.(tname{1});

        % cur chemical name
        [a,b,c]=fileparts(all.filepath{(curFile-1) });
        d=strsplit(b,'_'); chemName= d( length(d) );
        for curPlate=1:length(chem.file)

            curNSeconds=chem.file(curPlate).nseconds ;
            curNSpikes=[]; curAE=[];
            for ent = 1:60
                curNSpikes =[curNSpikes  chem.file(curPlate).channel(ent).nspikes ] ;
                % indicator of ae or not
                curAE=[curAE chem.file(curPlate).channel(ent).ae ] ;   

            end; %end of loop through all entities

            mfr_all=sum(curNSpikes)/(60*curNSeconds) ;
            curPerAE= sum( curAE )/60 ;
            mfr_AE = dot(curAE,curNSpikes)/dot(curAE, curNSeconds*ones(1,60) );

            %now do treatment
            curNSecondsTrt=trt.file(curPlate).nseconds ;
            curNSpikesTrt=[]; curAETrt=[];
            for ent = 1:60
                curNSpikesTrt =[curNSpikesTrt  trt.file(curPlate).channel(ent).nspikes ] ;
                % indicator of ae or not
                curAETrt=[curAETrt trt.file(curPlate).channel(ent).ae ] ;
            end; %end of loop through all entities

            mfr_all=sum(curNSpikesTrt)/(60*curNSecondsTrt) ;
            % BIC files mfr computed using AE from trt others mfr calculations
            % use pre-trt files
            if curFile==2
                %divide number of total spikes on AE in trt by total # seconds
               mfr_AE_trt = dot(curAETrt,curNSpikesTrt)/dot(curAETrt, curNSecondsTrt*ones(1,60) );
            else 
                mfr_AE_trt = dot(curAE,curNSpikesTrt)/dot(curAE, curNSecondsTrt*ones(1,60) );
            end;

            ratio_mfr=mfr_AE_trt/mfr_AE ;
            ratio_nAE=sum(curAETrt)/sum(curAE) ;


            [p,n1,x] = fileparts( chem.filepath{curPlate} ) ;
            title= strcat( n1,' ' );

            fileTitle{curPlate}=title ;

            temp=vertcat(ratio_mfr , ratio_nAE );
            oChunk=horzcat(oChunk, temp);


        end;%loop through plates in file
        celldata= cellstr( ['mfr_all' ;'mfr_AE ' ;'% AE   '] );

        %%%%%output
        xlswrite(xlsfilestr, chemName ,1,['B' num2str(startCell-2)] ) ;
        xlswrite(xlsfilestr, fileTitle ,1,['B' num2str(startCell-1)] ) ;
        xlswrite(xlsfilestr, 'mfr_pb_onBaselineAEChannels' ,1,['A' num2str(startCell)] ) ;
        xlswrite(xlsfilestr, 'pb_nAE' ,1,['A' num2str(startCell)] ) ;
        xlswrite(xlsfilestr, oChunk,1,['B' num2str(startCell)] ) ;

        startCell=startCell+5;
    end;  % end of loop through files

end; % if mac



%% ++++++++++++++++++++mac+++++++++++++++++++++++++++++++++++++++++++++++
mac=1;
xlsfilestr='/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/high_Freq/high_freq_03_01_2016_percentBL.csv';
title1={'High Frequency Information for All Files'};
title=title1{1};
fileID = fopen(xlsfilestr,'w');
fprintf( fileID,'%s\n',  title );


root='/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/';

file1=strcat(root, '01_27_2014_bic.mat' );
file2=strcat(root, '10_30_2014_biccon.mat');

file3=strcat(root, '01_02_2015_car5.mat');
file4=strcat(root, '01_02_2015_car5con.mat');

file5=strcat(root, '01_02_2015_car15.mat') ;
file6=strcat(root, '01_02_2015_car15con.mat') ;

file7=strcat(root, '01_02_2015_car30.mat') ;
file8=strcat(root, '01_02_2015_car30con.mat') ;

file9=strcat(root, '02_11_2014_per25.mat');
file10=strcat(root, '10_30_2014_per25con.mat');

file11=strcat(root, '01_02_2015_per50.mat' );
file12=strcat(root, '01_02_2015_conPER50.mat') ;

file13=strcat(root, '10_30_2014_da.mat');
file14=strcat(root, '10_30_2014_conDA.mat');

file15=strcat(root, '8_9_2014_tri.mat');
file16=strcat(root, '8_9_2014_contri.mat');

file17=strcat(root, '01_02_2015_ace.mat');
file18=strcat(root, '01_02_2015_conace.mat');

file19=strcat(root, '01_02_2015_h20.mat') ;
file20=strcat(root, '01_02_2015_conh20.mat' );

file21=strcat(root, '01_02_2015_gly.mat' );
file22=strcat(root, '01_02_2015_congly.mat' );

file23=strcat(root, '12_30_2014_dmso.mat');
file24=strcat(root, '12_30_2014_condmso.mat');

file25=strcat(root, '02_15_2016_lindane.mat');
file26=strcat(root, '02_15_2016_conlin.mat');

all.filepath = {file1, file2, file3, file4, file5,...
    file6,file7,file8,file9,file10,file11,file12,file13,file14,...
    file15,file16, file17, file18, file19, file20, file21,file22,file23,...
    file24, file25, file26}  ;
        


% 
for curFile=[2 4 6 8 10 12 14 16 18 20 22 24 26]
    %load control file
    s=load( all.filepath{curFile}  );
    name=fieldnames(s);
    
    chem = s.(name{1}) ;
    oChunk=[]; fileTitle=[];
    
    trtemp=load( all.filepath{curFile-1});
    tname=fieldnames(trtemp) ;
    trt=trtemp.(tname{1});
    
    % cur chemical name
    [a,b,c]=fileparts(all.filepath{(curFile-1) });
    d=strsplit(b,'_'); chemName= d( length(d) );
    for curPlate=1:length(chem.file)
 
        curNSeconds=chem.file(curPlate).nseconds ;
        curNSpikes=[]; curAE=[];
        for ent = 1:60
            curNSpikes =[curNSpikes  chem.file(curPlate).channel(ent).nspikes ] ;
            % indicator of ae or not
            curAE=[curAE chem.file(curPlate).channel(ent).ae ] ;   

        end; %end of loop through all entities

        mfr_all=sum(curNSpikes)/(60*curNSeconds) ;
        curPerAE= sum( curAE )/60 ;
        mfr_AE = dot(curAE,curNSpikes)/dot(curAE, curNSeconds*ones(1,60) );
        
        %now do treatment
        curNSecondsTrt=trt.file(curPlate).nseconds ;
        curNSpikesTrt=[]; curAETrt=[];
        for ent = 1:60
            curNSpikesTrt =[curNSpikesTrt  trt.file(curPlate).channel(ent).nspikes ] ;
            % indicator of ae or not
            curAETrt=[curAETrt trt.file(curPlate).channel(ent).ae ] ;
        end; %end of loop through all entities

        mfr_all=sum(curNSpikesTrt)/(60*curNSecondsTrt) ;
        % BIC files mfr computed using AE from trt others mfr calculations
        % use pre-trt files
        if curFile==2
            %divide number of total spikes on AE in trt by total # seconds
           mfr_AE_trt = dot(curAETrt,curNSpikesTrt)/dot(curAETrt, curNSecondsTrt*ones(1,60) );
        else 
            mfr_AE_trt = dot(curAE,curNSpikesTrt)/dot(curAE, curNSecondsTrt*ones(1,60) );
        end;
        
        ratio_mfr=mfr_AE_trt/mfr_AE ;
        ratio_nAE=sum(curAETrt)/sum(curAE) ;

        if strcmp( chemName,'lindane' )
            [p,n1,x] = fileparts( chem.filepath{curPlate} ) ;
        else
            temp=strsplit(chem.filepath{curPlate} , '\');
            n1=temp{end};
        end;
        title2= strcat( n1,' ' );
        
        fileTitle{curPlate}=title2;
        
        temp=vertcat(ratio_mfr , ratio_nAE );
        oChunk =horzcat(oChunk, temp);
                    
                
    end;%loop through plates in file
    celldata= cellstr( ['mfr_all' ;'mfr_AE ' ;'% AE   '] );
    
    %%%%%output
    % fileID = fopen(xlsfilestr,'w');
    temp=chemName;
    fprintf(fileID, '%s\n', temp{1}  );
    % file header
    temp= horzcat('__', fileTitle(1:(end-1)), strcat( fileTitle(end), '\n') ) ;
    formatSpec = strcat( repmat('%s ,', 1,(size(fileTitle,2)) ) , ' %s\n');
    fprintf(fileID, formatSpec, temp{1:end} );
    % mfr_ration
    formatSpec = strcat(repmat('%d ,', 1,(size(oChunk,2))) , ' %d\n');
    fprintf(fileID, '%s', 'ratio_mfr_onAEinTrt_over_mfr_onAEinBL');
    fprintf(fileID, formatSpec, horzcat(0 ,oChunk(1,:)) ) ;
    % nAE_ration 
    formatSpec = strcat(repmat('%d ,', 1,(size(oChunk,2))) , ' %d\n');
    fprintf(fileID, '%s', 'ratio_nAEinTrt_over_nAEinBL');
    fprintf(fileID, formatSpec, horzcat(0 ,oChunk(2,:)) ) ;
    
    fprintf(fileID, '\n');

    
end;  % end of loop through files
    
fclose(fileID) ;

% generate_derived_LFP_data.m
% diana hall
% november 5, 2013
%%%Use NeuroShare
%you need the FIND toolbox and the NeuroShare



xlsfilestr='F:\EXTRAP\matlab_Cina_Herr_Diana\high_Freq\high_freq_10-29-2014.xls';
title='High Frequency Information for All Files';
xlswrite(xlsfilestr,title, 1,'A1' ) ;


file1='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_bic.mat';
file2='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\02_11_2014_per.mat';
file3='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_29_2014_car.mat';
file4='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_28_2014_con.mat';

file5='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\03_17_2014_da.mat';
file6='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\03_17_2014_conDA.mat';

file7='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_tri.mat';
file8='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_contri.mat';

file9='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\07_30_2014_ace.mat';
file10='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\07_30_2014_conace.mat';

file11='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\07_7_2014_h20.mat';
file12='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\08_05_2014_conh20.mat';

file13='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_21_2014_per50.mat';
file14='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_21_2014_conPER50.mat';

file15='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_28_2014_gly.mat';
file16='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_28_2014_congly.mat';


all.filepath = {file1, file2, file3, file4, file5,...
    file6,file7,file8,file9,file10,file11,file12,file13,file14,...
    file15,file16}  ;
        
startCell=5;
for curFile=1:length(all.filepath)
    %load file
    s=load( all.filepath{curFile}  );
    name=fieldnames(s);
    chem = s.(name{1}) ;
    oChunk=[]; fileTitle=[];
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
        
        [p,n1,x] = fileparts( all.filepath{curFile} ) ;
        title= strcat( n1,' Plate ',num2str(curPlate) );
        
        fileTitle{curPlate}=title ;
        
        temp=vertcat(mfr_all, mfr_AE, curPerAE);
        oChunk=horzcat(oChunk, temp);
                    
                
    end;%loop through plates in file
    celldata= cellstr( ['mfr_all' ;'mfr_AE ' ;'% AE   '] );
    
    %%%%%output
    xlswrite(xlsfilestr, fileTitle ,1,['B' num2str(startCell-1)] ) ;
    xlswrite(xlsfilestr, fileTitle ,1,['B' num2str(startCell)] ) ;
    xlswrite(xlsfilestr, oChunk,1,['C' num2str(startCell)] ) ;
       
    startCell=startCell+5;
end;  % end of loop through files
    





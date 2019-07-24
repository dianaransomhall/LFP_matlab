%% RUN T-TEST on %Change in LFP at 5 Different Bands
% t-test_anova.m
% Diana Hall
% 01/27/2014
% code to do anova

% add directory
mac=1;
if ~mac
    addpath(genpath('F:\EXTRAP\matlab_Cina_Herr_Diana\Scripts' ) )
else
    addpath(genpath('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/Scripts' ) )
end;



%% ++++++++++++++++++++++ACE++++++++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_ace.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conace.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_ace.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_conace.mat');
end;

% make analysis data set
filepath.ace=ace.filepath;
filepath.conace=con.filepath;
chem=ace;  
cont=con; cont.file = con.file;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

ACEanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
ACE = do_t_anova(ACEanova, want ) ;
clear('ace');clear('con');



%% ++++++++++++++++++++++++++  BIC ++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
  load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_bic.mat');
  load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_biccon.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_27_2014_bic.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/10_30_2014_biccon.mat');
end;


%  t-test or anova
chem=bic; 
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

filepath.bic=bic.filepath;
filepath.conbic=con.filepath ;

BICanova = make_t_anova_data_AEOnly(chem, cont, nFiles, want) ;

BIC = do_t_anova(BICanova, want ) ;
clear('bic'); clear('con');



%% +++++++++++++++++++++++CAR 5 uM +++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car5.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car5con.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_car5.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_car5con.mat');
end;

% make analysis data set
chem=car;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;
% get filepaths for writing
filepath.car5 = car.filepath;
filepath.concar5=con.filepath ;

CAR5anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
CAR5 = do_t_anova(CAR5anova, want ) ;
%results: fail to reject at all bands 
clear('con'); clear('car');


%% ++++++++++++++++++++++++++++++++++++++++++++++++++CAR 15 uM +++++++++++++++
% load data
if ~mac
    load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car15.mat');
    load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car15con.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_car15.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_car15con.mat');
end;

% make analysis data set
chem=car;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;
% get filepaths for writing
filepath.car15 = car.filepath;
filepath.concar15=con.filepath ;
echo off;
CAR15anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
echo off;
CAR15 = do_t_anova(CAR15anova, want ) ;
%results: fail to reject at all bands 
for i=1:5
disp( strcat( ' band', num2str(i) ) )
disp(  CAR15{i}.t )
end;
clear('con'); clear('car');


%% +++++++++++++++++++++++++++++++++++++++++++++++++CAR 30 +++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30con.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_car30.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_car30con.mat');
end;

% make analysis data set
chem=car;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;
% get filepaths for writing
filepath.car30 = car.filepath;
filepath.concar30=con.filepath ;

CAR30anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
CAR30 = do_t_anova(CAR30anova, want ) ;
%results: fail to reject at all bands 
clear('car'); clear('con');




%% +++++++++++++DOMOIC ACID+++++++++++++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_da.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_conDA.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/10_30_2014_da.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/10_30_2014_conDA.mat');
end;

% make analysis data set
filepath.da=da.filepath;
filepath.conda=con.filepath;
chem=da;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

DAanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
DA = do_t_anova(DAanova, want ) ;
% DA{1}
clear('da');clear('con');


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++DMSO+++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_dmso.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_condmso.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/12_30_2014_dmso.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/12_30_2014_condmso.mat');
end;

% make analysis data set
filepath.dmso=dmso.filepath;
filepath.condmso=con.filepath;
chem=dmso;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

DMSOanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
DMSO = do_t_anova(DMSOanova, want ) ;
clear('dmso');clear('con');



%% +++++++++++++GLY+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_gly.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_congly.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_gly.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_congly.mat');
end;

% make analysis data set
filepath.gly=gly.filepath;
filepath.congly=con.filepath;
chem=gly;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

GLYanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
GLY = do_t_anova(GLYanova, want ) ;
clear('gly');clear('con');





%% +++++++++++++++++++++++++++H20+++++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_h20.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conh20.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_h20.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_conh20.mat');
end;

% make analysis data set
filepath.h20=h20.filepath;
filepath.conh20=con.filepath;
chem=h20;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

H20anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
H20 = do_t_anova(H20anova, want ) ;
clear('h20');clear('con');






%% +++++++++++++LINDANE+++++++++++++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
    load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\02_15_2016_lindane.mat');
    load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\02_15_2016_conlin.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/02_15_2016_lindane.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/02_15_2016_conlin.mat');
end;

% make analysis data set
filepath.lin=lin.filepath;
filepath.conlin=con.filepath;
chem=lin;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

LINanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
LIN = do_t_anova(LINanova, want ) ; % LIN{1}
clear('lin');clear('con');




%% ++++++++++++++++++++++++++++++PER25++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\02_11_2014_per25.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_per25con.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/02_11_2014_per25.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/10_30_2014_per25con.mat');
end;

% make analysis data set
chem=per;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

filepath.per25=per.filepath;
filepath.conper25=con.filepath ;

PER25anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
PER25 = do_t_anova(PER25anova, want ) ;
%results: fail to reject at bands = 1,2,3,4; sig. diff at band 5
clear('per');clear('con');



%% +++++++++++++++++++++PERM 50uM+++++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_per50.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conPER50.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_per50.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/01_02_2015_conPER50.mat');
end;

% make analysis data set
filepath.per50=per.filepath;
filepath.per50con=con.filepath;
chem=per;   %which files in con are PER baseline
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

PER50anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
PER50 = do_t_anova(PER50anova, want ) ;
%results: fail to reject at bands = 1,2,3,4; sig. diff at band 5
clear('per');clear('con');






%% +++++++++++++++++++++++++++++++++TRI+++++++++++++++++++++++++++++++++++++
% load data
if ~mac
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_tri.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_contri.mat');
else
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/8_9_2014_tri.mat');
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/8_9_2014_contri.mat');
end;

% make analysis data set
filepath.tri=tri.filepath;
filepath.conh20=con.filepath;
chem=tri;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0;want.diff=0;want.pchange=1;

TRIanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
TRI = do_t_anova(TRIanova, want ) ;
clear('tri');clear('con');









%% Save results

 if ~mac
    filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_03_2015_pChangeanova.mat' ;
 else 
     filename='/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/03_01_2016_pChangeanova.mat' ;
 end;
 save(filename, 'BIC','BICanova','ACE', 'ACEanova', ...
     'CAR5','CAR5anova','CAR15','CAR15anova','CAR30','CAR30anova',...
     'DA','DAanova','DMSO','DMSOanova','LIN', 'LINanova',...
     'GLY', 'GLYanova','H20','H20anova',...
     'PER25', 'PER25anova', 'PER50', 'PER50anova','TRI','TRIanova', ...
     'filepath') ;



 
 
 
 

%% T-test vs. baseline: Write Results
% ++++++++++++++++++++++++++++++++++++++
% write results with t-test H_0: t=1, H_a: t~=1
mac=1;
if ~mac
    addpath(genpath('F:\EXTRAP\matlab_Cina_Herr_Diana\Scripts' ) )
else
    addpath(genpath('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/Scripts' ) )
end;

wantData=1;
if wantData
    load('/Volumes/BACKUP/EXTRAP/matlab_Cina_Herr_Diana/finshed_data/03_01_2016_pChangeanova.mat');
end;



    xlsfilestr='/Volumes/BACKUP/EXTRAP/documents/Results/lfp_results_t-testdiff1_03_01_2016.csv';
    title1={'Results: % Baseline chem vs. % baseline control, mean all Electrodes, t-test: eq. var'};
    title=title{1};
    %da,gly,ace,tri are H20 controls
    %bic, car, per30, per50 are dmso controls
    chemName= {'ace','bic', 'da','dmso','car5','car15', 'car30','gly', 'h20',...
        'lin','per25','per50','tri' };

    fileID = fopen(xlsfilestr,'w');
    fprintf( fileID,'%s\n',  title );
    startCell=3;
    for k=1:size(chemName,2)
    % write files

        if strcmp( chemName{1,k},'per25' )
            anovaData = PER25anova ; test = PER25;
            for i=1:length(filepath.per25)
                if mac
                    temp=strsplit(filepath.per25{1}, '\');
                    n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.per25{i} );
                end;
                n{i}=n1 ;
            end;
            % 2 sample t-test to control
            for curB=1:5
                [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(PER25anova{curB}.pchange, DMSOanova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.per25);

        elseif strcmp( chemName{1,k},'per50')
            anovaData = PER50anova ; test = PER50
            for i=1:length(filepath.per50)
                if mac
                    temp=strsplit( filepath.per50{i}, '\');
                    n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.per50{i} );
                end;
                
                n{i}=n1 ;
            end;
            for curB=1:5
                [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(PER50anova{curB}.pchange, DMSOanova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.per50);


        elseif strcmp( chemName{1,k},'car5')
            anovaData = CAR5anova ; test = CAR5
            for i=1:length(filepath.car5)
                 if mac
                    temp=strsplit(  filepath.car5{i}, '\');
                    n1=temp{end};
                else
                     [p,n1,x] = fileparts( filepath.car5{i} );
                end;
               
                n{i}=n1 ;
            end;
            for curB=1:5
                [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(anovaData{curB}.pchange, DMSOanova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.car5);

        elseif strcmp( chemName{1,k},'car15')
            anovaData = CAR15anova ; test = CAR15 ;
            for i=1:length(filepath.car15 )
                if mac
                    temp=strsplit(  filepath.car15{i}, '\');
                    n1=temp{end};
                else
                     [p,n1,x] = fileparts( filepath.car15{i} );
                end;

                n{i}=n1 ;
            end;
            for curB=1:5
                [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(anovaData{curB}.pchange, DMSOanova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.car15);

        elseif strcmp( chemName{1,k},'car30')
            anovaData = CAR30anova ; test = CAR30 ;
            for i=1:length(filepath.car30)
                if mac
                    temp=strsplit( filepath.car30{i} , '\');
                    n1=temp{end};
                else
                     [p,n1,x] = fileparts( filepath.car30{i} );
                end;
               
                n{i}=n1 ;
            end;
            for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(anovaData{curB}.pchange, DMSOanova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.car30);

        elseif strcmp( chemName{1,k},'bic' )
            anovaData = BICanova ; test = BIC;
            for i=1:length(filepath.bic)
                 if mac
                    temp=strsplit( filepath.bic{i}  , '\');
                    n1=temp{end};
                else
                     [p,n1,x] = fileparts( filepath.bic{i} );
                end;
               
                n{i}=n1 ;
            end;
            for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(BICanova{curB}.pchange, DMSOanova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.bic);

        elseif strcmp( chemName{1,k},'da')
            anovaData = DAanova ; test = DA
            for i=1:length(filepath.da)
                 if mac
                    temp=strsplit(filepath.da{i} , '\');
                    n1=temp{end};
                else
                     [p,n1,x] = fileparts( filepath.da{i} );
                end;
               
                n{i}=n1 ;
            end;
            for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(DAanova{curB}.pchange, H20anova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.da);


        elseif strcmp( chemName{1,k},'ace')
            anovaData = ACEanova ; test = ACE
            for i=1:length(filepath.ace )
                if mac
                    temp=strsplit( filepath.ace{i} , '\');
                    n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.ace{i} );
                end;
                
                n{i}=n1 ;
            end;
            for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(ACEanova{curB}.pchange, H20anova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.ace);

        elseif strcmp( chemName{1,k},'h20')
            anovaData = H20anova ; test = H20
            for i=1:length(filepath.h20)
                if mac
                    temp=strsplit(filepath.h20{i} , '\');
                    n1=temp{end};
                else
                   [p,n1,x] = fileparts( filepath.h20{i} );
                end;
                
                n{i}=n1 ;
            end;
            nlines=length(filepath.h20);

        elseif strcmp( chemName{1,k},'dmso')
            anovaData = DMSOanova ; test = DMSO
            for i=1:length(filepath.dmso)

                if mac
                    temp=strsplit(filepath.dmso{i}, '\');
                    n1=temp{end};
                else
                   [p,n1,x] = fileparts( filepath.dmso{i} );
                end;
 
                n{i}=n1 ;
            end;
            nlines=length(filepath.dmso);


        elseif strcmp( chemName{1,k},'tri')
            anovaData = TRIanova ; test = TRI;
            for i=1:length(filepath.tri)
                 if mac
                    temp=strsplit( filepath.tri{i} , '\');
                    n1=temp{end};
                else
                  [p,n1,x] = fileparts( filepath.tri{i} );
                end;
                
                n{i}=n1 ;
            end;
            for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(TRIanova{curB}.pchange, H20anova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.tri );

        elseif strcmp( chemName{1,k},'gly')
            anovaData = GLYanova ; test = GLY ;
            for i=1:length(filepath.gly)
                if mac
                    temp=strsplit(  filepath.gly{i} , '\');
                   n1=temp{end};
                else
                  [p,n1,x] = fileparts( filepath.gly{i} );
                end;
                
                n{i}=n1 ;
            end;
            for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(GLYanova{curB}.pchange, H20anova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.gly );
        elseif strcmp( chemName{1,k},'lin')
            anovaData = LINanova ; test = LIN ;
            for i=1:length(filepath.lin)

                  [p,n1,x] = fileparts( filepath.lin{i} );

                
                n{i}=n1 ;
            end;
            for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
                anova{curB}.ci, anova{curB}.stats ]=...
                ttest2(LINanova{curB}.pchange, DMSOanova{curB}.pchange ) %significant
            end;
            nlines=length(filepath.lin );

        end;

        ex{1,1}=' '; ex{2,1}='1-4Hz'; ex{3,1}='4-8Hz'; ex{4,1}='8-14Hz'; 
        ex{5,1}='14-30Hz'; ex{6,1}='30-50Hz';
        nFiles = length(anovaData{1}.pchange ) ;
        for curfile=1:nFiles
            ex{1,1+curfile}=n{curfile} ;

        end;

        for j=1:5
            for i=1:length(anovaData{1}.pchange)
                ex{1+j,i+1} = anovaData{j}.pchange(i);
            end;
        end;
        ex{1,2+nFiles}='mean';
        ex{1,3+nFiles}='se';
        ex{1,4+nFiles}='sd'; 
        ex{1,5+nFiles}='t-test diff from 1'; 
        ex{1,6+nFiles}='p-value'; 
        ex{1,7+nFiles}='lci'; 
        ex{1,8+nFiles}='uci';
        ex{1,9+nFiles}='n';
        for j=1:5
            ex{1+j,2+nFiles}= mean(anovaData{j}.pchange);
            ex{1+j,3+nFiles}= test{j}.stats.sd/length(anovaData{j}.pchange);
            ex{1+j,4+nFiles}= test{j}.stats.sd; 
            ex{1+j,5+nFiles}= test{j}.stats.tstat; 
            ex{1+j,6+nFiles}= test{j}.tp_value; 
            ex{1+j,7+nFiles}= test{j}.ci(1); 
            ex{1+j,8+nFiles}= test{j}.ci(2); 
            ex{1+j,9+nFiles}= test{j}.stats.df+1;
        end;
        
        [nrows,ncols] = size(ex);
        % fileID = fopen(xlsfilestr,'w');
        fprintf(fileID, strcat( chemName{1,k}, '\n')  );
        for row = 1:nrows
            if row==1
                formatSpec = strcat( repmat('%s ,', 1,(size(ex,2)-1)) , ' %s\n');
                fprintf(fileID,formatSpec, ex{row,:});
                
            else
                formatSpec = strcat('%s ,',repmat('%d ,',1,(size(ex,2)-2)), ' %d\n'  ); 
                fprintf(fileID,formatSpec,ex{row,:});
            end;
        end;
    
        clear('ex');

        %set startCell to empty lines
        startCell=startCell+10 ;
    end; % end of for loop

    % close the file
    fclose(fileID);






    
    
%% +++++++ OUTPUT T-Test vs Control ++++++++++++++++++++++++++++++++++++++
% write results with t-test H_0: t=1, H_a: t~=1
if ~mac
    xlsfilestr='F:\EXTRAP\documents\Results\lfp_results_t-testwControl_01_06_2015.csv';
else
    xlsfilestr='/Volumes/BACKUP/EXTRAP/documents/Results/lfp_results_t-testvControl_03_01_2016.csv';
end;
title1={'Results: % Baseline chem vs. % baseline control, mean all Electrodes, t-test uneq. var'};
title=title1{1};
fileID = fopen(xlsfilestr,'w');
fprintf( fileID,'%s\n',  title );

%da,gly,ace,tri are H20 controls
%bic, car, per30, per50 are dmso controls
chemName= {'ace','bic', 'da','dmso','car5','car15', 'car30','gly', 'h20',...
    'lin','per25','per50','tri' };

     
startCell=3;
for k=1:size(chemName,2)
% write files

    if strcmp( chemName{1,k},'per25' )
        anovaData = PER25anova ; test = PER25;
        for i=1:length(filepath.per25)
            if mac
                 temp=strsplit( filepath.per25{i} , '\');
                   n1=temp{end};
                else
                   [p,n1,x] = fileparts( filepath.per25{i} );
                end;

            n{i}=n1 ;
        end;
        % 2 sample t-test to control
        for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(PER25anova{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.per25);
        
    elseif strcmp( chemName{1,k},'per50')
        anovaData = PER50anova ; test = PER50
        for i=1:length(filepath.per50)
             if mac
                 temp=strsplit( filepath.per50{i} , '\');
                   n1=temp{end};
                else
                   [p,n1,x] = fileparts( filepath.per50{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(PER50anova{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.per50);
        
        
    elseif strcmp( chemName{1,k},'car5')
        anovaData = CAR5anova ; test = CAR5
        for i=1:length(filepath.car5)
            if mac
                 temp=strsplit(  filepath.car5{i}  , '\');
                   n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.car5{i} );
                end;
           
            n{i}=n1 ;
        end;
        for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(anovaData{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.car5);
      
    elseif strcmp( chemName{1,k},'car15')
        anovaData = CAR15anova ; test = CAR15 ;
        for i=1:length(filepath.car15 )
            if mac
                 temp=strsplit(  filepath.car15{i}  , '\');
                   n1=temp{end};
                else
                     [p,n1,x] = fileparts( filepath.car15{i} );
                end;
           
            n{i}=n1 ;
        end;
        for curB=1:5
            [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(anovaData{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.car15);
     
    elseif strcmp( chemName{1,k},'car30')
        anovaData = CAR30anova ; test = CAR30 ;
        for i=1:length(filepath.car30)
            if mac
                 temp=strsplit(  filepath.car30{i}  , '\');
                   n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.car30{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(anovaData{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.car30);
        
    elseif strcmp( chemName{1,k},'bic' )
        anovaData = BICanova ; test = BIC;
        for i=1:length(filepath.bic)
            if mac
                 temp=strsplit(filepath.bic{i}  , '\');
                   n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.bic{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(BICanova{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.bic);
        
    elseif strcmp( chemName{1,k},'da')
        anovaData = DAanova ; test = DA
        for i=1:length(filepath.da)
              if mac
                 temp=strsplit(filepath.da{i}  , '\');
                   n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.da{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(DAanova{curB}.pchange, H20anova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=H20anova;
        nlines=length(filepath.da);
        
    
    elseif strcmp( chemName{1,k},'ace')
        anovaData = ACEanova ; test = ACE
        for i=1:length(filepath.ace )
            if mac
                
                    temp=strsplit( filepath.ace{i}   , '\');
                    n1=temp{end};
                else
                   [p,n1,x] = fileparts( filepath.ace{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(ACEanova{curB}.pchange, H20anova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=H20anova;
        nlines=length(filepath.ace);
            
    elseif strcmp( chemName{1,k},'h20')
        anovaData = H20anova ; test = H20
        for i=1:length(filepath.h20)
            if mac
                    temp=strsplit( filepath.h20{i}   , '\');
                    n1=temp{end};
                else
                   [p,n1,x] = fileparts( filepath.h20{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(H20anova{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.h20);
        
    elseif strcmp( chemName{1,k},'dmso')
        anovaData = DMSOanova ; test = DMSO
        for i=1:length(filepath.dmso)
            if mac
                    temp=strsplit( filepath.dmso{i}   , '\');
                    n1=temp{end};
                else
                   [p,n1,x] = fileparts( filepath.dmso{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(DMSOanova{curB}.pchange, H20anova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=H20anova;
        nlines=length(filepath.dmso);
        
        
    elseif strcmp( chemName{1,k},'tri')
        anovaData = TRIanova ; test = TRI;
        for i=1:length(filepath.tri)
            if mac
                    temp=strsplit( filepath.tri{i}  , '\');
                    n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.tri{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(TRIanova{curB}.pchange, H20anova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=H20anova;
        nlines=length(filepath.tri );
    
    elseif strcmp( chemName{1,k},'gly')
        anovaData = GLYanova ; test = GLY ;
        for i=1:length(filepath.gly)
             if mac
                    temp=strsplit(filepath.gly{i} , '\');
                    n1=temp{end};
                else
                    [p,n1,x] = fileparts( filepath.gly{i} );
                end;
            
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(GLYanova{curB}.pchange, H20anova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=H20anova;
        nlines=length(filepath.gly );
        
   elseif strcmp( chemName{1,k},'lin')
        anovaData = LINanova ; test = LIN ;
        for i=1:length(filepath.lin)
            [p,n1,x] = fileparts( filepath.lin{i} );
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(LINanova{curB}.pchange, DMSOanova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=DMSOanova;
        nlines=length(filepath.lin );
     
    end;
    
    
    ex{1,1}=' '; ex{2,1}='1-4Hz'; ex{3,1}='4-8Hz'; ex{4,1}='8-14Hz'; 
    ex{5,1}='14-30Hz'; ex{6,1}='30-50Hz';
    nFiles = length(anovaData{1}.pchange ) ;
    for curfile=1:nFiles
        ex{1,1+curfile}=n{curfile} ;
    end;
    
    for j=1:5
        for i=1:length(anovaData{1}.pchange)
            ex{1+j,i+1} = anovaData{j}.pchange(i);
        end;
    end;
    ex{1,2+nFiles}='mean trt';
    ex{1,3+nFiles}='sd trt';
    ex{1,4+nFiles}='mean cont';
    ex{1,5+nFiles}='sd con'; 
    ex{1,6+nFiles}='t-test diff from contol'; 
    ex{1,7+nFiles}='p-value'; 
    ex{1,8+nFiles}='lci diff'; 
    ex{1,9+nFiles}='uci diff';
    ex{1,10+nFiles}='n trt';
    
    for j=1:5
        ex{1+j,2+nFiles}= mean(anovaData{j}.pchange);
        ex{1+j,3+nFiles}= anova{j}.stats.sd(1) ; 
        ex{1+j,4+nFiles}= mean(CONanova{j}.pchange);
        ex{1+j,5+nFiles}= anova{j}.stats.sd(2) ;
        ex{1+j,6+nFiles}= anova{j}.stats.tstat; 
        ex{1+j,7+nFiles}= anova{j}.tp_value; 
        ex{1+j,8+nFiles}= anova{j}.ci(1); 
        ex{1+j,9+nFiles}= anova{j}.ci(2); 
        ex{1+j,10+nFiles}= test{j}.stats.df+1;
    end;
    [nrows,ncols] = size(ex);
        % fileID = fopen(xlsfilestr,'w');
        fprintf(fileID, strcat( chemName{1,k}, '\n')  );
        for row = 1:nrows
            if row==1
                formatSpec = strcat( repmat('%s ,', 1,(size(ex,2)-1)) , ' %s\n');
                fprintf(fileID,formatSpec, ex{row,:});
                
            else
                formatSpec = strcat('%s ,',repmat('%d ,',1,(size(ex,2)-2)), ' %d\n'  ); 
                fprintf(fileID,formatSpec,ex{row,:});
            end;
        end;
    
        
    clear('ex');
    
    
    %set startCell to empty lines
    startCell=startCell+10 ;
end; 

% close the file
fclose(fileID);















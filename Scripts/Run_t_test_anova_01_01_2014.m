% t-test_anova.m
% Diana Hall
% 01/27/2014
% code to do anova

% add directory
mac=1;
if ~mac
    addpath(genpath('F:\EXTRAP\matlab_Cina_Herr_Diana\Scripts' ) )
end;


% load data
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_bic.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_biccon.mat');


%+++++++++ t-test or anova
%+++++++  BIC +++++++++++++++++++++++
chem=bic; 
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

filepath.bic=bic.filepath;
filepath.conbic=con.filepath ;

BICanova = make_t_anova_data_AEOnly(chem, cont, nFiles, want) ;

BIC = do_t_anova(BICanova, want ) ;
%clear('bic'); clear('con');



%%%+++++++++++++CAR5 uM +++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car5.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car5con.mat');
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


%%%+++++++++++++CAR15 uM +++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car15.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car15con.mat');
% make analysis data set
chem=car;  
cont=con; cont.file = con.file ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;
% get filepaths for writing
filepath.car15 = car.filepath;
filepath.concar15=con.filepath ;

CAR15anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
CAR15 = do_t_anova(CAR15anova, want ) ;
%results: fail to reject at all bands 
clear('con'); clear('car');





%%%+++++++++++++CAR 30 +++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_car30con.mat');
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



%%%+++++++++++++PER25 +++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\02_11_2014_per25.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_per25con.mat');
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



%%%+++++++++++++PERM 50 uM+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_per50.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conPER50.mat');
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





%%%+++++++++++++DA+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_da.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\10_30_2014_conDA.mat');
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



%%%+++++++++++++Ace+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_ace.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conace.mat');
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





%%%+++++++++++++H20+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_h20.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_conh20.mat');
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






%%%+++++++++++++TRI+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_tri.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_9_2014_contri.mat');
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





%%%+++++++++++++dmso+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_dmso.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_30_2014_condmso.mat');
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


%%%+++++++++++++gly+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_gly.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_02_2015_congly.mat');
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



 filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_03_2015_pChangeanova.mat' ;
 save(filename, 'PER25', 'PER25anova', 'PER50', 'PER50anova', ...
     'DA','DAanova','DMSO','DMSOanova', 'GLY', 'GLYanova',...
    'TRI', 'TRIanova', 'BIC','BICanova', ...
    'CAR5','CAR5anova','CAR15','CAR15anova','CAR30','CAR30anova',...
     'H20','H20anova', 'ACE','ACEanova', ...
     'filepath') ;





% ++++++++++++++++++++++++++++++++++++++
% write results with t-test H_0: t=1, H_a: t~=1
xlsfilestr='F:\EXTRAP\documents\Results\lfp_results_t-testdiff1_01_06_2015.xls';
title='Results of Mean % Control in Power averaged across all Electrodes';
xlswrite(xlsfilestr,title, 1,'A1' ) ;
%da,gly,ace,tri are H20 controls
%bic, car, per30, per50 are dmso controls
chemName= {'h20','dmso','bic','car5','car15', 'car30','per25','per50',...
    'da','ace', 'tri', 'gly' };
startCell=3;
for k=1:size(chemName,2)
% write files

    if strcmp( chemName{1,k},'per25' )
        anovaData = PER25anova ; test = PER25;
        for i=1:length(filepath.per25)
            [p,n1,x] = fileparts( filepath.per25{i} );
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
            [p,n1,x] = fileparts( filepath.per50{i} );
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
            [p,n1,x] = fileparts( filepath.car5{i} );
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
            [p,n1,x] = fileparts( filepath.car15{i} );
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
            [p,n1,x] = fileparts( filepath.car30{i} );
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
            [p,n1,x] = fileparts( filepath.bic{i} );
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
            [p,n1,x] = fileparts( filepath.da{i} );
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
            [p,n1,x] = fileparts( filepath.ace{i} );
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
            [p,n1,x] = fileparts( filepath.h20{i} );
            n{i}=n1 ;
        end;
        nlines=length(filepath.h20);
        
    elseif strcmp( chemName{1,k},'dmso')
        anovaData = DMSOanova ; test = DMSO
        for i=1:length(filepath.dmso)
            [p,n1,x] = fileparts( filepath.dmso{i} );
            n{i}=n1 ;
        end;
        nlines=length(filepath.dmso);
        
        
    elseif strcmp( chemName{1,k},'tri')
        anovaData = TRIanova ; test = TRI;
        for i=1:length(filepath.tri)
            [p,n1,x] = fileparts( filepath.tri{i} );
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
            [p,n1,x] = fileparts( filepath.gly{i} );
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(GLYanova{curB}.pchange, H20anova{curB}.pchange ) %significant
        end;
        nlines=length(filepath.gly );
     
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
    xlswrite(xlsfilestr, { chemName{1,k} },1,['B' num2str(startCell-1)] ) ;
    xlswrite(xlsfilestr, ex,1,['C' num2str(startCell)] ) ;
    clear('ex');
    
    %set startCell to empty lines
    startCell=startCell+10 ;
end; 








% ++++++++++++++++++++++++++++++++++++++
% write results with t-test H_0: t=1, H_a: t~=1
xlsfilestr='F:\EXTRAP\documents\Results\lfp_results_t-testwControl_01_06_2015.xls';
title='Results of Mean % Control in Power averaged across all Electrodes';
xlswrite(xlsfilestr,title, 1,'A1' ) ;
%da,gly,ace,tri are H20 controls
%bic, car, per30, per50 are dmso controls
chemName= {'h20','dmso','bic','car5','car15', 'car30','per25','per50',...
    'da','ace', 'tri', 'gly' };
startCell=3;
for k=1:size(chemName,2)
% write files

    if strcmp( chemName{1,k},'per25' )
        anovaData = PER25anova ; test = PER25;
        for i=1:length(filepath.per25)
            [p,n1,x] = fileparts( filepath.per25{i} );
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
            [p,n1,x] = fileparts( filepath.per50{i} );
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
            [p,n1,x] = fileparts( filepath.car5{i} );
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
            [p,n1,x] = fileparts( filepath.car15{i} );
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
            [p,n1,x] = fileparts( filepath.car30{i} );
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
            [p,n1,x] = fileparts( filepath.bic{i} );
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
            [p,n1,x] = fileparts( filepath.da{i} );
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
            [p,n1,x] = fileparts( filepath.ace{i} );
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
            [p,n1,x] = fileparts( filepath.h20{i} );
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
            [p,n1,x] = fileparts( filepath.dmso{i} );
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
            [p,n1,x] = fileparts( filepath.tri{i} );
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
            [p,n1,x] = fileparts( filepath.gly{i} );
            n{i}=n1 ;
        end;
        for curB=1:5
        [anova{curB}.t, anova{curB}.tp_value, ...
            anova{curB}.ci, anova{curB}.stats ]=...
            ttest2(GLYanova{curB}.pchange, H20anova{curB}.pchange,'vartype','unequal' ) %significant
        end;
        CONanova=H20anova;
        nlines=length(filepath.gly );
     
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
    xlswrite(xlsfilestr, { chemName{1,k} },1,['B' num2str(startCell-1)] ) ;
    xlswrite(xlsfilestr, ex,1,['C' num2str(startCell)] ) ;
    clear('ex');
    
    %set startCell to empty lines
    startCell=startCell+10 ;
end; 
















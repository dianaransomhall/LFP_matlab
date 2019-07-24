% t-test_anova.m
% Diana Hall
% 01/27/2014
% code to do anova

% add directory
addpath(genpath('F:\EXTRAP\matlab_Cina_Herr_Diana\Scripts' ) )


% load data
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_27_2014_bic.mat');
 load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_28_2014_con.mat');
% load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_all.mat')


%+++++++++ t-test or anova
%+++++++  BIC +++++++++++++++++++++++
chem=bic; conFIndex=[1:5]; %index of con files
cont=con; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

filepath.bic=bic.filepath;
filepath.conbic=con.filepath(conFIndex) ;

BICanova = make_t_anova_data_AEOnly(chem, cont, nFiles, want) ;

BIC = do_t_anova(BICanova, want ) ;
%clear('bic'); clear('con');




%%%+++++++++++++CAR+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_29_2014_car.mat');
% make analysis data set
chem=car;  conFIndex=[11:15]; %which files in con are CAR baseline
cont=con; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;
% get filepaths for writing
filepath.car=car.filepath;
filepath.concar=con.filepath(conFIndex);

CARanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
CAR = do_t_anova(CARanova, want ) ;
%results: fail to reject at all bands 




%%%+++++++++++++PERM +++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\02_11_2014_per.mat');
% make analysis data set
chem=per;  conFIndex=[6:10]; %which files in con are PER baseline
cont=con; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

filepath.per=per.filepath;
filepath.conper=con.filepath(conFIndex);

PERanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
PER = do_t_anova(PERanova, want ) ;
%results: fail to reject at bands = 1,2,3,4; sig. diff at band 5
clear('per');clear('con');



%%%+++++++++++++PERM 50 uM+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\08_26_2014_per50.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\08_26_2014_conPER50.mat');
% make analysis data set
filepath.per50=per.filepath;
filepath.per50con=con.filepath;
chem=per;  conFIndex=[1:3]; %which files in con are PER baseline
cont=con; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

PER50anova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
PER50 = do_t_anova(PER50anova, want ) ;
%results: fail to reject at bands = 1,2,3,4; sig. diff at band 5
clear('per');clear('con');





%%%+++++++++++++DA+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\03_17_2014_da.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\03_17_2014_conDA.mat');
% make analysis data set
filepath.da=da.filepath;
filepath.conda=con.filepath;
chem=da;  conFIndex=[1:5]; %which files in con are PER baseline
cont=con; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

DAanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
DA = do_t_anova(DAanova, want ) ;
% DA{1}
clear('da');clear('con');



%%%+++++++++++++Ace+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\07_30_2014_ace.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\07_30_2014_conace.mat');
% make analysis data set
filepath.ace=ace.filepath;
filepath.conace=con.filepath;
chem=ace;  conFIndex=[1:5]; %which files in con are PER baseline
cont=con; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0; want.diff=0;want.pchange=1;

ACEanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
ACE = do_t_anova(ACEanova, want ) ;
clear('ace');clear('con');





%%%+++++++++++++H20+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\07_7_2014_h20.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\08_05_2014_conh20.mat');
% make analysis data set
filepath.h20=h20.filepath;
filepath.conh20=con.filepath;
chem=h20;  conFIndex=[1:10]; %which files in con are PER baseline
cont=con; cont.file = con.file(conFIndex) ;
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
chem=tri;  conFIndex=[1:7]; %which files in con are PER baseline
cont=con; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 
want.sec=0; want.ch=0; want.well=0;want.diff=0;want.pchange=1;

TRIanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
TRI = do_t_anova(TRIanova, want ) ;
clear('tri');clear('con');

% filename='F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\8_26_2013_pChangeanova.mat';
 save(filename, 'PER', 'PERanova', 'DA','DAanova',...
    'TRI', 'TRIanova', 'BIC','BICanova', 'CAR','CARanova',...
    'PER50','PER50anova', 'H20','H20anova', 'ACE','ACEanova', 'filepath') ;





% ++++++++++++++++++++++++++++++++++++++
% write results
xlsfilestr='F:\EXTRAP\documents\Results\lfp_results_8-26-2014.xls';
title='Results of Mean % Control in Power averaged across All Electrodes';
xlswrite(xlsfilestr,title, 1,'A1' ) ;

chemName= {'bicuculline','carbaryl','permethrin','permethrin50',...
    'DemoicAcid','Acetaminophen','H20', 'TRI'}
for k=1:size(chemName,2)
% write files
    if strcmp( chemName{1,k},'permethrin' )
        startCell=3;
        anovaData = PERanova ; test = PER;
        [p,n1,x] = fileparts( filepath.per{1} );
        [p,n2,x] = fileparts( filepath.per{2} );
        [p,n3,x] = fileparts( filepath.per{3} );
        [p,n4,x] = fileparts( filepath.per{4} );
        [p,n5,x] = fileparts( filepath.per{5} );
    elseif strcmp( chemName{1,k},'carbaryl')
        startCell=11;
        anovaData = CARanova ; test = CAR
        [p,n1,x] = fileparts( filepath.car{1} );
        [p,n2,x] = fileparts( filepath.car{2} );
        [p,n3,x] = fileparts( filepath.car{3} );
        [p,n4,x] = fileparts( filepath.car{4} );
        [p,n5,x] = fileparts( filepath.car{5} );

    elseif strcmp( chemName{1,k},'bicuculline' )
        startCell=20;
        anovaData = BICanova ; test = BIC;
        [p,n1,x] = fileparts( filepath.bic{1} );
        [p,n2,x] = fileparts( filepath.bic{2} );
        [p,n3,x] = fileparts( filepath.bic{3} );
        [p,n4,x] = fileparts( filepath.bic{4} );
        [p,n5,x] = fileparts( filepath.bic{5} );
        
    elseif strcmp( chemName{1,k},'DemoicAcid')
        startCell=29;
        anovaData = DAanova ; test = DA
        [p,n1,x] = fileparts( filepath.da{1} );
        [p,n2,x] = fileparts( filepath.da{2} );
        [p,n3,x] = fileparts( filepath.da{3} );
        [p,n4,x] = fileparts( filepath.da{4} );
        [p,n5,x] = fileparts( filepath.da{5} );
            
    elseif strcmp( chemName{1,k},'permethrin50')
        startCell=38;
        anovaData = PER50anova ; test = PER50
        [p,n1,x] = fileparts( filepath.per50{1} );
        [p,n2,x] = fileparts( filepath.per50{2} );
        [p,n3,x] = fileparts( filepath.per50{3} );
        [p,n3,x] = fileparts( filepath.per50{4} );
        [p,n3,x] = fileparts( filepath.per50{5} );
        [p,n3,x] = fileparts( filepath.per50{6} );
        
    elseif strcmp( chemName{1,k},'Acetaminophen')
        startCell=47;
        anovaData = ACEanova ; test = ACE
        [p,n1,x] = fileparts( filepath.ace{1} );
        [p,n2,x] = fileparts( filepath.ace{2} );
        [p,n3,x] = fileparts( filepath.ace{3} );
        [p,n4,x] = fileparts( filepath.ace{4} );
        [p,n5,x] = fileparts( filepath.ace{5} );
            
    elseif strcmp( chemName{1,k},'H20')
        startCell=56;
        anovaData = H20anova ; test = H20
        [p,n1,x] = fileparts( filepath.h20{1} );
        [p,n2,x] = fileparts( filepath.h20{2} );
        [p,n3,x] = fileparts( filepath.h20{3} );
        [p,n4,x] = fileparts( filepath.h20{4} );
        [p,n5,x] = fileparts( filepath.h20{5} );
        [p,n6,x] = fileparts( filepath.h20{6} );
        [p,n7,x] = fileparts( filepath.h20{7} );
        [p,n8,x] = fileparts( filepath.h20{8} );
        [p,n9,x] = fileparts( filepath.h20{9} );
        [p,n10,x] = fileparts( filepath.h20{10} );
        
       elseif strcmp( chemName{1,k},'TRI')
        startCell=69;
        anovaData = TRIanova ; test = TRI;
        [p,n1,x] = fileparts( filepath.tri{1} );
        [p,n2,x] = fileparts( filepath.tri{2} );
        [p,n3,x] = fileparts( filepath.tri{3} );
        [p,n4,x] = fileparts( filepath.tri{4} );
        [p,n5,x] = fileparts( filepath.tri{5} );
        [p,n6,x] = fileparts( filepath.tri{6} );
        [p,n7,x] = fileparts( filepath.tri{7} );
       
    end;
    ex{1,1}=' '; ex{2,1}='1-4Hz'; ex{3,1}='4-8Hz'; ex{4,1}='8-14Hz'; 
    ex{5,1}='14-30Hz'; ex{6,1}='30-50Hz';
    nFiles = length(anovaData{1}.pchange ) ;
    for curfile=1:nFiles
        ex{1,1+curfile}=evalin( 'base', strcat('n', num2str(curfile) ) ) ;
        %ex{1,6}=n5;
    end;
    
    for j=1:5
        for i=1:length(anovaData{1}.pchange)
            ex{1+j,i+1} = anovaData{j}.pchange(i);
        end;
    end;
    ex{1,2+nFiles}='mean';
    ex{1,3+nFiles}='se';
    ex{1,4+nFiles}='sd'; 
    ex{1,5+nFiles}='t'; 
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
end; 












%%++  whole channel psd analysis
%%%+++++++++++++CAR+++++++++++++++
% load data
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_29_2014_car.mat');
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\01_29_2014_carcon.mat');
% make analysis data set
chem=car;  % conFIndex=[11:15]; %which files in con are CAR baseline
cont=carcon; cont.file = con.file(conFIndex) ;
nFiles=length(cont.file) ; 

% get AE vector
temp = [1:60] ;

index_want=find(4<=f&<=8);
chipDiff = [];
for curFile=1:5
    chDiff = [];
    ae = temp([con.file(10+curFile).channel([1:60]).ae] );
    for curCh=ae
        temp2 = ...
            mean(car.file(curFile).channel(curCh).psd(index_want)...
            - carcon.file(curFile).channel(curCh).psd(index_want) ) ;
        chDiff = [chDiff temp2];
    end;
    chipDiff = [chipDiff mean(chDiff) ];
end;

CARanova = make_t_anova_data(chem, cont, nFiles, want) ;
% run analysis
CAR = do_t_anova(CARanova, want ) ;
%results: fail to reject at all bands 







% gen_psd_electrode_fig.m
% diana hall
% 11/06/2013
% purpose: to make lattice plot of psd pre and post trt
%  

% load data
 load( 'F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_all.mat' )
 
% add to path code that allows multiple pdf appending
addpath('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_functions' ) ;
% fisrt arg is new pdf to be made, next are all the ones you want combines
% Example:
%append_pdfs( 'G:\data\papers\Hellooo.pdf', 'G:\data\papers\3.pdf',...
%    'G:\data\papers\6.pdf' )

% try for one electrode first
curFile=1; ent= 55 ;  sec = 1 ; band = 1 ;

% take mean of obs. at each freq within band across seconds in channel
figure() ;
for band=1:5

    plot(bic.fbands{band} , mean ( bic.file(curFile).channel(ent).allPw{band}, 1 ) )

end;





%++++ All PSD +++++++++++++++++++++++++++++++++++++

% ++++++ BIC ++++++++++++++++++++++++++++++++
% plot of all psd
% why didn't we just take psd of whole file, how would it compare to
% average over seconds
curFile=1 ; 
[pathstr,name,ext] = fileparts( files.bic{curFile} ) ;
for ent=1:4:60
    
    %filename_ind = strcat('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\',...
    % 'PSD_File_',num2str( name(1:6) ),'_E',num2str(ent),'.pdf'  ) ;
    filename_ind = strcat('C:\Users\Diana\Downloads\PSD_File_',num2str( name(1:6) ),'_E',num2str(ent),'.pdf'  )
    
    for inc=0:3
        
        subplot(2,2, inc+1 )
                
        plot(con.file(1).channel(1).f , ...
            mean ( con.file(curFile).channel(ent+inc).allPwallBands , 1 ),'b',...
            bic.file(1).channel(1).f,...
            mean ( bic.file(curFile).channel(ent+inc).allPwallBands , 1 ),'r')
        hold on;
        lhandle = legend('cont', 'bic' ) ;
        set(lhandle, 'Location', 'NorthEast' )
        
        [pathstr,name,ext] = fileparts( files.bic{curFile} ) ;
        title( {strcat('PSD File ',name(1:6) , ' Elect',num2str(ent+inc) ), '  '} )
                xlabel('Frequency (Hz)')
                ylabel('Magnitude')

        xlim([0,50])
        ylim([0, 2*10^-6] )
        hold off;
        
    end;
    %pdf will be in color
    print('-dpdf','-r300',filename_ind )
  
end;

filename_all = strcat('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\',...
    'PSD_BIC&Con_File_',name(1:6),name(8:12),'_01_27.pdf'  ) ;
append_pdfs( filename_all,...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E1.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E5.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E9.pdf',...    
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E13.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E17.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E21.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E25.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E29.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E33.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E37.pdf',...  
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E41.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E45.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E49.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E53.pdf',...    
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_090413_E57.pdf')






% ++++++ CAR 5 uM++++++++++++++++++++++++++++++++
% plot of all psd
% why didn't we just take psd of whole file, how would it compare to
% average over seconds
curFile=4 ; 
[pathstr,name,ext] = fileparts( files.car{curFile} ) ;
for ent=1:4:60
    
    filename_ind = strcat('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\',...
     'PSD_File_',name(1:6),'_E',num2str(ent),'.pdf'  ) ;
    for inc=0:3
        
        subplot(2,2, inc+1 )
                
        plot(con.file(1).channel(1).f , ...
            mean ( con.file(curFile+length(bic.file) ).channel(ent+inc).allPwallBands , 1 ),'b',...
            car.file(1).channel(1).f,...
            mean ( car.file(curFile).channel(ent+inc).allPwallBands , 1 ),'r')
        hold on;
        lhandle = legend('cont', 'car' ) ;
        set(lhandle, 'Location', 'NorthEast' )
        
        [pathstr,name,ext] = fileparts( files.car{curFile} ) ;
        title( {strcat('PSD File ',name(1:6) , ' Elect',num2str(ent+inc) ), '  '} )
                xlabel('Frequency (Hz)')
                ylabel('Magnitude')

        xlim([0,50])
        ylim([0, 2*10^-6] )
        hold off;
        
    end;
    %pdf will be in color
    print('-dpdf','-r300',filename_ind )
  
end;

filename_all = strcat('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\',...
    'PSD_CAR_5uM&Con_File_',name(1:6),name(8:12),'.pdf'  ) ;
append_pdfs( filename_all,...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E1.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E5.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E9.pdf',...    
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E13.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E17.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E21.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E25.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E29.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E33.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E37.pdf',...  
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E41.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E45.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E49.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E53.pdf',...    
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113 19855 CAR FPSPK0002_E57.pdf')









% ++++++ CAR 30 uM++++++++++++++++++++++++++++++++
% plot of all psd
% why didn't we just take psd of whole file, how would it compare to
% average over seconds
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_car.mat') ;
load('F:\EXTRAP\matlab_Cina_Herr_Diana\finshed_data\12_04_2013_con.mat') ;

 file1='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 19860 CAR\090513 19860 CAR FPSPK0001.mcd';
 file2='F:\EXTRAP\FP SPIKE\09-2013 BIC & CAR\September5, 2013 BIC and CAR\September 5,2013 Carbaryl\090513 18416 CAR\090513 18416 CAR FPSPK0001.mcd';
 file3='F:\EXTRAP\FP SPIKE\September 19, 2013 CAR\09-19-13 CAR  20613\091913 20613 CAR FPSPK0001.mcd' ;
 file4='F:\EXTRAP\FP SPIKE\Oct 1 2013 CAR TRI\October 1 2013 CAR 15157\100113 15157 CAR FPSPK0001.mcd';
 % 20596 is broken
 %file5='F:\EXTRAP\FP SPIKE\Oct 1 2013 CAR TRI\October 1 2013 CAR 20596\100113 20596 CAR FPSPK0002.mcd' ;
 file5='F:\EXTRAP\FP SPIKE\December 2013\December 4, 2013 CAR\120413 18397 CAR FPSPK0001.mcd' ;
 file6='F:\EXTRAP\FP SPIKE\December 2013\December 17, 2013 CAR\December 17, 2013 CAR 18412\121713 18412 CAR FPSPK0001.mcd';
 file7='F:\EXTRAP\FP SPIKE\December 2013\December 17, 2013 CAR\December 17, 2013 CAR 19858\121713 19858 CAR FPSPK0001.mcd';
 file8='F:\EXTRAP\FP SPIKE\December 2013\December 17, 2013 CAR\December 17, 2013 CAR 20599\121713 20599 CAR FPSPK0001.mcd';

files.car={file1, file2, file3, file4, file5, file6, file7, file8} ;

curFile=5 ;  numbic=4;
[pathstr,name,ext] = fileparts( files.car{curFile} ) ;
for ent=1:4:60
    
    filename_ind = strcat('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\',...
     'PSD_File_',name(1:6),'_E',num2str(ent),'.pdf'  ) ;
    for inc=0:3
        
        subplot(2,2, inc+1 )
                
        plot(con.file(1).channel(1).f , ...
            mean ( con.file(curFile+numbic ).channel(ent+inc).allPwallBands , 1 ),'b',...
            car.file(1).channel(1).f,...
            mean ( car.file(curFile).channel(ent+inc).allPwallBands , 1 ),'r')
        hold on;
        lhandle = legend('cont', 'car' ) ;
        set(lhandle, 'Location', 'NorthEast' )
        
        [pathstr,name,ext] = fileparts( files.car{curFile} ) ;
        title( {strcat('PSD File ',name(1:6) , ' Elect',num2str(ent+inc) ), '  '} )
                xlabel('Frequency (Hz)')
                ylabel('Magnitude')

        xlim([0,50])
        ylim([0, 3*10^-7] )
        hold off;
        
    end;
    %pdf will be in color
    print('-dpdf','-r300',filename_ind )
  
end;

filename_all = strcat('F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\',...
    'PSD_CAR_30uM&Con_File_',name(1:6),name(8:12),'.pdf'  ) ;
append_pdfs( filename_all,...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E1.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E5.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E9.pdf',...    
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E13.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E17.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E21.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E25.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E29.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E33.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E37.pdf',...  
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E41.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E45.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E49.pdf',...
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E53.pdf',...    
    'F:\EXTRAP\matlab_Cina_Herr_Diana\matlab_figures\PSD_File_112113_E57.pdf')







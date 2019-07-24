% make_t_anova_data.m
% Diana Hall
% 1-27-2014
% purpose: to create data and analyze 


function[CHEManova ] = make_t_anova_data(chem, cont, nFiles, want)
% this function makes a data set
% example call: 
% chem=bic; 
% cont=con; cont.file = con.file(conFIndex) ;
% nFiles=length(cont.file) ; conFIndex=1:5; 
% want.sec=0; want.ch=0; want.well=0; want.diff=1;
% CHEManova = make_t_anova_data(chem, cont, nFiles, want) ;


 % make ANOVA data
 %  by second
 % set your denomination
    if want.sec
        for curB=1:5
            CHEManova{curB}.sec=[];
            for curFile=1:nFiles
                for curEnt=1:60
                    temp = [cont.file(curFile).bandM(1:300,curEnt,curB)  chem.file(curFile).bandM(1:300,curEnt,curB)];
                    CHEManova{curB}.sec = vertcat(CHEManova{curB}.sec, temp) ;
                end;

            end;

        end;
    elseif want.ch
        for curB=1:5
            CHEManova{curB}.ch=[];
            for curFile=1:nFiles
                for curEnt=1:60
                    temp = [mean( cont.file(curFile).bandM(1:300,curEnt,curB))...
                        mean( chem.file(curFile).bandM(1:300,curEnt,curB)) ];
                    CHEManova{curB}.ch = vertcat(CHEManova{curB}.ch, temp) ;
                end;

            end;

        end;
    elseif want.well
        % by well
        for curB=1:5
            CHEManova{curB}.well=[];
            for curFile=1:nFiles
                temp_well=[];
                for curEnt=1:60
                    temp = [mean( cont.file(curFile).bandM(1:300,curEnt,curB))...
                        mean( chem.file(curFile).bandM(1:300,curEnt,curB)) ];
                    temp_well= vertcat(temp_well, temp) ;
                end;
                CHEManova{curB}.well = vertcat(CHEManova{curB}.well , mean(temp_well) ) ;
            end;

        end;

    elseif want.diff

        for curB=1:5
            CHEManova{curB}.diff=[];
            for curFile=1:nFiles
                
                % by well paired differences
                s=[];
                for i=1:60 s=vertcat(s,cont.file(curFile).channel(i).ae); end;
                seq=1:60;
                chWanted = transpose(s).*seq;
                chWanted = chWanted(chWanted~=0) ;
                temp_diff=[];
                % precaution for files containing less than 300 seconds
                nsecc = min(cont.file(curFile).nseconds, 300 );
                nsect = min(chem.file(curFile).nseconds, 300 );
                for curEnt=1:60
                    temp = mean( chem.file(curFile).bandM(1:nsect,curEnt,curB))-...
                       mean( cont.file(curFile).bandM(1:nsecc,curEnt,curB)) ;
                    temp_diff= vertcat(temp_diff, temp) ;
                end;
                CHEManova{curB}.diff = vertcat(CHEManova{curB}.diff, mean(temp_diff) ) ;
            end;

        end;

        elseif want.pchange

        for curB=1:5
            CHEManova{curB}.pchange=[];
            for curFile=1:nFiles
                
                % by well paired differences
                s=[];
                for i=1:60 s=vertcat(s,cont.file(curFile).channel(i).ae); end;
                seq=1:60;
                chWanted = transpose(s).*seq;
                chWanted = chWanted(chWanted~=0) ;
                temp_pchange=[];
                % precaution for files containing less than 300 seconds
                nsecc = min(cont.file(curFile).nseconds, 300 );
                nsect = min(chem.file(curFile).nseconds, 300 );
                for curEnt=1:60
                    temp = mean( chem.file(curFile).bandM(1:nsect,curEnt,curB))/...
                       mean( cont.file(curFile).bandM(1:nsecc,curEnt,curB)) ;
                    temp_pchange = vertcat(temp_pchange, temp) ;
                end;
                CHEManova{curB}.pchange = vertcat( CHEManova{curB}.pchange, mean(temp_pchange) ) ;
            end; %through curFiles

        end; %through curBand
    end; %end of if data making

end %end of make_data


% do_t_anova.m
% Diana Hall
% 1-27-2014
% purpose: to take data output from "make_t_anova_data.m
%  and test anova or paired t-test model of differences


function[ anova ] = do_t_anova( CHEManova, want )
    if want.diff 
        % diff t-test
        [anova{1}.t, anova{1}.tp_value, anova{1}.ci, anova{1}.stats ]=...
            ttest(CHEManova{1}.diff) %significant
        [anova{2}.t, anova{2}.tp_value, anova{2}.ci, anova{2}.stats ]=...
            ttest(CHEManova{2}.diff ) %significant
        [anova{3}.t, anova{3}.tp_value, anova{3}.ci, anova{3}.stats ]=...
            ttest(CHEManova{3}.diff) %significant
        [anova{4}.t, anova{4}.tp_value, anova{4}.ci, anova{4}.stats ]=...
            ttest(CHEManova{4}.diff) %significant
        [anova{5}.t, anova{5}.tp_value, anova{5}.ci, anova{5}.stats ]=...
            ttest(CHEManova{5}.diff) %significant

    elseif want.well

        % well level
        anova{1}.wanova = anova1(CHEManova{1}.well) %significant
        anova{2}.wanova = anova1(CHEManova{2}.well) %significant
        anova{3}.wanova = anova1(CHEManova{3}.well) %significant
        anova{4}.wanova = anova1(CHEManova{4}.well) %significant
        anova{5}.wanova = anova1(CHEManova{5}.well) %not significant 
    elseif want.ch
        % channellevel
        anova{1}.chanova=anova1(CHEManova{1}.ch) %significant
        anova{2}.chanova=anova1(CHEManova{2}.ch) %significant
        anova{3}.chanova=anova1(CHEManova{3}.ch) %significant
        anova{4}.chanova=anova1(CHEManova{4}.ch) %significant
        anova{5}.chanova=anova1(CHEManova{5}.ch) %not significant 
    elseif want.pchange
               % diff t-test
        [anova{1}.t, anova{1}.tp_value, anova{1}.ci, anova{1}.stats ]=...
            ttest(CHEManova{1}.pchange,  1) 
        [anova{2}.t, anova{2}.tp_value, anova{2}.ci, anova{2}.stats ]=...
            ttest(CHEManova{2}.pchange, 1 ) %significant
        [anova{3}.t, anova{3}.tp_value, anova{3}.ci, anova{3}.stats ]=...
            ttest(CHEManova{3}.pchange, 1 ) %significant
        [anova{4}.t, anova{4}.tp_value, anova{4}.ci, anova{4}.stats ]=...
            ttest(CHEManova{4}.pchange, 1 ) %significant
        [anova{5}.t, anova{5}.tp_value, anova{5}.ci, anova{5}.stats ]=...
            ttest(CHEManova{5}.pchange, 1 ) %significant 

    end; 
    
end % function










function [res,fit_param_factor_chosen,i_chosen, fit_param_sorted] = GreedyLSQ(Coefficients,param)
close all; clc;

% format long e
NumFac=size(Coefficients,2);
options = optimset('MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9,'TolX',1e-9);

%GENERAL 

%first factor
i_chosen=0;
for i=1:NumFac 
    if min((i ~= i_chosen)) == 1
        param_0 = ones(2,1);
        Chosen_factor=Coefficients(:,i);
        [fit_param_factors ,res_fit_param] = lsqcurvefit(@fit_param_model,param_0,[Chosen_factor],param,[],[],options) ;
        res_fit_param_norm(i)=res_fit_param ;
        fit_param_factors_all(:,i)=fit_param_factors; 
    end
end

[res(1), i_t]=min(res_fit_param_norm);
fit_param_factor_chosen=fit_param_factors_all(:,i_t);
i_chosen=i_t;

for j=1:NumFac-1 
    clear fit_param_factors_all
    for i=1:NumFac   
        if min((i ~= i_chosen)) == 1
            param_0 = [fit_param_factor_chosen;1];
            Chosen_factor=Coefficients(:,i);
            [fit_param_factors ,res_fit_param] = lsqcurvefit(@fit_param_model,param_0,[Coefficients(:,i_chosen), Chosen_factor],param,[],[],options) ;
            res_fit_param_norm(i)=res_fit_param ;
            fit_param_factors_all(:,i)=fit_param_factors; 
        end
    end

    [res(j+1), i_t]=min(res_fit_param_norm);
    fit_param_factor_chosen=fit_param_factors_all(:,i_t);
    i_chosen=[i_chosen i_t];

end

res=res';


for i=1:NumFac 
    %sorts the factor coeff based on the factor_tag 
    fit_param_sorted(i_chosen(i),1)=fit_param_factor_chosen(i+1); 
    
end


return

%%
function fit_k_estimate = fit_param_model(m,factordata)

x=m(2:end);

fit_k_estimate = m(1) + factordata*x ;

return
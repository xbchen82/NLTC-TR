function [RMSE, MAE, R2,RELErr1,RELErr2] = compute_measure(v_pred,v_real)
% m1 = mean(v_pred);
% m2 = mean(v_real);
% s1 = std(v_pred,1);
% s2 = std(v_pred,1);
% PA = (v_pred-m1)'*(v_real-m2)/(length(v_pred)-1)/s1/s2;
R2 = corr(v_pred,v_real)^2;
MAE = mean(abs(v_pred-v_real));
RMSE = sqrt(mean((v_pred-v_real).^2));
RELErr1 = mean(abs(v_pred-v_real)./v_real);
RELErr2 = norm(v_pred-v_real)/norm(v_real)*100;
% fprintf(1,[method_name ' RMSE: %.4f MAE: %.4f R2: %.4f \n'],RMSE,MAE,R2);
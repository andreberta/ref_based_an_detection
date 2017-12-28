%#ok<*AGROW>
%#ok<*SAGROW>
%% explanation
%in this script the data generated from the train is loaded and manipulated
%to be visualized to study the results


%% load data
load('thresh.mat')
load('thresh.mat')
%% plot thresh(sigma)

figure,plot([sigma{:}],[thres_fov_sigma{:}],'xk')
title('Thres(sigma) fov')
xlabel('sigma')
ylabel('thresh(sigma)')

figure,plot([sigma{:}],[thres_NL_sigma{:}],'xk')
title('Thres(sigma) NL')
xlabel('sigma')
ylabel('thresh(sigma)')


%% plot relation AUC - PSNR

marker_colour = {'xm','xc','xr','xb','xk'};

%NL
figure
for ii = 1:sigma_lenght
    plot([AUC_NL{ii,:}],[PSNR_NL{ii,:}],marker_colour{ii})
    hold on
end
title('NL AUC-PSNR')
xlabel('AUC')
ylabel('PSNR')
legend('10','20','30','50','70');

%fov
figure
for ii = 1:sigma_lenght
    plot([AUC_fov{ii,:}],[PSNR_fov{ii,:}],marker_colour{ii})
    hold on
end
title('fov AUC-PSNR')
xlabel('AUC')
ylabel('PSNR')
legend('10','20','30','50','70');


%% plot mean relation AUC - PSNR for every value of sigma
%NL
figure
for ii = 1:sigma_lenght
    plot(mean([AUC_NL{ii,:}]),mean([PSNR_NL{ii,:}]),marker_colour{ii})
    hold on
end
title('NL AUC-PSNR')
xlabel('AUC')
ylabel('PSNR')
legend('10','20','30','50','70');

%fov
figure
for ii = 1:sigma_lenght
    plot(mean([AUC_fov{ii,:}]),mean([PSNR_fov{ii,:}]),marker_colour{ii})
    hold on
end
title('fov AUC-PSNR')
xlabel('AUC')
ylabel('PSNR')
legend('10','20','30','50','70');
%% mean ROC
line_colour = {'m','c','r','b','k'};
figure;
for jj = 1:sigma_lenght
    whole_fpr_NL = [];
    whole_tpr_NL = [];
    whole_cut_off_NL = [];
    for ii = 1:images_length
        [~,fpr_NL,tpr_NL] = ROC( anomaly_NL{jj,ii} , compare , 0 );
        
        whole_fpr_NL = [whole_fpr_NL ; fpr_NL]; 
        whole_tpr_NL = [whole_tpr_NL ; tpr_NL];
    end
    
    [m,n] = size(whole_fpr_NL);
    
    for kk = 1:n
        final_fpr_NL(kk) = mean(whole_fpr_NL(:,kk));
        final_tpr_NL(kk) = mean(whole_tpr_NL(:,kk));
    end
    
    
    plot(final_fpr_NL,final_tpr_NL,line_colour{jj}), hold on ,
    title('ROC NL')

    trapz(final_fpr_NL,final_tpr_NL);
end
xlabel('false positive rate') % x-axis label
ylabel('true positive rate') % y-axis label
legend('10','20','30','50','70')

figure;
for jj = 1:sigma_lenght
    whole_fpr_fov = [];
    whole_tpr_fov = [];
    whole_cut_off_fov = [];
    for ii = 1:images_length
        [~,fpr_fov,tpr_fov] = ROC( anomaly_fov{jj,ii} , compare , 0 );

        whole_fpr_fov = [whole_fpr_fov ; fpr_fov];
        whole_tpr_fov = [whole_tpr_fov ; tpr_fov];  
    end
        
    [m,n] = size(whole_fpr_fov);
    
    for kk = 1:n
        final_fpr_fov(kk) = mean(whole_fpr_fov(:,kk));
        final_tpr_fov(kk) = mean(whole_tpr_fov(:,kk)); 
    end
    

    plot(final_fpr_fov,final_tpr_fov,line_colour{jj}), hold on ,
    title('ROC fov')

    trapz( final_fpr_fov,final_tpr_fov);
end
xlabel('false positive rate') % x-axis label
ylabel('true positive rate') % y-axis label
legend('10','20','30','50','70')






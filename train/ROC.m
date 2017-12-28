function [AUC,fpr,tpr] = ROC( anomaly , A , plot_ROC )

k = min(anomaly(:)):abs(min(anomaly(:))-max(anomaly(:)))/100:max(anomaly(:));
k_length = numel(k);
false_positive_rate_a = zeros([1,k_length]);
true_positive_rate_a = zeros([1,k_length]);

for ii=1:k_length
    an = anomaly > k(ii);
    true_positive = sum(sum(an(A)));
    false_positive = sum(sum(an(~A)));
    true_negative = sum(sum(~an(~A)));
    false_negative = sum(sum(~an(A)));
    false_positive_rate_a(ii) = false_positive / (false_positive + true_negative);
    true_positive_rate_a(ii) = true_positive/ (true_positive + false_negative);
end
fpr = [ 0, false_positive_rate_a(length(false_positive_rate_a):-1:1), 1];
tpr = [0,true_positive_rate_a(length(true_positive_rate_a):-1:1),1];
AUC =  trapz( fpr,tpr );

if(plot_ROC == 1)
    figure;
    plot(fpr,tpr), hold on ,
    title('ROC')
    xlabel('false positive rate') % x-axis label
    ylabel('true positive rate') % y-axis label
    legend(strcat('AUC = ',num2str(AUC)))
end


end


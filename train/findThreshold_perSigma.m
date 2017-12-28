%#ok<*AGROW>

se_dilate = strel('square',4);
selector = imdilate(compare,se_dilate);

thres_NL_sigma = cell(sigma_lenght,1);
thres_fov_sigma = cell(sigma_lenght,1);

for jj=1:sigma_lenght
    whole_anomaly_NL = [];
    whole_anomaly_fov = [];
    for ii=1:images_length
        tmp_NL = anomaly_NL{jj,ii};
        tmp_fov = anomaly_fov{jj,ii};
        whole_anomaly_NL = [whole_anomaly_NL; tmp_NL(~selector);];
        whole_anomaly_fov = [whole_anomaly_fov; tmp_fov(~selector)];
    end
    thres_NL_sigma{jj} = prctile(whole_anomaly_NL,99);
    thres_fov_sigma{jj} = prctile(whole_anomaly_fov,99);
end



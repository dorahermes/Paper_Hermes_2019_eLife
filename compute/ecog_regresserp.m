function data_epoch = ecog_regresserp(data_epoch,t_base,stims)
% function to regress the ERP from some data
% USAGE:
% data_epoch = ecog_regresserp(data_epoch,t_base,stims)
% 
% data_epoch % electrodes X epochs X time
% t_base % indices of the baseline period
% stims % different code for different conditions

% baseline correct
for k=1:size(data_epoch,1)%channels
    for m=1:size(data_epoch,2)%epochs
        x=squeeze(data_epoch(k,m,:));
        x=x-mean(x(t_base));
        data_epoch(k,m,:)=x;
    end
end

% regress erp out
for k=1:size(data_epoch,1)%channels
    disp(['regressing erp from channel ' int2str(k)])
    for m=1:size(data_epoch,2)%epochs
        x=squeeze(data_epoch(k,m,:));
        % regress ERP out
        s=stims(m);
        av_erp=squeeze(mean(data_epoch(k,stims==s,:),2));
        [~,~,reg_R] = regress(x,av_erp);
        data_epoch(k,m,:)=reg_R;
    end
end

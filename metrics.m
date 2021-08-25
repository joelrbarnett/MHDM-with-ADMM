%Function metrics: returns metrics for regular multiscale decomposition
% F_orig is original image
% F_data is noisy image
% xkArray is the array of multiscale restored images
function [xk_f_norm2,rmse_final,stopCrit,snr]= metrics(F_orig,F_data,xkArray,numScales,tightFlag)
[m,n]=size(F_orig);
delta_2= sum(sum(F_data./F_orig + log(F_orig) - log(F_data) -ones(size(F_orig))));
rmse_final=zeros(numScales,1);
xk_f_norm2=zeros(numScales,1);
stopCrit=zeros(numScales,1);
snr=zeros(numScales,1);
for k=1:numScales
    %L2 error and RMSE error at each iteration
    xk_f_norm2(k)=norm(F_orig-xkArray(:,:,k),'fro');% difference from true image
    rmse_final(k)=xk_f_norm2(k)/sqrt(m*n);   %RMSE error from true image

    %Bregman distance from F_data to Txk (restored, blurred image)
    D_f_data_Txk=sum(sum(F_data./xkArray(:,:,k) + log(xkArray(:,:,k)) - log(F_data)-ones(size(F_orig))  ));
    if tightFlag(1)==0 %stopping criteria for regular 
        stopCrit(k)= (D_f_data_Txk)/(delta_2);%ratio of bregman distances
    elseif tightFlag(1)==1 %stopping criteria for tight and refined versions
        alp = tightFlag(2)/(k)^(3/2);
        stopCrit(k)= (D_f_data_Txk+ alp*TV(log(xkArray(:,:,k))))/(delta_2);%ratio of bregman distances
    end
    %SNR at each iteration
    snr(k)=SNR(xkArray(:,:,k),F_orig);
end

%to compute SNR 
function snr=SNR(im,ref)
    snr= 20.*log10(norm(ref,'fro')/norm(im-ref,'fro'));
end

%function to compute the non-regularized total variation of image A
function tvA=TV(A)
    A=padarray(A,[1,1],'replicate');
    DXF=A(3:end,2:end-1)-A(2:end-1,2:end-1);
    DYF=A(2:end-1,3:end)-A(2:end-1,2:end-1);
    tvA=sum(sum(DXF.^2+DYF.^2));
end


end
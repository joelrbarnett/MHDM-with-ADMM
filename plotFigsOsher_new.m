
function plotFigsOsher_new(F_orig, F_data, xkArray,params,filePrefix,figPrefix,saveFlag,tightFlag)

[m,n]=size(F_orig);
numScales = length(xkArray(1,1,1,:));
%tightFlag=[0,0]; %tightFlag(1)=indicates tight or not, tightFlag(2)= value of alp0
if length(params)==6
    paramsCell=num2cell(params);
    [maxIters, dt, epsilon, lambda0, q,alp0]=paramsCell{:};
else 
    paramsCell=num2cell(params);
    [maxIters, dt, epsilon, lambda0,q]=paramsCell{:};
end
[xk_f_norm2,rmse_final,stopCrit,snr]= metrics(F_orig,F_data,squeeze(xkArray),numScales,tightFlag);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find min error and stopping criteria, original RMSE
%Min error and index of that error. Divide by (m*n) to get RMSE
[minVal,mink]=min(rmse_final);
%mink=mink-1;
%k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
k_star=min(find((stopCrit<=1)==1));
if k_star==1
    k_star=1;
elseif ~isempty(k_star)&&k_star>1
    k_star=k_star-1;
end
noisyRMSE=norm(F_orig-F_data,'fro')/sqrt(m*n); %original RMSE error
noisySNR=20.*log(norm(F_orig,'fro')/norm(F_orig-F_data,'fro'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maxIters=params(1); %time iterations in solving for wk
% dt=params(2);
% epsilon=params(3); %for regularizing TV
% lambda0=params(4); 
% q = params(5);
% alp0=params(6); %for tight/refined


%Plot the original, degraded, and recovered images
figure()
image(F_orig); axis image; axis off; colormap(gray(256));
title('Original Image','FontSize',16)
% 
figure()
image(F_data); axis image; axis off; colormap(gray(256));
title(['Noisy Image. RMSE=',num2str(noisyRMSE),', SNR=',num2str(noisySNR)],'FontSize',16)
if saveFlag==1
    figName=filePrefix+figPrefix+"noisy.fig";
    savefig(figName)
end
%Restored image
figure()
image(xkArray(:,:,1,mink)); axis image; axis off; colormap(gray(256));
subtitle="\lambda_k ="+lambda0+ "*"+q +"^k";
titleStr="Restored image: k="+num2str(mink-1)+" RMSE="+num2str(minVal)+", SNR="+num2str(snr(mink));
title([titleStr,subtitle],'FontSize',16)
if saveFlag==1
    figName=filePrefix+figPrefix+"restored.fig";
    savefig(figName)
end
%k_star image
if (~isempty(k_star)) && (k_star~=mink)
    figure()
    image(xkArray(:,:,1,k_star)); axis image; axis off; colormap(gray(256));
    subtitle="\lambda_k ="+lambda0+ "*"+q +"^k";
 titleStr="Restored image: k^*="+num2str(k_star-1)+" RMSE="+num2str(rmse_final(k_star))+", SNR="+num2str(snr(k_star));
title([titleStr,subtitle],'FontSize',16)
 if saveFlag==1
        figName=filePrefix+figPrefix+"restored_kstar.fig";
        savefig(figName)
    end
end


%Plot montage of multiscales near optimal
% figure()
% numPlots =min(9,length(xkArray(1,1,1,:)));
% montage(xkArray(:,:,:,end-numPlots+1:end),gray(256))
% 
% titleStr = "Multiscales k = "+num2str(numScales -numPlots+1)+" through k= "+num2str(numScales);
% title(titleStr,'FontSize',16)
% 
% if saveFlag==1
%     figName=filePrefix+figPrefix+"multiscales.png";
%     saveas(gcf,figName)
% end
% 
% %Plot montage of residuals near optimal
% figure()
% numPlots =min(9,length(xkArray(1,1,1,:)));
% residuals = xkArray;
% residuals(:,:,1,:) = abs(F_orig-residuals(:,:,1,:))+125;
% montage(residuals(:,:,:,end-numPlots+1:end),gray(256))
% 
% titleStr = "Residuals k = "+num2str(numScales -numPlots+1)+" through k= "+num2str(numScales);
% title(titleStr,'FontSize',16)
% 
% if saveFlag==1
%     figName=filePrefix+figPrefix+"residuals.png";
%     saveas(gcf,figName)
% end

%Plot the RMSE and SNR, as well as the stopping criterion plot
figure('position',[100,100,1150,400])
subplot(1,2,1)
yyaxis left
plot(0:numScales-1,rmse_final)
xlabel('Multiscales: k','FontSize',16)
ylabel('RMSE','FontSize',16)
title(['RMSE & SNR vs multiscale-decompositions, k_{min}=',num2str(mink-1)],'FontSize',16)
yyaxis right
plot(0:numScales-1, snr)
ylabel('SNR','FontSize',16)

subplot(1,2,2)
semilogy(0:numScales-1, stopCrit,0:numScales-1,ones(numScales,1))
xlabel('Multiscales: k','FontSize',16)
if tightFlag(1)==1
    title(['(D(F_{data},Tx_k)+\alpha_kTV(log(x_k)))/D(F_{data},Tu), k^*=',num2str(k_star-1)],'FontSize',16)    
else
    title(['D(F_{data},Tx_k)/D(F_{data},Tu), k^*=',num2str(k_star-1)],'FontSize',16)
end

if saveFlag==1
    figName=filePrefix+figPrefix+"metrics.fig";
    savefig(figName)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              \
%                       GET REGIONAL STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram for capillary density per region
% myDataProcessed = i,j,k,q,r= (groups, patient, CBV Scan, region, [mean, std, numVoxels, mode, GOF])
clearvars myDataProcessed ModeVals myDataHist myDataRaw

resolution=0.005;   % resolution of what??
filtsize=3;
NumOfProblemRegions=0;

mkdir([resultPathName DASH 'PeakFit_figs']);
mkdir([resultPathName DASH 'PeakFit_figsNew']);

fig = figure;
set(fig, 'position', [100, 100, 2048, 1440]);
myDataProcessed=NaN(NumGroups,max([ MyData(:).GroupNumPatients ]),max([MyData(:).NumUTEscans]),RegionNums,5);
for i=3:NumGroups % Groups
    fprintf('Group %d, ',i)
    load([resultPathName DASH GroupFolders(i).name '_myDataRaw'])
    names = arrayfun(@(Patient) {[Patient.Subject, Patient.Day1]}, MyData(i).PatientDesignation);
    myDataHist=cell(MyData(i).GroupNumPatients,max([MyData(i).NumUTEscans]),RegionNums); % [j, k, q]
    for j=1:MyData(i).GroupNumPatients
        fprintf('Patient %s, Image...\n',names{j});

        kq_myDataProcessed=zeros(MyData(i).NumUTEscans(j),RegionNums,5);
        kq_myDataHist=cell(MyData(i).NumUTEscans(j),RegionNums);
        for k=PostconStart:MyData(i).NumUTEscans(j)    % 1:MyData(i).NumUTEscans-2 % Notice i put ute-2 since 4 postcons
        % parfor
            fprintf(' %2d, region:   ',k);
            if isempty(MyData(i).name.NII{j,k}) == 0
                q_myDataProcessed=zeros(RegionNums,5);
                q_myDataHist=cell(1,RegionNums);
                for q=1:RegionNums
                    % myDataProcessed
                    fprintf('\b\b\b%3d',q);
                    % myDataProcessed: mean,mode,size
                    raw_temp=cell2mat(myDataRaw(j,k,q));

                    if ~isempty(raw_temp)
                    q_myDataProcessed(q,1)=mean(raw_temp);
                    q_myDataProcessed(q,2)=std(raw_temp);
                    [~,q_myDataProcessed(q,3)]=size(raw_temp);

                    maxVal=max(raw_temp);
                    minVal=min(raw_temp);
                    nbins=round((maxVal-minVal)/resolution); %200;
                    [histCounts,histBins]=hist(reshape(raw_temp,[],1),nbins);

                    % myDataProcessed: ModeVals
                    %                 clearvars X temp_hist
                    q_myDataHist(q)={[histCounts(:) histBins(:)]};
                    temp_hist=cell2mat(q_myDataHist(q));
                    %                     if isempty(temp_hist)
                    %                         NumOfProblemRegions=1+NumOfProblemRegions;
                    %                         ProblemRegions=[NumOfProblemRegions, q];
                    %                     else
                    X=[];
                    X(:,1)=temp_hist(:,2);
                    X(:,2)=temp_hist(:,1);

                    [M,I]=max(X(:,2));
                    signal=X;
                    center=X(I,1);
                    window=2.0;
                    NumPeaks=1;
                    peakshape=0;
                    extra=0;
                    NumTrials=5;

                    [FitResults, GOF,~,~,~,xi,yi]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials);
                    %                     FitInfo(j,k,q)={[FitResults]}; % COULD SAVE
                    %                     RSqrOfModeGaussian(j,k,q)=GOF(1,2); % COULD SAVE
                    pause(1)
                    eval(['print -depsc ' resultPathName DASH 'PeakFit_figs' DASH 'fit_' names{j} sprintf('-Image%02d-reg%03d.eps',k,q) ]);

                    [M,I]=max(yi);
                    q_myDataProcessed(q,4)=xi(I); % Mode finally!

                    % Now we quantify the gaussian-ness of the whole
                    % distribution
                    center=0;
                    window=0;
                    NumPeaks=1;
                    peakshape=0;
                    extra=0;
                    NumTrials=5;
                    [~,GOF,~,~,~,~,~]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials);
                    q_myDataProcessed(q,5)=GOF(1,2); % COULD SAVE


% % new fit, to hijack the old fit
%                     dists = { ...
%                         'Beta', 'Burr', 'Gamma', ...
%                         'gev', 'Nakagami', 'Rayleigh', ...
%                         'Rician', 'Weibull', 'Kernel'};
% 
%                     m = 3; n = 3;
%                     s = raw_temp; % / max(raw_temp(:));
%                     % nbins = 300;
%                     count = 1;
%                     sgtitle('Histogram fit to various distributions. ');
%                     for im = 1 : m
%                         for in = 1 : n
%                             if count > length(dists)
%                                 continue;
%                             end
%                             % fprintf('\b\b\b%3d', count);
%                             subplot(m, n, count);
%                             try
%                                 if im == m && in == n
%                                     disp('Kernel fit;');
%                                 end
%                                 %[histCounts,histBins] = histogram(s, nbins);
%                                 h = histfit(s, 300, dists{count});
%                                 % from handle (h), we have x1,y1 being the histogram; x2,y2 is
%                                 % the fitting.
%                                 y21 = interp1(h(1).XData, h(1).YData, h(2).XData, "linear", "extrap");  % signal
%                                 err = h(2).YData - y21;
%                                 hold on; plot(h(2).XData, 1.2 * max(h(1).YData) + err); hold off; % plot error
% 
%                                 SSE = sum(err.^2);  % Sum of Squared Errors
%                                 MSE = mean(err.^2);  % Mean Squared Error
%                                 RMSE = sqrt(MSE);  % Root Mean Squared Error
%                                 MAE = mean(abs(err));  % Mean Absolute Error
% 
%                                 % Compute the total sum of squares (SST), which is the variance of the signal
%                                 SST = sum((y21 - mean(y21)).^2);
% 
%                                 % Compute R-squared (Coefficient of Determination)
%                                 R_squared = 1 - SSE/SST;
% 
%                                 % p-value
%                                 mdl = fitlm(y21, h(2).YData);
%                                 pv = mdl.ModelFitVsNullModel.Pvalue;
%                                 title(sprintf('%s fitting with \\sigma^2 = %.2f, R^2 = %.2f, p-value = %.3f', dists{count}, RMSE*RMSE, R_squared, pv));
% 
%                             catch exeception
%                                 fprintf('Error encountered: %s\n', exeception.message);
%                             end
% 
%                             if count == 9 % select the fit;
%                                 maxVal = max(y21);
%                                 sigma2 =   RMSE*RMSE;
%                                 R2 = R_squared;
%                             end
% 
%                             count = count + 1;
%                         end
%                     end
%                     q_myDataProcessed(q,4) = maxVal;   % Mode
%                     q_myDataProcessed(q,5) = R2;       % R square
%                     pause(1)
%                     eval(['print -depsc ' resultPathName DASH 'PeakFit_figsNew' DASH 'fit_' names{j} sprintf('-Image%02d-reg%03d.eps',k,q) ]);
% 
% % new fit end                    
                    end
                end
                kq_myDataProcessed(k,:,:)=q_myDataProcessed(:,:);
                kq_myDataHist(k,:)=q_myDataHist(:);
            end
            fprintf('\n');
        end
        myDataProcessed(i,j,1:MyData(i).NumUTEscans(j),:,:)=kq_myDataProcessed;
        myDataHist(j,:,:)=kq_myDataHist(:,:);
    end
    save([resultPathName DASH GroupFolders(i).name '_myDataHist'],'myDataHist')
end
save([resultPathName DASH 'myDataProcessed'],'myDataProcessed')
COMPLETED_MyDataProcessed=1
toc % 14479 secodns or about 4 hours
%
%PLOT CBV WHOLE BRAIN
myDataProcessed(myDataProcessed == 0) = NaN



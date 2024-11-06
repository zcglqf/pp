%FIXED AUTOMATED BLOODMAP PROCESS3



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Imaginostics QUTE-CE MRI Image Processing
% This is the code used on cohort 2
% Current Version: 20/09/2019, Imaginostics
% Author: Codi Gharagouzloo, PhD
%
% All steps after pre-processing (no B1+ here)
% Preprocessing = Reconstruction, B1- Correction, Motion Correction
%
% Version: 1.0
% Application: ApoE4 HFD Rats (n=4 WT, n=5 ApoE4 HFD)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% SETUP AND USE INFOMRATION:
% User selects folder containing groups.
% The group name should be the folder name. Inside each group folder there
% should be three sub-folders: NORM_UTE | ATLAS | BLOODMAP
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clc
close all



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % %   INFORMATION      %  % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Select root folder, where Subject folders e.g. F12, M31, etc, are stored
pathname = uigetdir('/media/czhao/WD_BLACK/Imaginostics/');
if pathname == 0
    disp('Folder selection canceled.');
    return;
else
    disp(['Folder selected: ', pathname]);
end

% Find database file: a file that records all the necessary data for each
% study (a set of pre- and post- contrast scan data, Blood map, Atlas map, etc)
scanDataBase = fullfile(pathname, 'RatScans.xlsx');
if ~exist(scanDataBase, "file")
    disp('Scan Database file not found!');
    return;
end

atlasRegionLabels = fullfile(pathname, 'ATLAS_MODS.xlsx');
if ~exist(atlasRegionLabels, "file")
    disp('Atlas file not found!');
    return;
end

scanTable = readtable(scanDataBase, 'Sheet', 1);
labelTable = readtable(atlasRegionLabels);

% use {'Var1', 'Var7', 'Var8'} for small table
% use {'Var1', 'Var5', 'Var6'} for medium table
% use {'Var1', 'Var3', 'Var4'} for large table
% As such, the first column contains the integer label of each atlas region
% (totally 174 if using NEU atlas); the second column has the integer
% labels of the clustered atalas, and the third column stores the text name of each region.  
labelTable = labelTable(:, {'Var1', 'Var7', 'Var8'});

segbins = unique(labelTable.Var7(~isnan(labelTable.Var7))); % Sum(~isnan(labelTable{:, 2}));
segbins = [segbins; segbins + length(segbins)];             % Left right split

RegionNums = length(segbins)+1;


PostconStart=2;
Co2Start=NaN;




% CREATE RESULTS PATH
resultPathName = fullfile(pathname, 'RESULTS');
% GET NUMBER OF GROUPS
NumGroups = length(GroupFolders); % Number of Groups



%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % %   PROCESSING      %  % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%
% CREATE DIRECTORY FOR SAVING RESULTS
mkdir(resultPathName)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: CALCULATE BLOOD ROI VALUES, assuming one blood map per set of post UTE Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for rowIdx = 1 : height(scanTable)
    if ~scanTable.Process(rowIdx)
        fprintf('Skip row %d\n', rowIdx);
        continue;
    end
    
    subjectName = scanTable.Subject{rowIdx};
    subjectFolder = fullfile(pathname, subjectName);
    if ~exist(subjectFolder, "dir")
        fprintf('Subject folder %s does not exist', subjectFolder);
    end
    dateLabel = scanTable.DateLabel{rowIdx};
    studyNo = scanTable.StudyNo(rowIdx);
    fprintf('Start to process dynmaic scan data for subject %s, %s, study no %d', subjectName, dateLabel, studyNo);

    % create results folder for current subject, current study
    resultPathName = fullfile(subjectFolder, 'RESULTS', [dateLabel, '_Study', num2str(studyNo)]);
    if exist(resultPathName, "dir")
        prompt = sprintf('Result folder %s already exist, overwrite or skip?', resultPathName);
        choice = questdlg(prompt, 'Confirmation', 'Overwrite', 'Skip', 'Skip');

        switch choice
            case 'Overwrite'
                fprintf('User chooses to overwrite folder %s!\n', resultPathName);
            case 'Skip'
                fprintf('User chooses to skip processing current study %s, %s, %d\n', subjectName, dateLabel, studyNo);
                continue;
            otherwise
                fprintf('Folder %s will be overwritten.\n', resultPathName);
        end
    end
    mkdir(resultPathName);

    % data folder path: e.g. .../F12/d0
    dataFolder = fullfile(subjectFolder, dateLabel);
    if ~exist(dataFolder, "dir")
        fprintf('Data folder %s does not exist', dataFolder);
    end

    % Get all nii files from under the data folder
    allFiles = dir(fullfile(dataFolder, '**', '*'));
    allFiles = allFiles(~[allFiles.isdir]);
    niiFiles = allFiles(endsWith({allFiles.name}, '.nii') | endsWith({allFiles.name}, '.nii.gz'));

    % Get pre-contrast file
    preconFileName = scanTable.Pre_contrast{rowIdx};
    preconFile = niiFiles(strcmp({niiFiles.name}, preconFileName));
    if isempty(preconFile)
        fprintf('Pre-contrast file %s not found.', preconFileName);
        continue;
    end
    if numel(preconFile) > 1
        fprintf('More than one pre-contrast file found, check consistency.');
        continue;
    end
    fprintf('Found pre-contrast file: %s.\n', preconFile.name);
    
    % Get Post-contrast files
    postconFileRegExp = [scanTable.Post_contrast{rowIdx}, '_\d{2}\.(?:nii|nii\.gz)']; % (?:...) means no capture.
    postconFiles = niiFiles(~cellfun('isempty', regexp({niiFiles.name}, postconFileRegExp)));
    fprintf('Found post-contrast files:\n%s.\n', strjoin({postconFiles.name}, newline));

    % Get blood map
    bloodMapFileName = scanTable.BloodMap{rowIdx};
    bloodMapFile = niiFiles(strcmp({niiFiles.name}, bloodMapFileName));
    if isempty(bloodMapFile)
        fprintf('Blood map file %s not found.', bloodMapFile);
        continue;
    end
    if numel(bloodMapFile) > 1
        fprintf('More than one blood map file found, check consistency.');
        continue;
    end    

     % Get Atlas map
     bloodMapFileName = scanTable.BloodMap{rowIdx};
     bloodMapFile = niiFiles(strcmp({niiFiles.name}, bloodMapFileName));
     if isempty(bloodMapFile)
        fprintf('Blood map file %s not found.', bloodMapFile);
        continue;
    end
    if numel(bloodMapFile) > 1
        fprintf('More than one blood map file found, check consistency.');
        continue;
    end    

end
mkdir(fullfile(resultPathName, 'BLOODMAPS'));
IB_CALC = NaN(NumGroups, max([ MyData.GroupNumPatients ]), max([MyData.NumUTEscans])); % (NumGroups, Num Animals / group, numUTEscans)

v_IB_CALC = NaN(NumGroups, max([ MyData.GroupNumPatients ]), max([MyData.NumUTEscans])- 1);

for i = 1 : NumGroups % i GROUPS
    for j = 1 : MyData(i).GroupNumPatients % j PATIENTS

        TempBloodNII = special_load_nii(MyData(i).name.Bloodmap{j});
        TempBloodNII.img(TempBloodNII.img > 1) = 1;
        mask_all = zeros(size(TempBloodNII.img));

        for k = 1 : PostconStart-1
            MyData(i).ScanTime{j,k}=NaN; % Relative Scan Time
        end
        MyData(i).ScanTime{j, PostconStart}=0; % Relative Scan Time

        for k = PostconStart : MyData(i).NumUTEscans(j) % for postcon UTE
            % fprintf('k: %d ',k);

            TempUTENII = special_load_nii(MyData(i).name.NII{j,k});

            % myData.ScanTime  & myData.absScanTime
            descrip = getfield(TempUTENII, 'hdr', 'hist', 'descrip');
            tmp = findstr( descrip, ' :: ' );
            MyData(i).absScanTime{j,k} = datetime( descrip(1:tmp) ); % Abs Scan Time
            MyData(i).absScanTime{j,k}.Format = 's';  % set the date-time object to output seconds by default
            if k >= PostconStart+1
                MyData(i).ScanTime{j,k}= seconds( MyData(i).absScanTime{j,k} - MyData(i).absScanTime{j,PostconStart} ); % Relative Scan Time [sec] to first postcon
                % myData(i).ScanTime{j,k}=etime( datevec(datenum(myData(i).absScanTime{j,k})), datevec(datenum(myData(i).absScanTime{j,PostconStart})) ); % Relative Scan Time [sec] to first postcon
            end

            % Rough ROI Mask
            TempCalculation = TempBloodNII.img .* TempUTENII.img; 

            % Create Mask (one voxel / slice)
            mask = zeros(size(TempCalculation)); % Starts as zeros
            for gg = 1:size(TempCalculation,2)
                tempy2 = TempCalculation(:,gg,:); % The Multiplication
                tempy1 = zeros(size(TempCalculation(:,gg,:)));
                if max(max((tempy2)))~=0
                    tempy1(tempy2==max(max((tempy2))))=1; % Max of Multiplication =1
                end
                mask(:,gg,:)=tempy1;
            end
            mask_all = mask_all+mask; % Keep track of voxels used

            maskedBloodNII = mask .* TempUTENII.img; % Values along Z for IB

            % Save Mask
            nii = make_nii(mask);
            save_nii(nii, fullfile(resultPathName, 'BLOODMAPS', [MyData(i).PatientDesignation(j).Subject, MyData(i).PatientDesignation(j).Day1, '_BloodMap', sprintf('%d',k) '.nii']));

            % IB
            IB(i,j,k,1) = mean(nonzeros(maskedBloodNII)); % The mean
            IB(i,j,k,2) = std(nonzeros(maskedBloodNII)); % The STD
            IB(i,j,k,3) = length(nonzeros(maskedBloodNII)); % The Number of Voxels

        end

        % Save Mask from all voxels used
        mask_all(mask_all>1)=1;
        nii = make_nii(mask_all);
        save_nii(nii, fullfile(resultPathName, 'BLOODMAPS', [MyData(i).PatientDesignation(j).Subject, MyData(i).PatientDesignation(j).Day1, '_BloodMap_all.nii']));

        % Get Precon Blood values
        for k=1:PostconStart-1
            TempUTENII=special_load_nii(MyData(i).name.NII{j,k});
            maskedBloodNII = mask_all.*TempUTENII.img; % Values along Z for IB precon
            IB_precon(i,j,k,1) = mean(nonzeros(maskedBloodNII));
            IB_precon(i,j,k,2) = std(nonzeros(maskedBloodNII));
            [IB_precon(i,j,k,3),~] = size(nonzeros(maskedBloodNII));

        end

    end

end

% % Save scantmimes
% scantimes = MyData.ScanTime
% for x=1:size(scantimes,1)
%     for y=1:size(scantimes,2)
%         if isempty(scantimes{x,y});
%             scantimes{x,y}=NaN
%         end
%         %          if isnan(scantimes{x,y});
%         %             scantimes{x,y}=0
%         %         end
%     end
% end
% % scantimes=cell2mat(scantimes)
% save([resultPathName DASH 'scantimes.mat'],'scantimes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   PLOT IB PRECON VALUES %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(fullfile(resultPathName, 'PRECON_GRAPHS'));
%
% fontSize = 16
for i=1:NumGroups

    for j=1:MyData(i).GroupNumPatients % j PATIENTS

        figure

        for k=1:PostconStart-1
            IB_CALC(i,j,k)=IB_precon(i,j,k,1);  % Assign Precon Val in IB_Calc

            errorbar(IB_precon(i,j,k,1),IB_precon(i,j,k,2),'-o')
            title('IB PRECON VALUES')
            names = {[MyData(i).PatientDesignation(j).Subject, MyData(i).PatientDesignation(j).Day1]};
            set(gca,'xtick', 1 : MyData(i).GroupNumPatients, 'xticklabel', names)
            hold on
        end
        saveas(gcf, fullfile(resultPathName, 'PRECON_GRAPHS', ['group' sprintf('%d',i) '_subject' sprintf('%d',j) '.png'])); % Save Figures

    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE BLOOD ROI VALUES  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(fullfile(resultPathName, 'IB_FITS'));
% This the Time Axis for Postcon UTE Scans
% i is the group, j is the animal and k is the postcon scan
close all
fontSize=16;
for i=1:NumGroups
    for j=1:MyData(i).GroupNumPatients
        figure
        x=cell2mat(MyData(i).ScanTime(j,PostconStart:end)); % All postcon scan times measured from first
        y=nonzeros(squeeze(IB(i,j,PostconStart:end,1))) ;% Y AXIS for postcons
        error=nonzeros(squeeze(IB(i,j,PostconStart:end,2))) ;% Y AXIS for postcons
        x=x';
        format long
        errorbar(x,y,error,'o')

        p = polyfit(x, y, 1);
        slope=p(1);
        intercept=p(2);
        yfit=p(1)*x+p(2);
        hold on;
        plot(x,yfit)   %plots trendline with calculated IB

        % ylim([0 2600]) %sets Y axis range
        % set(gca,'FontSize',fontSize)  %sets font size
        pos = get(gca,'position');
        set(gca,'position',pos + [ 0 +0.1 0 -0.1 ]);

        xlabel('Scan Time (s)','fontsize',fontSize)
        ylabel('Intensity','fontsize',fontSize)
        title(['Linear Regression for ' MyData(i).PatientDesignation(j).Subject, MyData(i).PatientDesignation(j).Day1],'fontsize',fontSize)
        grid on
        
        Rsq=1-sum((y-yfit).^2)/sum((y-mean(y)).^2); %finds r^2 value; %uses yfit, IS THIS also not needed now?

        % text(3000,2200,sprintf('Rsq=%f',Rsq),'FontSize', 14) %add Rsq value to chart

        % place equation on graph
        xl = xlim;
        yl = ylim;
        xt = 0.05 * (xl(2)-xl(1)) + xl(1);
        yt = 0.90 * (yl(2)-yl(1)) + yl(1);
        caption = sprintf('y = %4.2e * x + %4.2e', p(1), p(2));
        % text(xt, yt, [ caption '  ' sprintf('Rsq=%f',Rsq) ], 'FontSize', fontSize)

        subplot('position',[ pos(1) pos(2)-0.1 pos(3)  0.1])
        text(0,0.5,[caption '  ' sprintf('Rsq=%f',Rsq) ])
        axis off;
        pause(0.5);
        
        eval([ 'print -depsc ' resultPathName DASH 'IB_FITS' DASH MyData(i).PatientDesignation(j).Subject MyData(i).PatientDesignation(j).Day1  '.eps']);
        
        saveas(gcf, fullfile(resultPathName, ['IB_FITS', MyData(i).PatientDesignation(j).Subject, MyData(i).PatientDesignation(j).Day1, '.fig'])); % Save Figures with Linear Regression per animal for postcon

        v_IB_CALC(i,j,1:size(yfit,1))=yfit;
        IB_CALC(i,j,PostconStart:size(yfit,1)+PostconStart-1)=v_IB_CALC(i,j,1:size(yfit,1)); %CHANGE, added i
    end

end

close all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     SAVE IB CALCULTED MATRIX INCLUDING FIT     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([resultPathName DASH 'IB.mat'],'IB')
save([resultPathName DASH 'IB_CALC.mat'],'IB_CALC')
save([resultPathName DASH 'IB_precon.mat'],'IB_precon')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW PART JUST TO SEE IF WE USE AVERAGE OF POSTCONS WHAT THE RESULTS LOOK LIKE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:NumGroups
%     for j=1:myData(i).GroupNumPatients
%        avg(i,j,:)=mean(IB_CALC(i,j,PostconStart:end))
%        for k=postConStart:end
%         IB_CALC_avg(i,j,k)=avg(i,j)
%        end
%
%     end
% end
% save([resultPathName DASH 'AVERAGE_IB_POSTCON.mat'],'IB_CALC_avg')

%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CALCULATE CBV VALUES           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CBV OF REGULAR POSTCON
for i=1:NumGroups
    % Get all patient names from current group
    names = arrayfun(@(Patient) {[Patient.Subject, Patient.Day1]}, MyData(i).PatientDesignation);
    fprintf('Group: %d\n',i)
    if ~isnan(Co2Start)
        for j=1:MyData(i).GroupNumPatients
            fprintf('Patient: %s ',names{j});
            TempUTENII_pre=special_load_nii(MyData(i).name.NII{j,PostconStart-1}); %scan 2 first normal precon, scan 1 was co2 precon
            TempUTENII_pre.img(TempUTENII_pre.img==0)=NaN;
            for k=PostconStart:1:Co2Start-1 %THIS IS THE FIRST POSTCON IMG WITHOUT CO2
                fprintf('Image: ');
                if isempty(MyData(i).name.NII{j,k}) == 0
                    fprintf(' %d ',k);
                    TempUTENII=special_load_nii(MyData(i).name.NII{j,k});
                    TempUTENII.img(TempUTENII.img==0)=NaN;
                    %             v_CBV(j,:,:,:)=(TempUTENII.img-TempUTENII_pre.img)./(IB_CALC(i,j,k)-IB_CALC(i,j,2));
                    CBV(:,:,:,i,j,k)=(TempUTENII.img-TempUTENII_pre.img)./(IB_CALC(i,j,k)-IB_CALC(i,j,PostconStart-1));
                end
            end
            fprintf('\n');
        end
    end
    

    if isnan(Co2Start)
        for j=1:MyData(i).GroupNumPatients
            fprintf('Patient: %s ',names{j});
            TempUTENII_pre=special_load_nii(MyData(i).name.NII{j,PostconStart-1}); %scan 2 first normal precon, scan 1 was co2 precon
            TempUTENII_pre.img(TempUTENII_pre.img==0)=NaN;
            fprintf('Image: %d ');
            for k=PostconStart:1:MyData(i).NumUTEscans(j) %THIS IS THE FIRST POSTCON IMG WITHOUT CO2
                fprintf(' %d ',k);
                if isempty(MyData(i).name.NII{j,k}) == 0
                    TempUTENII=special_load_nii(MyData(i).name.NII{j,k});
                    TempUTENII.img(TempUTENII.img==0)=NaN;
                    %             v_CBV(j,:,:,:)=(TempUTENII.img-TempUTENII_pre.img)./(IB_CALC(i,j,k)-IB_CALC(i,j,2));
                    CBV(:,:,:,i,j,k)=(TempUTENII.img-TempUTENII_pre.img)./(IB_CALC(i,j,k)-IB_CALC(i,j,PostconStart-1));
                end
            end
            fprintf('\n');
        end
    end


end

%CBV OF CO2 POSTCON

for i=1:NumGroups
    if ~isnan(Co2Start)
        FB1 = CBV;                                                     %size is (:,:,:,i,j,3 & 5) %is this still k-1?  Because the uploaded form might just be k lowkey
        FT1= 1-FB1;                                                    %size is (:,:,:,2,6,3 & 5)
        IB1 = IB_CALC;                                                 %specify post con scans (:,:,k); where k=3,5
        IB2 = IB_CALC;                                                 %specify co2 scans(:,:,k+1);    where k=4,6
        IB0 = IB_CALC;                                                 %specify precon scan(:,;,1); 2
        for j=1:MyData(i).GroupNumPatients
            fprintf('Patient %d, ',j);
            IM0=special_load_nii(MyData(i).name.NII{j,PostconStart-1}); IM0.img(IM0.img==0)=NaN;    %LOAD PRECON, SCAN 2 IS PRECON HERE
            IM1=special_load_nii(MyData(i).name.NII{j,PostconStart}); IM1.img(IM1.img==0)=NaN;    %SCANS 3 is first postcon

            IT = (IM1.img-(FB1(:,:,:,i,j,PostconStart).*IB1(i,j,PostconStart)))./(FT1(:,:,:,i,j,PostconStart)); %(180 180 180) Intensity of tissue of one post con scan
            for k=Co2Start:MyData(i).NumUTEscans % Should be same for all groups, but ...
                k
                if isempty(MyData(i).name.NII{j,k}) == 0
                    IM2 = special_load_nii(MyData(i).name.NII{j,k});IM2.img(IM2.img==0)=NaN;              %SCANS 4,6 HAVE POST CON CO2 CHALLENGEZ

                    CBV(:,:,:,i,j,k)=(IM2.img-IM1.img+(FB1(:,:,:,i,j,PostconStart).*(IB1(i,j,PostconStart)-IT)))./(IB2(i,j,k)-IT); %EQUATION=CBV=(IM2.img-IM1.img+(FB1.*(IB1-IT)))./(IB2-IT);
                end
            end
        end
    end
end

%
save([resultPathName DASH 'CBV.mat'],'CBV','-v7.3')
% CBV=single(CBV);
COMPLETED_CBV=1
toc % 102 seconds
%%
%
tic
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        GET REGIONAL VOXELS          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
erorrs_accum=ones(NumGroups,MyData(i).GroupNumPatients,RegionNums); % check for region outside FoV: [Group, Patient, Region]
%
for i=1:NumGroups %Groups
    names = arrayfun(@(Patient) {[Patient.Subject, Patient.Day1]}, MyData(i).PatientDesignation);
    fprintf('Group %d: \n',i);
    myDataRaw=cell(MyData(i).GroupNumPatients,max(MyData(i).NumUTEscans(:)),RegionNums); % Patients, CBV scan, Region
    for j=1:MyData(i).GroupNumPatients %Animals
    % parfor
        fprintf('%s \n',names{j}); % Patient: 
        currentAtlas=special_load_nii(MyData(i).name.ATLAS{j});
        [~, loc] = ismember(currentAtlas.img, [0; labelTable.Var1; labelTable.Var1+174]); % returns location index to Var1
        newMapTable = [0; labelTable.Var7; labelTable.Var7+length(segbins)];   % map first column to 7-th column: [0; labelTable.Var7];          
        currentAtlas.img = newMapTable(loc);  % mappedAtlas
        kq_myDataRaw=cell(MyData(i).NumUTEscans(j),RegionNums);
        for k=PostconStart:MyData(i).NumUTEscans(j) % post con scans
            fprintf('Image %d ',k);
            if isempty(MyData(i).name.NII{j,k}) == 0
                q_myDataRaw=cell(1,RegionNums);
                fprintf('atlas region:    ');
                for q=1: RegionNums   % (length(segbins)+1) 
                    % fprintf('\b\b\b%3d',q);
                    ATLAS= imresize3( currentAtlas.img, size(CBV,1) / size(currentAtlas.img,1), 'nearest' );
                    if q==(length(segbins)+1) %  RegionNums
                        ATLAS(ATLAS>1)=1;    %grab whole brain for last "region"
                        ATLAS(ATLAS~=1)=NaN;
                    else
                        ATLAS( ATLAS ~= segbins(q) )=NaN;
                        ATLAS( ATLAS == segbins(q) )=1;
                    end
                    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3d, voxels %3d',q, nnz(ATLAS ~= 0 & ~isnan(ATLAS)));
                    % CBVreg
                    CBVreg=squeeze(CBV(:,:,:,i,j,k)).*ATLAS;

                    % Compare #voxels in brain region vs. atlas
                    if erorrs_accum(i,j,q)==1
                        IMAGE=ATLAS.*CBVreg;
                        [division1,~]=size(~isnan(ATLAS));
                        [division2,~]=size(~isnan(IMAGE));
                        division=division2/division1;
                        if division<.98
                            erorrs_accum(i,j,q)=NaN;
                        end
                    end

                    % CBVreg save myDataRaw
                    CBVreg=reshape(CBVreg,1,[]);
                    CBVreg(isnan(CBVreg))=[];
                    q_myDataRaw(q)={CBVreg};
                end
                kq_myDataRaw(k,:)=q_myDataRaw(:);
            end
            fprintf('\n');
        end
        myDataRaw(j,:,:)=kq_myDataRaw(:,:);
    end
    fprintf('\n');
    save([resultPathName DASH GroupFolders(i).name '_myDataRaw'],'myDataRaw')
end

save([resultPathName DASH 'erorrs_accum'],'erorrs_accum')
COMPLETED_RegionalVoxels=1
toc % 545 seconds (~9min)
MyData
%%
tic
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
                    window=0.5;
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

% new fit end                    
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






%%
% -------------------------------------------------------------------------
% - Sub function to return files filtered.                                -
% -------------------------------------------------------------------------
% input: 
% filesIn, file struct array return by command dir
% filter, regular express to which the file names should match
% output:
% filesOut, filtered file struct array
function filesOut = filterFilesByName(filesIn, filter)

if ~all(isfield(filesIn, {'name', 'folder', 'date', 'bytes', 'isdir', 'datenum'}))
    error('Input file list should be a string array returned by command dir.');
end

if ~(ischar(filter) || isstring(filter))
    error('Regular expression must be a string.');
end

names = {filesIn.name};
matches = regexp(names, filter, 'match');
matchingIndices = ~cellfun('isempty', matches); 
filesOut = filesIn(matchingIndices);

end

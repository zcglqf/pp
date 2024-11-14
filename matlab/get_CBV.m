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
scanTable = scanTable(scanTable.StudyNo == 2, :);

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


fontSize = 16;

% CREATE RESULTS PATH
prefixResultFolder = fullfile(pathname, 'RESULTS');
% GET NUMBER OF GROUPS




%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % %   PROCESSING      %  % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%
% CREATE DIRECTORY FOR SAVING RESULTS
mkdir(prefixResultFolder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: CALCULATE BLOOD ROI VALUES, assuming one blood map per set of post UTE Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxNumDynamics = 100;
templRegion = struct('mean', 0, 'std', 0, 'mode', 0, 'length', 0);

templSubject = struct( ...
    'Subject', '', ...
    'Date', '', ...
    'PrefixArchive', '', ...
    'StudyID', '', ...
    'PostConDynamics', 0, ...
    'IB_precon', struct( ...                    % pre-contrast IB
         'mean', 0, ...
         'std', 0, ...
         'length', 0), ...
    'IB', struct( ...
         'mean', zeros(maxNumDynamics,1), ...   % post-contrast IB
         'std', zeros(maxNumDynamics,1), ...
         'length', zeros(maxNumDynamics,1)), ...
    'IB_FIT', zeros(maxNumDynamics,1), ...
    'TimeStamp', repmat(datetime(0,1,1),maxNumDynamics, 1), ... % time stamp of the post-contrast scans.
    'ScanTime', zeros(maxNumDynamics,1), ...                    % relevant time of the post-contrast scans in s.
    'PreconTimeStamp', datetime(0,1,1), ...
    'PreconScanTime', 0, ...
    'Regions', repmat(templRegion, RegionNums, 1));

results = repmat(templSubject, height(scanTable), 1);

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
    fprintf('Start to process dynmaic scan data for subject %s, %s, study no %d\n', subjectName, dateLabel, studyNo);

    % create results folder for current subject, current study
    % resultPathName = fullfile(subjectFolder, 'RESULTS', [dateLabel, '_Study', num2str(studyNo)]);
    % if exist(resultPathName, "dir")
    %     prompt = sprintf('Result folder %s already exist, overwrite or skip?', resultPathName);
    %     choice = questdlg(prompt, 'Confirmation', 'Overwrite', 'Skip', 'Skip');
    % 
    %     switch choice
    %         case 'Overwrite'
    %             fprintf('User chooses to overwrite folder %s!\n', resultPathName);
    %         case 'Skip'
    %             fprintf('User chooses to skip processing current study %s, %s, %d\n', subjectName, dateLabel, studyNo);
    %             continue;
    %         otherwise
    %             fprintf('Folder %s will be overwritten.\n', resultPathName);
    %     end
    % end
    % mkdir(resultPathName);
    % results(rowIdx).PrefixArchive = resultPathName;

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
    preconFile = fullfile(preconFile.folder, preconFile.name);
    fprintf('Found pre-contrast file:\n%s.\n', preconFile);
    
    % Get Post-contrast files
    postconFileRegExp = [scanTable.Post_contrast{rowIdx}, '_\d{2}\.(?:nii|nii\.gz)']; % (?:...) means no capture.
    % postconFiles = niiFiles(~cellfun('isempty', regexp({niiFiles.name}, postconFileRegExp)));
    postconFiles = filterFilesByName(niiFiles, postconFileRegExp);
    postconFiles = arrayfun(@(x) fullfile(x.folder, x.name), postconFiles, 'UniformOutput', false);
    fprintf('Found post-contrast files:\n%s.\n', strjoin(postconFiles, newline));

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
    bloodMapFile = fullfile(bloodMapFile.folder, bloodMapFile.name);
    fprintf('Found blood map files %s.\n', bloodMapFile);

    % Get Atlas map
    atlasMapFileName = scanTable.AtlasLabelRegistered{rowIdx};
    atlasMapFile = niiFiles(strcmp({niiFiles.name}, atlasMapFileName));
    if isempty(atlasMapFile)
        fprintf('Atlas map file %s not found.', atlasMapFile);
        continue;
    end
    if numel(atlasMapFile) > 1
        fprintf('More than one atlas map file found, check consistency.');
        continue;
    end
    atlasMapFile = fullfile(atlasMapFile.folder, atlasMapFile.name);
    fprintf('Found atlas map files %s.\n', atlasMapFile);

    % Get T2W
    t2wFileName = scanTable.T2W{rowIdx};
    t2wFile = niiFiles(strcmp({niiFiles.name}, t2wFileName));
    t2wFile = fullfile(t2wFile.folder, t2wFile.name);
    assert(exist(t2wFile, 'file'));

    % Visually Check the registration of UTE scan, blood map, atlas map
    fprintf('Execute the following command in a terminal if you would like to check the registration across the data:\n');
    fprintf('~/fsl/bin/fsleyes %s %s %s %s %s\n', t2wFile, atlasMapFile, bloodMapFile, preconFile, postconFiles{1});
    
    % load blood map data
    niiBlood = special_load_nii(bloodMapFile);
    niiBlood.img(niiBlood.img > 1) = 1;

    % load atlas map data
    niiAtlas = special_load_nii(atlasMapFile);

    % load pre-contrast UTE data
    niiPreconUTE = special_load_nii(preconFile);

    % fill the strut holding the study info and result
    results(rowIdx).PrefixArchive = prefixResultFolder;
    results(rowIdx).StudyID = [subjectName, dateLabel, '_Study', num2str(studyNo)]; % File name for resultant plots.
    results(rowIdx).Subject = subjectName;
    results(rowIdx).Date = dateLabel;

    % load post-contrast UTE data
    mask_all = zeros(size(niiBlood.img));
    numDynamics = length(postconFiles);
    results(rowIdx).PostConDynamics = numDynamics;
    for k = 1 : numDynamics
        niiPostconUTE = special_load_nii(postconFiles{k});

        % Get time stamp and derive relative scan time.
        descrip = getfield(niiPostconUTE, 'hdr', 'hist', 'descrip');
        tmp = strfind( descrip, ' :: ' );
        absScanTime = datetime( descrip(1:tmp) ); % Abs Scan Time
        results(rowIdx).TimeStamp(k) = absScanTime;
        if k > 1
            results(rowIdx).ScanTime(k) = seconds(absScanTime - results(rowIdx).TimeStamp(1));
        end

        % Rough ROI Mask
        TempCalculation = niiBlood.img .* niiPostconUTE.img;

        % Create Mask (one voxel / slice)
        mask = zeros(size(TempCalculation));  % Starts as zeros
        for gg = 1 : size(TempCalculation,2)  % slice index
            [maxValue, linearIndex] = max(squeeze(TempCalculation(:,gg,:)), [], 'all');
            [rowIndex, colIndex] = ind2sub(size(squeeze(TempCalculation(:,gg,:))), linearIndex);
            if maxValue ~= 0
                mask(rowIndex,gg,colIndex) = 1;
            end
        end
        mask_all = mask_all + mask; % Keep track of voxels used

        maskedBloodNII = mask .* niiPostconUTE.img; % Values along Z for IB
        
        % IB
        results(rowIdx).IB.mean(k) = mean(nonzeros(maskedBloodNII)); % The mean
        results(rowIdx).IB.std(k) = std(nonzeros(maskedBloodNII)); % The STD
        results(rowIdx).IB.length(k) = length(nonzeros(maskedBloodNII)); % The Number of Voxels

    end

    % Get time stamp of pre-contrast scan
    descrip = getfield(niiPreconUTE, 'hdr', 'hist', 'descrip');
    tmp = strfind( descrip, ' :: ' );
    absScanTime = datetime( descrip(1:tmp) ); % Abs Scan Time
    results(rowIdx).PreconTimeStamp = absScanTime;
    results(rowIdx).PreconScanTime = seconds(absScanTime - results(rowIdx).TimeStamp(1));


    % Get Precon Blood values
    mask_all(mask_all > 1) = 1;
    maskedBloodNII = mask_all .* niiPreconUTE.img; % Values along Z for IB precon
    results(rowIdx).IB_precon.mean = mean(nonzeros(maskedBloodNII));
    results(rowIdx).IB_precon.std =  std(nonzeros(maskedBloodNII));
    results(rowIdx).IB_precon.length = length(nonzeros(maskedBloodNII));

    % plot IB dynamics from pre-contrast to post-contrast scans
    results(rowIdx) = plotIBAndArchievePlots(results(rowIdx));

end 








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

%%
% -------------------------------------------------------------------------
% - Sub function to plot IB_precon and IB and save figs                               -
% -------------------------------------------------------------------------
% input: 
% result, a struct of type 'templSubject'
% output: Plots and saved plots to given folder (IB_FIT)
function output = plotIBAndArchievePlots(input)
    % assert(isequal(fieldnames(input), fieldnames(templSubject)), 'Input must be a struct with same type as templSubject');
    output = input;
    figure;
    fontSize = 16;
    % plot pre-contrast errorbar first
    errorbar(output.PreconScanTime, output.IB_precon.mean, output.IB_precon.std, '-o');
       
            % title('IB PRECON VALUES')
    name = [output.Subject, output.Date];
            % set(gca,'xtick', 1 : MyData(i).GroupNumPatients, 'xticklabel', names)
    hold on
    % saveas(gcf, fullfile(resultPathName, 'PRECON_GRAPHS', ['group' sprintf('%d',i) '_subject' sprintf('%d',j) '.png'])); % Save Figures

    numDynamics = output.PostConDynamics;
    x = output.ScanTime(1:numDynamics); % All postcon scan times measured from first
    y = output.IB.mean(1:numDynamics);  % Y AXIS for postcons
    error = output.IB.std(1:numDynamics);
    format long;
    errorbar(x, y, error, 'o');
    p = polyfit(x, y, 1);
    slope = p(1);
    intercept = p(2);
    yfit = p(1)*x+p(2);
    hold on;
    plot(x, yfit)   %plots trendline with calculated IB

    pos = get(gca,'position');
    set(gca,'position',pos + [ 0 +0.1 0 -0.1 ]);

    xlabel('Scan Time (s)','fontsize',fontSize)
    ylabel('Intensity','fontsize',fontSize)
    title(['Linear Regression for ' output.Subject, output.Date], 'fontsize', fontSize)
    grid on

    Rsq=1-sum((y-yfit).^2)/sum((y-mean(y)).^2); % finds r^2 value; %uses yfit, IS THIS also not needed now?

    % place equation on graph
    xl = xlim;
    yl = ylim;
    xt = 0.05 * (xl(2)-xl(1)) + xl(1);
    yt = 0.90 * (yl(2)-yl(1)) + yl(1);
    caption = sprintf('y = %4.2e * x + %4.2e', p(1), p(2));

    subplot('position',[ pos(1) pos(2)-0.1 pos(3)  0.1])
    text(0,0.5,[caption '  ' sprintf('Rsq=%f',Rsq) ])
    axis off;
    pause(0.5);

    output.IB_FIT(1:numDynamics) = yfit;

    targetFolder = fullfile(output.PrefixArchive, 'IB_FIT');
    mkdir(targetFolder);

    targetFile = fullfile(targetFolder, [output.StudyID '.png']);
    saveas(gcf, targetFile);
end

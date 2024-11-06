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
% Variables that are indexed along with Postcon CBV
% myData(i).ScanTime
%

% myData.AnimalNum=cell(1,AnimalTypes);
% myData.AnimalNum(1,1)={[1 2 5 6 9 10]}; % 6 WT
% myData.AnimalNum(1,2)={[3 4 7 8 11 12]}; % 6 WT+H
% myData.AnimalNum(1,3)={[13 17 18 19 20]}; % 5 APOE
% myData.AnimalNum(1,4)={[15 21 22 23 24]}; % 5 APOE+H



clear
clc
close all

% SET DASHES
% Compatible with PC and MAC/Linux
DASH = filesep;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % %   INFORMATION      %  % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% USER SELECTS FILE WITH GROUPS
% pathname='/scratch/gharagouzloo/imaginostics/TBI/DATA_extracted/PatientData/R10patientdata';
if  ~exist('pathname','var')
    [pathname] = uigetfile_n_dir();
end
pathname=char(pathname);


% GET NUMBER OF REGIONS
f = dir([pathname DASH '*Atlas*.txt']);
fid = fopen([f.folder DASH f.name]);
segbins = [];
while ~feof(fid)
    fln = fgetl(fid);
    if ~strcmp(fln(1),'#');  % skip the comments
        tmp = sscanf( fln, '%d', 1 );
        if ( tmp ~= 0 ), segbins = [segbins; tmp ]; end;
    end
end
fclose(fid);
RegionNums = length(segbins)+1;
% RegionNums=234+1; % SIGMA Anatomical InVivo Atlas

tic
%
% % % % % % % % % % %       NEW   % % % % % % % % % % % % % % % % % %
% prompt = {'Type first post con #','Type first Co2 scan #','Choose Cohort/Animal name','What scans would you like to compare in ttest analysis?'};
% dims = [1 25; 1 25; 1 50; 1 50]; definput = {'2','NaN',' ',''};
% Information = (inputdlg(prompt,'STUDY INFORMATION',dims,definput));

% Comp=Information(4,1); Comp=str2num(Comp{1})
% PostconStart=Information(1,1); PostconStart=str2num(PostconStart{1})
% Co2Start=Information(2,1); Co2Start=str2num(Co2Start{1})
% COHORT=cell2mat(Information(3,1));
PostconStart=2;
Co2Start=NaN;
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%

% GET DIRECTORY INFO
dirData = dir(pathname); dirData=dirData(3:end); % Get the data for the current directory, always two extra folders '.' and '..'
dirIndex = [dirData.isdir];  % Index for directories
[temp,~]=size(dirData);
%
for i=1:temp
    if strcmp(dirData(i).name,'RESULTS') % Exclude Results directory
        dirIndex(i)=0;
    end
end
%
dirData=dirData(dirIndex==1); % Remove non-directory indicies


% CREATE RESULTS PATH
resultPathName=[dirData(1).folder DASH 'RESULTS'];
% GET NUMBER OF GROUPS
[NumGroups,~]=size(dirData); % Number of Groups
%
% STORE GROUP SIZE AND PATIENT DESIGNATIONS
for i=1:NumGroups
    TEMP_dir=dir([dirData(1).folder DASH dirData(i).name DASH 'ATLAS']); % Create dir for group BLOODMAP

    [TEMP_size,~]=size(TEMP_dir);% Remove non-ATLAS indicies
    TEMP_param=zeros(TEMP_size,1);
    for j=1:TEMP_size
        ParamPos=strfind(TEMP_dir(j).name,'ATLAS');
        if ~isempty(ParamPos)
            ParamPos=strfind(TEMP_dir(j).name,'._');
            if isempty(ParamPos)
                TEMP_param(j,1)=1;
            end
        end
    end
    TEMP_dir=TEMP_dir(TEMP_param==1);
    [TEMP_size,~]=size(TEMP_dir);% Find Patient Designations
    myData(i).GroupNumPatients=TEMP_size;
    for j=1:TEMP_size
        ParamPos=strfind(TEMP_dir(j).name,'_ATLAS');
        TEMP_designation=TEMP_dir(j).name;
        TEMP_designation(ParamPos:end)=[];
        myData(i).PatientDesignation{j}=TEMP_designation; % Store Designations
    end

end
%
% STORE NUMBER OF UTE SCANS

for i=1:NumGroups
    myData(i).NumUTEscans=5;
end

%
% STORE NUMBER OF BLOODMAPS
for i=1:NumGroups
    myData(i).NumBloodmaps=1;
end


%
% QUICK NAMES: UTE, BLOODMAP, ATLAS
for i=1:NumGroups % i GROUPS
    for j=1:myData(i).GroupNumPatients % j PATIENTS
        %myData(i).name.BLOODMAP{j}=[dirData(1).folder DASH dirData(i).name DASH 'BLOODMAP' DASH myData(i).PatientDesignation{j} '_BLOODMAP.nii']; % One for Each Scan and Animal
        myData(i).name.ATLAS{j}=[dirData(1).folder DASH dirData(i).name DASH 'ATLAS' DASH myData(i).PatientDesignation{j} '_ATLAS.nii.gz']; % One for Each Group and Subject

        for k=1:myData(i).NumUTEscans % k UTE Scans
            myData(i).name.NII{j,k}=[dirData(1).folder DASH dirData(i).name DASH 'NORM_UTE' DASH myData(i).PatientDesignation{j} '_UTE' sprintf('%d',k) '_NORM.nii']; % Precons and Postcons
        end
        
        for k=1:myData(i).NumBloodmaps % k bloodmaps per animal
            myData(i).name.Bloodmap{j,k}=[dirData(1).folder DASH dirData(i).name DASH 'BLOODMAP' DASH myData(i).PatientDesignation{j} '_BLOODMAP.nii.gz']; % Postcon bloodmaps
        end

    end
end

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
mkdir([resultPathName DASH 'BLOODMAPS'])
IB_CALC=NaN(NumGroups,max([ myData(:).GroupNumPatients ]),myData(i).NumUTEscans); % (NumGroups, Num Animals / group, numUTEscans)

v_IB_CALC=NaN(NumGroups,max([ myData(:).GroupNumPatients ]),myData(i).NumUTEscans-1);

for i=1:NumGroups % i GROUPS
    for j=1:myData(i).GroupNumPatients % j PATIENTS

        TempBloodNII=special_load_nii(myData(i).name.Bloodmap{j});
        TempBloodNII.img(TempBloodNII.img>1)=1;
        mask_all=zeros(size(TempBloodNII.img));

        for k=1:PostconStart-1
            myData(i).ScanTime{j,k}=NaN; % Relative Scan Time
        end
        myData(i).ScanTime{j,PostconStart}=0; % Relative Scan Time

        for k=PostconStart:myData(i).NumUTEscans % for postcon UTE
            % fprintf('k: %d ',k);

            TempUTENII=special_load_nii(myData(i).name.NII{j,k});

            % myData.ScanTime  & myData.absScanTime
            descrip = getfield(TempUTENII,'hdr','hist','descrip');
            tmp = findstr( descrip, ' :: ' );
            myData(i).absScanTime{j,k} = datetime( descrip(1:tmp) ); % Abs Scan Time
            myData(i).absScanTime{j,k}.Format = 's';  % set the date-time object to output seconds by default
            if k>=PostconStart+1
                myData(i).ScanTime{j,k}= seconds( myData(i).absScanTime{j,k} - myData(i).absScanTime{j,PostconStart} ); % Relative Scan Time [sec] to first postcon
                % myData(i).ScanTime{j,k}=etime( datevec(datenum(myData(i).absScanTime{j,k})), datevec(datenum(myData(i).absScanTime{j,PostconStart})) ); % Relative Scan Time [sec] to first postcon
            end

            % Rough ROI Mask
            TempCalculation=TempBloodNII.img.*TempUTENII.img; 

            % Create Mask (one voxel / slice)
            mask=zeros(size(TempCalculation)); % Starts as zeros
            for gg=1:size(TempCalculation,2)
                tempy2=TempCalculation(:,gg,:); % The Multiplication
                tempy1=zeros(size(TempCalculation(:,gg,:)));
                if max(max((tempy2)))~=0
                    tempy1(tempy2==max(max((tempy2))))=1; % Max of Multiplication =1
                end
                mask(:,gg,:)=tempy1;
            end
            mask_all=mask_all+mask; % Keep track of voxels used

            maskedBloodNII=mask.*TempUTENII.img; % Values along Z for IB

            % Save Mask
            nii=make_nii(mask);
            save_nii(nii,[resultPathName DASH 'BLOODMAPS' DASH myData(i).PatientDesignation{j} '_BloodMap' sprintf('%d',k) '.nii'])

            % IB
            IB(i,j,k,1) = mean(nonzeros(maskedBloodNII)); % The mean
            IB(i,j,k,2) = std(nonzeros(maskedBloodNII)); % The STD
            IB(i,j,k,3) = length(nonzeros(maskedBloodNII)); % The Number of Voxels

        end

        % Save Mask from all voxels used
        mask_all(mask_all>1)=1;
        nii=make_nii(mask_all);
        save_nii(nii,[resultPathName DASH 'BLOODMAPS' DASH myData(i).PatientDesignation{j} '_BloodMap_all.nii'])

        % Get Precon Blood values
        for k=1:PostconStart-1
            TempUTENII=special_load_nii(myData(i).name.NII{j,k});
            maskedBloodNII=mask_all.*TempUTENII.img; % Values along Z for IB precon
            IB_precon(i,j,k,1)=mean(nonzeros(maskedBloodNII));
            IB_precon(i,j,k,2)=std(nonzeros(maskedBloodNII));
            [IB_precon(i,j,k,3),~]=size(nonzeros(maskedBloodNII));

        end

    end

end

% Save scantmimes
scantimes = myData.ScanTime
for x=1:size(scantimes,1)
    for y=1:size(scantimes,2)
        if isempty(scantimes{x,y});
            scantimes{x,y}=NaN
        end
        %          if isnan(scantimes{x,y});
        %             scantimes{x,y}=0
        %         end
    end
end
% scantimes=cell2mat(scantimes)
save([resultPathName DASH 'scantimes.mat'],'scantimes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   PLOT IB PRECON VALUES %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir([resultPathName DASH 'PRECON_GRAPHS'])
%
fontSize=16
for i=1:NumGroups

    for j=1:myData(i).GroupNumPatients % j PATIENTS

        figure

        for k=1:PostconStart-1
            IB_CALC(i,j,k)=IB_precon(i,j,k,1);  % Assign Precon Val in IB_Calc

            errorbar(IB_precon(i,j,k,1),IB_precon(i,j,k,2),'-o')
            title('IB PRECON VALUES')
            names=(myData(i).PatientDesignation)'
            set(gca,'xtick',[1:myData(i).GroupNumPatients],'xticklabel',names)
            hold on
        end
        saveas(gcf,[resultPathName DASH 'PRECON_GRAPHS' DASH 'group' sprintf('%d',i) '_subject' sprintf('%d',j) '.png']); % Save Figures

    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE BLOOD ROI VALUES  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir([resultPathName DASH 'IB_FITS'])
% This the Time Axis for Postcon UTE Scans
% i is the group, j is the animal and k is the postcon scan
close all
fontSize=16;
for i=1:NumGroups
    for j=1:myData(i).GroupNumPatients
        figure
        x=cell2mat(myData(i).ScanTime(j,PostconStart:end)); % All postcon scan times measured from first
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
        title(['Linear Regression for ' myData(i).PatientDesignation{j}],'fontsize',fontSize)
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
        
        eval([ 'print -depsc ' resultPathName DASH 'IB_FITS' DASH myData(i).PatientDesignation{j} '.eps']);
        
        saveas(gcf,[resultPathName DASH 'IB_FITS' DASH myData(i).PatientDesignation{j} '.fig']); % Save Figures with Linear Regression per animal for postcon


        v_IB_CALC(i,j,1:size(yfit,1))=yfit;
        IB_CALC(i,j,PostconStart:size(yfit,1)+PostconStart-1)=v_IB_CALC(i,j,1:size(yfit,1));%CHANGE, added i
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
    fprintf('Group: %d\n',i)
    if ~isnan(Co2Start)
        for j=1:myData(i).GroupNumPatients
            fprintf('Patient: %s ',names{j});
            TempUTENII_pre=special_load_nii(myData(i).name.NII{j,PostconStart-1}); %scan 2 first normal precon, scan 1 was co2 precon
            TempUTENII_pre.img(TempUTENII_pre.img==0)=NaN;
            for k=PostconStart:1:Co2Start-1 %THIS IS THE FIRST POSTCON IMG WITHOUT CO2
                fprintf('Image: ');
                if isempty(myData(i).name.NII{j,k}) == 0
                    fprintf(' %d ',k);
                    TempUTENII=special_load_nii(myData(i).name.NII{j,k});
                    TempUTENII.img(TempUTENII.img==0)=NaN;
                    %             v_CBV(j,:,:,:)=(TempUTENII.img-TempUTENII_pre.img)./(IB_CALC(i,j,k)-IB_CALC(i,j,2));
                    CBV(:,:,:,i,j,k)=(TempUTENII.img-TempUTENII_pre.img)./(IB_CALC(i,j,k)-IB_CALC(i,j,PostconStart-1));
                end
            end
            fprintf('\n');
        end
    end
    

    if isnan(Co2Start)
        for j=1:myData(i).GroupNumPatients
            fprintf('Patient: %s ',names{j});
            TempUTENII_pre=special_load_nii(myData(i).name.NII{j,PostconStart-1}); %scan 2 first normal precon, scan 1 was co2 precon
            TempUTENII_pre.img(TempUTENII_pre.img==0)=NaN;
            fprintf('Image: %d ');
            for k=PostconStart:1:myData(1).NumUTEscans %THIS IS THE FIRST POSTCON IMG WITHOUT CO2
                fprintf(' %d ',k);
                if isempty(myData(i).name.NII{j,k}) == 0
                    TempUTENII=special_load_nii(myData(i).name.NII{j,k});
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
        for j=1:myData(i).GroupNumPatients
            fprintf('Patient %d, ',j);
            IM0=special_load_nii(myData(i).name.NII{j,PostconStart-1}); IM0.img(IM0.img==0)=NaN;    %LOAD PRECON, SCAN 2 IS PRECON HERE
            IM1=special_load_nii(myData(i).name.NII{j,PostconStart}); IM1.img(IM1.img==0)=NaN;    %SCANS 3 is first postcon

            IT = (IM1.img-(FB1(:,:,:,i,j,PostconStart).*IB1(i,j,PostconStart)))./(FT1(:,:,:,i,j,PostconStart)); %(180 180 180) Intensity of tissue of one post con scan
            for k=Co2Start:myData(i).NumUTEscans % Should be same for all groups, but ...
                k
                if isempty(myData(i).name.NII{j,k}) == 0
                    IM2 = special_load_nii(myData(i).name.NII{j,k});IM2.img(IM2.img==0)=NaN;              %SCANS 4,6 HAVE POST CON CO2 CHALLENGEZ

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
erorrs_accum=ones(NumGroups,myData(i).GroupNumPatients,RegionNums); % check for region outside FoV: [Group, Patient, Region]
%
for i=1:NumGroups %Groups
    fprintf('Group %d: ',i);
    myDataRaw=cell(myData(i).GroupNumPatients,myData(i).NumUTEscans,RegionNums); % Patients, CBV scan, Region
    for j=1:myData(i).GroupNumPatients %Animals
    % parfor
        fprintf('%s ',names{j}); % Patient: 
        currentAtlas=special_load_nii(myData(i).name.ATLAS{j});
        
        kq_myDataRaw=cell(myData(i).NumUTEscans,RegionNums);
        for k=PostconStart:myData(i).NumUTEscans % post con scans
            fprintf('Image %d ',k);
            if isempty(myData(i).name.NII{j,k}) == 0
                q_myDataRaw=cell(1,RegionNums);
                fprintf('atlas region:\n');
                for q=1:(length(segbins)+1) % RegionNums
                    fprintf('\r%d',q);
                    ATLAS= imresize3( currentAtlas.img, size(CBV,1) / size(currentAtlas.img,1), 'nearest' );
                    if q==(length(segbins)+1) %  RegionNums
                        ATLAS(ATLAS>1)=1;    %grab whole brain for last "region"
                        ATLAS(ATLAS~=1)=NaN;
                    else
                        ATLAS( ATLAS ~= segbins(q) )=NaN;
                        ATLAS( ATLAS == segbins(q) )=1;
                    end

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
    save([resultPathName DASH dirData(i).name '_myDataRaw'],'myDataRaw')
end

save([resultPathName DASH 'erorrs_accum'],'erorrs_accum')
COMPLETED_RegionalVoxels=1
toc % 545 seconds (~9min)

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
         
myDataProcessed=NaN(NumGroups,max([ myData(:).GroupNumPatients ]),myData(1).NumUTEscans,RegionNums,5);
for i=1:NumGroups % Groups
    fprintf('Group %d, ',i)
    load([resultPathName DASH dirData(i).name '_myDataRaw'])

    myDataHist=cell(myData(i).GroupNumPatients,myData(i).NumUTEscans,RegionNums); % [j, k, q]
    for j=1:myData(i).GroupNumPatients
        fprintf('Patient %s, Image...\n',names{j});

        kq_myDataProcessed=zeros(myData(i).NumUTEscans,RegionNums,5);
        kq_myDataHist=cell(myData(i).NumUTEscans,RegionNums);
        for k=PostconStart:myData(i).NumUTEscans    % 1:myData(i).NumUTEscans-2 % Notice i put ute-2 since 4 postcons
        % parfor
            fprintf(' %2d, region:   ',k);
            if isempty(myData(i).name.NII{j,k}) == 0
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
                    nbins=200; % round((maxVal-minVal)/resolution);
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
                    window=0.1;
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
                    end
                end
                kq_myDataProcessed(k,:,:)=q_myDataProcessed(:,:);
                kq_myDataHist(k,:)=q_myDataHist(:);
            end
            fprintf('\n');
        end
        myDataProcessed(i,j,:,:,:)=kq_myDataProcessed;
        myDataHist(j,:,:)=kq_myDataHist(:,:);
    end
    save([resultPathName DASH dirData(i).name '_myDataHist'],'myDataHist')
end
save([resultPathName DASH 'myDataProcessed'],'myDataProcessed')
COMPLETED_MyDataProcessed=1
toc % 14479 secodns or about 4 hours
%
%PLOT CBV WHOLE BRAIN
myDataProcessed(myDataProcessed == 0) = NaN

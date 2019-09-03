function [nBadTs,nBadElecs,nVisBadT,meanBadT,meanVisBadT,nBadEyeTs,numTs,slopeValsVsFreq] = findBadTrialsWithEEG_now(subjectParams,badTrialThreshold)
%     clear;clc;
%     folderIn = 'Z:\ADProject\ADGammaProject_EEG_rawData';
%     subjectName = '163MR_F1';
%     expDate = '280918';
%     capType  = 'actiCap64';
%     subjectFolder = fullfile('Z:\Data\human\ADGammaProject\Data',subjectName,'EEG');
%     protocolName = 'GAV_0003';    
%     channelsNotForConsideration = [65 66]; % TODO: [33 34]
%     tapersPSD = 1;

    %% Initializations

    folderIn = subjectParams.folderIn;
    subjectName = subjectParams.subjectName;
    expDate = subjectParams.expDate;
    capType = subjectParams.capType;
    subjectFolder = subjectParams.subjectFolder;
    protocolName = subjectParams.protocolName;
    channelsNotForConsideration = subjectParams.channelsNotForConsideration;
    tapersPSD= subjectParams.tapersPSD;

    highPassCutOff = 1.6; % Hz
    FsEye = 500; % Hz
    checkPeriod = [-0.500 0.750]; % s
    impedanceTag = 'Start'; % Start or End
    ImpedanceCutOff = 25; % KOhm
    time_threshold  = 6;
    psd_threshold = 3; % also 3
    badTrialThreshold = 30; % Percentage
%     badElecThreshold = 10; % Percentage
%     coocurence_flag = 1;
    
    slopeRange = {[55 85]}; %Hz % slope range used to compute slopes
    freqsToAvoid = {[0 0] [8 12] [46 54] [96 104]}; %Hz

    folderName = fullfile(subjectFolder,expDate,protocolName);
    folderSegment = fullfile(folderName,'segmentedData');
    folderLFP = fullfile(folderSegment,'LFP');
    lfpInfo = load(fullfile(folderLFP,'lfpInfo'),'analogChannelsStored','timeVals');
    analogChannelsStored = lfpInfo.analogChannelsStored;
    timeVals = lfpInfo.timeVals;
    Fs = 1/(timeVals(2) - timeVals(1)); %Hz

    checkBaselinePeriod = [-0.5 0-1/Fs]; %s  %For computing slopes for artifact rejection
    
    if exist('channelsNotForConsideration','var') && ~isempty(channelsNotForConsideration)
        analogChannelsStored(ismember(analogChannelsStored,channelsNotForConsideration))=[];
    end
    numElectrodes = length(analogChannelsStored);

    % Defining filter
    if exist('highPassCutOff','var') || ~isempty(highPassCutOff)
        d1 = designfilt('highpassiir','FilterOrder',8, ...
            'PassbandFrequency',highPassCutOff,'PassbandRipple',0.2, ...
            'SampleRate',Fs);
    end

    % Setting MT parameters
    params.tapers   = [3 5];
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 200];
    params.trialave = 0;

    %% 1. Get bad trial from eye trials

    if exist('FsEye','var') && ~isempty(FsEye)
        badEyeTrials = findBadTrialsFromEyeData(subjectFolder,expDate,protocolName,FsEye,checkPeriod)'; % added by MD 10-09-2017
    else
        badEyeTrials = [];
    end

    %% 2. Get electrode impedances for rejecting noisy electrodes (impedance > 25k)

    filenameAtStart = fullfile(folderIn,[subjectName expDate],[subjectName expDate 'GAV_ImpedanceAt' impedanceTag '.txt']);
    [elecNames,elecImpedance] = getImpedanceEEG(filenameAtStart);
    EEGelectrodeLabels = load('EEGelectrodeLabels.mat','EEGelectrodeLabels');
    EEGelectrodeLabels = EEGelectrodeLabels.EEGelectrodeLabels;
    electrodeLabelsList = EEGelectrodeLabels(2:end,(ismember(EEGelectrodeLabels(1,:),capType)));
    if(strcmp(capType,'actiCap31Posterior')); electrodeLabelsList(32:end)=[]; end

    clear elecInds
    for iML = 1:length(electrodeLabelsList)
        temp = find(strcmp(electrodeLabelsList(iML),elecNames));
        if(isempty(temp)); elecInds(iML)=0;
        else elecInds(iML) = temp; end
    end
    unRegisteredElecs = ~elecInds;
    if(sum(unRegisteredElecs))
        elecImpedance(end+1)=NaN; % inserting NaN at the end
        elecInds(unRegisteredElecs)=length(elecImpedance); % so that unregistered ones also get NaNs
    end
    elecImpedance = elecImpedance(elecInds); % Remap the electrodes according to the standard montage workspace
    GoodElec_Z = elecImpedance<ImpedanceCutOff; % (1) along with this, (2) NaN condition is handled somehow!

    % "nBadElecs{1}" - Removing bad impedance electrodes & '--' electrodes
    nBadElecs{1} = (~GoodElec_Z) | unRegisteredElecs; % (3) Unregistered electrodes condition


    %% 3. Analysis for each trial and each electrode

    allBadTrials = cell(1,numElectrodes);
    hW1 = waitbar(0,'Processing electrodes...');

    for iElec=1:numElectrodes

        if ~GoodElec_Z(iElec); allBadTrials{iElec} = NaN; continue; end % Analyzing only those electrodes with impedance < 25k

        clear analogData
        electrodeNum=analogChannelsStored(iElec);
        waitbar((iElec-1)/numElectrodes,hW1,['Processing electrode: ' num2str(electrodeNum)]);
        analogData = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeNum) '.mat']),'analogData');
        analogData = analogData.analogData;
        analogDataAllElecs(iElec,:,:) = analogData; % to use later for PSDslopes
        analogData(badEyeTrials,:) = [];

        if exist('highPassCutOff','var') || ~isempty(highPassCutOff)    % high pass filter
            clear analogDataSegment; analogDataSegment = analogData;
            clear analogData; analogData = filtfilt(d1,analogDataSegment')';
            clear analogDataSegment
        end

        % determine indices corresponding to the check period
        checkPeriodIndices = timeVals>=checkPeriod(1) & timeVals<=checkPeriod(2);
        analogData = analogData(:,checkPeriodIndices);

        % subtract dc
        analogData = analogData - repmat(mean(analogData,2),1,size(analogData,2));

        % Check time-domain waveforms
        numTrials = size(analogData,1);                            % excluding bad eye trials
        meanTrialData = nanmean(analogData,1);                     % mean trial trace
        stdTrialData = nanstd(analogData,[],1);                    % std across trials

        tDplus = (meanTrialData + (time_threshold)*stdTrialData);    % upper boundary/criterion
        tDminus = (meanTrialData - (time_threshold)*stdTrialData);   % lower boundary/criterion

        tBoolTrials = sum((analogData > ones(numTrials,1)*tDplus) | (analogData < ones(numTrials,1)*tDminus),2);

        clear tmpBadTrialsTime;
        tmpBadTrialsTime = find(tBoolTrials>0);

        % Check PSD
        clear powerVsFreq;
        [powerVsFreq,~] = mtspectrumc(analogData',params);
        powerVsFreq = powerVsFreq';

        clear meanTrialData stdTrialData tDplus
        meanTrialData = nanmean(powerVsFreq,1);                     % mean trial trace
        stdTrialData = nanstd(powerVsFreq,[],1);                    % std across trials

        tDplus = (meanTrialData + (psd_threshold)*stdTrialData);    % upper boundary/criterion
        clear tBoolTrials; tBoolTrials = sum((powerVsFreq > ones(numTrials,1)*tDplus),2);
        clear tmpBadTrialsPSD; tmpBadTrialsPSD = find(tBoolTrials>0);

        tmpBadTrialsAll = unique([tmpBadTrialsTime;tmpBadTrialsPSD]);

        % Remap bad trial indices to original indices
        originalTrialInds = 1:size(analogDataAllElecs,2);
        originalTrialInds(badEyeTrials) = [];
        allBadTrials{iElec} = originalTrialInds(tmpBadTrialsAll);

    end
    close(hW1);

    %% 4. x% bad trials threshold - "nBadElecs{2}" - Remove electrodes containing more than x% bad trials

    badTrialUL = (badTrialThreshold/100)*numTrials;
    badTrialLength=cellfun(@length,allBadTrials);
    badTrialLength(nBadElecs{1})=0; % Removing the bad impedance electrodes
    nBadElecs{2} = logical(badTrialLength>badTrialUL)';

    %% 5. Trimming bad Trials

    badTrials = trimBadTrials(allBadTrials,nBadElecs);

    %% 6. PSD Slope calculation across baseline period
    elec_specific = 0; %SWITCH
    checkPeriodIndicesPSD = timeVals>=checkBaselinePeriod(1) & timeVals<=checkBaselinePeriod(2);    
%     analogDataAllElecs_good = analogDataAllElecs(:,setdiff(originalTrialInds,badTrials),checkPeriodIndicesPSD);
    params.tapers   = [(tapersPSD+1)/2 tapersPSD];
    for iElec=1:numElectrodes
        if(elec_specific)
            analogDataAllElecs_good{iElec} = analogDataAllElecs(iElec,setdiff(originalTrialInds,allBadTrials{iElec}),checkPeriodIndicesPSD);
        else
            analogDataAllElecs_good{iElec} = analogDataAllElecs(iElec,setdiff(originalTrialInds,badTrials),checkPeriodIndicesPSD);
        end
        % Computing slopes
%         analogDataPSD = squeeze(analogDataAllElecs_good(iElec,:,:));
        analogDataPSD = squeeze(analogDataAllElecs_good{iElec});
        analogDataPSD = analogDataPSD - repmat(mean(analogDataPSD,2),1,size(analogDataPSD,2));
        clear powerVsFreq;
        [powerVsFreq,freqVals] = mtspectrumc(analogDataPSD',params);
        powerVsFreq = powerVsFreq';
%         powerVsFreqAllElecs_nomean(iElec,:,:)=powerVsFreq;
        powerVsFreqAllElecs_cell{iElec}=powerVsFreq;
        powerVsFreqAllElecs(iElec,:)= mean(powerVsFreq,1);
        slopeValsVsFreq{iElec} = getSlopesPSDBaseline_v2(log10(mean(powerVsFreq,1)),freqVals,slopeRange,[],freqsToAvoid);
    end
    [goodSlopeFlag] = getGoodSlopeFlag_v3(slopeValsVsFreq,capType);
    nBadElecs{3}=~goodSlopeFlag; % we need bad electrode labels

    %% PSD fit (MSE) figure
    plot_psds_slope = 0;
    if(plot_psds_slope)
        temp_freqStr = load('freq_ft_str_temp.mat');
        freqBaseline = temp_freqStr.freq_ft;
        freqBaseline.powspctrm = log10(powerVsFreqAllElecs);
        freqBaseline.freq = freqVals;
        
        freqBaselineFit = temp_freqStr.freq_ft;
        freqBaselineFit.powspctrm = log10(psdSlopefit_v2(slopeValsVsFreq,freqVals));
        freqBaselineFit.freq = freqVals;
        
        fPos = union(find(freqVals<slopeRange{1}(1)),find(freqVals>slopeRange{1}(2)));
        %     freqBaselineFit.powspctrm(:,fPos) = 0;
        
        [~, ftdir] = ft_version;
        x=load(fullfile(ftdir,'template','layout','acticap-64ch-standard2.mat'));
        layout = x.lay;
        
        cfg = [];
        cfg.layout        = layout;
        cfg.showlabels    = 'yes';
        cfg.showoutline   = 'yes';
        cfg.ylim = [-2 1];
        figure; ft_multiplotER(cfg,freqBaseline,freqBaselineFit);
        
        %%
        figure;
        slopes=cellfun(@get2ndElement,slopeValsVsFreq);
        topoplot(slopes,'D:\Santosh\data\elec_locs\valid_locs.sph','maplimits',[0 6],'emarker',{'.','k',20,1});
        colorbar;
        title([subjectName ' ' protocolName]);
    end
    %% 7. Update final bad Trial
    % badTrials should be redefined after rejecting electrodes with bad PSD slopes
    % which hasn't been performed now.

    %% Assigning outputs

    nBadTs = length(badTrials)/numTrials; % excluding bad eye trials
    if strcmp(capType,'actiCap64');elecs = [24 57 58 26 61 62 63 29 30 31]; end
    if strcmp(capType,'actiCap31Posterior'); elecs = [15 16 18 19 23 24 25 28 29 30]; end
    badTrials_vis= trimBadTrials(allBadTrials(elecs),nBadElecs);
    nVisBadT = length(badTrials_vis)/numTrials;
    BTcountOverelecs = cellfun(@length,allBadTrials);
    meanBadT = mean(BTcountOverelecs)/numTrials;
    meanVisBadT = mean(BTcountOverelecs(elecs))/numTrials;
    numTs = numTrials + length(badEyeTrials); % excluding bad eye trials
    nBadEyeTs = length(badEyeTrials)/numTs;
   
    
end

function [newBadTrials] =  trimBadTrials(allBadTrials,nBadElecs)
    badElecThreshold = 10; % Percentage
    % a. Taking union across bad electrodes for conditions 1 and 2
    newBadTrials=[];
    numElectrodes = length(allBadTrials);
    for iElec=1:numElectrodes
        if ~(nBadElecs{2}(iElec)|| nBadElecs{1}(iElec)); newBadTrials=union(newBadTrials,allBadTrials{iElec}); end
    end

    % b. Co-occurence condition - Counting the trials which occurs in more than x% of the electrodes
    for iTrial = 1:length(newBadTrials)
        badTrialElecs(iTrial)=0;
        for iElec = 1:numElectrodes
            if isnan(allBadTrials{1,iElec}); continue; end % Discarding the electrodes where the bad trials are NaN
            % because of this NaN entries in badTrials have zero in 'badTrialElecs'
            if find(newBadTrials(iTrial)==allBadTrials{1,iElec})
                badTrialElecs(iTrial) = badTrialElecs(iTrial)+1;
            end
        end
    end
    newBadTrials(badTrialElecs<(badElecThreshold/100.*numElectrodes))=[];
end
% Used this function to debug, not needed now.
function out = get2ndElement(a)
out=a{2};
end
function PlayWavNP_LED(ratID, params)
    %% Basic parameters
    % Time of one trial (in seconds)
    load gong;
    ttrial = 60*30; % 30 min
    
    Fs = 192000; % Wav write and playback sampling rate
    totaldurSec = 1; % length of a stimulus in seconds
    totaldur = Fs * totaldurSec; % samples
    vol = 1; % volume of sound
    maxleverpresses = 200; %ceil(ttrial/totaldurSec);
    
    vSampleRate = 44; % nosepoke voltage sampling rate
        
    stimtype = 'soundNP_LED';
    
    %% Arduino
    voltage1 = zeros(ttrial*vSampleRate,1);
    voltage2 = zeros(ttrial*vSampleRate,1);
    disp('Connecting to arduino 1 (COM3)...');
    ard1 = arduino('COM3','uno');
    disp('Connected to arduino 1');
    
    disp('Connecting to arduino 2 (COM4)...');
    ard2 = arduino('COM4','uno');
    disp('Connected to arduino 2');

    %% Set parameters from args
    if nargin < 1 || isempty(ratID)
        % Rat ID
        ratID = 'test';
    end
    if nargin < 2 || isempty(params)
        % Other parameters
        s1 = 'gong';
        s2 = 'none';
        % TODO: Make sure to counterbalance var1 and var2 for jitter
        % TODO: Make sure to counterbalance rarity of stims for unigrams
    else
        % Read parameters from arg params
        % params = struct('N',N,'W',W,'C0',C0,'G0',G0,'U',U,'C',C,'T',T,'D',D,'THETA',THETA);
%        avgint = params.avgint;  % ms
        s1 = params.s1;
        s2 = params.s2;
    end

    %% Directory for saving
    if IsWin
        savedir = 'C:\Users\OIST\ownCloud\Thesis\BehExp\MusicLounge\';
        sep = '\';
    else
        savedir = '/Users/yoriko/ownCloud/Thesis/BehExp/MusicLounge/';
        sep='/';
    end

    %% Read sounds
    disp('Reading wav files');
    if strcmp(s1,'pinga')
        wavedata1 = audioread('pinga1_5.wav')';
        wavedata1 = wavedata1(1:totaldur);
    elseif strcmp(s1,'none')
        wavedata1 = zeros(1,totaldur);
    elseif strcmp(s1,'gong')
        wavedata1 = repmat(y(1:floor(length(y)*0.6)),8,1);
        wavedata1 = wavedata1(1:totaldur)';
    end

    if strcmp(s2,'pinga')
        wavedata2 = audioread('pinga1_5.wav')';
        wavedata2 = wavedata2(1:totaldur);
    elseif strcmp(s2,'none')
        wavedata2 = zeros(1,totaldur);
    elseif strcmp(s2,'gong')
        wavedata2 = repmat(y(1:floor(length(y)*0.6)),8,1);
        wavedata2 = wavedata2(1:totaldur)';
    end
    disp('Loaded wav files');

    %% Keyboard
    KbName('UnifyKeyNames');
    escape = KbName('Escape');
    enter = KbName('Return');

    %% Initialize PsychSound
    % Select one channel mono playback:
    nrchannels = 1;

    % Open sound device 'pahandle' with specified frequency and number of audio
    % channels for playback in timing precision mode on the default audio
    % device:
    InitializePsychSound(1);
    pamaster = PsychPortAudio('Open', [], 1+8, 1, Fs, nrchannels,[],[]); %Somehow this works!
    PsychPortAudio('Start', pamaster, 0, 0, 1);
    %pahandle = PsychPortAudio('Open', [], [], 0, Fs, nrchannels);
    PsychPortAudio('Volume', pamaster, vol);
    
    pahandle1 = PsychPortAudio('OpenSlave', pamaster, 1);
    PsychPortAudio('FillBuffer', pahandle1, wavedata1);
        
    pahandle2 = PsychPortAudio('OpenSlave', pamaster, 1);
    PsychPortAudio('FillBuffer', pahandle2, wavedata2);
    
    disp('Initialized PsychSound.');
    
    %% Variables for recording
    %baselineNP = zeros(50,2);
    %bnpind = 1;
    choiceseq = zeros(maxleverpresses,1);
    started = 0;
    aborted = 0;
    trials = 0;
    nosepoke = 0;
    triggerTimes = zeros(maxleverpresses,2);
    
    %% Loop to record nosepokes
    disp('Hit any key to start recording nosepokes');
    KbStrokeWait;
    disp('Beginning trials...');
    beginTimeRecord = GetSecs;
    beginDateRecord = datetime('now');
    
    vindex = 1;    
    disp('Hit Enter to start trial, Esc to abort.');
    j = 1; % count lever presses
    
    while vindex <= length(voltage1)
        voltage1(vindex) = readVoltage(ard1,'A0')+1; % read photo transistor value from arduino
        voltage2(vindex) = readVoltage(ard2,'A0')+1; % read photo transistor value from arduino
            %disp(voltage(vindex));
            %plot(voltage); % plot voltage value
            %drawnow;
            
        [ keyIsDown, keyTime, keyCodeTemp ] = KbCheck;
        if keyIsDown
            fprintf('\n About %d seconds remaining...\n',floor((length(voltage1)-vindex)/vSampleRate)); % press any key to show remaining time
            keyCode = keyCodeTemp;
            KbReleaseWait;

            if keyCode(escape) % To abort
                abortTime = GetSecs;
                fprintf('Aborting recording at %d sec\n', abortTime - beginTimeRecord);
                aborted = 1;
                break;
            elseif keyCode(enter) % To begin trial
                started = 1;
                tStartSound = uint64(0);
                j_trial = 1;
                triggerTimes = zeros(maxleverpresses,2);
                beginDateTrial = datetime('now');
                savefilenameStart = strcat(stimtype,sep,ratID,'_',datestr(beginDateTrial,'mm-dd-yyyy_HH-MM-SS-FFF'));
                beginTimeTrial = GetSecs;
                tic; % tic/toc pair 1
                fprintf(strcat('Beginning trial %d at\t', datestr(beginDateTrial,'HH:MM:SS')),trials+1);
            end
            % Reset keyCode array
            keyCode = zeros(size(keyCode));
        end
        
        %% Play sound in response to nosepokes
        if started % if trial has started
            if (toc < ttrial) && (j <= maxleverpresses)% tic/toc pair 1
                if (voltage1(vindex) < 4) && (toc(tStartSound) > totaldurSec) % tic/toc pair 2
                    choiceseq(j) = 1;
                    triggerTimes(j_trial,1) = GetSecs;
                    if ~strcmp(s1,'none')
                        writeDigitalPin(ard1,'D12',1);
                        disp('Turned on LED 1');
                        PsychPortAudio('Start', pahandle1, [], triggerTimes(j_trial,1));
                        disp('Playing stimulus 1...');
                    else
                        disp('!Not playing anything');
                    end
                    
                    tStartSound = tic; % tic/toc pair 2

                    if ~strcmp(s1,'none')
                        PsychPortAudio('Stop', pahandle1, 3);
                        writeDigitalPin(ard1,'D12',0);
                        disp('Turned off LED 1');
                    end
                    
                    voltage1(vindex+1:vindex+totaldurSec*vSampleRate) = ones(totaldurSec*vSampleRate,1);
                    voltage2(vindex+1:vindex+totaldurSec*vSampleRate) = 10*ones(totaldurSec*vSampleRate,1);
                    vindex = vindex + totaldurSec*vSampleRate; % 
                    j = j+1;
                    j_trial = j_trial+1;
                elseif (voltage2(vindex) < 4) && (toc(tStartSound) > totaldurSec)
                    choiceseq(j) = 2;
                    triggerTimes(j_trial,2) = GetSecs;
                    if ~strcmp(s2,'none')
                        writeDigitalPin(ard2,'D12',1);
                        disp('Turned on LED 2');
                        PsychPortAudio('Start', pahandle2, [], triggerTimes(j_trial,2));
                        disp('Playing stimulus 2...');
                    else
                        disp('!Not playing anything');
                    end
                    tStartSound = tic; % tic/toc pair 2
                    if ~strcmp(s2,'none')
                        PsychPortAudio('Stop', pahandle2, 3);
                        writeDigitalPin(ard2,'D12',0);
                        disp('Turned off LED 2');
                    end

                    voltage1(vindex+1:vindex+totaldurSec*vSampleRate) = 10*ones(totaldurSec*vSampleRate,1);
                    voltage2(vindex+1:vindex+totaldurSec*vSampleRate) = ones(totaldurSec*vSampleRate,1);
                    vindex = vindex + totaldurSec*vSampleRate; % 
                    j = j+1;
                    j_trial = j_trial+1;
                end
            else % close trial
                started = 0;
                fprintf('Ending trial %d at %f sec \n', trials, GetSecs-beginTimeRecord);
                trials = trials + 1;
                voltage1(vindex+1) = -1;
                voltage2(vindex+1) = -1;
                vindex = vindex + 1;
                % Save nosepoke times to CSV file
                csvwrite(strcat(savedir,savefilenameStart,'_NPtime.csv'),triggerTimes-beginTimeTrial);
                disp('Hit Enter to start trial, Esc to abort.');
            end
        else
             if (voltage1(vindex) < 3) && nosepoke == 0
%                 baselineNP(bnpind,2) = GetSecs;
                 disp('Nosepoke 1!');
                 nosepoke = 1;
%                 baselineNP(bnpind,1) = 1;
%                 bnpind = bnpind + 1;
             elseif (voltage2(vindex) < 3) && nosepoke == 0
% %                 baselineNP(bnpind,2) = GetSecs;
                 disp('Nosepoke 2!');
                 nosepoke = 2;
% %                 baselineNP(bnpind,1) = 1;
             else
                 nosepoke = 0;
             end
        end
        vindex = vindex + 1;
    end

    % End of session, close down driver & arduino:
    endTimeRecord = GetSecs;
    fprintf('Duration of recording %d seconds \n',endTimeRecord-beginTimeRecord);
    fprintf('Nosepoke samples %d \n',nnz(voltage1));

    PsychPortAudio('Close');
    writeDigitalPin(ard1,'D13',0);
    clear ard1 ard2;
    endTime = datetime('now');
    disp('Ending session...');
    %BaselineNP(:,2) = BaselineNP(:,2) - beginTimeRecord;
    
    %% Save summaries
    savefilenameEnd = strcat(stimtype,sep,ratID,'_',datestr(endTime,'mm-dd-yyyy_HH-MM-SS-FFF'));
    
    % Save baseline nosepokes to CSV file
    %csvwrite(strcat(savedir,savefilenameStart,'_BLNP.csv'),BaselineNP);
    
    % Save nosepoke voltage trace to CSV file
    csvwrite(strcat(savedir,savefilenameEnd,'_voltage1.csv'),nonzeros(voltage1));
    csvwrite(strcat(savedir,savefilenameEnd,'_voltage2.csv'),nonzeros(voltage2));
    
    % Save summary to TXT file
    fh = fopen(strcat(savedir,savefilenameEnd,'.txt'),'w');

    fprintf(fh,'------------Experiment summary------------\n');
    fprintf(fh,strcat('BEGIN SESSION:\t',datestr(beginDateRecord),'\n'));
    fprintf(fh,strcat('END SESSION:\t',datestr(endTime),'\n'));
    fprintf(fh,strcat('Number of completed trials:\t',num2str(trials),'\n'));
    fprintf(fh,strcat('Assumed voltage sampling rate:\t',num2str(vSampleRate),'\n'));
    
    fprintf(fh,strcat('Stim 1 (COM3):\t',s1,'\t Stim 2 (COM4):\t',s2,'\n'));
    
    fprintf(fh,strcat('Rat ID:\t\t',ratID,'\n'));

    fprintf(fh,strcat('Lever Press times:\t',mat2str(round(10*(triggerTimes-beginTimeRecord))/10),'\n'));
    fprintf(fh,strcat('Choices:\t',mat2str(choiceseq),'\n'));

    fprintf(fh,strcat('Saved NP voltage trace in:\t',savefilenameEnd,'_voltage1.csv,',savefilenameEnd,'_voltage2.csv\n'));
    
    if aborted
        fprintf(fh,'Aborted session at %d sec \n',abortTime - beginTimeRecord);
    end

    fprintf(fh,'------------------------------------------\n');
    fclose(fh);

    disp('Saved summaries.');
    
    %% Plot nosepoke voltage time series
    plot(nonzeros(voltage1),'b');
    hold on;
    plot(nonzeros(voltage2),'r');
    ylim([0,6]);
    legend(strcat(s1,'(COM3)'), strcat(s2,'(COM4)'));
    savefig(strcat(savedir,savefilenameEnd,'_voltage.fig'));
end
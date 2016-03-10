function PlayOddballNP_LED(ratID, session, arduinoON, params)
    %% Basic parameters
    % Time of one trial (in seconds)
    load gong;
    ttrial = 60*25; % 25 min
    
    Fs = 192000; % Wav write and playback sampling rate
    totaldurSec = 1; % length of a stimulus in seconds
    totaldur = Fs * totaldurSec; % samples
    vol = 1; % volume of sound
    maxnosepokes = 50; %ceil(ttrial/totaldurSec);
    
    vSampleRate = 46; % nosepoke voltage sampling rate
        
    stimtype = 'oddballNP_LED';
    
    voltage1 = zeros(1.5*ttrial*vSampleRate,1);
    voltage2 = zeros(1.5*ttrial*vSampleRate,1);
    
    %% Set parameters from args
    if nargin < 1 || isempty(ratID)
        % Rat ID
        ratID = 'test';
    end
    if nargin < 2 || isempty(session)
        session = 0;
    end
    if nargin < 3 || isempty(arduinoON)
        arduinoON = 0; % test mode
    end
    if nargin < 4 || isempty(params)
        % Other parameters
        s1 = 'gong';
        s2 = 'none';
        nodd = 0;
        % TODO: Make sure to counterbalance var1 and var2 for jitter
        % TODO: Make sure to counterbalance rarity of stims for unigrams
    else
        % Read parameters from arg params
        % params = struct('N',N,'W',W,'C0',C0,'G0',G0,'U',U,'C',C,'T',T,'D',D,'THETA',THETA);
%        avgint = params.avgint;  % ms
        s1 = params.s1;
        s2 = params.s2;
        nodd = params.nodd;
    end

    %% Arduino
    if arduinoON
        disp('Connecting to arduino 1 (COM3)...');
        ard1 = arduino('COM3','uno');
        disp('Connected to arduino 1');

        disp('Connecting to arduino 2 (COM4)...');
        ard2 = arduino('COM4','uno');
        disp('Connected to arduino 2');
    else
        v1Sample = repmat(csvread('voltage1.csv'),2,1);
        v2Sample = repmat(csvread('voltage2.csv'),2,1);
    end
    
    fprintf('Number of oddballs:\t %d \n',nodd);
    
    %% Stimulus files
    filename1 = strcat(s1,'_oddball_seq_a_');
    filename2 = strcat(s2,'_oddball_seq_b_');
    
    %% Create directory for saving
    if IsWin
        savedir = 'C:\Users\OIST\ownCloud\Thesis\BehExp\behavioral';
        sep = '\';
    else
        savedir = '/Users/yoriko/ownCloud/Thesis/BehExp/behavioral';
        sep='/';
    end
    savedir = strcat(savedir,sep,ratID,int2str(session));
    mkdir(savedir);
    disp(strcat('Working in directory ',savedir));

%     %% Read sounds
%     disp('Reading wav files');
%     if strcmp(s1,'pinga')
%         wavedata1 = audioread('pinga1_5.wav')';
%         wavedata1 = wavedata1(1:totaldur);
%     elseif strcmp(s1,'none')
%         wavedata1 = zeros(1,totaldur);
%     elseif strcmp(s1,'gong')
%         wavedata1 = repmat(y(1:floor(length(y)*0.6)),20,1);
%         wavedata1 = wavedata1(1:totaldur)';
%     end
% 
%     if strcmp(s2,'pinga')
%         wavedata2 = audioread('pinga1_5.wav')';
%         wavedata2 = wavedata2(1:totaldur);
%     elseif strcmp(s2,'none')
%         wavedata2 = zeros(1,totaldur);
%     elseif strcmp(s2,'gong')
%         wavedata2 = repmat(y(1:floor(length(y)*0.6)),20,1);
%         wavedata2 = wavedata2(1:totaldur)';
%     end
%     disp('Loaded wav files');

    %% Generate ntrial*2 stimulus files for ntrial sequences each    
    x = audioread('pinga1_5.wav');
    beep = zeros(2,floor(length(y)*0.6));
    
    beep(1,:) = y(1:floor(length(y)*0.6))';
    %beep(2,:) = [zeros(1,floor(length(y)*0.6)-floor(length(y)*0.3)),y(1:floor(length(y)*0.3))'];
    beep(2,:) = x(length(x)-floor(length(y)*0.6)+1:length(x))';
    
    nbeeps = 8;
    [intvs1,stims1] = gen_odd(nodd,nbeeps,maxnosepokes,0);
    [intvs2,stims2] = gen_odd(nodd,nbeeps,maxnosepokes,0);
    
    for k=1:maxnosepokes
        %% Generate sequence of beeps 1
        beeps1 = beep(stims1(k,1),:);
        
        for j = 2:nbeeps
            beeps1=[beeps1,beep(stims1(k,j),:)];
        end
        % don't truncate at totaldur
        % beeps1 = beeps1(1:totaldur);
        
        if ~strcmp(s1,'none')
            audiowrite(strcat(savedir,sep,filename1,int2str(k),'.wav'),beeps1,Fs);
        end
            
        %% Generate sequence of beeps 2
        beeps2 = beep(stims2(k,1),:);
        
        for j = 2:nbeeps
            beeps2=[beeps2,beep(stims2(k,j),:)];
        end
        % don't truncate at totaldur
        % beeps2 = beeps2(1:totaldur);

        if ~strcmp(s2,'none')
            audiowrite(strcat(savedir,sep,filename2,int2str(k),'.wav'),beeps2,Fs);
        end
    end
    
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
    
    pahandle1 = zeros(maxnosepokes,1);
    pahandle2 = zeros(maxnosepokes,1);
    
    for ii = 1:maxnosepokes
        if ~strcmp(s1,'none')
            pahandle1(ii) = PsychPortAudio('OpenSlave', pamaster, 1);
            wavedata1 = audioread(strcat(savedir,sep,filename1,int2str(ii),'.wav'))';
            PsychPortAudio('FillBuffer', pahandle1(ii), wavedata1);
            disp(strcat('Loaded ',filename1,int2str(ii),'.wav'));
        end
        if ~strcmp(s2,'none') 
            pahandle2(ii) = PsychPortAudio('OpenSlave', pamaster, 1);
            wavedata2 = audioread(strcat(savedir,sep,filename2,int2str(ii),'.wav'))';
            PsychPortAudio('FillBuffer', pahandle2(ii), wavedata2);
            disp(strcat('Loaded ',filename2,int2str(ii),'.wav'));
        end
    end
    
    %pahandle1 = PsychPortAudio('OpenSlave', pamaster, 1);
    %PsychPortAudio('FillBuffer', pahandle1, wavedata1);
        
    %pahandle2 = PsychPortAudio('OpenSlave', pamaster, 1);
    %PsychPortAudio('FillBuffer', pahandle2, wavedata2);
    
    disp('Initialized PsychSound.');
    
    %% Variables for recording
    choiceseq = zeros(maxnosepokes,1);
    started = 0;
    aborted = 0;
    nosepoke = 0;
    wait1 = 0;
    wait2 = 0;
    triggerTimes = zeros(maxnosepokes,2);
    
    %% Loop to record nosepokes
    disp('Hit any key to start recording nosepokes');
    KbStrokeWait;
    beginTimeRecord = GetSecs;
    beginDateRecord = datetime('now');
    fprintf(strcat('Beginning trials at \t',datestr(beginDateRecord,'HH:MM:SS'),'...\n'));
    
    vindex = 1;    
    disp('Hit Enter to start trial, Esc to abort.');
    j = 1; % count lever presses
    
    while vindex <= length(voltage1)
        if arduinoON
            voltage1(vindex) = readVoltage(ard1,'A0')+1; % read photo transistor value from arduino
            voltage2(vindex) = readVoltage(ard2,'A0')+1; % read photo transistor value from arduino
        else
            voltage1(vindex) = v1Sample(vindex); % for testing code
            voltage2(vindex) = v2Sample(vindex); % for testing code
        end
            
        [ keyIsDown, keyTime, keyCodeTemp ] = KbCheck;
        if keyIsDown
            fprintf('\n About %d seconds remaining...\n',floor((length(voltage1)-vindex)/vSampleRate)); % press any key to show remaining time
            keyCode = keyCodeTemp;
            KbReleaseWait;

            if keyCode(escape) % To abort
                abortTime = GetSecs;                    
                fprintf('Aborting recording at %d sec\n', abortTime - beginTimeRecord);
                voltage1(vindex+1) = -1;
                voltage2(vindex+1) = -1;
                vindex = vindex + 1;
                aborted = 1;
                break;
            elseif keyCode(enter) % To begin trial
                started = 1;
                tStartSound1 = uint64(0);
                tStartSound2 = uint64(0);
                j_trial = 1;
                beginDateTrial = datetime('now');
                beginTimeTrial = GetSecs;
                voltage1(vindex+1) = 6;
                voltage2(vindex+1) = 6;
                vindex = vindex + 1;
                tic; % tic/toc pair 1
                fprintf(strcat('Beginning trial at\t', datestr(beginDateTrial,'HH:MM:SS'),'...\n'));
            end
            % Reset keyCode array
            keyCode = zeros(size(keyCode));
        end
        
        %% Play sound in response to nosepokes
        if started % if trial has started
            if (toc < ttrial) && (j <= maxnosepokes)% tic/toc pair 1
                
                if ~wait1 % If voltage has been up
                    if (voltage1(vindex) < 4) && (toc(tStartSound1) > totaldurSec) && (toc(tStartSound2) > totaldurSec) % tic/toc pair 2
                        choiceseq(j) = 1;
                        triggerTimes(j_trial,1) = GetSecs;
                        if ~strcmp(s1,'none')
                            if arduinoON
                                writeDigitalPin(ard1,'D12',1);
                            end
                            %disp('Turned on LED 1');
                            PsychPortAudio('Start', pahandle1(j), [], triggerTimes(j_trial,1));
                        else
                            disp('!Not playing anything');
                        end
                        tStartSound1 = tic; % tic/toc pair 2
                        wait1 = 1; % waiting for voltage1 to go back up
                        
                        if ~strcmp(s1,'none')
                            fprintf(strcat('Playing\t',filename1,int2str(j),'.wav...\n'));
                            PsychPortAudio('Stop', pahandle1(j), 3);
                            if arduinoON
                                writeDigitalPin(ard1,'D12',0);
                            end
                            %disp('Turned off LED 1');
                        end

                        voltage1(vindex+1:vindex+totaldurSec*vSampleRate) = 0.5*ones(totaldurSec*vSampleRate,1);
                        voltage2(vindex+1:vindex+totaldurSec*vSampleRate) = 5.45*ones(totaldurSec*vSampleRate,1);
                        vindex = vindex + totaldurSec*vSampleRate; % 
                        j = j+1;
                        j_trial = j_trial+1;
                    elseif voltage1(vindex)<4
                        disp('---still playing sound 1');
                    end
                else % If voltage has not gone up
                    if mean(voltage1(vindex-vSampleRate/2+1:vindex))>4
                        wait1 = 0;
                    else
                        %disp('Waiting for voltage 1 to go back up');
                    end
                end
                
                if ~wait2 % If voltage has been up
                    if (voltage2(vindex) < 4) && (toc(tStartSound2) > totaldurSec) && (toc(tStartSound1) > totaldurSec) % tic/toc pair 2
                        choiceseq(j) = 2;
                        triggerTimes(j_trial,2) = GetSecs;
                        if ~strcmp(s2,'none')
                            if arduinoON
                                writeDigitalPin(ard2,'D12',1);
                            end
                            %disp('Turned on LED 2');
                            PsychPortAudio('Start', pahandle2(j), [], triggerTimes(j_trial,2));
                        else
                            disp('!Not playing anything');
                        end
                        tStartSound2 = tic; % tic/toc pair 2
                        wait2 = 1; % waiting for voltage2 to go back up
                        
                        if ~strcmp(s2,'none')
                            fprintf(strcat('Playing\t',filename2,int2str(j),'.wav...\n'));
                            PsychPortAudio('Stop', pahandle2(j), 3);
                            if arduinoON
                                writeDigitalPin(ard2,'D12',0);
                            end
                            %disp('Turned off LED 2');
                        end

                        voltage1(vindex+1:vindex+totaldurSec*vSampleRate) = 5.35*ones(totaldurSec*vSampleRate,1);
                        voltage2(vindex+1:vindex+totaldurSec*vSampleRate) = 0.5*ones(totaldurSec*vSampleRate,1);
                        vindex = vindex + totaldurSec*vSampleRate; % 
                        j = j+1;
                        j_trial = j_trial+1;
                    elseif voltage2(vindex)<4
                        disp('---still playing sound 2');
                    end
                else % If voltage has not gone up
                    if mean(voltage2(vindex-vSampleRate/2+1:vindex))>4
                        wait2 = 0;
                    else
                        %disp('Waiting for voltage 2 to go back up');
                    end
                end
                
            else % close trial
                fprintf('Ending trial at %f sec \n', GetSecs-beginTimeRecord);
                voltage1(vindex+1) = -1;
                voltage2(vindex+1) = -1;
                vindex = vindex + 1;
                
                disp('Hit Esc to abort recording.');
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

    %% End of session, close down driver & arduino:
    endTimeRecord = GetSecs;
    fprintf('Duration of recording %d seconds \n',endTimeRecord-beginTimeRecord);
    fprintf('Nosepoke samples %d \n',nnz(voltage1));

    PsychPortAudio('Close');
    if arduinoON
        writeDigitalPin(ard1,'D12',0);
        writeDigitalPin(ard2,'D12',0);
    
        clear ard1 ard2;
    end
    endDateRecord = datetime('now');
    fprintf(strcat('Ending session at\t',datestr(endDateRecord,'HH:MM:SS'),'...\n'));
    
    %% Save nosepoke time data
    if ~isempty(nonzeros(triggerTimes(:,1)))
        csvwrite(strcat(savedir,sep,'NPtime1.csv'),nonzeros(triggerTimes(:,1))-beginTimeTrial);
    else
        disp('No nosepokes in COM3');
    end
    if ~isempty(nonzeros(triggerTimes(:,2)))
        csvwrite(strcat(savedir,sep,'NPtime1.csv'),nonzeros(triggerTimes(:,1))-beginTimeTrial);
    else
        disp('No nosepokes in COM4');
    end

    %% Save voltage data
    csvwrite(strcat(savedir,sep,'voltage1.csv'),nonzeros(voltage1));
    csvwrite(strcat(savedir,sep,'voltage2.csv'),nonzeros(voltage2));
    
    %% Save stimuli in CSV
    if ~strcmp(s1,'none')
        csvwrite(strcat(savedir,sep,filename1,'.csv'),stims1');
    elseif ~strcmp(s2,'none')
        csvwrite(strcat(savedir,sep,filename2,'.csv'),stims2');
    end
    
    %% Save choice sequence in CSV
    csvwrite(strcat(savedir,sep,'choices.csv'),choiceseq);
    
   %% Save actually used stimuli to CSV file
    used_stims = zeros(nbeeps,maxnosepokes);
%     used_intvs = zeros(nbeeps,maxleverpresses);
    for lp = 1:maxnosepokes
        if choiceseq(lp) == 1
            used_stims(:,lp) = stims1(lp,:)';
%             used_intvs(:,lp) = intvs1(:,lp);
        elseif choiceseq(lp) == 2
            used_stims(:,lp) = stims2(lp,:)';
%             used_intvs(:,lp) = intvs2(:,lp);
        end
    end
    csvwrite(strcat(savedir,sep,'stim_used.csv'),used_stims);    
    
    %% Save session data to MAT file
<<<<<<< HEAD
    if ~aborted
        save(strcat(savedir,sep,'sessiondata.mat'),'ratID','arduino','s1','s2','session','ttrial','nodd','beginTimeRecord','beginDateRecord','beginTimeTrial','beginDateTrial','endTimeRecord','endDateRecord');
    else
        save(strcat(savedir,sep,'sessiondata.mat'),'ratID','arduino','s1','s2','session','ttrial','nodd','beginTimeRecord','beginDateRecord','beginTimeTrial','beginDateTrial','endTimeRecord','endDateRecord','abortTime');
=======
    if aborted && started
        save(strcat(savedir,sep,'sessiondata.mat'),'ratID','arduinoON','s1','s2','session','maxnosepokes','nodd','ttrial','Fs','totaldurSec','vSampleRate','started','aborted','beginTimeRecord','beginDateRecord','beginTimeTrial','beginDateTrial','endTimeRecord','endDateRecord','abortTime');
    elseif started
        save(strcat(savedir,sep,'sessiondata.mat'),'ratID','arduinoON','s1','s2','session','maxnosepokes','nodd','ttrial','Fs','totaldurSec','vSampleRate','started','aborted','beginTimeRecord','beginDateRecord','beginTimeTrial','beginDateTrial','endTimeRecord','endDateRecord');
    elseif aborted
        save(strcat(savedir,sep,'sessiondata.mat'),'ratID','arduinoON','s1','s2','session','maxnosepokes','nodd','ttrial','Fs','totaldurSec','vSampleRate','started','aborted','beginTimeRecord','beginDateRecord','endTimeRecord','endDateRecord','abortTime');
    else
        save(strcat(savedir,sep,'sessiondata.mat'),'ratID','arduinoON','s1','s2','session','maxnosepokes','nodd','ttrial','Fs','totaldurSec','vSampleRate','started','aborted','beginTimeRecord','beginDateRecord','endTimeRecord','endDateRecord');
>>>>>>> refs/remotes/origin/new-NP-counting
    end
    
    
    disp('Saved data.');
    
    %% Plot nosepoke voltage time series
    plot(nonzeros(voltage1),'b');
    hold on;
    plot(nonzeros(voltage2),'r');
    ylim([0,6]);
    legend(strcat(s1,'(COM3)'), strcat(s2,'(COM4)'));
    savefig(strcat(savedir,sep,'voltage.fig'));
end
vSampleRate = 46;

names = cellstr(['fuku  ';'usagi ';'kuro  ';'hidari';'tsuki ';'migi  ';'sen   ';'mohi1 ']);
condition = [1,1,1,1,2,2,2,2];
activeside = [1,1,2,2,1,1,2,2];

if exist('nosepokes_total.csv', 'file') == 2
    delete('nosepokes_total.csv'); % delete old file if it exists
end

for j = 1:8
    files = dir(strcat(names{j},'*'));
    for k = 1:length(files)
        if isdir(files(k).name)
            cd(files(k).name);
            trial = str2double(files(k).name(length(names{j})+1));
            fprintf('Reading trial:\t%d\n',trial);

            %% Read and get nosepokes
            v1=csvread('voltage1.csv');
            v2=csvread('voltage2.csv');

            v1bin = (v1 == 0.5);
            v2bin = (v2 == 0.5);

            %% Flag consecutive nosepokes as bites
            flag1 = zeros(size(v1bin));
            flag2 = zeros(size(v2bin));
            %sum1 = zeros(size(v1bin));
            %sum2 = zeros(size(v2bin));
            mean1 = zeros(size(v1bin));
            mean2 = zeros(size(v2bin));

    %         for i=1:length(v1bin)-1000
    %             temp1 = v1bin(i:i+1000);
    %             temp2 = v2bin(i:i+1000);
    %             sum1(i) = sum(temp1);
    %             sum2(i) = sum(temp2);
    %             if sum1(i) >= 3*vSampleRate
    %                 flag1(i:i+1000) = 1;
    %             end
    %             if sum2(i) >= 3*vSampleRate
    %                 flag2(i:i+1000) = 1;
    %             end
    %         end

            for i=1:length(v1)-vSampleRate*5
                temp1 = v1(i:i+vSampleRate*5);
                temp2 = v2(i:i+vSampleRate*5);
                mean1(i) = mean(temp1);
                mean2(i) = mean(temp2);
                if mean1(i) <= 2.5
                    flag1(i:i+vSampleRate*5) = 1;
                end
                if mean2(i) <= 2.5
                    flag2(i:i+vSampleRate*5) = 1;
                end
            end

            %% Compute number of bites and subtract from pokes
            bite1 = flag1'*v1bin;
            bite2 = flag2'*v2bin;

            fprintf('number of bite samples in COM3:\t%d\n',bite1);
            fprintf('number of bite samples in COM4:\t%d\n',bite2);
            fprintf('number of bites in COM3:\t%d\n',bite1/vSampleRate);
            fprintf('number of bites in COM4:\t%d\n',bite2/vSampleRate);

            np1 = (nnz(v1bin) - bite1)/vSampleRate + 1;
            np2 = (nnz(v2bin) - bite2)/vSampleRate + 1;

            fprintf('number of nosepokes in COM3:\t%d\n',np1);
            fprintf('number of nosepokes in COM4:\t%d\n',np2);

            %% Record nosepokes as active or passive
            load('sessiondata.mat');

            if strcmp(s1,'gong')
                newline = [j,condition(j),activeside(j),trial,np1,np2]; % write active in fifth column
            else
                newline = [j,condition(j),activeside(j),trial,np2,np1]; % write active in fifth column
            end

            cd('../'); % go back to behavioral directory
            dlmwrite('nosepokes_total.csv',newline,'delimiter',',','-append');

            disp(strcat('Wrote data from:',files(k).name));

            %% Plot bite overlap
%             cd('../');
%             clf;
%             if strcmp(s1,'gong')
%                 subplot(2,1,1); % plot active on top
%             else
%                 subplot(2,1,2);
%             end
%             hold on;
%             scatter(1:length(v1bin),v1bin,'b');
%             plot(flag1,'r');
%             plot(mean1/max(mean1),'Color',[0.4,0.4,0.4]);
%             ylim([0,1]);
%             ax = gca;
%             ax.YTick = [0,0.25,0.5,0.75,1];
%             ax.YTickLabel = cellstr(num2str([0;0.25;0.5;0.75;1]*max(mean1)))';
%             title(strcat(files(k).name,'(hole 1:COM3)'));
% 
%             if strcmp(s1,'gong')
%                 subplot(2,1,2);
%             else
%                 subplot(2,1,1);
%             end
%             hold on;
%             scatter(1:length(v2bin),v2bin,'b');
%             plot(flag2,'g');        
%             plot(mean2/max(mean2),'Color',[0.4,0.4,0.4]);
%             ylim([0,1]);
%             ax = gca;
%             ax.YTick = [0,0.25,0.5,0.75,1];
%             ax.YTickLabel = cellstr(num2str([0;0.25;0.5;0.75;1]*max(mean2)))';
%             title(strcat(files(k).name,'(hole 2:COM4)'));
% 
%             savefig(strcat(files(k).name,'/bites.fig'));
%             saveas(gcf,strcat(files(k).name,'_bites.png'));
%             disp(strcat('Saved figure in:',files(k).name));
        end
    end
end


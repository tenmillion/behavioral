% Script for plotting summary of data

%% Old data
% fuku = [7,7;12,9;12,14];
% usagi = [12,6;11,7;15,9;14,7];
% kuro = [17,9;11,7;9,5;9,8];
% hidari = [20,17;4,4;8,7;9,6];
% 
% sen = [19,12;15,12;11,7];
% tsuki = [10,11;13,7;7,1];
% migi = [17,12;19,6;24,11];
% mohi1 = [4,3;0,0];
olddata1 = csvread('../MusicLounge/soundNP_LED/nosepokes_total.csv');
olddata2 = csvread('../MusicLounge/oddballNP_LED/nosepokes_total.csv');

%% New data
% readvoltage.m creates a count of (nosepokes - bites) in nosepokes_total.csv
newdata = csvread('nosepokes_total.csv');

alldata = [olddata1;olddata2;newdata];
allLEdata = [olddata1(olddata1(:,1)~=8,:);olddata2(olddata2(:,1)~=8,:);newdata(newdata(:,1)~=8,:)]; % All data minus Mohi 1 (young long evans)

%names = cellstr(['fuku  ';'usagi ';'kuro  ';'hidari';'tsuki ';'migi  ';'sen   ';'mohi1 ']);
%newline = [j,condition(j),activeside(j),trial,np1,np2]; % write active in fifth column

%% make plots
fuku = alldata(alldata(:,1)==1,:);
usagi = alldata(alldata(:,1)==2,:);
kuro = alldata(alldata(:,1)==3,:);
hidari = alldata(alldata(:,1)==4,:);
tsuki = alldata(alldata(:,1)==5,:);
migi = alldata(alldata(:,1)==6,:);
sen = alldata(alldata(:,1)==7,:);
mohi1 = alldata(alldata(:,1)==8,:);

figure;
subplot(2,4,1); bar(fuku(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 1 (no oddball)');%legend('active','passive');
subplot(2,4,2); bar(usagi(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 2 (no oddball)');%legend('active','passive');
subplot(2,4,3); bar(kuro(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 3 (no oddball)');%legend('active','passive');
subplot(2,4,4); bar(hidari(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 4 (no oddball)');%legend('active','passive');
subplot(2,4,5); bar(tsuki(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 5 (2 oddballs)');%legend('active','passive');
subplot(2,4,6); bar(migi(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 6 (2 oddballs)');%legend('active','passive');
subplot(2,4,7); bar(sen(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 7 (2 oddballs)');%legend('active','passive');
subplot(2,4,8); bar(mohi1(:,5:6)); xlim([0,5.5]); ylim([0,30]); title('rat 8 (2 oddballs)');legend('active','inactive');

%% plot overall summary across trials
% column 5 is active hole, column 6 is inactive hole
figure;
hold on;
bar([mean(allLEdata(:,5)),mean(allLEdata(:,6))]); ylim([0,30]);
errorbar([mean(allLEdata(:,5)),mean(allLEdata(:,6))],[std(allLEdata(:,5)),std(allLEdata(:,6))],'.')
title('Excluding old Wistar');

figure;
hold on;
bar([mean(alldata(:,5)),mean(alldata(:,6))]); ylim([0,30]);
errorbar([mean(alldata(:,5)),mean(alldata(:,6))],[std(alldata(:,5)),std(alldata(:,6))],'.')
title('All rats');

c1 = alldata(alldata(:,2)==1,:);
%c1 = c1(c1(:,4)<=5,:); % only using till trial 5
c2 = alldata(alldata(:,2)==2,:);
c2LE = c2(c2(:,1)~=8,:); % leaving out mohi1

figure;
boxplot([alldata(:,5),alldata(:,6)]); ylim([0,30]);
title('All rats box plot');

figure;
boxplot([allLEdata(:,5),allLEdata(:,6)]); ylim([0,30]);
title('No old Wistar box plot');

%% plot by condition
% figure;
% boxplot([c1(:,5),c2(:,5)]); ylim([0,30]);
% title('No oddball vs. 2 oddballs (all rats)');

% c1_sorted = sortrows(c1,1); % sort based on ID
% c2_sorted = sortrows(c2LE,1); % sort based on ID
% c1_active_byrat = reshape(c1_sorted(:,5),[20/4,1*4]);
% c2_active_byrat = reshape(c2_sorted(:,5),[15/3,1*3]);
% c1_passive_byrat = reshape(c1_sorted(:,6),[20/4,1*4]);
% c2_passive_byrat = reshape(c2_sorted(:,6),[15/3,1*3]);

% c1a_p = [prctile(c1_active_byrat',90)-mean(c1_active_byrat');mean(c1_active_byrat')-prctile(c1_active_byrat',10)];
% c2a_p = [prctile(c2_active_byrat',90)-mean(c2_active_byrat');mean(c2_active_byrat')-prctile(c2_active_byrat',10)];
% c1p_p = [prctile(c1_passive_byrat',90)-mean(c1_passive_byrat');mean(c1_passive_byrat')-prctile(c1_passive_byrat',10)];
% c2p_p = [prctile(c2_passive_byrat',90)-mean(c2_passive_byrat');mean(c2_passive_byrat')-prctile(c2_passive_byrat',10)];

c1_sorted = sortrows(c1,1); % sort based on ID
c2_sorted = sortrows(c2,1); % sort based on ID
c1_active_byrat = reshape(c1_sorted(:,5),[24/4,1*4]);
c2_active_byrat = reshape(c2_sorted(:,5),[20/4,1*4]);
c1_passive_byrat = reshape(c1_sorted(:,6),[24/4,1*4]);
c2_passive_byrat = reshape(c2_sorted(:,6),[20/4,1*4]);

c1a_p = std(c1_active_byrat');
c2a_p = std(c2_active_byrat');
c1p_p = std(c1_passive_byrat');
c2p_p = std(c2_passive_byrat');

% figure;
% subplot(2,2,1);
% hold on;
% shadedErrorBar(1:5,mean(c1_active_byrat,2),c1a_p,'k-',1);
% plot(c1_active_byrat,'x--');
% ylim([0,30]);
% title('no oddball active');
% 
% subplot(2,2,3);
% hold on;
% shadedErrorBar(1:5,mean(c1_passive_byrat,2),c1p_p,'k-',1);
% plot(c1_passive_byrat,'x--');
% ylim([0,30]);
% title('no oddball inactive');

subplot(2,2,1);
hold on;
shadedErrorBar(1:6,mean(c1_active_byrat,2),c1a_p,'k-',1);
plot(c1_active_byrat,'x--');
ylim([0,30]);
title('no oddball active');

subplot(2,2,3);
hold on;
shadedErrorBar(1:6,mean(c1_passive_byrat,2),c1p_p,'k-',1);
plot(c1_passive_byrat,'x--');
ylim([0,30]);
title('no oddball inactive');

subplot(2,2,2);
hold on;
shadedErrorBar(1:5,mean(c2_active_byrat,2),c2a_p,'k-',1);
plot(c2_active_byrat,'x--');
ylim([0,30]);
title('oddball active');

subplot(2,2,4);
hold on;
shadedErrorBar(1:5,mean(c2_passive_byrat,2),c2p_p,'k-',1);
plot(c2_passive_byrat,'x--');
ylim([0,30]);
title('oddball inactive');

%% hypothesis testing
% [h1,p1,t1]=ttest2(c1(:,5),c1(:,6))
% 
% [h2,p2,t2]=ttest2(c2(:,5),c2(:,6))
% [h3,p3,t3]=ttest2(c1(:,5),c2(:,5))
% 
% [h4,p4,t4]=ttest2(c2LE(:,5),c2LE(:,6))
% [h5,p5,t5]=ttest2(c1(:,5),c2LE(:,5))

% ratio
fukuR = fuku(:,5)./fuku(:,6);
usagiR = usagi(:,5)./usagi(:,6);
kuroR = kuro(:,5)./kuro(:,6);
hidariR = hidari(:,5)./hidari(:,6);

senR = sen(:,5)./sen(:,6);
tsukiR = tsuki(:,5)./tsuki(:,6);
migiR = migi(:,5)./migi(:,6);
mohiR = mohi1(:,5)./mohi1(:,6);

figure;
subplot(2,4,1); bar(fukuR); xlim([0,5.5]); ylim([0,5]); title('rat 1 (no oddball)');%legend('active','passive');
subplot(2,4,2); bar(usagiR); xlim([0,5.5]); ylim([0,5]); title('rat 2 (no oddball)');%legend('active','passive');
subplot(2,4,3); bar(kuroR); xlim([0,5.5]); ylim([0,5]); title('rat 3 (no oddball)');%legend('active','passive');
subplot(2,4,4); bar(hidariR); xlim([0,5.5]); ylim([0,5]); title('rat 4 (no oddball)');%legend('active','passive');
subplot(2,4,5); bar(tsukiR); xlim([0,5.5]); ylim([0,5]); title('rat 5 (2 oddballs)');%legend('active','passive');
subplot(2,4,6); bar(migiR); xlim([0,5.5]); ylim([0,5]); title('rat 6 (2 oddballs)');%legend('active','passive');
subplot(2,4,7); bar(senR); xlim([0,5.5]); ylim([0,5]); title('rat 7 (2 oddballs)');%legend('active','passive');
subplot(2,4,8); bar(mohiR); xlim([0,5.5]); ylim([0,5]); title('rat 8 (2 oddballs)');legend('active:inactive');

c1R = [fukuR,usagiR,kuroR,hidariR];
c2R = [senR,tsukiR,migiR,mohiR];

figure;
boxplot([reshape(c1R(1:5,:),[1,20])',reshape(c2R,[1,20])']);
title('no oddball vs. 2 oddballs (active:inactive ratio)');
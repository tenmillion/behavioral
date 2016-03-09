fuku = [7,7;12,9;12,14];
usagi = [12,6;11,7;15,9;14,7];
kuro = [17,9;11,7;9,5;9,8];
hidari = [20,17;4,4;8,7;9,6];

sen = [19,12;15,12;11,7];
tsuki = [10,11;13,7;7,1];
migi = [17,12;19,6;24,11];
mohi1 = [4,3;0,0];

% make plots
c1 = [fuku;usagi;kuro;hidari];
c2 = [sen;tsuki;migi];

% column 1 is active hole, column 2 is inactive hole
hold on;
bar([mean(c1),mean(c2)])
errorbar([mean(c1),mean(c2)],[std(c1),std(c2)],'.')

figure;
boxplot([c1(1:9,:),c2]);
[h,p]=ttest2(c1(:,1),c1(:,2))
[h,p]=ttest2(c2(:,1),c2(:,2))
[h,p]=ttest2(c1(:,1),c2(:,1))

% ratio
fukuR = fuku(:,1)./fuku(:,2);
usagiR = usagi(:,1)./usagi(:,2);
kuroR = kuro(:,1)./kuro(:,2);
hidariR = hidari(:,1)./hidari(:,2);

senR = sen(:,1)./sen(:,2);
tsukiR = tsuki(:,1)./tsuki(:,2);
migiR = migi(:,1)./migi(:,2);

c1R = [fukuR,usagiR(2:4),kuroR(2:4),hidariR(2:4)];
c2R = [senR,tsukiR,migiR];

boxplot([reshape(c1R(:,2:4),[1,9])',reshape(c2R,[1,9])'])
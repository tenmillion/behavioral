% script to analyze data
clear;
myfolder = dir('*_NPtime.csv');

for i = 1:length(myfolder)
    data=csvread(myfolder(i).name);
    data2=data(data(:,2)>0,2);
    data1=data(data(:,1)>0,1);
    data2=data2(data2<30*60); % cutoff at 30 minutes
    data1=data1(data1<30*60); % cutoff at 30 minutes
    figure;
    hold on;
    %scatter(data1,ones(size(data1)),'b.');
    %scatter(data2,2*ones(size(data2)),'r.');
    histogram(data1,'BinWidth',10);
    histogram(data2,'BinWidth',10);
    legend('Stim 1','Stim 2');
    title(myfolder(i).name);
    ylim([0,15]);
    xlim([0,1800]);
    disp(myfolder(i).name);
    disp(size(data1));
    disp(size(data2));
    
    header = myfolder(i).name(1:findstr(myfolder(i).name,'2016')+3);
    txtfile = dir(strcat(header,'*.txt'));
    fid = fopen(txtfile.name);
    disp(txtfile.name);
    for j = 1:6
        tline = fgetl(fid);
    end
    disp(tline);
    fclose(fid);
    
    disp('-----');
    waitforbuttonpress;
end
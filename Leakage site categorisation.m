function leakage_June12(data2)
close all;clc;
%% File where leaks are calulated and the leakage points are plotted. 
% First "load" the .mat (Input_trachea/Ear/Diaphragm) files to the
% workspace using load "Input_X.mat". You have to be in the same folder as
% the files are in. 
% Input: 
% Run the above code from the command line by substituing the "data2" value with 
% (a) Ear (b) trachea or (c) dia 
% the code has explanations where possible. But please contact me if there
% is a bug or something is not clear. 
% Output: 
% (a) ROI_10 (b) ROI_20 and (c) Images. 
% (a) and (b) are stored in the Results.mat file created in the same folder
% while (c) is Nmmber of figures == mouse number ; Subplots within each
% figure == ROIS/Sheets for that mouse. 
% the labelling of the figures would be different for the diaphragm where
% there are no "mouse" excels. 
for i=1:length(data2)
    animal=data2{i,1};
    %% function for charectarizing leakage at 10pixel window. Returns the ROI1 structure that the details of leakage.  
    ROI_10(i,:)=leakage_10_pixels(animal);
    %% Using the ROI_10 details from above as input look at adjacent bins i.e. look at 20pixel distance to quantify leakage sites. 
    ROI_20(i,:)=leakage_20_pixels(ROI_10,i);
    %% Function for plotting the peaks at the leakage sites 
    plotting_peaks(ROI_20);
end
save('Results.mat',"ROI_10","ROI_20")
end

function ROI_10=leakage_10_pixels(animal)
ROI_10=struct();
for j=1:size(animal,1)
    window=[];leakage=[];
    data=animal{j,1};                            % mouse number
    roi=animal{j,2};                            % sheet number
    mouse=animal{j,3};
    a=1;i=1;
    window_size=10;                              % create a window of 10 pixels 
    while i<=(length(data))
        green_count=0;red_count=0;
        if (length(data)-i)< window_size        % gather the data in the window 
            window=data(i:end,:);
        else
            window=data(i:i+window_size,:);
        end
        green_count=sum(window(:,2)>0);         % see if there are leaks in green channel within window
        if (green_count<2)
            green_count=0;
        elseif (green_count>=2)                 % if there are 2 or more  green leaks classify this as a green count
            green_count=1;
        end
        red_count=sum(window(:,3)>0);           % same process as above for the red channel 
        if (red_count<2)
            red_count=0;
        elseif (red_count>=2)
            red_count=1;
        end
        leakage(a,1)=data(i);                   % collect the first point of the window
        leakage(a,2)=green_count;               % collect the green count within 10pixel 
        leakage(a,3)=red_count;                 % collect red count within 10pixel 
        a=a+1;                                  % increase counter to the matrix
        i=i+window_size;                        % move the window by one pixel 
    end                                         % repeat till the whole sheet is completed. 
    
    %% Structure saving all the data from 10pixel size 
    
    ROI_10.stats{j,:}=leakage;                  % leakage details 
    ROI_10.sheet{j,1}=roi;                      % sheet number  
    ROI_10.mouse{j,1}=mouse;                    % mouse number 
end
end


%% Once the 10pixel window is done (above) use that data and look at the adjacent windows i.e. 20pixels. 
function ROI_20=leakage_20_pixels(ROI_10,i)
no=i;
ROI_20=struct();
sheets=length(ROI_10(no).sheet);                % unwrapping the ROI_10 structure to get sheet info 
for j=1:sheets
    data=ROI_10(no).stats{j,1};                 % unwrap the ROI_10 struct. to get 10pixel count info. 
    count_matrix=[];
    a=1;i=1;
    leakage2=[];green_positons=[];red_positons=[];yellow_positons=[];
    green_count=0;red_count=0;yellow_count=0;green_high=[];red_high=[];yellow_high=[];
    window_size=2;                             % the window size = 2 to look at only adjacent bins (each 10pixel long).
    while i<=length(data)
        green_high=0;red_high=0;yellow_high=0;
        if (length(data)-i)< window_size
            window=data(i:end,:);
        else
            window=data(i:i+window_size-1,:);
        end
        green_high=sum(window(:,2)>0);          % get the count value of greens in both bins. 
        red_high=sum(window(:,3)>0);             % get count of red from both bins 
        if (green_high>=1) && (red_high>=1)      % if both green AND red counts are >=1 then call this a yellow leak 
            yellow_high=1;
            green_high=0;
            red_high=0;
            leakage2(a,1)=data(i,1);
            leakage2(a,2)=green_high;
            leakage2(a,3)=red_high;
            leakage2(a,4)=yellow_high;
            leakage2(a,5)=0;                    % am not sure why I added these last 3 columns now :) 
            leakage2(a,6)=0;
            leakage2(a,7)=1;
            i=i+window_size;                    % move to the next 20 bins. 
       elseif (green_high>=1) && (red_high==0)  % if only green count>=1 then call it a green leak 
           yellow_high=0;
           green_high=1;
           red_high=0;
           leakage2(a,1)=data(i,1);
           leakage2(a,2)=green_high;
           leakage2(a,3)=red_high;
           leakage2(a,4)=yellow_high;
           leakage2(a,5)=1;
           leakage2(a,6)=0;
           leakage2(a,7)=0;
           i=i+window_size;
        elseif (green_high==0) && (red_high>=1) % if only red count >=1 then call it a red leak.     
           yellow_high=0;
           green_high=0;
           red_high=1;
           leakage2(a,1)=data(i,1);
           leakage2(a,2)=green_high;
           leakage2(a,3)=red_high;
           leakage2(a,4)=yellow_high;
           leakage2(a,5)=0;
           leakage2(a,6)=1;
           leakage2(a,7)=0;
           i=i+window_size;
        elseif (green_high==0) && (red_high==0) % last option is nobody leaks. 
            yellow_high=0;
            green_high=0;
            red_high=0;
            leakage2(a,1)=data(i,1);
            leakage2(a,2)=green_high;
            leakage2(a,3)=red_high;
            leakage2(a,4)=yellow_high;
            leakage2(a,5)=0;
            leakage2(a,6)=0;
            leakage2(a,7)=0;
            i=i+1;
        end    
        if i>=length(data)                  % break and quit the loop when you reach the end. 
            break
        end
        a=a+1;
     end
    green_count=sum(leakage2(:,2));         % Total count of green leaks 
    red_count=sum(leakage2(:,3));           % Totla red count 
    yellow_count=sum(leakage2(:,4));        % Total yellow count
    count_matrix=[green_count,red_count,yellow_count];
    ROI_20.count{j,1}=count_matrix;         
    ROI_20.leakage{j,1}=leakage2;
    total_length=leakage2(end,1);
    leak_length=(count_matrix./total_length)*100;    
    ROI_20.pervessel{j,:}=leak_length;      % Calculate average leak similar to the last row in excel.
%    figure(no),subplot(2,ceil(sheets/2),j),stem(ROI_20.leakage{j, 1}(:,[1,5:end])), hold on
 end
end


function plotting_peaks(ROI_20)

for i=1:size(ROI_20,1)                   % number of mice
    yellow_leaks=[];red_leaks=[];green_leaks=[];
    for j=1:size(ROI_20(i).leakage,1)            % number of ROIS per mice 
        n=2;m=round(size(ROI_20(i).leakage,1)/2);
        yellow_leaks(:,2)=ROI_20(i).leakage{1,1}((find(ROI_20(i).leakage{1,1}(:,4)==1)),1);
        yellow_leaks(:,1)=1;
        red_leaks(:,2)=ROI_20(1).leakage{1,1}((find(ROI_20(i).leakage{1,1}(:,3)==1)),1);
        red_leaks(:,1)=1;
        green_leaks(:,2)=ROI_20(1).leakage{1,1}((find(ROI_20(i).leakage{1,1}(:,2)==1)),1);
        green_leaks(:,1)=1;
        figure(i);
        subplot(n,m,j);stem(yellow_leaks(:,2),yellow_leaks(:,1),'Color','y');hold on
        subplot(n,m,j);stem(red_leaks(:,2),red_leaks(:,1),'Color','r');hold on
        subplot(n,m,j);stem(green_leaks(:,2),green_leaks(:,1),'Color','g');hold on
    end
end

end










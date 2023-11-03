close all

%Updated Script which now runs on a folder of images and develops radius
%calculations for all of the droplets in those folders and presents a
%histogram of all the radi, it also saves images showing the circles it has
%found in the CircleCheck Folder, the raw data it uses to create the
%histogram is seperated by image in teh RawData folder, and is concatenated
%in the AllRawData file, Individual histograms for each image are included
%in the individual histogram folder
%Alex


%This script was written for the use of droplet & circle identification
%feel free to use it as needed, if you do have questions though, MATLAB
%help might be the best option to start with, search "imfindcircles" as a
%starting point, or let me know if there are issues and we can def improve
%upon it
%SCLP

%%
SRadius = 25;   % The minimum droplet/microsphere radius (in pixels) the script will seek out (and recognize)
LRadius = 35;   % The maximum droplet/microsphere radius (in pixels) the script will seek out (and recognize)
%Generally you want to know relative pixel size, if your radius range is
%too big accuracy for the script goes down. However, if you do not know the
%approximate pixel size, try a run from 10 to 1000 and identify approximate
%size and then update accordingly
% 
Counts = [];%array which keeps track of # of beads for later use

SensitivityFactor = 0.98;   % The sensitivity of the circle-hunting command.  Lower values will 
                            % lead the script to be more selective, while higher values will 
                            % lead to it being less selective.


ImgsToAnalyze = '20231003_DropletGen_DP'; %Replace with the name of the folder with the gel trail photos to be analyzed

FolderInfo = dir(ImgsToAnalyze);
RawDataFolder = append(ImgsToAnalyze,'RawData');
CircleCheckFolder = append(ImgsToAnalyze,'CircleCheck');
HistogramFolder = append(ImgsToAnalyze,'IndividualHistograms');

%if(!)
mkdir(RawDataFolder);
mkdir(CircleCheckFolder);
mkdir(HistogramFolder);

A = zeros([3 length(FolderInfo)-2]);
MinuteCount = zeros(length(FolderInfo)-2);

for c = 3:length(FolderInfo) %loop portion
%c = 3; %tmp c to run on one image

stripFileExtention = FolderInfo(c).name(1:end-4); %creates a file name without the extension to be used late
% to save respective files



RawImage = imread(append(ImgsToAnalyze,'/',FolderInfo(c).name));  % Change value to the filename you want to analyze, also if not in the same folder index appropriately
% this will require you to write out the path directory to the specific subfolder and file

%will crop image if extra on top/bottom
[height,width] = size(RawImage);
xcrop=0;
ycrop=30;
RawImage = imcrop(RawImage,[xcrop ycrop width-2*xcrop height-2*ycrop]);

% figure, imshow(RawImage);% Will show you your image, and contains data we will analyze, You can technically comment this out, or else you will have redundant image

B = RawImage; % This is all your information for analysis
LofI = length(B); % Determining the length of the matrix, to know dimensions (normally all images are squares so one dimension is fine)
%%
%Isolating the data for your analysis later on

red = B(:,:,1); % Red channel data
green = B(:,:,2); % Green channel data
blue = B(:,:,3); % Blue channel data
a = zeros(size(B, 1), size(B, 2)); %Makes an all zero matrix to blot out channels that are not of interest, i.e they become 0
%%
%Producing matrices that are only composed of a single channel's worth of data

just_red = cat(3, red, a, a);     %Creates a new set of data only composed of red channel data
just_green = cat(3, a, green, a); %Creates a new set of data only composed of green channel data
just_blue = cat(3, a, a, blue);  %Creates a new set of data only composed of blue channel data
back_to_original_img = cat(3, red, green, blue); %Back to your original image such that you can evaluate
%%
%This portion below is only necessary to visualize, the different versions
%of the image you may have created, names should be descriptive enough to use
%I would suggest to comment most or all of them out if you are not using them

% figure, imshow(B), title('Original image')
% figure, imshow(just_red), title('Red channel')
% figure, imshow(just_green), title('Green channel')
% figure, imshow(just_blue), title('Blue channel')
% figure, imshow(back_to_original_img), title('Back to original image')
%%




%Preparing image analysis and output

L = bwlabel(red); %labels connected components in a 2D array/
s = regionprops(L, 'PixelIdxList', 'PixelList', 'Area', 'Centroid', 'FilledArea'); %measures properties of image (i.e. property definitions) should not need to change

f1=figure('visible','off'); 
imshow(RawImage); %Produces an image 
[centersBright, radiiBright] = imfindcircles(just_red,[SRadius LRadius],'ObjectPolarity','bright','Sensitivity',SensitivityFactor);%Algorithm that identifies the circles
%pause(1);
viscircles(centersBright, radiiBright,'EdgeColor','r') %places circles around your objects, should help visualize the specific beads or droplets identified

set(f1,'visible','on');
savefig(append(CircleCheckFolder,'/',stripFileExtention,'CircleCheck.fig'))
close(f1)

%%
%If using the Celestron microscope images make sure to put in the
%appropriate pixel to micron conversion, same goes for nikon, etc... 
%multipliers are in micron/pixel

%Celestron 
%For 4X use 1.6447 as multiplier
%For 10X use 0.6964 as multiplier

%Nikon
%For 10x use 1.24 as multiplier
%For 20x use 0.62 as multiplier
%

%EVOS LP Lab
%2X objective = 1302 pixels per 5618um
%4X objective = 1216 pixels per 2620um
%10X objective = 1724 pixels per 1500um
%10X multiplier would be 1500um/1724pix
%20X objective = 1793 pixels per 800um\

%Stereoscope LP Lab
% 1100um per 940pixels = 1.16 multiplier
% 1100um per 570pixels = ___ multiplier

DiaBrightAdj = radiiBright*(2000/(height-60))*2; %Gives you diameter of the beads/droplets... NEED TO MAKE SURE CORRECT MULTIPLIER IS USED  

edges = linspace(0,150,151);%creates the bins from 0um to 150um of size 1um
f2=figure('visible','off');histogram(DiaBrightAdj,'BinEdges',edges); % if want auto bins replace w/ figure,histogram(\\DiaBirghtAdj)

title('Droplet size distribution','FontSize', 14) % title of your plot
xlabel('Diameter (um)', 'FontSize', 14) % x-axis label of your plot
ylabel('Number of beads', 'FontSize', 14) % y-axis label of your plot
countDia = size(DiaBrightAdj,1); %Data that may be useful to have present
DiaMean = mean(DiaBrightAdj);    %Data that may be useful to have present
DiaSTD = std(DiaBrightAdj);      %Data that may be useful to have presen
%%Information that will be placed on the histogram, useful if you want to
%%visualize output, if you want to eliminate comment out
%A =[countDia; DiaMean; DiaSTD];




Counts = [Counts; countDia];%keeps track of #s of beads

A(1,c-2) = countDia;% sets an array w the count, mean, and std
A(2,c-2) = DiaMean;
A(3,c-2) = DiaSTD;

grid on;%adds grid lines to histogram
xticks(0:10:150)%makes the x axis tick marks every 10um


DataInfo= {'Samples Counted: %.3f';...
           'Mean Diameter: %.3f';...
           'St.Dev: %.3f'};
str = compose(DataInfo, A(1:3,c-2));
dim = [0.55 0.6 0.3 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');%formats and displays count,mean,std in box



set(f2,'visible','on');
savefig(append(HistogramFolder,'/',stripFileExtention,'histogram.fig'))%saves the histogram to the individualhistogram folder
close(f2)


writematrix(DiaBrightAdj,append(RawDataFolder,'/',stripFileExtention,'RawData'))%writes the raw data to the raw data folder

end
%%


%StartTime = 60*MinuteCount(1,1)+MinuteCount(1,2);
%for c = 3:length(FolderInfo) %loop portion
%    MinuteCount(c-2,3) = 60*MinuteCount(c-2,1)+MinuteCount(c-2,2)-StartTime;
%end
%x = MinuteCount(:,3);

%x = A(2,1:length(FolderInfo)-2);
%y = A(1,1:length(FolderInfo)-2);

%y = A(2,1:length(FolderInfo)-2);
%scatter(x,y)


%scatter(x,y)

%sd_vct = [A(3,:)];
%errorbar(x,y,sd_vct)



DataInfo = dir(RawDataFolder);%gets all the names of the raw data files
M = [];%matrix of all raw data in list
length(Counts)
disp(max(Counts))
TotDataCols =zeros(max(Counts),length(Counts));
for d = 3:length(DataInfo)
    TotDataCols(1:Counts(d-2),d-2) = readmatrix(append(RawDataFolder,'/',DataInfo(d).name));
    M = [M; readmatrix(append(RawDataFolder,'/',DataInfo(d).name))];

end

writematrix(M,append(ImgsToAnalyze,'_AllRawDataList'));
writematrix(TotDataCols,append(ImgsToAnalyze,'AllRawDataColSeperated'))

edges = linspace(0,150,151);%creates the bins from 0 to 150 of size 1
figure,histogram(M,'BinEdges',edges)% if want auto bins replace w/ figure,histogram(M)

title('Droplet size distribution','FontSize', 14) % title of your plot
xlabel('Diameter (um)', 'FontSize', 14) % x-axis label of your plot
ylabel('Number of beads', 'FontSize', 14) % y-axis label of your plot

countDia = size(M,1); %Data that may be useful to have present
DiaMean = mean(M);    %Data that may be useful to have present
DiaSTD = std(M);     %Data that may be useful to have present

grid on;
xticks(0:10:150)

A(1,c-2) = countDia;
A(2,c-2) = DiaMean;
A(3,c-2) = DiaSTD;

DataInfo= {'Samples Counted: %.0f';...
           'Mean Diameter: %.3f';...
           'St.Dev: %.3f'};
str = compose(DataInfo, A(1:3,c-2));
dim = [0.55 0.6 0.3 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

savefig(append(ImgsToAnalyze,'_Histogram.fig'))










%% Just starting point for myself
% fwrite(radiiBright)%
% 
% SensitivityFactor = 0.92;   % The sensitivity of the circle-hunting command.  Lower values will 
%                             % lead the script to be more selective, while higher values will 
%                             % lead to it being less selective.
% 
% RawImage = imread('Droplets_TEMP\4min-3.jpg');  % Change value to the filename you want to analyze
% 
% figure, imshow(RawImage)
% 
% B = RawImage;
% 
% L = bwlabel(B);
% s = regionprops(L, 'PixelIdxList', 'PixelList', 'Area', 'Centroid', 'FilledArea');
% 
% figure, imshow(RawImage)
% [centersBright, radiiBright] = imfindcircles(B,[SRadius LRadius],'ObjectPolarity','dark','Sensitivity',SensitivityFactor);
% viscircles(centersBright, radiiBright,'EdgeColor','b');

close all


%This script was written for the use of droplet & circle identification
%feel free to use it as needed, if you do have questions though, MATLAB
%help might be the best option to start with, search "imfindcircles" as a
%starting point, or let me know if there are issues and we can def improve
%upon it
%SCLP

%%
SRadius = 15;   % The minimum droplet/microsphere radius (in pixels) the script will seek out (and recognize)
LRadius = 30;   % The maximum droplet/microsphere radius (in pixels) the script will seek out (and recognize)
%Generally you want to know relative pixel size, if your radius range is
%too big accuracy for the script goes down. However, if you do not know the
%approximate pixel size, try a run from 10 to 1000 and identify approximate
%size and then update accordingly
% 

SensitivityFactor = 0.98;   % The sensitivity of the circle-hunting command.  Lower values will 
                            % lead the script to be more selective, while higher values will 
                            % lead to it being less selective.


ImgsToAnalyze = 'ImageToProcess'; %Replace with the name of the folder with the gel trail photos to be analyzed

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

%for c = 3:length(FolderInfo) %loop portion
c = 3; %tmp c to run on one image

stripFileExtention = FolderInfo(c).name(1:end-4); %creates a file name without the extension to be used late
% to save respective files

%disp(append(ImgsToAnalyze,'/',FolderInfo(c).name))

RawImage = imread(append(ImgsToAnalyze,'/',FolderInfo(c).name));  % Change value to the filename you want to analyze, also if not in the same folder index appropriately
% this will require you to write out the path directory to the specific subfolder and file

% figure, imshow(RawImage);% Will show you your image, and contains data we will analyze, You can technically comment this out, or else you will have redundant image

%finds the lines in the image to crop out the edges
BWImage = rgb2gray(RawImage);
BW = edge(BWImage, 'sobel');
[H,theta,rho] = hough(BWImage);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,theta,rho,P,'FillGap',40,'MinLength',500);
disp(lines);
disp(theta);

figure, imshow(BW), hold on
max_len = 0;
for k = 1:length(lines)
  xy = [lines(k).point1; lines(k).point2];
  if abs(xy(2,2) - xy(1,2)) <=1
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

    % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    
    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
     if ( len > max_len)
        max_len = len;
        xy_long = xy;
     end
   end
end

% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');

%{
%this will crop the RawImage to remove borders around droplet channel

CroppedImage = imcrop(RawImage,[400 400 3240 1560]);
%CroppedImage = rgb2gray(CroppedImage);

b = imsharpen(CroppedImage);
imshow(b)

B = CroppedImage; % This is all your information for analysis
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

% This line is for something I doing, don't worry about it % imwrite(RawImage, 'test.tif'); % C = imfinfo('test.tif');


%Preparing image analysis and output

L = bwlabel(red); %labels connected components in a 2D array/
s = regionprops(L, 'PixelIdxList', 'PixelList', 'Area', 'Centroid', 'FilledArea'); %measures properties of image (i.e. property definitions) should not need to change

f=figure('visible','off'); 
imshow(CroppedImage); %Produces an image 
[centersBright, radiiBright] = imfindcircles(just_red,[SRadius LRadius],'ObjectPolarity','bright','Sensitivity',SensitivityFactor);%Algorithm that identifies the circles
viscircles(centersBright, radiiBright,'EdgeColor','r') %places circles around your objects, should help visualize the specific beads or droplets identified

set(f,'visible','on');
savefig(append(CircleCheckFolder,'/',stripFileExtention,'CircleCheck.fig'))
close(f)

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

DiaBrightAdj = radiiBright*(250/230)*2; %Gives you diameter of the beads/droplets... NEED TO MAKE SURE CORRECT MULTIPLIER IS USED  
f=figure('visible','off');histogram(DiaBrightAdj); %Gives you a histogram of the data

title('Droplet size distribution','FontSize', 14) % title of your plot
xlabel('Diameter (um)', 'FontSize', 14) % x-axis label of your plot
ylabel('Number of beads', 'FontSize', 14) % y-axis label of your plot
countDia = size(DiaBrightAdj,1); %Data that may be useful to have present
DiaMean = mean(DiaBrightAdj);    %Data that may be useful to have present
DiaSTD = std(DiaBrightAdj);      %Data that may be useful to have presen
%%Information that will be placed on the histogram, useful if you want to
%%visualize output, if you want to eliminate comment out
%A =[countDia; DiaMean; DiaSTD];

A(1,c-2) = countDia;
A(2,c-2) = DiaMean;
A(3,c-2) = DiaSTD;

%MinuteCount(c-2,1) = str2num(FolderInfo(c).name(10:end-8));
%MinuteCount(c-2,2) = str2num(FolderInfo(c).name(12:end-6));

DataInfo= {'Samples Counted: %.3f';...
           'Mean Diameter: %.3f';...
           'St.Dev: %.3f'};
str = compose(DataInfo, A(1:3,c-2));
dim = [0.55 0.6 0.3 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');



set(f,'visible','on');
savefig(append(HistogramFolder,'/',stripFileExtention,'histogram.fig'))
close(f)

writematrix(DiaBrightAdj,append(RawDataFolder,'/',stripFileExtention,'RawData'))

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



%DataInfo = dir(RawDataFolder);
%M = [];
%for d = 3:length(DataInfo)
%    M = [M; readmatrix(append(RawDataFolder,'/',DataInfo(d).name))];

%end

%writematrix(M,append(ImgsToAnalyze,'AllRawData'));
%figure,histogram(M)

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
%}
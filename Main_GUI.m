function varargout = Main_GUI(varargin)
% MAIN_GUI MATLAB code for Main_GUI.fig
%      MAIN_GUI, by itself, creates a new MAIN_GUI or raises the existing
%      singleton*.
%
%      H = XSQMAIN_GUI returns the handle to a new MAIN_GUI or the handle to
%      the existing singleton*.
%
%      MAIN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_GUI.M with the given input arguments.
%
%      MAIN_GUI('Property','Value',...) creates a new MAIN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main_GUI

% Last Modified by GUIDE v2.5 14-Feb-2020 14:45:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Main_GUI is made visible.
function Main_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main_GUI (see VARARGIN)

% Choose default command line output for Main_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Main_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

addpath('subfn\');

% --- Outputs from this function are returned to the command line.
function varargout = Main_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % -- Performance Measures -- %
    
    global Labeled_Image Merged_Regions 
    % -- Performance Measures -- %
    
    Actual = Labeled_Image;
    Predicted = Merged_Regions;
    
    Performance = perf(Actual,Predicted);
    
    Precision = Performance.precision;
    Recall = Performance.recall;
    F_Measure = Performance.Fmeasure;
    
    rnames = {};
    cnames = {'Precision','Recall','F-Measure'};
%     f=figure('Name','Performance Measures','NumberTitle','off');

set(handles.uitable1,'data',[Precision,Recall,F_Measure],'ColumnName',cnames,'RowName',rnames);
%     t = uitable('Parent',f,'Data',[Precision,Recall,F_Measure],'ColumnName',cnames,'RowName',rnames);


axes(handles.axes5);
    
    bar([Precision*100,Recall*100,F_Measure*100],0.5);
    hold on;
    plot([Precision*100,Recall*100,F_Measure*100],'-r*');
    hold off;
set(gca,'XTickLabel',{'Prec','Recall','F-Measure'});
title('Performance Graph');
grid on;

% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% -- Box Feature Matching -- %
global SIFTimg row col input
global Merged_Regions Labeled_Image

% -- Box Feature Matching -- %

count = 1;
for i = 1:30:row
    for j = 1:30:col
        SIFTimage{count} = SIFTimg(i:i+29,j:j+29,1:3);
        count = count + 1;
    end
end
[size1,size2,size3] = size(SIFTimage{1});
for i = 2:(count-1)
    
    tSIFTimage1 = SIFTimage{i};
    tSIFTimage2 = SIFTimage{i-1};
    cnt1 = 1;
    cnt2 = 1;
    for ii = 1:size1
        for jj = 1:size2
            if tSIFTimage1(ii,jj,1:3) == 255;
                xkeypoint1(cnt1) = ii;
                ykeypoint1(cnt1) = jj;
                cnt1 = cnt1+1;
            end
            if tSIFTimage2(ii,jj,1:3) == 255;
                xkeypoint2(cnt2) = ii;
                ykeypoint2(cnt2) = jj;
                cnt2 = cnt2+1;
            end
        end
    end
    
end

for i = 1:length(xkeypoint1)
    xa = xkeypoint1(i);
    xb = xkeypoint2(i);
    ya = ykeypoint1(i);
    yb = ykeypoint2(i);
    d_fa_fb(i) = ( (xa - xb) + (ya - yb) )^(1/2);
end


% -- Labeled Feature Points -- %

dataset = dir('Dataset');
for itr = 3:length(dataset)
    filename  = dataset(itr).name;
    sampinput = imread(['Dataset\',filename]);
    sampinput = imresize(sampinput,[300,300]);
    Train(itr-2,:) = mean(sampinput(:));
end

Test = mean(input(:));

for i = 1:length(Train)
    xa = mean(Train(i,:));
    xb = mean(Test);
    d_fafb(i) = ( abs(xa - xb) )^(1/2);
end

Label(1:4) = 1;
Label(5:8) = 2;
% Label(23:33) = 3;
% Label(34:44) = 4;
% Label(45:55) = 5;
% Label(56:66) = 6;
% Label(67:77) = 7;
% Label(78:88) = 8;
% Label(89:99) = 9;
% Label(100:110) = 10;

Distance = d_fafb;
[val1,loc1] = find(Distance == 0);
class = Label(loc1(1));

Original_Image = imread(['samp\',num2str(class),'.jpg']);
Original_Image = imresize(Original_Image,[300,300]);

R1 = input(:,:,1);
G1 = input(:,:,2);
B1 = input(:,:,3);

R2 = Original_Image(:,:,1);
G2 = Original_Image(:,:,2);
B2 = Original_Image(:,:,3);

count = 1;
for i = 1:10:row
    for j = 1:10:col
        image1{count} = input(i:i+9,j:j+9,1:3);
        image2{count} = input(i:i+9,j:j+9,1:3);
        count = count + 1;
    end
end

Theta = [45,90,135,180,225,270,315,360];
TRsim = 15;

for i = 1:(count-1)
    for j = 1:(count-1)
        
        diffval_1 = image1{i};
        diffval_2 = image2{j};
        sdiffval = diffval_1 - diffval_2;
        diffval = rgb2gray(sdiffval);
        differeceval = diffval<=TRsim;
        
        if unique(diffval(:)) == 0
            detectregion = ones(10,10);
        else
            detectregion = diffval;
        end
        
        finaldetectregion{i} = detectregion;
        
    end
end

findetectregion = zeros(row,col);

cnt = 1;

for i = 1:10:row
    for j = 1:10:col
        findetectregion(i:i+9,j:j+9) = finaldetectregion{cnt};
        cnt = cnt + 1;
    end
end


tLabeled_Image = rgb2gray(input) - rgb2gray(Original_Image);
Labeled_Image = tLabeled_Image > 0;

figure('Name','Labeled Feature Points','NumberTitle','Off');
imshow(Labeled_Image);
axis off;
title('Labeled Feature Points','fontname','Times New Roman','fontsize',12);

if unique(Labeled_Image(:)) == 0
    
    msgbox('The Input Image Is Original');
    
else
    
    
    % -- Merge Regions -- %
    
    se = strel('disk',5);
    Merged_Regions = imopen(Labeled_Image,se);
    figure('Name','Merged Regions','NumberTitle','Off');
    imshow(Merged_Regions);
    axis off;
    title('Merged Regions','fontname','Times New Roman','fontsize',12);
    
    boundaries = bwboundaries(Merged_Regions);
    numberOfBoundaries = length(boundaries);
    figure('Name','Detected Forged Region','NumberTitle','Off');
    imshow(input);axis off;hold on;
    for i1 = 1:numberOfBoundaries
        thisBoundary = boundaries{i1};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'color' , 'r' , 'LineWidth', 2);
    end
    title('Detected Forged Region','fontname','Times New Roman','fontsize',12);
    
    for ii = 1:size(input,1);
        for jj = 1:size(input,2)
            if Merged_Regions(ii,jj) == 0
                Detect_Img(ii,jj,1) = input(ii,jj,1);
                Detect_Img(ii,jj,2) = input(ii,jj,2);
                Detect_Img(ii,jj,3) = input(ii,jj,3);
            elseif Merged_Regions(ii,jj) == 1
                Detect_Img(ii,jj,1) = 200;
                Detect_Img(ii,jj,2) = 100;
                Detect_Img(ii,jj,3) = 50;
            end
        end
        
    end
    figure('Name','Detected Image');
    imshow(uint8(Detect_Img));
    title('Tampered Detected')
end

% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Merged_Regions input

    boundaries = bwboundaries(Merged_Regions);
    numberOfBoundaries = length(boundaries);
    
    axes(handles.axes3);
    imshow(input);
    axis off;
    hold on;
    
    for i1 = 1:numberOfBoundaries
        thisBoundary = boundaries{i1};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'color' , 'r' , 'LineWidth', 2);
    end
    title('Detected Forged Region','fontname',...
        'Times New Roman','fontsize',12);
    
    for ii = 1:size(input,1);
        for jj = 1:size(input,2)
            
            if Merged_Regions(ii,jj) == 0
                
                Detect_Img(ii,jj,1) = input(ii,jj,1);
                Detect_Img(ii,jj,2) = input(ii,jj,2);
                Detect_Img(ii,jj,3) = input(ii,jj,3);
                
            elseif Merged_Regions(ii,jj) == 1
                
                Detect_Img(ii,jj,1) = 250;
                Detect_Img(ii,jj,2) = 60;
                Detect_Img(ii,jj,3) = 50;
                
            end
        end
        
    end
    
    axes(handles.axes4);

    imshow(uint8(Detect_Img));
    title('Tampered Detected')
    

% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global input row col SIFTimg SIFTpoints SLICIMG

[row,col,cha] = size(input);

    % -- Box Segmentation -- %
    
    [CA{1}, CH{1}, CV{1}, CD{1}] = dwt2(input,'db1');
    [CA{2}, CH{2}, CV{2}, CD{2}] = dwt2(input,'db1');
    [CA{3}, CH{3}, CV{3}, CD{3}] = dwt2(input,'db1');
    [CA{4}, CH{4}, CV{4}, CD{4}] = dwt2(input,'db1');
    CA4 = CA{4};
    ELF = sum(sum(abs(CA4(:))));
    EHF = [];
    for i = 1:length(CA)
        CDi = CD{i};
        CHi = CH{i};
        CVi = CV{i};
        tEHF = sum(sum(CDi(:))) + sum(sum(CHi(:))) + sum(sum(CVi(:)));
    end
    EHF = sum(tEHF);
    PLF = ( ELF/(ELF+EHF) ) * 100;
    M = row;
    N = col;
    if PLF > 50
        S = (0.02*M*N)^(1/2);
    elseif PLF <= 50
        S = (0.01*M*N)^(1/2);
    end
    
    [l, Am, Sp, d] = slic(input, round(S) , 10, 1, 'median');
    
    SLICIMG  = drawregionboundaries(l, input, [255 0 0]);
axes(handles.axes3);
imshow(SLICIMG);axis off;   % Showing SLIC Image
    title('SLIC Image','fontname','Times New Roman','fontsize',12);

    % -- Box Feature Extraction -- %
% SIFT

HarrisThresh = 10;
k = 0.055;NOfWindows = 4;
BorderDistance = 4*NOfWindows;
NBHOOD = ones(3);
sigma_nmbr = 9;
[x y z]=size(SLICIMG);
if z==3
    img1gr = rgb2gray(SLICIMG);
else
    img1gr = SLICIMG;
end

SIFTpoints = SIFT( img1gr, NBHOOD, BorderDistance, HarrisThresh, k, sigma_nmbr );

SIFTxcoord = SIFTpoints(:,1);
SIFTycoord = SIFTpoints(:,2);
axes(handles.axes4);
imshow(SLICIMG);axis off;hold on;
plot(SIFTycoord,SIFTxcoord,'or')
title('SIFT Features Extracted','fontname','Times New Roman','fontsize',12);

SIFTimg = SLICIMG;
for i = 1:length(SIFTxcoord)
    SIFTimg(SIFTxcoord(i),SIFTycoord(i),1:3) = 255;
end


% -- Box Feature Extraction -- %
% SIFT

HarrisThresh = 10;
k = 0.055;NOfWindows = 4;
BorderDistance = 4*NOfWindows;
NBHOOD = ones(3);
sigma_nmbr = 9;
[x y z]=size(SLICIMG);
if z==3
    img1gr = rgb2gray(SLICIMG);
else
    img1gr = SLICIMG;
end

SIFTpoints = SIFT( img1gr, NBHOOD, BorderDistance, HarrisThresh, k, sigma_nmbr );

SIFTxcoord = SIFTpoints(:,1);
SIFTycoord = SIFTpoints(:,2);

% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SIFTpoints SLICIMG LBPval

SEG = rgb2gray(SLICIMG);

filtDims = [2 3];
% inImg = rgb2gray(SEG);
inImg = SEG;
imgSize=size(inImg);
filtDims=filtDims+1-mod(filtDims,2);
filt=zeros(filtDims);
nNeigh=numel(filt)-1;

iHelix=snailMatIndex(filtDims);
filtCenter=ceil((nNeigh+1)/2);
iNeight=iHelix(iHelix~=filtCenter);
filt(filtCenter)=1;
filt(iNeight(1))=-1;
sumLBP=zeros(imgSize);

for i=1:length(iNeight)
    currNieghDiff=filter2(filt, inImg, 'same');
    sumLBP=sumLBP+2^(i-1)*(currNieghDiff>0);
    
    if i<length(iNeight)
        
        filt( iNeight(i) )=0;
        filt( iNeight(i+1) )=-1;
        
    end
    
end

LBPimg=sumLBP;

filtDimsR=floor(filtDims/2);
iNeight(iNeight>filtCenter)=iNeight(iNeight>filtCenter)-1;

zeroPadRows=zeros(filtDimsR(1), imgSize(2));
zeroPadCols=zeros(imgSize(1)+2*filtDimsR(1), filtDimsR(2));

inImg=cat(1, zeroPadRows, inImg, zeroPadRows);
inImg=cat(2, zeroPadCols, inImg, zeroPadCols);

imgSize=size(inImg);

neighMat=true(filtDims);

neighMat( floor(nNeigh/2)+1 )=false;
weightVec= (2.^( (1:nNeigh)-1 ));

LBPimg=zeros(imgSize);

for iRow=( filtDimsR(1)+1 ):( imgSize(1)-filtDimsR(1) )
    for iCol=( filtDimsR(2)+1 ):( imgSize(2)-filtDimsR(2) )
        
        subImg=inImg(iRow+(-filtDimsR(1):filtDimsR(1)), iCol+(-filtDimsR(2):filtDimsR(2)));
        
        diffVec=repmat(inImg(iRow, iCol), [nNeigh, 1])-subImg(neighMat);
        LBPimg(iRow, iCol)= weightVec*(diffVec(iNeight)>0);
        
    end
end


LBPimg = LBPimg(( filtDimsR(1)+1 ):( end-filtDimsR(1) ),( filtDimsR(2)+1 ):( end-filtDimsR(2) ));

LBPfeature = mean(LBPimg);

LBPval = LBPfeature(1,:);

% figure,
% 
% td = uitable('data',LBPval);

set(handles.uitable1,'data',LBPval);

SIFTxcoord = SIFTpoints(:,1);
SIFTycoord = SIFTpoints(:,2);

axes(handles.axes2);
figure('Name','SIFT Features Extracted','NumberTitle','Off');
imshow(SLICIMG);axis off;hold on;
plot(SIFTycoord,SIFTxcoord,'or')
title('SIFT Features Extracted','fontname','Times New Roman','fontsize',12);

SIFTimg = SLICIMG;
for i = 1:length(SIFTxcoord)
    SIFTimg(SIFTxcoord(i),SIFTycoord(i),1:3) = 255;
end

% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global input 

input = imresize(input,[300,300]);

axes(handles.axes2);
imshow(input);
axis off;
title('Resized Image','fontsize',12,...
    'fontname','Times New Roman','color','Black');

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename pathname input

[ filename , pathname ]=uigetfile('Dataset\*.jpg','Select an Image');
input = imread([pathname,filename]);
input = imresize(input,[300,300]);
axes(handles.axes1);
imshow(input);
axis off;
title('Input Image','fontsize',12,'fontname','Times New Roman','color','Black');

[row,col,cha] = size(input);



% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear all
close all
Main_GUI


% --------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global LBPval

load Trainfea

Testfea = LBPval;

Labels = [1 2 2 2 1 2 2 2];

Class = knnclassify(Testfea,Trainfea,Labels);

if Class == 1
    msgbox('Non-Tampered Image');
    
elseif Class == 2
    msgbox('Tampered Image');
    
end

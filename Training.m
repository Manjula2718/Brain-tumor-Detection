for ijk = 1:8
    
    input = imread(['Dataset\',num2str(ijk),'.jpg']);
    input = imresize(input,[300,300]);
    % axes(handles.axes1);
    % imshow(input);
    % axis off;
    % title('Input Image','fontsize',12,'fontname','Times New Roman','color','Black');
    
    [row,col,cha] = size(input);
    
    
    input = imresize(input,[300,300]);
    
    % axes(handles.axes2);
    % imshow(input);
    % axis off;
    % title('Resized Image','fontsize',12,...
    %     'fontname','Times New Roman','color','Black');
    
    
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
    % axes(handles.axes3);
    % imshow(SLICIMG);axis off;   % Showing SLIC Image
    % title('SLIC Image','fontname','Times New Roman','fontsize',12);
    
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
    
    % SIFTpoints = SIFT( img1gr, NBHOOD, BorderDistance, HarrisThresh, k, sigma_nmbr );
    %
    % SIFTxcoord = SIFTpoints(:,1);
    % SIFTycoord = SIFTpoints(:,2);
    % axes(handles.axes4);
    % imshow(SLICIMG);axis off;hold on;
    % plot(SIFTycoord,SIFTxcoord,'or')
    % title('SIFT Features Extracted','fontname','Times New Roman','fontsize',12);
    
    % SIFTimg = SLICIMG;
    % for i = 1:length(SIFTxcoord)
    %     SIFTimg(SIFTxcoord(i),SIFTycoord(i),1:3) = 255;
    % end
    %
    %
    % % -- Box Feature Extraction -- %
    % % SIFT
    %
    % HarrisThresh = 10;
    % k = 0.055;NOfWindows = 4;
    % BorderDistance = 4*NOfWindows;
    % NBHOOD = ones(3);
    % sigma_nmbr = 9;
    % [x y z]=size(SLICIMG);
    % if z==3
    %     img1gr = rgb2gray(SLICIMG);
    % else
    %     img1gr = SLICIMG;
    % end
    %
    % SIFTpoints = SIFT( img1gr, NBHOOD, BorderDistance, HarrisThresh, k, sigma_nmbr );
    %
    % SIFTxcoord = SIFTpoints(:,1);
    % SIFTycoord = SIFTpoints(:,2);
    
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
    
    LBPimg = sumLBP;
    
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
    
    LBPimg = zeros(imgSize);
    
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
    
    Trainfea(ijk,:) = LBPval;
    
    
    
    close all
    clc
    ijk
end

save Trainfea Trainfea
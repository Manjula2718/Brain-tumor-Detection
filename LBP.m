function LBP= LBP(inImg, filtDims, isEfficent)
if nargin<3
   isEfficent=true;
   if nargin<2
      filtDims=[3,3];
      if nargin==0
         error('Input image matrix/file name is missing.')
      end
   end
end

if ischar(inImg) && exist(inImg, 'file')==2 % In case of file name input- read graphical file
   inImg=imread(inImg);
end

if size(inImg, 3)==3
   inImg=rgb2gray(inImg);
end

inImgType=class(inImg);
isDoubleInput=strcmpi(inImgType, 'double');
if ~isDoubleInput
   inImg=double(inImg);
end
imgSize=size(inImg);

% verifiy filter dimentions are odd, so a middle element always exists
filtDims=filtDims+1-mod(filtDims,2); 

filt=zeros(filtDims, 'double');
nNeigh=numel(filt)-1;

if nNeigh<=8
  outClass='uint8';
elseif nNeigh>8 && nNeigh<=16
   outClass='uint16';
elseif nNeigh>16 && nNeigh<=32
   outClass='uint32';
elseif nNeigh>32 && nNeigh<=64
   outClass='uint64';   
else
   outClass='double';
end

iHelix=snailMatIndex(filtDims);
filtCenter=ceil((nNeigh+1)/2);
iNeight=iHelix(iHelix~=filtCenter);


if isEfficent
   %% Better filtering/concolution based attitude

   filt(filtCenter)=1;
   filt(iNeight(1))=-1;
   sumLBP=zeros(imgSize);
   for i=1:length(iNeight)
      currNieghDiff=filter2(filt, inImg, 'same');
      sumLBP=sumLBP+2^(i-1)*(currNieghDiff>0); % Thanks goes to Chris Forne for the bug fix

      if i<length(iNeight)
         filt( iNeight(i) )=0;
         filt( iNeight(i+1) )=-1;
      end
   end   
   if strcmpi(outClass, 'double')
      LBP=sumLBP;
   else
      LBP=cast(sumLBP, outClass);
   end
else % if isEfficent
   
   %% Primitive pixelwise solution
   filtDimsR=floor(filtDims/2); % Filter Radius
   iNeight(iNeight>filtCenter)=iNeight(iNeight>filtCenter)-1; % update index values.
   
   % Padding image with zeroes, to deal with the edges
   zeroPadRows=zeros(filtDimsR(1), imgSize(2));
   zeroPadCols=zeros(imgSize(1)+2*filtDimsR(1), filtDimsR(2));

   inImg=cat(1, zeroPadRows, inImg, zeroPadRows);
   inImg=cat(2, zeroPadCols, inImg, zeroPadCols);
   imgSize=size(inImg);

   neighMat=true(filtDims);

   neighMat( floor(nNeigh/2)+1 )=false;
   weightVec= (2.^( (1:nNeigh)-1 ));
   LBP=zeros(imgSize, outClass);
   for iRow=( filtDimsR(1)+1 ):( imgSize(1)-filtDimsR(1) )
      for iCol=( filtDimsR(2)+1 ):( imgSize(2)-filtDimsR(2) )
         subImg=inImg(iRow+(-filtDimsR(1):filtDimsR(1)), iCol+(-filtDimsR(2):filtDimsR(2)));
         % find differences between current pixel, and it's neighours
         diffVec=repmat(inImg(iRow, iCol), [nNeigh, 1])-subImg(neighMat);  
         LBP(iRow, iCol)=cast( weightVec*(diffVec(iNeight)>0),  outClass);   % convert to decimal. 
      end % for iCol=(1+filtDimsR(2)):(imgSize(2)-filtDimsR(2))
   end % for iRow=(1+filtDimsR(1)):(imgSize(1)-filtDimsR(1))
   
   % crop the margins resulting from zero padding
   LBP=LBP(( filtDimsR(1)+1 ):( end-filtDimsR(1) ),...
      ( filtDimsR(2)+1 ):( end-filtDimsR(2) ));
end % if isEfficent
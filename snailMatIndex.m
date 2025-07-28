function iHelix=snailMatIndex(inMatDims)
if length(inMatDims)==1
   inMatDims=repmat(inMatDims, 1, 2);
end

nElemems=inMatDims(1)*inMatDims(2);
indMat =reshape(1:nElemems, inMatDims);
iHelix=[];

while ~isempty(indMat)
   if size(indMat, 2)==1
      iHelix=cat( 2, iHelix, transpose(indMat(:, 1)) );
      indMat=[];
   else
      iHelix=cat( 2, iHelix, indMat(1, :) );
      indMat(1, :)=[];        % remove the current top row
      if ~isempty(indMat)
         indMat=rot90(indMat);   % rotate index matrix 90° clock-wise
      end % if ~isempty(indMat)
   end % if size(indMat, 2)==1
end % while ~isempty(indMat)
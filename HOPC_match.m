function [CP_Ref,CP_Sen] = HOPC_match(im_Ref,im_Sen,CP_Check_file,disthre,tranFlag, templateSize, searchRad,SFlag);

% detect the tie points using HOPC  (histogram of oriented phase
% congruency) by the template matching strategy.As the matlab cannot handle
% the georeference information of images, we use some control points ('CP_Check_file')
% to determine the search region

% If you use this implementation please cite:
%Yuanxin Ye et al.: "HOPC: A NOVEL SIMILARITY METRIC BASED ON GEOMETRIC STRUCTURAL 
%PROPERTIES FOR MULTI-MODEL REMOTE SENSING IMAGE MATCHING"
% ISPRS XXIII congress (2016)

% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------

%
% Contact: yeyuanxin110@163.com
%input Parameters: 
%                 im_Ref: the reference image
%                 im_Sen: the sensed image
%                 CP_Check_file: the file path of check points. they are used to determine the searchRegion, and judge if the matching points are correct
%                 disthre: the threshold of the errors of matching point pair, if the error of the point pair is less than dis, they are regrad as the correct match                 
%                 tranFlag:  the geometric transformation between two images. 0:affine, 1: projective, 2: Quadratic polynomial,3: cubic polynomial,the default is 3
%                 templateSize: the template size, its minimum value must be more than 20, the default is 100
%                 searchRad: the radius of searchRegion, the default is 10.the maxinum should be less than 20
%                 sFlag: 1: compute the NCC of the hopc descriptors; 0:compute the SSD of the hopc descriptors (more fast).the default is 1
%                 
%                 
%                 
%                
%                
%                 

%return vaules:                
%                 CP_Ref: the coordinates ([x,y]) in the reference image
%                 CP_Sen: the coordinates ([x,y]) in the sensed image


if nargin < 3
    disp('the number input parameters must be >= 3 ');
    return;
end

if nargin < 4
    disthre = 1.5;        % % the threshod of match errors,here we do not use that
end
if nargin < 5
    tranFlag = 3;            % the type of geometric transformation between images
end
if nargin < 6
    templateSize = 100;% the template size
end
if nargin < 7
    searchRad = 10;    % the radius of search region
end
if nargin < 7
    SFlag = 0;    % compute the SSD of hopc descriptors
end

% tranfer the rgb to gray
[k1,k2,k3] = size(im_Ref);
if k3 == 3
    im_Ref = rgb2gray(im_Ref);
end
im_Ref = double(im_Ref);

[k1,k2,k3] = size(im_Sen);
if k3 == 3
    im_Sen = rgb2gray(im_Sen);
end
im_Sen = double(im_Sen);

[im_RefH,im_RefW] = size(im_Ref);
[im_SenH,im_SenW] = size(im_Sen);

templateRad = round(templateSize/2);    %the template radius
marg=templateRad+searchRad+2;           %the boundary. we don't detect tie points out of the boundary

C = 0;;%the number of correct match 
CM = 0 ;%the number of total match 
C_e = 0;%the number of mismatch
e = 0.0000001; % avoid divided by zero


%extract the interest points using blocked-harris
im1 = im_Ref(marg:im_RefH-marg,marg:im_RefW-marg);% remove the pixel near the boundary
Value = harrisValue(im1);                        % harris intensity value
[r,c,rsubp,cubp] = nonmaxsupptsgrid(Value,3,0.3,5,8); % non-maxima suppression in regular
                                                       % here is 5*5 grid
                                                       % 8 points in each grid, in total 200 interet points 

points1 =[r,c] + marg - 1;
pNum = size(points1,1); % the number of interest points



%caculate the dense block-HOPC descriptor for each pixel

cellsize=3; % the pixel number (cellsize*cellsize) in a cell
orbin = 8;  % the number of orientation bin
nCells =3;  % the cell number (nCell*nCell) of in a block 
blockoverlay = 0.5;% the overlay degree between adjacent block
interval =round(cellsize*nCells*(1-blockoverlay));% the pixel interval between adjact block
if interval == 0
    interval =1;
end
%interval =1;
%caculate the dense block-HOPC descriptor 
blockHOPC_Ref = denseBlockHOPC(single(im_Ref),cellsize,orbin,nCells);
blockHOPC_Sen = denseBlockHOPC(single(im_Sen),cellsize,orbin,nCells);

%read check points from file;
checkPt = textread(CP_Check_file);
refpt = [checkPt(:,1),checkPt(:,2)]; %the check points in the referencing image
senpt = [checkPt(:,3),checkPt(:,4)]; %the check points in the sensed image

% solve the geometric tranformation parameter
% tran 0:affine, 1: projective, 2: Quadratic polynomial,3: cubic polynomial,the default is 3
tform = [];
if tranFlag == 0
    tform = cp2tform(refpt,senpt,'affine'); 
    T = tform.tdata.T;
elseif tranFlag == 1
    tform = cp2tform(refpt,senpt,'projective');
    T = tform.tdata.T;
    else
    T = solvePoly(refpt,senpt,tranFlag);
end
H = T';%the geometric transformation parameters from im_Ref to im_Sen

%detect the tie points by the template matching strategy for each
%interestin pionts
for n = 1: pNum
    
    %the x and y coordinates in the reference image
    X_Ref=points1(n,2);
    Y_Ref=points1(n,1);
    
    % get the HOPC descriptor of the template window centered on (X_Ref,Y_Ref)
    HOPC_Ref = single(getDesc(blockHOPC_Ref,Y_Ref,X_Ref,templateRad,interval));

    
    %transform the (x,y) of reference image to sensed image by the geometric relationship of check points 
    %to determine the search region
    tempCo = [X_Ref,Y_Ref];
    tempCo1 = transferTo(tform,tempCo,H,tranFlag);
    
    %tranformed coordinate (X_Sen_c, Y_Sen_c)
    X_Sen_c = tempCo1(1);
    Y_Sen_c = tempCo1(2);
    X_Sen_c1=round(tempCo1(1));
    Y_Sen_c1 =round(tempCo1(2)); 

    %judge whether the transformed points are out the boundary of right image.

    if (X_Sen_c1 < marg+1 | X_Sen_c1 > size(im_Sen,2)-marg | Y_Sen_c1<marg+1 | Y_Sen_c1 > size(im_Sen,1)-marg)
        %if out the boundary, this produre enter the next cycle
        continue;
    end
        
    corr = zeros(2*searchRad + 1); % the NCC of HOPC descriptor
     
    % caculate the NCC of HOPC for the search region 
    for i = -searchRad:searchRad
        for j = -searchRad:searchRad
            
            Y_Sen_c2 = Y_Sen_c1+i;
            X_Sen_c2 = X_Sen_c1+j;
            
            % get the HOPC descriptor of the template window centered on (X_Sen,Y_Sen)
            HOPC_Sen = single(getDesc(blockHOPC_Sen,Y_Sen_c2,X_Sen_c2,templateRad,interval));
           
           if SFlag 
               %calculate the NCC between two HOPC descriptors 
               temp=corrcoef(HOPC_Ref,HOPC_Sen);
               corr(i + searchRad +1,j + searchRad + 1)=temp(1,2);
           else
               %calculate the SSD between two HOPC descriptors 
               temp = -norm(HOPC_Ref(:)-HOPC_Sen(:));
               corr(i + searchRad +1,j + searchRad + 1)=temp;
           end
        end
    end
    %judge if corr is nan
     nan = isnan(corr);
     if nan(1) > 0
         continue;
     end
    
    %get coordinates with the maximum NCC in the search region 
     maxCorr = max(max(corr));
     max_index = find(corr == maxCorr);
      if(size(max_index,1) > 1);
            % if two maxumal appear, it go to nexe cycle;
         continue;
      end
       
      [max_i,max_j] = ind2sub(size(corr),max_index);
      
      %the (matchY,matchX) coordinates of match
      Y_match = Y_Sen_c1-searchRad + max_i-1;
      X_match = X_Sen_c1-searchRad + max_j-1;
     % calculate the match errors      
       diffY = abs(Y_match-Y_Sen_c);
       diffX = abs(X_match-X_Sen_c);
      diff = sqrt(diffX.^2+diffY.^2);
      
      % calculate the numbers of correct match, mismatch and total match

          C = C+1; % the number of correct matches
          corrp(C,:)=[X_Ref,Y_Ref,X_match,Y_match,diff];% the coordinates of correct matches

end

    CP_Ref = corrp(:,1:2);
    CP_Sen = corrp(:,3:4);
end


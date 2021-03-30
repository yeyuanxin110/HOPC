
clear;
% optical and sar matching
 im_Ref = imread('.\data\optical_ref.png');
 im_Sen = imread('.\data\SAR_sen.png');
 CP_Check_file = '.\data\OpticaltoSAR_CP.txt';


%lidar depth and optical matching
% im_Ref = imread('.\data\LidarDepth_ref.png');
% im_Sen = imread('.\data\optical1_sen.png');
% CP_Check_file = '.\data\LidartoOptical_CP.txt';


% % %visible and infrared image matching
%im_Ref = imread('.\data\visible_ref.tif');
%im_Sen = imread('.\data\infrared_sen.tif');
%CP_Check_file = '.\data\VisibletoInfrared_CP.txt';

%im_Ref = imread('.\data\Stest3_ref.tif');
%im_Sen = imread('.\data\Stest3_sen.tif');
%CP_Check_file = '.\data\Stest3gcp.pts';

errorthre = 1.5;        % the threshod for error detection, the deflaut is 1.5. for
                      % high resolution image covering urban areas, we
                      % should set it to a larger threshod (such as 2.0).
                      % This is beccause that the geometric distortions between such images
                      % is very complicated, and the transfrom model used
                      % by us (such as projective model) can only prefit
                      % the geometric distortion. Therefore, a larger threshod
                      % have the more flexibility  

% template matching using HOPC
tic;
[CP_Ref,CP_Sen] = HOPC_match(im_Ref,im_Sen,CP_Check_file);
fprintf('the total matching time is %fs\n',toc);


%detect the error
[corrRefPt,corrSenPt] = ErrorDect(CP_Ref,CP_Sen,0,errorthre);

%wirite the point in the envi format
corrPT = [corrRefPt,corrSenPt];%correct match

%[H1, inlier_ransac] = ransacfithomography(CP_Ref', CP_Sen', 0.001);

%cleanedPoints11 = matchedPoints1(inliersIndex, :);
%cleanedPoints22 = matchedPoints2(inliersIndex, :);
%corrRefPt1 = CP_Ref(inlier_ransac, :);
%corrSenPt1 =  CP_Ref(inlier_ransac, :);


path = 'D:\test1.pts'
Wenvifile1(corrPT',path);

figure;
imshow(im_Ref),hold on;
plot(corrRefPt(:,1),corrRefPt(:,2),'yo','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',5);hold on;
title('reference image');

figure;
imshow(im_Sen),hold on;
plot(corrSenPt(:,1),corrSenPt(:,2),'yo','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',5);hold on;
title('sensed image');



function denseDesc = getDesc(desc,y,x,templateRad,interval);
% This function collects the block-HOPC descriptor at interval to form the
% HOPC descriptor of a template window

%inpute parameters:
%              desc     - the block-HOPC decriptors for a image
%              y        - the y coordinate of the input interest point
%              x        - the x coordinate of the input interest point
%        templateRad    - the radius of the template size
%         interval      - the interval for collecting the block-HOPC
%                         descriptors

%return value:
%              denseDesc- the HOPC descriptor of the template window


[h,w,d]=size(desc);

upY1 = y - templateRad;
downY2 = y + templateRad-interval;

leftX1 = x - templateRad;
rightX2 = x + templateRad-interval;

Ybegin = upY1;
Yend = downY2 ;

Xbegin = leftX1;
Xend = rightX2;

if Ybegin < 1
    Ybegin = 1;
end
if Xbegin < 1
    Xbegin= 1;
end
if Yend > h
    Yend = h;
end
if Xend > w
    Xend = w;
end

denseDesc = desc(Ybegin:interval:Yend,Xbegin:interval:Xend,:);
denseDesc = permute(denseDesc,[3,2,1]);

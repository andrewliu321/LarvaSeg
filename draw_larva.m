title('Input outline of Larva')
[x2,y2]=getpts(f1);
%round the points
x2=round(x2);
y2=round(y2);
%set zero to collect pixel list
xpixel2 = [];
ypixel2 = [];
    
for count = 1:size(x2,1)

    %if last point, connect to the first point
    if count == size(x2,1)
        count1=1;
    else %else, compare to the next point
        count1=count+1;
    end
    
    %x or y has bigger difference?
    dif_x2=abs(x2(count)-x2(count1));
    dif_y2=abs(y2(count)-y2(count1));

    if dif_x2 > dif_y2
        %which have bigger difference, find all points in between
        if x2(count)<x2(count1)
            xpixel2_temp=(x2(count):x2(count1));
        else
            xpixel2_temp=(x2(count):-1:x2(count1));
        end
        ypixel2_temp=round(interp1([x2(count);x2(count1)],[y2(count);y2(count1)],xpixel2_temp));
    else
        if y2(count)<y2(count1)
            ypixel2_temp=(y2(count):y2(count1));
        else
            ypixel2_temp=(y2(count):-1:y2(count1));
        end
        xpixel2_temp=round(interp1([y2(count);y2(count1)],[x2(count);x2(count1)],ypixel2_temp));
    end

    %add each segment to the previous one
    xpixel2 = [xpixel2 xpixel2_temp];
    ypixel2 = [ypixel2 ypixel2_temp]; 
end
%%
%draw the boundaries of the larva
larva=zeros(size(im(:,:,1)),'double');
indexes2 = sub2ind(size(im(:,:,1)), ypixel2, xpixel2);
larva(indexes2) = 1;
larva = logical(larva);

larva=larva(any(larva,2),any(larva,1));
larva = padarray(larva,[100 100],0,'both');

%fill fill holes to get pouch
larva2 = imfill(larva,'holes');

larva_record{position} = larva;

% area_larva(position) = bwarea(larva2); % calculated outside this script

%% show skeleton

larva3 = bwmorph(larva2,'skel',inf);

larva0 = larva3;
larva0 = uint8(larva0);
larva0(:,:,3)= larva2*100;
larva0(:,:,1) = larva3*254;

f4 = figure(4);
close(f4)
f4 = figure(4)
imshow(larva0)
title({'Computed skeleton of image','Input Midline'})
%% input midline manually


%show whole project stack
[x1,y1]=getpts(f4);
%round the points
x1=round(x1);
y1=round(y1);
%set zero to collect pixel list
xpixel1 = [];
ypixel1 = [];
    
for count = 1:size(x1,1)-1

    %x or y has bigger difference?
    dif_x1=abs(x1(count)-x1(count+1));
    dif_y1=abs(y1(count)-y1(count+1));

    if dif_x1 > dif_y1
        %which have bigger difference, find all points in between
        if x1(count)<x1(count+1)
            xpixel1_temp=(x1(count):x1(count+1));
        else
            xpixel1_temp=(x1(count):-1:x1(count+1));
        end
        ypixel1_temp=round(interp1(x1(count:count+1),y1(count:count+1),xpixel1_temp));
    else
        if y1(count)<y1(count+1)
            ypixel1_temp=(y1(count):y1(count+1));
        else
            ypixel1_temp=(y1(count):-1:y1(count+1));
        end
        xpixel1_temp=round(interp1(y1(count:count+1),x1(count:count+1),ypixel1_temp));
    end

    %add each segment to the previous one
    xpixel1 = [xpixel1 xpixel1_temp];
    ypixel1 = [ypixel1 ypixel1_temp]; 
end

%plot line
hold on
plot(x1,y1,'b')
plot(xpixel1,ypixel1,'c')%,'LineWidth',4)

%make im of line
indexes1 = sub2ind(size(larva2), ypixel1, xpixel1);
mid = zeros(size(larva2),'uint8');
mid(indexes1)=255;
% outside the boundary, erase the ap line
mid(larva2==0)=0;

larva0(:,:,1) = mid;

f4 = figure(4);
imshow(larva0)

% larva_length(position) = length(indexes1);

%% calculate perpendicular lines and shortest intersect w/ larva1:
% which is the boundary of larva segmentation

point11=[];
point22=[];


[mid_index,~] = bwboundaries(mid,8,'noholes');

%get endpoints of the line
mid_end = bwmorph(mid,'endpoints');
%get the pixel list of the two end points
[pts(:,1), pts(:,2)] = find(mid_end == 1);
%find where these poitns lie within the pixel list of the lines
[~,points] = ismember(pts,mid_index{1},'rows');
%remake index with only that list
mid_index{1} = mid_index{1}(min(points):max(points),:);


f4 = figure(4);
plot(mid_index{1}(:,2),mid_index{1}(:,1),'*')

%record_mid_index = new points identifying the midline
record_mid_index{position} = mid_index{1};

parfor count1 = 1:size(record_mid_index{position},1)

    % for each pixel, define the center:
    center=[record_mid_index{position}(count1,2) record_mid_index{position}(count1,1)];
    
   
    if count1 < 4 
        %if 1st and 2nd pixel, then use itself as left bound
        left = center; 
    else
        %or else, left pixel = count1-2 (2 pixels to the left)
        left = [record_mid_index{position}(count1-3,2) record_mid_index{position}(count1-3,1)];
    end
    
    if count1 > size(record_mid_index{position},1) - 3
        %if last or 2nd til last pixel, then use itself as right bound
        right = center;
    else
        %or else right pixel = count+2
        right = [record_mid_index{position}(count1+3,2) record_mid_index{position}(count1+3,1)];
    end

    
    %calculate slope - inverse negative slope = normal
    rise = (right(2)-left(2));
    run = (right(1)-left(1));
%     slope=rise/run;
    slope_inv = -run/rise;

    % if slope is horizontal, then just define the end points
    if rise==0
        solx=[center(1) center(1)];
        soly=[center(2)+100 center(2)-100];
    else
        %if slope isn't horizontal, then solve for end points using slope and 
        %distance equations
        x = sym('x');
        y = sym('y');
        eqns = [(y-center(2))/(x-center(1))==slope_inv, 100^2 == (x-center(1))^2 + (y-center(2))^2];
        vars = [x y];
        [solx, soly] = solve(eqns,vars);
    end

    %turn solution into round numbers
    solx=round(double(solx));
    soly=round(double(soly));
    
    %define the points based on solutions...
    point11(count1,:) = [(solx(1)) (soly(1))];
    point22(count1,:) = [(solx(2)) (soly(2))];

    
end



%% plot the perpendicular

figure(4)
hold off
imshow(larva)
hold on

clear c
clear cx
clear cy
clear d


for i = 1:size(point11,1)
    perpen(1,:) = point11(i,:);
    perpen(2,:) = point22(i,:);
    plot(perpen(:,1), perpen(:,2))

    %find intersection of the perpendicular line & larva outline
    [cx,cy,c] = improfile(larva,perpen(:,1),perpen(:,2));

    %if it intersects in two points:
    if sum(c)==2
    
        %find the pixel list of the two points
        X(:,1)=cx(c==1);
        X(:,2)=cy(c==1);
        %calculate distance and record it
        d(i) = pdist(X,'euclidean');
        
    end

end

%calculate average width by ignoring the 0's (ones where the perpendicular
%does not intersect 2 points)

% mean(nonzeros(d))  > calcualted outside this script












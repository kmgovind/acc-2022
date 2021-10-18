function [goalLat,goalLong] = goalSelection(environment,vehicle)
%goalSelection selects the endpoint for the next path generation
%   Current edition, find minimum points on map, choose the farthest one
%   from current location
%   Inputs: map - matrix representing current cost map
%           vehicle - struct representing attributes of vehicle
%   Outputs: goalLong - longitude of goal location
%            goalLat - latitude of goal location

% find lowest points of coverage and their locations

dataCell = struct2cell(environment.coverageMap);
dataMatrix = cell2mat(dataCell);
lats = squeeze(dataMatrix(1,:,:));
longs = squeeze(dataMatrix(2,:,:));
coverages = squeeze(dataMatrix(3,:,:));

% find boat starting location linear index
lats = round(lats, 2);
longs = round(longs,2);
tempLat = round(vehicle.latitude, 2);
tempLong = round(vehicle.longitude, 2);
latIndex = find(lats == tempLat);
longIndex = find(longs == tempLong);
start = latIndex(find(ismember(latIndex, longIndex))); % this is the starting point on the grid (boat location)

[boatRow, boatCol] = ind2sub(size(coverages), start); %boat x,y

tempCoverages = coverages; %matrix to calculate sums of coverage from boat location to gridpoint

% figure;
% contourf(1 - coverages);
% title('Coverages');

% calculate distances from current location to minimum locations
distMatrix = zeros(size(environment.coverageMap));
for distCheck = 1:numel(distMatrix)
    distMatrix(distCheck) = distance('gc', vehicle.latitude, vehicle.longitude, environment.coverageMap(distCheck).latitude, environment.coverageMap(distCheck).longitude);
end

figure;
subplot(2,2,1);
contourf(distMatrix);
title('DistMatrix');

% Calculate estimated gain in coverage from boat location to every point in
% domain

m = 0.01/max(vehicle.motorSpeed, 3);
% r = km2deg(6);

for index = 1:numel(tempCoverages)
    [goalRow, goalCol] = ind2sub(size(coverages), index);
    [xVals, yVals] = linEst(boatRow, boatCol, goalRow, goalCol);
    indices = sub2ind(size(coverages), yVals, xVals);
    tempCoverages(index) = sum(1-coverages(indices) - m.*distMatrix(indices).*(1 - index/numel(indices)));
%     tempCoverages(index) = sum(r.*(1-coverages(indices)) - (0.01)*coverages(indices));
end

subplot(2,2,2);
contourf(tempCoverages);
title('Temp Coverages');


timeMatrix = distMatrix./(vehicle.motorSpeed);%finds times to get to each point

subplot(2,2,3);
contourf(timeMatrix);
title('timeMatrix');

% TODO: Increase quality of estimation by interpolating more points
% (not just endpoints)
% r = 6; %6 nm sensing radius
% distMatrix = r.*distMatrix.*(1 - tempCoverage); %solve for coverage gain with a weight towards further locations

selCoverages = tempCoverages./max(timeMatrix,0.05); % add constant k
% selCoverages = tempCoverages./timeMatrix;

subplot(2,2,4);
contourf(selCoverages);
title('Goal Selection');

% find which of the coverageGains is the largest
maxDistIndex= find(selCoverages == max(max(selCoverages)));

% [~,maxDistIndex] = max(selCoverages, 'all');


%     maxDistIndex = maxDistIndex(1);

close all;

% assign goal as the first of the largest distances (in case there are
% multiple distances of the same value that are the farthest)
goalLat = environment.coverageMap(maxDistIndex).latitude;
goalLong = environment.coverageMap(maxDistIndex).longitude;
end

function [xVals, yVals] = linEst(boaty, boatx, goaly, goalx)

slope = (goaly - boaty)/(goalx - boatx); %slope of linear estimation
xVals = 1;
yVals = 1;

if abs(goaly - boaty) > abs(goalx - boatx) % if there is more room in the y-points sample using y-values
    if boaty < goaly
        yVals = boaty:1:goaly;
        xVals = floor(((yVals-boaty)./slope) + boatx);
        xVals(xVals == 0) = 1;
    elseif boaty > goaly
        yVals = boaty:-1:goaly;
        xVals = floor(((yVals-boaty)./slope) + boatx);
        xVals(xVals == 0) = 1;
    elseif boaty == goaly
        if boatx > goalx
            xVals = boatx:-1:goalx;
        else
            xVals = boatx:1:goalx;
        end
        xVals(xVals == 0) = 1;
        yVals = zeros(1,numel(xVals));
        yVals(yVals == 0) = boaty;
    end
else %sample using x-values
    if boatx < goalx
        xVals = boatx:1:goalx;
        yVals = floor(slope.*(xVals - boatx) + boaty);
        yVals(yVals == 0) = 1;
    elseif boatx > goalx
        xVals = boatx:-1:goalx;
        yVals = floor(slope.*(xVals - boatx) + boaty);
        yVals(yVals == 0) = 1;
    elseif boatx == goalx
        if boaty > goaly
            yVals = boaty:-1:goaly;
        else
            yVals = boaty:1:goaly;
        end
        yVals(yVals == 0) = 1;
        xVals = zeros(1,numel(yVals));
        xVals(xVals == 0) = boatx;
    end
    
end

end

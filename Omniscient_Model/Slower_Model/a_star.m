function final = a_star(environment,vehicle,target)
%UNTITLED Summary of this function goes here
%   final - Returns a list of coordinates that represent path in terms of
%   latitude and longitude

% Default return - empty for failure (which shouldn't happen)
global currentTime;
final = [];

% Create simplified coverage map for a* to operate on and establish
% boat in a grid point on the map
dataCell = struct2cell(environment.coverageMap);
dataMatrix = cell2mat(dataCell);
lats = squeeze(dataMatrix(1,:,:));
longs = squeeze(dataMatrix(2,:,:));
coverages = squeeze(dataMatrix(3,:,:));
lats = round(lats, 2);
longs = round(longs, 2);

% find boat starting location linear index
tempLat = round(vehicle.latitude, 2);
tempLong = round(vehicle.longitude, 2);
latIndex = find(lats == tempLat);
longIndex = find(longs == tempLong);
start = latIndex(find(ismember(latIndex, longIndex))); % this is the starting point on the grid (boat location)

% find goal location linear index

tempLat = round(target.lat, 2);
tempLong = round(target.long, 2);
latIndex = find(lats == tempLat);
longIndex = find(longs == tempLong);
goal = latIndex(find(ismember(latIndex, longIndex))); % this is the goal point

% Create cost barrier to edges of domain
% left barrier
boundind = find(longs(1,:) == round(environment.longitudeRange(1) + km2deg(5),2));
coverages(1:boundind,:) = 2;

% right barrier
boundind = find(longs(1,:) == round(environment.longitudeRange(end) - km2deg(5),2));
coverages(boundind:end,:) = 2;

% bottom barrier
boundind = find(lats(:,1) == round(environment.latitudeRange(1) + km2deg(5),2));
coverages(:,1:boundind) = 2;

% top barrier
boundind = find(lats(:,1) == round(environment.latitudeRange(end) - km2deg(5),2));
coverages(:,boundind:end) = 2;

% Rewrite flowmap in terms of cost (pointing towards goal is low cost)
% find heading to goal from current location
overallLongs = -79:0.01:-74;
overallLats = 32:0.01:37;

latitudeRange = environment.latitudeRange;
longitudeRange = environment.longitudeRange;

latindStart = find(overallLats == latitudeRange(1));
latindEnd = find(overallLats == latitudeRange(end));
longindStart = find(overallLongs == longitudeRange(1));
longindEnd = find(overallLongs == longitudeRange(end));


goalHeading = vehicle.headingCalc(target.lat, target.long);

currentDay = round(days(minutes(currentTime)));
flowus = environment.enviroData.u(latindStart:latindEnd,longindStart:longindEnd, currentDay);
flowvs = environment.enviroData.v(latindStart:latindEnd,longindStart:longindEnd, currentDay);

[flowVels, flowHeads] = vehicle.flowHeading(flowus, flowvs);


flowHeads = (-cosd(abs(flowHeads - goalHeading)) + 1)./2;
% flowHeads = -cosd(abs(flowHeads - goalHeading));
flowMap = flowVels.*flowHeads;


% Setup map with coverage values
map = coverages + flowMap;
% map = coverages;

% if currentTime >= 9000 && currentTime <= 12000
%     hold on;
%     contourf(map);
%     keyboard
% end


% Initialize open set with start
mapSize = size(map);
mapNumEl = numel(mapSize);

%Initialize open set, starting with start
openSet = false(mapSize);
openSet(start) = true;

% Initialize close set (close set is for visited locations)
closedSet = false(mapSize);
cameFrom = zeros(1, mapNumEl);

% G-score, cost from start to current
gScore = inf(mapSize);
gScore(start) = 0;

%linear index -> row, col for the goal
[goalRow, goalCol] = ind2sub(mapSize, goal);

% F-score, g + h
fScore = inf(mapSize);
fScore(start) = costCalculation(map, start, goal, vehicle, flowus, flowvs);
% fScore(start) = costCalculation(map, start, goal);

S2 = sqrt(2);
% While the open set is not empty
while any(openSet(:) > 0)
    
    % Find the minimum fScore within the open set
    [~, current] = min(fScore(:));
    
    % If we've reached the goal
    if current == goal
        % Get the full path and return it
        final = get_path(cameFrom, current);
        return
    end
    
    % Linear index -> row, col subscripts
    rc = rem(current - 1, mapSize(1)) + 1;
    cc = (current - rc) / mapSize(1) + 1;
    
    % Remove CURRENT from openSet
    openSet(rc, cc) = false;
    % Place CURRENT in closedSet
    closedSet(rc, cc) = true;
    
    fScore(rc, cc) = inf;
    gScoreCurrent = gScore(rc, cc) + map(rc, cc);
    
    % Get all neighbors of CURRENT. Neighbors are adjacent indices on
    %   the map, including diagonals.
    % Col 1 = Row, Col 2 = Col, Col 3 = Distance to the neighbor
    n_ss = [ ...
        rc + 1, cc + 1, S2 ; ...
        rc + 1, cc + 0, 1 ; ...
        rc + 1, cc - 1, S2 ; ...
        rc + 0, cc - 1, 1 ; ...
        rc - 1, cc - 1, S2 ; ...
        rc - 1, cc - 0, 1 ; ...
        rc - 1, cc + 1, S2 ; ...
        rc - 0, cc + 1, 1 ; ...
        ];
    
    % keep valid indices only
    valid_row = n_ss(:,1) >= 1 & n_ss(:,1) <= mapSize(1);
    valid_col = n_ss(:,2) >= 1 & n_ss(:,2) <= mapSize(2);
    n_ss = n_ss(valid_row & valid_col, :);
    % subscripts -> linear indices
    neighbors = n_ss(:,1) + (n_ss(:,2) - 1) .* mapSize(1);
    % only keep neighbors in the map and not in the closed set
    ixInMap = map(neighbors) & ~closedSet(neighbors);
    neighbors = neighbors(ixInMap);
    % distance to each kept neighbor
    dists = n_ss(ixInMap, 3);
    
    % Add each neighbor to the open set
    openSet(neighbors) = true;
    
    % TENTATIVE_GSCORE is the score from START to NEIGHBOR.
    tentative_gscores = gScoreCurrent + (map(neighbors).*0.6) .* (dists.*0.4);
    
    % IXBETTER indicates where a better path was found
    ixBetter = tentative_gscores < gScore(neighbors);
    bestNeighbors = neighbors(ixBetter);
    
    % For the better paths, update scores
    cameFrom(bestNeighbors) = current;
    gScore(bestNeighbors) = tentative_gscores(ixBetter);
    fScore(bestNeighbors) = gScore(bestNeighbors) + costCalculation(map, bestNeighbors, goal, vehicle, flowus, flowvs);
%     fScore(bestNeighbors) = gScore(bestNeighbors) + costCalculation(map, bestNeighbors, goal);
    
end % while

end

function cost = costCalculation(map,current,goal, vehicle, flowu, flowv)
%Calculate the cost associated with traveling from the current spot to the
%goal. This includes the current cells coverage as part of the weighting
%and the euclidean distance to the goal

% Weights
coverageWeight = 0.35;
flowWeight = 0.60;
distToGoWeight = 0.05;

%linear index to row,col
mapSize = size(map);
[currentRow, currentCol] = ind2sub(mapSize, current);
[goalRow, goalCol] = ind2sub(mapSize, goal);

% Values
distToGoal = sqrt((currentRow - goalRow).^2 + (currentCol - goalCol).^2)./(10*sqrt(2));
coverageAtLoc = map(current);

% Vehicle Stuff
uvel = vehicle.motorSpeed_u + flowu(current);
vvel = vehicle.motorSpeed_v + flowv(current);
[vtrue, ~] = vehicle.flowHeading(uvel, vvel);
flowAtLoc = distToGoal./abs(vtrue);
if isnan(flowAtLoc)
    flowAtLoc = 0;
end

cost = coverageWeight.*coverageAtLoc + distToGoWeight.*distToGoal + flowWeight.*flowAtLoc;
% cost = coverageWeight.*coverageAtLoc + distToGoWeight.*distToGoal;
if isnan(cost)
    cost = 1;
end
end

function p = get_path(cameFrom, current)
% Returns the path. This function is only called once and therefore
%   does not need to be extraordinarily efficient
inds = find(cameFrom);
p = nan(1, length(inds));
p(1) = current;
next = 1;
while any(current == inds)
    current = cameFrom(current);
    next = next + 1;
    p(next) = current;
end
p(isnan(p)) = [];
end

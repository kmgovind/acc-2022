clear;clc; close all;
% Driver Code

filenamedtp = 'boatMotionD2P.gif';
filenameastar = 'boatMotionAstar.gif';
filenametransect = 'boatMotionTransect.gif';
global currentTime
speedConstant = convvel(3, 'kts', 'm/s');
wpThresh = 0.05;

%% D2P implementation
% % Initialize Environment and Vehicles
% domaindtp = Environment; % Mission Domain
% boatdtp = Vehicle(33.5,-76.1); % initialize boat at start lat/long
% 
% totalCoveraged2p = zeros(size(domaindtp.startTime:minutes(domaindtp.timeStep):domaindtp.endTime));
% 
% % Track boattransect path
% latListDTP = nan(1,domaindtp.endTime - domaindtp.startTime);
% longListDTP = nan(1,domaindtp.endTime - domaindtp.startTime);
% 
% % Initialize first goal
% % tic
% % List of upcoming waypoints
% [goal.lat, goal.long] = goalSelection(domaindtp, boatdtp);
% % toc
% 
% 
% for currentTime = domaindtp.startTime:minutes(domaindtp.timeStep):domaindtp.endTime
%     %     tic
%     %     fprintf('\nDist2Goal: %d', distance('gc',goal.lat, goal.long, boatdtp.latitude, boatdtp.longitude));
%     %     fprintf('\nDist2Goal nm: %d\n', deg2nm(distance('gc',goal.lat, goal.long, boatdtp.latitude, boatdtp.longitude)));
%     %     toc
%     
%     %     tic
%     % Find new goal if sufficiently close to previous waypoint
%     if deg2nm(distance('gc',goal.lat, goal.long, boatdtp.latitude, boatdtp.longitude)) <= wpThresh
%         [goal.lat, goal.long] = goalSelection(domaindtp, boatdtp);
%     end
%     
%     % Move boatdtp towards goal
%     boatdtp = boatdtp.moveBoat(domaindtp, speedConstant, goal.lat, goal.long);
%     %     fprintf('\nTime: %d\tLat: %d\tLong:%d', currentTime-domaindtp.startTime+1, boatdtp.latitude, boatdtp.longitude);
%     latListDTP(currentTime - domaindtp.startTime + 1) = boatdtp.latitude;
%     longListDTP(currentTime - domaindtp.startTime + 1) = boatdtp.longitude;
%     
%     % Update Coverage
%     [domaindtp, totalCoveraged2p(currentTime - domaindtp.startTime + 1)] = domaindtp.updateCoverage(boatdtp.latitude, boatdtp.longitude);
%     %     fprintf('\nTotalCoverage: %d\n', totalCoveraged2p(currentTime - domaindtp.startTime + 1));
%     
%     % PLOT UPDATE
%     hold on;
%     
%     % pull coverage map values out of data structure
%     dataCell = struct2cell(domaindtp.coverageMap);
%     dataMatrix = cell2mat(dataCell);
%     lats = squeeze(dataMatrix(1,:,:));
%     longs = squeeze(dataMatrix(2,:,:));
%     coverages = squeeze(dataMatrix(3,:,:));
%     
%     % plot coverage map
%     contourf(longs, lats, coverages);
%     
%     if currentTime == domaindtp.endTime
%         saveas(gcf, '32dtpPathOverCoverage.fig');
%     elseif currentTime == minutes(days(15))
%         saveas(gcf, '15dtpPathOverCoverage.fig');
%     end
%     
%     % plot boattransect location, goal location, and path traveled
%     plot(goal.long, goal.lat, 'r*', 'MarkerSize', 20) % Path goal
%     plot(boatdtp.longitude, boatdtp.latitude, 'md', 'MarkerSize', 10) % Current vehicle location
%     plot(longListDTP, latListDTP, 'r-', 'LineWidth', 2); % Path to follow
%     % set(gca, 'ydir', 'reverse')
%     axis([-76.1,-75.1,33.5,34.5]);
%     
%     
%     % Add each timestep as image in GIF
% %     drawnow
% %     frame = getframe(1);
% %     im = frame2im(frame);
% %     [imind,cm] = rgb2ind(im,256);
% %     if currentTime - domaindtp.startTime + 1 == 1
% %         imwrite(imind,cm,filenamedtp,'gif', 'DelayTime',0.1, 'Loopcount',inf);
% %     else
% %         imwrite(imind,cm,filenamedtp,'gif','WriteMode','append', 'DelayTime',0.1);
% %     end
%     
%     
%     
%     clf;
%     
% end
% % toc
% 
% % Path over time
% y = latListDTP(latListDTP ~= 0);
% x = longListDTP(longListDTP ~= 0);
% z = zeros(size(x));
% col = 1:length(x);  % This is the color, vary with x in this case.
% surface([x;x],[y;y],[z;z],[col;col],...
%     'facecol','no',...
%     'edgecol','interp',...
%     'linew',2);
% axis([-76.1,-75.1,33.5,34.5]);
% title('D2P Path Over Time');
% saveas(gcf, 'd2ppath.fig');
% 
% % Convert total coverage to percentage and plot
% figure;
% plot(totalCoveraged2p./numel(domaindtp.coverageMap));
% saveas(gcf, 'd2pcoverage.fig');
% 
% %% A* implementation
% % Initialize Environment and Vehicles
% domainastar = Environment; % Mission Domain
% boatastar = Vehicle(33.5,-76.1); % initialize boatastar at start lat/long
% 
% totalCoverageastar = zeros(size(domainastar.startTime:minutes(domainastar.timeStep):domainastar.endTime));
% 
% %Initialize first goal
% % tic
% %List of upcoming waypoints
% [goal.lat, goal.long] = goalSelection(domainastar, boatastar);
% overall.lat = goal.lat;
% overall.long = goal.long;
% % toc
% 
% latListAstar = nan(1,domainastar.endTime - domainastar.startTime);
% longListAstar = nan(1,domainastar.endTime - domainastar.startTime);
% 
% 
% paths = a_star(domainastar, boatastar, goal);
% paths = fliplr(paths);
% goal.lat = domainastar.coverageMap(paths(1)).latitude;
% goal.long = domainastar.coverageMap(paths(1)).longitude;
% 
% for currentTime = domainastar.startTime:minutes(domainastar.timeStep):domainastar.endTime
%     %     tic
%     %     fprintf('\nDist2Goal: %d', distance('gc',goal.lat, goal.long, boatastar.latitude, boatastar.longitude));
%     %     fprintf('\nDist2Goal nm: %d\n', deg2nm(distance('gc',goal.lat, goal.long, boatastar.latitude, boatastar.longitude)));
%     %     toc
%     
%     %     tic
%     %Continue navigating to goal while sufficiently far away
%     if deg2nm(distance('gc',goal.lat, goal.long, boatastar.latitude, boatastar.longitude)) <= wpThresh
%         
%         %Chop off previous goal
%         paths = paths(2:end);
%         
%         %If we have exhausted paths, select new goal and a* to it
%         if numel(paths) == 0
%             [goal.lat, goal.long] = goalSelection(domainastar, boatastar);
%             paths = a_star(domainastar, boatastar, goal);
%             paths = fliplr(paths);
%             %             if currentTime >= 9000 && currentTime <= 12000
%             %                 hold on;
%             %                 [tempCols, tempRows] = ind2sub(size(domainastar.coverageMap),paths);
%             %                 plot(tempRows, tempCols, 'r-', 'LineWidth', 2);
%             %                 keyboard
%             %                 hold off;
%             %             end
%         end
%         
%         goal.lat = domainastar.coverageMap(paths(1)).latitude;
%         goal.long = domainastar.coverageMap(paths(1)).longitude;
%     end
%     %Move boatastar towards goal
%     boatastar = boatastar.moveBoat(domainastar, speedConstant, goal.lat, goal.long);
%     %     fprintf('\nTime: %d\tLat: %d\tLong:%d', currentTime-domainastar.startTime+1, boatastar.latitude, boatastar.longitude);
%     latListAstar(currentTime - domainastar.startTime + 1) = boatastar.latitude;
%     longListAstar(currentTime - domainastar.startTime + 1) = boatastar.longitude;
%     
%     
%     %Update Coverage
%     [domainastar, totalCoverageastar(currentTime - domainastar.startTime + 1)] = domainastar.updateCoverage(boatastar.latitude, boatastar.longitude);
%     %     fprintf('\nTotalCoverage: %d\n', totalCoverageastar(currentTime - domainastar.startTime + 1));
%     
%     %PLOT UPDATE
%     hold on;
%     
%     %pull coverage map values out of data structure
%     dataCell = struct2cell(domainastar.coverageMap);
%     dataMatrix = cell2mat(dataCell);
%     lats = squeeze(dataMatrix(1,:,:));
%     longs = squeeze(dataMatrix(2,:,:));
%     coverages = squeeze(dataMatrix(3,:,:));
%     
%     %plot coverage map
%     contourf(longs, lats, coverages);
%     
%     if currentTime == domainastar.endTime
%         saveas(gcf, '32astarPathOverCoverage.fig');
%     elseif currentTime == minutes(days(15))
%         saveas(gcf, '15astarPathOverCoverage.fig');
%     end
%     
%     %plot boatastar location, goal location, and path traveled
%     plot(domainastar.coverageMap(paths(end)).longitude, domainastar.coverageMap(paths(end)).latitude, 'g*', 'MarkerSize', 20) % Path goal
%     plot(boatastar.longitude, boatastar.latitude, 'md', 'MarkerSize', 10) % Current vehicle location
%     plot(longListAstar, latListAstar, 'g-', 'LineWidth', 2); % Path to follow
%     plot([domainastar.coverageMap(paths).longitude], [domainastar.coverageMap(paths).latitude], 'r-', 'LineWidth', 2);
%     axis([-76.1,-75.1,33.5,34.5]);
%     
%     %Add each timestep as image in GIF
% %     drawnow
% %     frame = getframe(1);
% %     im = frame2im(frame);
% %     [imind,cm] = rgb2ind(im,256);
% %     if currentTime - domainastar.startTime + 1 == 1
% %         imwrite(imind,cm,filenameastar,'gif', 'DelayTime',0.1, 'Loopcount',inf);
% %     else
% %         imwrite(imind,cm,filenameastar,'gif','WriteMode','append', 'DelayTime',0.1);
% %     end
%     
%     
%     clf;
% end
% 
% % Path over time
% y = latListAstar(latListAstar ~= 0);
% x = longListAstar(longListAstar ~= 0);
% z = zeros(size(x));
% col = 1:length(x);  % This is the color, vary with x in this case.
% surface([x;x],[y;y],[z;z],[col;col],...
%     'facecol','no',...
%     'edgecol','interp',...
%     'linew',2);
% axis([-76.1,-75.1,33.5,34.5]);
% title('A* Path Over Time');
% saveas(gcf, 'astarpath.fig');
% 
% %Convert total coverage to percentage and plot
% figure;
% totalCoverageastar(totalCoverageastar == 0) = [];
% plot(totalCoverageastar./numel(domainastar.coverageMap));
% saveas(gcf, 'astarcoverage.fig');

%% Transect implementation
%Initialize Environment and Vehicles
domaintransect = Environment; % Mission Domain
boattransect = Vehicle(33.5,-76.1); % initialize boattransect at start lat/long
boattransect.latitude = boattransect.latitude + km2deg(5);
boattransect.longitude = boattransect.longitude + km2deg(5);
totalCoverageTransect = zeros(size(domaintransect.startTime:minutes(domaintransect.timeStep):domaintransect.endTime));

% Track boattransect path
latListTransect = nan(1,domaintransect.endTime - domaintransect.startTime);
longListTransect = nan(1,domaintransect.endTime - domaintransect.startTime);

% Initialize overall goal
% tic
% % List of upcoming waypoints
% [goal.lat, goal.long] = goalSelection(domaintransect, boattransect);
% overall.lat = goal.lat;
% overall.long = goal.long;
% toc

% Direction Leg used to determine part of transect strategy
% 0 - right; 1 - up; 2 - left; 3 - up
legCount = 0;
[goal.lat, goal.long] = transectWaypoint(legCount, boattransect, domaintransect);

for currentTime = domaintransect.startTime:minutes(domaintransect.timeStep):domaintransect.endTime
    %     tic
    %     fprintf('\nDist2Goal: %d', distance('gc',goal.lat, goal.long, boattransect.latitude, boattransect.longitude));
    %     fprintf('\nDist2Goal nm: %d\n', deg2nm(distance('gc',goal.lat, goal.long, boattransect.latitude, boattransect.longitude)));
    %     toc
    
    % Find new goal if sufficiently close to previous waypoint
    if deg2nm(distance('gc',goal.lat, goal.long, boattransect.latitude, boattransect.longitude)) <= wpThresh
        legCount = mod(legCount + 1, 4);
        [goal.lat, goal.long] = transectWaypoint(legCount, boattransect, domaintransect);
    end
    
    % Move boattransect towards goal
    boattransect = boattransect.moveBoat(domaintransect, speedConstant, goal.lat, goal.long);
    %     fprintf('\nTime: %d\tLat: %d\tLong:%d', currentTime-domaintransect.startTime+1, boattransect.latitude, boattransect.longitude);
    latListTransect(currentTime - domaintransect.startTime + 1) = boattransect.latitude;
    longListTransect(currentTime - domaintransect.startTime + 1) = boattransect.longitude;
    
    % Update Coverage
    [domaintransect, totalCoverageTransect(currentTime - domaintransect.startTime + 1)] = domaintransect.updateCoverage(boattransect.latitude, boattransect.longitude);
    %     fprintf('\ntotalCoverageTransect: %d\n', totalCoverageTransect(currentTime - domaintransect.startTime + 1));
    
    % PLOT UPDATE
    hold on;
    
    % pull coverage map values out of data structure
    dataCell = struct2cell(domaintransect.coverageMap);
    dataMatrix = cell2mat(dataCell);
    lats = squeeze(dataMatrix(1,:,:));
    longs = squeeze(dataMatrix(2,:,:));
    coverages = squeeze(dataMatrix(3,:,:));
    
    % plot coverage map
    contourf(longs, lats, coverages);
    
    if currentTime == domaintransect.endTime
        saveas(gcf, '32transectPathOverCoverage.fig');
    elseif currentTime == minutes(days(15))
        saveas(gcf, '15transectPathOverCoverage.fig');
    end
    
    % plot boattransect location, goal location, and path traveled
    plot(goal.long, goal.lat, 'r*', 'MarkerSize', 20) % Path goal
    plot(boattransect.longitude, boattransect.latitude, 'md', 'MarkerSize', 10) % Current vehicle location
    plot(longListTransect, latListTransect, 'r-', 'LineWidth', 2); % Path to follow
    % set(gca, 'ydir', 'reverse')
    axis([-76.1,-75.1,33.5,34.5]);
    
    % Add each timestep as image in GIF
%     drawnow
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if currentTime - domaintransect.startTime + 1 == 1
%         imwrite(imind,cm,filenametransect,'gif', 'DelayTime',0.1, 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filenametransect,'gif','WriteMode','append', 'DelayTime',0.1);
%     end
    
    
    clf;
    
end

% Path over time
y = latListTransect(latListTransect ~= 0);
x = longListTransect(longListTransect ~= 0);
z = zeros(size(x));
col = 1:length(x);  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2);
axis([-76.1,-75.1,33.5,34.5]);
title('Transect Path Over Time');
saveas(gcf, 'transectpath.fig');

% Convert total coverage to percentage and plot
figure;
plot(totalCoverageTransect./numel(domaintransect.coverageMap));
saveas(gcf, 'transectcoverage.fig');

%% Overall plots
clf;
hold on;
plot(totalCoveraged2p./numel(domaindtp.coverageMap), 'r');
plot(totalCoverageastar./numel(domainastar.coverageMap), 'g');
plot(totalCoverageTransect./numel(domaintransect.coverageMap), 'b');
ylabel('% total coverage');
xlabel('Time (min)');
legend('d2p', 'a*', 'transect');
saveas(gcf, 'totalcoverage.fig');

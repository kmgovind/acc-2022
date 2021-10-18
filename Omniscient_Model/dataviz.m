% close all;

%% Plot the net flow
% days = [1,10,20,30];
% 
% for day = days
%     day
%     currentDay = day;
%     flowus = domaindtp.enviroData.u(:,:, currentDay);
%     flowvs = domaindtp.enviroData.v(:,:, currentDay);
%     
%     
%     
%     [flowVels, flowHeads] = boatdtp.flowHeading(flowus, flowvs);
%     latitudeRange = domaindtp.latitudeRange;
%     longitudeRange = domaindtp.longitudeRange;
%     
%     overallLongs = -79:0.01:-74;
%     overallLats = 32:0.01:37;
%     
%     figure;
%     contourf(overallLongs, overallLats, flowVels)
%     title(['Flow Resource Day ' num2str(currentDay)]);
%     xlabel('Longitude');
%     ylabel('Latitude');
%     
%     hold on;
%     wall = zeros(1,numel(latitudeRange));
%     wall(:) = longitudeRange(1);
%     % plot(latitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
%     plot(wall, latitudeRange, 'Linewidth', 2, 'Color', 'r');
%     wall(:) = longitudeRange(end);
%     % plot(latitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
%     plot(wall, latitudeRange, 'Linewidth', 2, 'Color', 'r');
%     wall(:) = latitudeRange(1);
%     % plot(wall, longitudeRange, 'Linewidth', 2, 'Color', 'r');
%     plot(longitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
%     wall(:) = latitudeRange(end);
%     plot(longitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
%     
%     
% end
% 
% 
% 
% latindStart = find(overallLats == latitudeRange(1));
% latindEnd = find(overallLats == latitudeRange(end));
% longindStart = find(overallLongs == longitudeRange(1));
% longindEnd = find(overallLongs == longitudeRange(end));
% 
% imptu = flowus(latindStart:latindEnd, longindStart:longindEnd);
% imptv = flowvs(latindStart:latindEnd, longindStart:longindEnd);
% 
% imptVels = flowVels(latindStart:latindEnd, longindStart:longindEnd);
% imptHeadings = flowHeads(latindStart:latindEnd, longindStart:longindEnd);
% 
% % hold on;
% % contourf(imptVels);
% % quiver(imptu, imptv, 'm');
% % axis([0,100,0,100]);
% %
% % hold off;
% % figure;
% % contourf(flowVels);
% % axis([latitudeRange(1), latitudeRange(end), longitudeRange(1), longitudeRange(end)]);
% 
% %% Map of chosen area
% figure;
% % geobasemap satellite;
% geobasemap colorterrain;
% 
% hold on;
% wall = zeros(1,numel(latitudeRange));
% wall(:) = longitudeRange(1);
% geoplot(latitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
% % geoplot(wall, latitudeRange, 'Linewidth', 2, 'Color', 'r');
% wall(:) = longitudeRange(end);
% geoplot(latitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
% % geoplot(wall, latitudeRange, 'Linewidth', 2, 'Color', 'r');
% wall(:) = latitudeRange(1);
% geoplot(wall, longitudeRange, 'Linewidth', 2, 'Color', 'r');
% % geoplot(longitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
% wall(:) = latitudeRange(end);
% geoplot(wall, longitudeRange, 'Linewidth', 2, 'Color', 'r');
% % geoplot(longitudeRange, wall, 'Linewidth', 2, 'Color', 'r');
% geolimits([overallLats(1), overallLats(end)],[overallLongs(1), overallLongs(end)])
% 
% geoplot(34.22573,-77.94471,'r*');
% text(34.22573,-77.94471,'Wilmington');
% 
% geoplot(35.224622, -75.5301509,'r*');
% text(35.224622, -75.5301509 , 'Cape Hatteras');
% title('Mission Domain');
% 
% 
% % hold on;
% % plot(totalCoveraged2p./numel(domaindtp.coverageMap), 'r');
% % plot(totalCoverageastar./numel(domainastar.coverageMap), 'g');
% % plot(totalCoverageTransect./numel(domaintransect.coverageMap), 'b');
% % ylabel('% total coverage');
% % xlabel('Time (min)');
% % legend('d2p', 'a*', 'transect');
% 
% %% Test Plotting path vs time
% 
% % figure;
% %     y = latList(latList ~= 0);
% %     x = longList(longList ~= 0);
% %     z = zeros(size(x));
% %     col = 1:length(x);  % This is the color, vary with x in this case.
% %     surface([x;x],[y;y],[z;z],[col;col],...
% %         'facecol','no',...
% %         'edgecol','interp',...
% %         'linew',2);
% %     axis([-78,-77,32.5,33.5]);
% 
% 

%data is in 3 hour increments, so 8 datapoints per day
figure;
sunPowers = zeros(1,248);
xVals = 1:1:248;
for i = xVals
    sunPowers(i) = sun_ssr(25,25,i);
end
xVals = xVals/8;
plot(xVals,sunPowers);
title('Solar Irradiance vs. Time');
xlabel('Days');
ylabel('Solar Irradiance (W/m^2)');

figure;
avgPower = zeros(1,numel(xVals)/8);
for i = 1:numel(xVals)/8
    dates = 8*i - 7:1:8*i;
    avgPower(i) = mean(sun_ssr(25,25,dates));
end
stairs(1:numel(xVals)/8, avgPower);
title('Average Solar Irradiance per Day');
xlabel('Days');
ylabel('Average Solar Irradiance (W/m^2)');


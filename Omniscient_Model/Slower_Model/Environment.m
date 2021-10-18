classdef Environment
    % Properties of environment
    
    properties
        enviroData; % flow and sunlight data
        coverageMap; % initial coverage map
    end
    
    properties (Constant)
        % Time
        timeStep = minutes(1);
        startTime = minutes(days(1)); % Start at 1 Day
        endTime = minutes(days(32)); % End at 32 Day
%         endTime = minutes(days(3));
%         latitudeRange = 32:0.01:37;
%         longitudeRange = -79:0.01:-74;
%         latitudeRange = 32.5:0.01:33.5;
%         longitudeRange = -78:0.01:-77;

        latitudeRange = 33.5:0.01:34.5;
        longitudeRange = -76.1:0.01:-75.1;
    end
    
    %% Methods
    methods (Static)
        function obj = Environment()
            obj.enviroData = load('flowSunData3.mat');
            obj.coverageMap = repmat(struct('latitude', 0, 'longitude', 0, 'coverage', 0.001), length(obj.latitudeRange), length(obj.longitudeRange));
%             obj.flow = repmat(struct('latitude', 0, 'longitude', 0, 'speed', 0, 'heading', 0), length(
            
            
            for lats = 1:length(obj.latitudeRange)
                for longs = 1: length(obj.longitudeRange)
                    obj.coverageMap(lats, longs).latitude = obj.latitudeRange(lats);
                    obj.coverageMap(lats,longs).longitude = obj.longitudeRange(longs);
                end
            end
        end
    end
    
    methods
        function [flow_u, flow_v] = flowComponents(obj, latitude, longitude)
            
            % Interpolate data at required time
            global currentTime
            flowTime = days(minutes(currentTime));
            flow_u = interp3(obj.enviroData.flowlat, obj.enviroData.flowlon, obj.enviroData.flowt, obj.enviroData.u, latitude, longitude, flowTime);
            flow_v = interp3(obj.enviroData.flowlat, obj.enviroData.flowlon, obj.enviroData.flowt, obj.enviroData.v, latitude, longitude, flowTime);
            
%             % Convert from m/s to kts
%             flow_u = convvel(flow_u, 'm/s', 'kts');
%             flow_v = convvel(flow_v, 'm/s', 'kts');
        end
        
        function [obj, coverageSum] = updateCoverage(obj, boatLat, boatLong)
            
            % Initialize coverageSum
            coverageSum = 0;
            
            % Update coverage at each gridpoint
            dataCell = struct2cell(obj.coverageMap);
            dataMatrix = cell2mat(dataCell);
            lats = squeeze(dataMatrix(1,:,:));
            longs = squeeze(dataMatrix(2,:,:));
            coverages = squeeze(dataMatrix(3,:,:));
            
            
            % Coverage Loss
            % Currently 0.01% of current coverage
            lossConstant = 0.01/100; % 0.01 percent
            coverages = (1-lossConstant)*coverages;
            
            
            %             offset = deg2nm(distance('gc',gridLat, gridLong, boatLat, boatLong));
            
            nmPerLong = 60*cosd(34.5);
            nmPerLat = 60;
            offset = sqrt((nmPerLat*lats-nmPerLat*boatLat).^2 + (nmPerLong*longs - nmPerLong*boatLong).^2); %dist in nm
            
            
           % Update Coverages
%             coverages = max(coverages,1./(offset.^2+1));

%             exp(-dist^2/2*10^2)
            lengthScale = km2nm(10); %10 km length scale converted to nm
%             lengthScale = 10;
            coverages = max(coverages, exp( (-offset.^2)./(2*lengthScale^2)));
            
            
            % Update Coverage Sum
            % Compute coverage sum
            coverageSum = sum(coverages,'all');
            
            T = num2cell(coverages);
            [obj.coverageMap(:,:).coverage] = T{:};
            

            
            
            
        end
    end
end
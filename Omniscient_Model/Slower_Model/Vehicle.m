classdef Vehicle
    %% Properties of boat
    properties
        % Boat Position
        latitude; % Vehicle's current latitude
        longitude; % Vehicle's current longitude
        
        % Boat Net Speeds in kts
        speed; % Vehicle net speed (including flow)
        speed_u; % Vehicle net speed u component
        speed_v; % Vehicle net speed v component
        heading; % Vehicle's current heading relation to due north in degrees
        
        % Boat Motor Speeds in kts
        motorSpeed; % Vehicle's magnitude of speed from only motor
        motorSpeed_u; % Vehicle speed from motor in u direction
        motorSpeed_v; % Vehicle speed from motor in v direction
        
        % Battery Details
        charge; % Vehicle's current battery charge
        batteryCapacity = 6.75e3; % 6.75 kWh battery capacity
        propulsionDraw = 500; % 500W to run motor
    end
    
    
    %% Methods
    methods
        function obj = Vehicle(lat, long)
            % Assign Lat-Long
            obj.latitude = lat;
            obj.longitude = long;
            
            % Instantiate Net Speeds
            obj.heading = 0;
            obj.speed = 0;
            obj.speed_u = 0;
            obj.speed_v = 0;
            
            % Instantiate Motor Speeds
            obj.motorSpeed = 0;
            obj.motorSpeed_u = 0;
            obj.motorSpeed_v = 0;
            
            % Instantiate Battery
            obj.charge = obj.batteryCapacity;
            
        end
        
        function obj = moveBoat(obj, environment, speed, goalLat, goalLong)
            %obj is the current boat we are moving
            
            % Collect Boat assigned speed and heading towards goal
            obj.motorSpeed = speed;
            goalHeading = obj.headingCalc(goalLat, goalLong);
            
            % Pull flow components
            [flow_u, flow_v] = environment.flowComponents(obj.latitude, obj.longitude); % Pull flow components from environmental data
            if isnan(flow_u)
                flow_u = 0;
            end
            if isnan(flow_v)
                flow_v = 0;
            end
            [flowspeed, flowheading] = obj.flowHeading(flow_u, flow_v);
            
            % Compute heading that compensates for flow and points boat
            % towards goal
            obj.heading = goalHeading + asind(-(flowspeed/obj.motorSpeed) * sind(mod(flowheading - goalHeading, 360)));
            % TODO: Optimize lower-level strategy
            
            % Identify velocity components and add to flow components
            obj.motorSpeed_u = obj.motorSpeed * sind(obj.heading);
            obj.motorSpeed_v = obj.motorSpeed * cosd(obj.heading);

            obj.speed_u = obj.motorSpeed_u + flow_u;
            obj.speed_v = obj.motorSpeed_v + flow_v;
            
%             if obj.charge >= obj.propulsionDraw * (environment.timeStep/60)
                % Calculate resulting position only if battery can handle it
                % u moves longitude; v moves latitude
                obj.longitude = obj.longitude + nm2deg((obj.speed_u * hours(environment.timeStep)));
                obj.latitude = obj.latitude + nm2deg((obj.speed_v * hours(environment.timeStep)));
                
                % Update battery state
%                 obj.charge = useCharge(obj, environment);
%             end
            
        end
        function updatedCharge = useCharge(obj, environment)
            updatedCharge = obj.charge - obj.propulsionDraw * (hours(environment.timeStep));
            if updatedCharge <= 0
                updatedCharge = 0;
            end
        end
        
        function heading = headingCalc(obj, goalLat, goalLong)
            latDif = goalLat - obj.latitude;
            longDif = goalLong - obj.longitude;
            heading = mod(180 * atan2(longDif, latDif) / pi, 360);
        end
        
        function [speed, heading] = flowHeading(obj, flowu, flowv)
            speed = sqrt(flowu.^2 + flowv.^2);
            heading = mod(180 .* atan2(flowu, flowv)./pi, 360);
        end
        
    end
end
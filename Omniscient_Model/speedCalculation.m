% clc;
% % Winter it’s near 100 W/m^2 and summer it’s 220-240. Yearly average is ~175
% % Solar Irradiance integrated over time yields Energy J/m^2
% % Multiply Irradiance by solar panel area
% 
% %    1 kW/m2 × (24 h/day) = (24 kWh/m2)/day
% %     (24 kWh/m2)/day × (365 days/year) = (8760 kWh/m2)/year.
% 
% % 750 W @ 1 kW/m2 irradiance
% 
% percentage = 175/1000;
% wattGain = 750 * percentage;
% 
% fprintf('Wattage Gained from Solar Panels Assuming 175 W/m^2 irradiance: %f\n', wattGain);
% 
% % 4.5 kts at 500 W
% 
% %power = force * velocity
% 
% vel = convvel(4.5, 'kts', 'm/s');
% power = 500;
% 
% pushingForce = power/vel;
% 
% boatVel = wattGain/pushingForce
% convvel(boatVel, 'm/s', 'kts')


as = 0.75;
it = 175;
pe = 500;
aw = 8.688;
ns = 18;
cd = 0.0030;
v = 1.543;

nm = (2*(ns*as*it) - pe)/(aw*cd*v^3)


function [dHm, dm, xyz]=lldistkm2(latlon1,latlon2)
% format: [d1km d2km]=lldistkm(latlon1,latlon2)
% Distance:
% d1km: distance in km based on Haversine formula
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
% d2km: distance in km based on Pythagoras’ theorem
% (see: http://en.wikipedia.org/wiki/Pythagorean_theorem)
% After:
% http://www.movable-type.co.uk/scripts/latlong.html
%
% --Inputs:
%   latlon1: latlon of origin point [lat lon]
%   latlon2: latlon of destination point [lat lon]
%
% --Outputs:
%   d1km: distance calculated by Haversine formula
%   d2km: distance calculated based on Pythagoran theorem
%
% --Example 1, short distance:
%   latlon1=[-43 172];
%   latlon2=[-44  171];
%   [d1km d2km]=distance(latlon1,latlon2)
%   d1km =
%           137.365669065197 (km)
%   d2km =
%           137.368179013869 (km)
%   %d1km approximately equal to d2km
%
% --Example 2, longer distance:
%   latlon1=[-43 172];
%   latlon2=[20  -108];
%   [d1km d2km]=distance(latlon1,latlon2)
%   d1km =
%           10734.8931427602 (km)
%   d2km =
%           31303.4535270825 (km)
%   d1km is significantly different from d2km (d2km is not able to work
%   for longer distances).
%
% First version: 15 Jan 2012
% Updated: 17 June 2012
%--------------------------------------------------------------------------

r=6371000;
lat1=latlon1(:,1)*pi/180;
lat2=latlon2(:,1)*pi/180;
lon1=latlon1(:,2)*pi/180;
lon2=latlon2(:,2)*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;

x=r*deltaLon.*cos((lat1+lat2)/2);
y=r*deltaLat;
dm=sqrt(x.*x + y.*y); %Pythagoran distance
z=latlon2(:,3)-latlon1(:,3);
dHm=sqrt(x.^2+y.^2+z.^2);

xyz = [r*lon1.*cos(lat1) r*lat1 latlon1(:,3)];
end
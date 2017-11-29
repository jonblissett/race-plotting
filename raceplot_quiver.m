% bike-sim
% Copyright (C) 2017  Jonathan Blissett
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Contact, jonathan@blissett.me.uk

if (exist ('OCTAVE_VERSION', 'builtin') > 0)
    pkg load signal  % for octave compatibility
end

% ###### Import lap data and convert to cartesian ######
load Sun-Oct--2-15-20-15-2016   % Portimao circuit recording

%index1 = 12910:15390;
%index2 = 10516:12910;
index3 = (8050:10520)+1500;
%index = 5648:8123;
index = index3;
startline = [37.232093 -8.630912; 37.232136 -8.630722]*pi/180;
finishline = [37.229776 -8.630075; 37.229821 -8.629884]*pi/180;

r=6371000;  % earth

startline(:,2)=r*startline(:,2).*cos(startline(:,1));
startline(:,1)=r*startline(:,1);
finishline(:,2)=r*finishline(:,2).*cos(finishline(:,1));
finishline(:,1)=r*finishline(:,1);

nonz = lat ~= 0;

GPStime = GPStime(nonz);
lat = lat(nonz);
lon = lon(nonz);
N_sat = N_sat(nonz);
hDop = hDop(nonz);
h = h(nonz);
RTC = RTC(nonz);

latD = floor(lat/100);
latM = lat - 100*latD;
lat = latD+latM/60;
lonD = floor(lon/100);
lonM = lon - 100*lonD;
lon = lonD+lonM/60;

latlon1 = [lat lon h];
latlon2 = circshift(latlon1,1);

figure(1)
[dHm, dm, xyz]=lldistkm2(latlon1,latlon2);
plot(-xyz(index,1),xyz(index,2),startline(:,2),startline(:,1),finishline(:,2),finishline(:,1))
axis equal
rotate3d on
figure(2)
plot(dm(index)*20*2.23)
ylim([0 170])
xlabel('Distance (m)')
ylabel('Speed (mph)')

% ###### Make course 'map' ######

Portimao_map.dDist = dHm(index);
Portimao_map.dist = cumtrapz(dHm(index));
Portimao_map.h = h(index);
Portimao_map.x = xyz(index,1);
Portimao_map.y = xyz(index,2);
Portimao_map.z = xyz(index,3);

%figure(3);
%plot(Portimao_map.dist/1609.34,[Portimao_map.h Portimao_map.h(1)+cumtrapz(gradient(Portimao_map.h,Portimao_map.dist))]);
figure(3);
plot3(Portimao_map.x,Portimao_map.y,Portimao_map.z,Portimao_map.x,Portimao_map.y,zeros(size(Portimao_map.h)),'--')
xlabel('Position (m)')
ylabel('Position (m)')
zlabel('Elevation (m)')



d_spacing = 1/10*length(Portimao_map.dist);%66;%round(Portimao_map.dist(end)/245);
d_even = linspace(Portimao_map.dist(1),Portimao_map.dist(end),d_spacing)';
h_even = interp1(Portimao_map.dist,Portimao_map.h,d_even,'pchip');

%grad = diff([Portimao_map.h' Portimao_map.h(1)]')./Portimao_map.dDist;
grad = gradient(Portimao_map.h,Portimao_map.dist);
angle = 180/pi*atan(grad);

%grad_even = diff([h_even' h_even(1)]')/d_spacing;
grad_even = gradient(h_even,d_even);%gradient(filtfilt(b,a,h_even));
angle_even = 180/pi*atan(grad_even);

[b, a]=butter(3, 0.1/4);                  
filtangle = filtfilt(b,a,angle);
filtgrad = filtfilt(b,a,grad_even);%tan(angle_even*pi/180);

figure(4);
plot(Portimao_map.dist,angle,d_even,angle_even,Portimao_map.dist,filtangle)
xlabel('Distance (m)')
ylabel('Elevation change (m)')

Portimao_map.gradient = interp1(d_even,grad_even,Portimao_map.dist,'pchip');

figure(5);
plot(Portimao_map.dist,[Portimao_map.h Portimao_map.h(1)+cumtrapz(Portimao_map.dist,Portimao_map.gradient)])
xlabel('Distance (m)')
ylabel('Elevation (m)')


% ###### Smoothing ######

fGPS = 20;
tGPS = 1/fGPS*(0:(length(dm(index))-1))+1.15;
%vGPS = dm(index)*fGPS;
dm(isnan(dm)) = 0;      % NASTY!!!!!!
fvGPS = filtfilt(b,a,dm*fGPS);

% ###### Calculate lean angle ######
% Approximation based on resolving lateral and vertical acceleration
% ref, https://www.reddit.com/r/EngineeringStudents/comments/3i7v3g/calculate_lateral_acceleration_from_gps_data/

heading = atan2(gradient(xyz(:,1)),gradient(xyz(:,2)));

heading1 = heading;
heading1(heading1 > +pi/2) = heading1(heading1 > +pi/2) - 2*pi;
%heading1(heading1 < -pi/2) = heading1(heading1 < -pi/2) + pi;

dheading = gradient(heading);
dheading(dheading > +pi/2) = dheading(dheading > +pi/2) - pi;
dheading(dheading < -pi/2) = dheading(dheading < -pi/2) + pi;
w = dheading * fGPS; % dtheta/dt
w(isnan(w)) = 0;      % NASTY!!!, but matlab can't handle NaN in filt()
GPSFLat = fvGPS .* filtfilt(b,a,w);

%plot(Portimao_map.dist,GPSFLat)
lean = atan(GPSFLat/9.81);

lean = lean(index);
fvGPS = fvGPS(index);
dheading = dheading(index);
w = w(index);

Portimao_map.lean = lean;
Portimao_map.heading = cumtrapz(filtfilt(b,a,dheading));
Portimao_map.w = filtfilt(b,a,w);

figure(6)
plot(Portimao_map.dist,lean*180/pi, Portimao_map.dist, atan(filtfilt(b,a,Portimao_map.dDist*fGPS).*filtfilt(b,a,[0 diff(Portimao_map.heading)'.*fGPS/9.81]'))*180/pi)
xlabel('Distance (m)')
ylabel('Calculated lean angle (deg)')

figure(7)
plot(Portimao_map.dist,Portimao_map.heading)
xlabel('Distance (m)')
ylabel('Heading (rad)')

%figure(8)
%plot(Portimao_map.dist,fvGPS)

figure(8)
plot(Portimao_map.dist,filtfilt(b,a,w))
xlabel('Distance (m)')
ylabel('Rate of change of heading (rad/s)')
%
if (exist ('OCTAVE_VERSION', 'builtin') > 0)
    %clear -x Portimao_map RaceDM2 locsMin % housekeeping
    %save'Portimao_map.mat'
    %break
end

circ_correct = (2.16-2*pi*0.095*(1-cos(lean)));

figure(9)
plot(circ_correct/2/pi)
xlabel('Distance (m)')
ylabel('Tyre effective radius duo to lean (m)')

% ###### Quiver plotting, 3D ######

% input data
[b, a]=butter(3, 0.1);                  
x = Portimao_map.x+startline(4);
y = Portimao_map.y-startline(1);
z = Portimao_map.z;
lean_ang = lean;
fv = fvGPS;

% downsample, else arrows are too close
divsc = 5;
x = -decimate(x,divsc);
y = decimate(y,divsc);
z = decimate(z,divsc);
fv = decimate(fv,divsc);
lean_ang = decimate(lean_ang,divsc);

on = ones(length(x),1);
zo = zeros(length(x),1);
u = circshift(x,1)-x;
v = circshift(y,1)-y;
w = circshift(z,1)-z;
mag = sqrt(u.^2+v.^2+w.^2);
u = u./mag;
v = v./mag;
w = w./mag;
k = [u v w];    % unit vector of vehicle direction of motion

% define unit vector [0 0 1] for each point, and rotate about k to show lean angle
vec = [zo zo on];
cosTh = [cos(lean_ang) cos(lean_ang) cos(lean_ang)];
sinTh = [sin(lean_ang) sin(lean_ang) sin(lean_ang)];
%vecr = vec.*cosTh+cross(k,vec).*sinTh+k.*dot(k,vec).*(ones(size(cosTh))-cosTh);
vecr = vec.*cosTh+cross(k,vec).*sinTh+k.*repmat(dot(k,vec),length(k),1).*(ones(size(cosTh))-cosTh);

uc = vecr(:,1);
vc = vecr(:,2);
wc = vecr(:,3); 
figure(10)
quiver3(x,y,z,uc,vc,wc)
hold on
plot3(x,y,z,'r')
set(gca,'DataAspectRatio',[1 1 1])
rotate3d on
xlabel('Position (m)')
ylabel('Position (m)')
zlabel('Elevation (m)')

% ###### Quiver plotting, 2D with speed ######
fontsize = 16;
figure(11)
hold on
plot(-y,x,'k')
axis equal
hh = quiverwcolorbar2_big(-y,x,-vc,uc,50,fv);
set(gca,'linewidth', 2, 'fontsize', fontsize)
cbh = findobj( gcf(), 'tag', 'colorbar');   % octave only?
ti = round(get(cbh,'yTick'));
set(cbh, 'linewidth', 2,'yticklabel', ti, 'FontSize', 16)
%'tickdir', 'out',
%'ylabel','v (m/s)',
%'yLim', [15, 58],
%set(cbh, 'ytick', 0:10:60 , 'yticklabel', 0:10:60)

xlabel('Position (meters)', 'fontsize', fontsize)
ylabel('Position (meters)', 'fontsize', fontsize)
yp = get(gca(),'ylim');
xp = get(gca(),'xlim');
text(xp(2) + diff(xp)*0.1,yp(1) - diff(yp)*0.125,'v (m/s)','fontsize',16)

if (exist ('OCTAVE_VERSION', 'builtin') > 0)
    print -deps -color Portimao_color_quiver_R3.eps
    print('-dpng','-color','Portimao_color_quiver_R3.png')
end
% ###### Not-really quiver plotting, 2D with e.g. lean angle ######
figure(12)
hold on
axis equal
h = quiverwcolorbar2(y,x,v,u,20,abs(lean_ang)*180/pi);
xlabel('Distance (meters)')
ylabel('Distance (meters)')
yp = get(gca(),'ylim');
xp = get(gca(),'xlim');
text(xp(2) + diff(xp)*0.05,yp(1) - diff(yp)*0.125,'lean (�)')

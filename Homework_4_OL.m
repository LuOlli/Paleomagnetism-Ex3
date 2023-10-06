%% 1
clear all
% A12A 63.55 26.57 346.7 32.1 1.4
% A12B 63.55 26.57 345.1 33.0 1.4
% A12C 63.55 26.57 342.1 30.0 1.0
% A12D 63.55 26.57 344.3 45.4 4.5
% A12E 63.55 26.57 342.4 32.8 0.9
% A12F 63.55 26.57 339.6 35.0 0.6

slat = 63.55;
slon = 26.57;
D = [346.7 345.1 342.1 344.3 342.4 339.6];
I = [32.1 33.0 30.0 45.4 32.8 35.0];
alpha95 = [1.4 1.4 1.0 4.5 0.9 0.6];
li = cosd(I).*cosd(D);
mi = cosd(I).*sind(D);
ni = sind(I);
R = sqrt(sum(li)+sum(mi).^2+sum(ni).^2);%The resultant vector R
l = sum(li)/R;
m = sum(mi)/R;
n = sum(ni)/R;
% Mean declination:
Dm = atan2(m,l)*(180/pi)+360;
% Mean inclination
Im = asind(n);
N= length(D);
% Precision parameter
k = (N-1)/(N-R);
% 95% confidence limit
p=0.05;
%Alpha95 = acosd(1-(N-R)/R*(((1/p)^(1/(N-1))-1)));
A95 = 140/sqrt(k*N);
disp('     D          I          R        k        a95')
disp([Dm Im R k A95])
%% 2
VGP_mean = [4.7 166.5];
A95 = 6.3;
plat = -26;
%VGP coordinates
Plon = [161.8 175 166.5 174.8 180.9 183.2 197.1 197.1 161.5 157.4 164.5 ...
    186.5 183.5 133.7 158.3 153.2 158.6 149.6 152.6 141 149.1 173 159.2 ...
    174 158.5 160.7 174 181.4]';%degs
Plat = [9.4 6.4 9.4 10 10.9 11.5 13.9 2.8 8.1 3.5 -5.7 35.8 22.6 -3.5 7.1 -7.8 3.9...
    -10.1 1.3 -16.5 -2.4 5.2 -1.7 4.4 -1.3 0.5 1.7 8.6]';%degs
Amin = 82*28^(-0.63)
Amax = 12*28^(-0.40)
cosDelta = cosd(Plat).*cosd(Plon).*cosd(VGP_mean(1)).*cosd(VGP_mean(2))+...
    cosd(Plat).*sind(Plon).*cosd(VGP_mean(1)).*sind(VGP_mean(2))+...
    sind(Plat).*sind(VGP_mean(1));
Deltai = acosd(cosDelta);
N = length(Plon)
S_T  =sqrt(1/(N-1)*sum(Deltai.^2))
%% 3 Euler rotations
latE = 47.5;%°𝑁
lonE = 1.5;%°𝐸
angE = 49.0;%° (CCW)
age  = 1452e6;%a
latS = 15.2;%°𝑁
lonS = 177.1;%°𝐸
Pp = [cosd(lonS)*cosd(latS) sind(lonS)*cosd(latS) sind(latS)]';% Site location in cartesian
PE = [cosd(lonE)*cosd(latE) sind(lonE)*cosd(latE) sind(latE)]';% Pole location in cartesian
%Rotation matrix R 
R11 = (PE(1)*PE(1))*(1-cosd(angE))+cosd(angE);
R12 = (PE(1)*PE(2))*(1-cosd(angE))-PE(3)*sind(angE);
R13 = (PE(1)*PE(3))*(1-cosd(angE))+PE(2)*sind(angE);
R21 = (PE(2)*PE(1))*(1-cosd(angE))+PE(3)*sind(angE);
R22 = (PE(2)*PE(2))*(1-cosd(angE))+cosd(angE);
R23 = (PE(2)*PE(3))*(1-cosd(angE))-PE(1)*sind(angE);
R31 = (PE(3)*PE(1))*(1-cosd(angE))-PE(2)*sind(angE);
R32 = (PE(3)*PE(2))*(1-cosd(angE))+PE(1)*sind(angE);
R33 = (PE(3)*PE(3))*(1-cosd(angE))+cosd(angE);
R = [R11 R12 R13;
     R21 R22 R23;
     R31 R32 R33];
Ppr = R*Pp;%Rotated coordinates in cartesian
latSr = asind(Ppr(3));
lonSr = atan2(Ppr(2),Ppr(1))*(180/pi);
if lonSr<0
    lonSr = 360+lonSr;
end
disp('Latitude  Longitude')
disp([latSr lonSr])
%% 4
% R = 6370e3;%m 
% age1 = 1385e6;%a
% lat1 = -19.2;%°𝑁
% % lon1 = 
% 
% lon_p=-114.11
% lon_n=-61.85
% lat_p = 47.09
% lat_n = 56.53
% % 47.09,-114.11 pilcher
% % 56.53,-61.85 Nain
% cosDelta = cosd(lat_p).*cosd(lon_p).*cosd(lat_p).*cosd(lon_n)+...
%     cosd(lat_p).*sind(lon_p).*cosd(lat_p).*sind(lon_n)+...
%     sind(lat_p).*sind(lat_n);
% Deltai = acosd(cosDelta);
% d  =(Deltai/360)*2*pi*R

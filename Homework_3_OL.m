%MATLAB script for calculating the problems in Paleomagnetism 2023 EX3

%1.a)
disp('1')
D = [278 33 35 108]';%degrees, Declinations for i to iv
I = [17 -41 81 -43]';%degrees, inclinations for i to iv
Lat = [57 55 -15 -69]';%Latitudes
Lon = [-62 15 -59 78]';%longitudes
%Conversion to cartesian:
x1= cosd(D).*cosd(I);
x2= sind(D).*cosd(I);
x3= sind(I);
%Need to calculate these first:
Id = atand(2*tand(Lat));%Expected inclination from GAD field
% inclination of the paleofield vector projected onto the N-S plane:
alpha = atan2(x3,x1)*(180/pi);%this function always outputs radians, need to fix
%Rotation to (x1',x2',x3')
x1prime = sind(Id-alpha).*(x1.^2+x3.^2).^(0.5);
x2prime = x2;
x3prime = cosd(Id-alpha).*(x1.^2+x3.^2).^(0.5);
%Conversion
Bprime = sqrt(x1prime.^2+x2prime.^2+x3prime.^2)%Funny how this became an unit vector
Dprime = atan2(x2prime,x1prime)*(180/pi)
Iprime = atan2(x3prime,Bprime)*(180/pi)

%b)
%conventions:
phi_s = [360-62 15 360-59 78]';%longitudes in 0-360 format
th_s = 90-Lat; %observation site co-latitude
%!!!!Declinations between 180∘ and 360∘ are equivalent to D - 360 ∘ which
%are counter-clockwise with respect to North!!!
DD = [360-278 33 35 108]';%"conventional" declinations (I don't understand why)
%magnetic co-latitude
th_m=acotd(0.5*tand(I));
th_p = acosd((cosd(th_s).*cosd(th_m)+sind(th_s).*sind(th_m).*cosd(DD)));
%VGP latitude
lambda_p = 90-th_p
%angular difference between the pole and site longitude
Deltaphi =asind(sind(th_m).*(sind(DD)./sind(th_p)))
%If cosθm ≥ cosθs cosθp, then ϕp = ϕs + Δϕ. 
%If cosθm < cosθs cosθp then ϕp = ϕs + 180 − Δϕ. 
con = cosd(th_m)>=cosd(th_s).*cosd(th_p);
phi_p = con.*(phi_s+Deltaphi)+~con.*(phi_s+180-Deltaphi)
pos = [lambda_p phi_p]

%% 2.a)
clear all %get rid of p1 values
Age = [1535 1709 1800 1258 1885 1855 1741 1382 1397 1769]';%Ma
slat = [2.2 -14.3 61.3 62.4 16.4 -25.2 67.5 81.5 70.4 39.0]';%degs N
slon = [298.9 133.2 30.0 18.1 78.9 30.6 242.0 315.3 104.6 114.0]';%degs E
D = [132.2 100.4 347.0 43.0 129.1 126.9 136.5 265.2 223.6 12.5]';
I = [35.4 42.9 40.6 -40.1 4.2 -58.7 57.4 21.5 -20.9 -3.1]';
%procedure just copied from p1.
th_s = 90-slat;%observation site co-latitude 
D(D>180) = D(D>180)-360;%Declinations to -180 - 180 format
%magnetic co-latitude
th_m=acotd(0.5*tand(I));
th_p = acosd((cosd(th_s).*cosd(th_m)+sind(th_s).*sind(th_m).*cosd(D)));
%VGP latitude
lambda_p = 90-th_p;
%angular difference between the pole and site longitude
Deltaphi =asind(sind(th_m).*(sind(D)./sind(th_p)));
%If cosθm ≥ cosθs cosθp, then ϕp = ϕs + Δϕ. 
%If cosθm < cosθs cosθp then ϕp = ϕs + 180 − Δϕ. 
con = cosd(th_m)>=cosd(th_s).*cosd(th_p);
phi_p = con.*(slon+Deltaphi)+~con.*(slon+180-Deltaphi);
phi_p(phi_p>360) = phi_p(phi_p>360)-360;
pos = [lambda_p phi_p];
disp('2a): Virtual geomagnetic poles with GAD:')
disp('Latitude    Longitude')
disp(pos)
%2b)
G2 = 0;
G3 = 0.23;
%Easier to type this in parts:
a=G2*((9/2)*(cosd(th_m).^2-3/2));
b=G3*(10*cosd(th_m).^3-6*cosd(th_m));
c=G2*(3*cosd(th_m).*sind(th_m));
d=G3*((15/2)*cosd(th_m).^2.*sind(th_m)-3*sind(th_m));
tanI_obs = (2*cosd(th_m)+a+b)./(sind(th_m)+c+d);
I_obs = atan2((2*cosd(th_m)+a+b),(sind(th_m)+c+d))*(180/pi);%new inclinations
I_obs(I_obs>90) = I_obs(I_obs>90)-180;%inclinations to -90 - +90
%VGP calculation:
th_m2=acotd(0.5*tand(I_obs));
th_p2 = acosd((cosd(th_s).*cosd(th_m2)+sind(th_s).*sind(th_m2).*cosd(D)));
%VGP latitude
lambda_p2 = 90-th_p2;
Deltaphi2 =asind(sind(th_m2).*(sind(D)./sind(th_p2)));
clear con
con = cosd(th_m2)>=cosd(th_s).*cosd(th_p2);
phi_p2 = con.*(slon+Deltaphi2)+~con.*(slon+180-Deltaphi2);
phi_p2(phi_p2>360) = phi_p2(phi_p2>360)-360;
pos2 = [lambda_p2 phi_p2];
disp('2b): Virtual geomagnetic poles with octupole model:')
disp('Latitude    Longitude')
disp(pos2)
%c)
cosDelta = cosd(lambda_p).*cosd(phi_p).*cosd(lambda_p2).*cosd(phi_p2)+...
    cosd(lambda_p).*sind(phi_p).*cosd(lambda_p2).*sind(phi_p2)+...
    sind(lambda_p).*sind(lambda_p2);
Delta = acosd(cosDelta);
disp('2c): Arc distances:')
disp(Delta)

%% 4
disp('4')
clear all
lam_p = 23.7;%degs
phi_p = 191.4;%degs
lam_s = 60.0;%degs
phi_s = 20.0;%degs
%colatitudes
th_p=90-lam_p;
th_s=90-lam_s;
th_m = acosd(cosd(lam_p).*cosd(phi_p).*cosd(lam_s).*cosd(phi_s)+...
    cosd(lam_p).*sind(phi_p).*cosd(lam_s).*sind(phi_s)+...
    sind(lam_p).*sind(lam_s));
I = atand(2*cotd(th_m))
cosD = (cosd(th_p)-cosd(th_s)*cosd(th_m))/(sind(th_s)*sind(th_m));
D = acosd(cosD)
















function [ b,Mb ] = accel_error(f_input,phi,h,d_phi,d_theta)
%ACCEL_ERROR Function
%   Die accelerrorfunktion berechnet die BIAS-Fehler, den Skalefaktor und
%   die Misaligmentmatrix.
%   Uebergabeparameter ist eine Matrix, bestehend aus 6
%   symmetrischen Positionen. Die Positionen haben folgende Reihenfolge:
%
%       Pos1: (0,0,0)       (x,y,z)
%       Pos2: (0,180,90)    (x,y,z)
%       Pos3: (90,0,90)     (x,y,z)
%       Pos4: (-90,0,180)   (x,y,z)
%       Pos5: (-90,-90,0)   (x,y,z)
%       Pos6: (180,90,0)    (x,y,z)
%
% Die Parameter muessen bereits gemittelt sein, damit eine Berechnung
% durchgefuehrt werden kann.
%
%[b,Mb]=accel_error(F_in,phi,h,d_phi,d_theta)
%
%   F_in=[fpose1x;fpose1y;fpose1z,fpose2x;fpose2y;fpose2z, ...]
%   phi => Breitengrad  (Rotationstisch)                Default=53[°]
%   h   => Hoehe        (Rotationstisch)
%   d_phi=> Stationaere Abweichung vom Rotationstisch
%   d_theta=>Stationaere Abweichung vom Rotationstisch
switch nargin
    case 4
        d_theta=0;
    case 3
        d_theta=0;
        d_phi=0;
    case 2
        d_theta=0;
        d_phi=0;
        h=0;
    case 1
        d_theta=0;
        d_phi=0;
        h=0;
        phi=53;
end
a=6378137;
f=1/298.57223563;
e=sqrt(f*(2-f));
Rn=a*(1-e^2)/(1-e^2*(sind(phi)^2))^(3/2);
Re=a*(1/(sqrt(1-e^2*(sind(phi)^2))));
R0=sqrt(Re*Rn);
g0=9.780318;
g1=5.3024e-3;
g2=5.9e-6;
g=g0*(1+g1*(sind(phi)^2+g2*(sind(2*phi)^2)*(R0/(R0-h))^2));

%Positionen generieren
pose.pose1=rotxyz(0+d_phi,0+d_theta,0);
pose.pose2=rotxyz(0+d_phi,180+d_theta,90);
pose.pose3=rotxyz(90+d_phi,0+d_theta,90);
pose.pose4=rotxyz(-90+d_phi,0+d_theta,180);
pose.pose5=rotxyz(-90+d_phi,-90+d_theta,0);
pose.pose6=rotxyz(180+d_phi,90+d_theta,0);
%Gemessene Wertepare zusammenfassen

fb1=[f_input(1,1);f_input(2,1);f_input(3,1)];
fb2=[f_input(1,2);f_input(2,2);f_input(3,2)];
fb3=[f_input(1,3);f_input(2,3);f_input(3,3)];
fb4=[f_input(1,4);f_input(2,4);f_input(3,4)];
fb5=[f_input(1,5);f_input(2,5);f_input(3,5)];
fb6=[f_input(1,6);f_input(2,6);f_input(3,6)];
%Referenzwerte Berechnen

A1=[eye(3,3),kron(fb1',eye(3,3))];
A2=[eye(3,3),kron(fb2',eye(3,3))];
A3=[eye(3,3),kron(fb3',eye(3,3))];
A4=[eye(3,3),kron(fb4',eye(3,3))];
A5=[eye(3,3),kron(fb5',eye(3,3))];
A6=[eye(3,3),kron(fb6',eye(3,3))];
Y1=pose.pose1*[0;0;g]-fb1;
Y2=pose.pose2*[0;0;g]-fb2;
Y3=pose.pose3*[0;0;g]-fb3;
Y4=pose.pose4*[0;0;g]-fb4;
Y5=pose.pose5*[0;0;g]-fb5;
Y6=pose.pose6*[0;0;g]-fb6;
A=[A1;A2;A3;A4;A5;A6];
Y=[Y1;Y2;Y3;Y4;Y5;Y6];

eta1=inv(A'*A)*A'*Y;            %Eta berechnen A*eta1=Y
b=[eta1(1);eta1(2);eta1(3)]*3;    %Biaswerte extrahieren
Mb=(eta1(4:12));                %Mb Matrix extrahieren
Mb=reshape(Mb,[3 3]);           %Mb Matrix zuruecktransponieren

end


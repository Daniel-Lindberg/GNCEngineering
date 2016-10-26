function [ b,Mb ] = gyro_error( f_omega,phi,rot,d_phi,d_theta)
%GYRO_ERROR gyro_error
%   Die Gyroscopeerrorfunktion berechnet die BIAS-Fehler, den Skalefaktor und
%   die Misaligmentmatrix.
%   Uebergabeparameter ist eine Matrix, bestehend aus 6
%   symmetrischen Rotationspositionen. Die Positionen haben folgende Reihenfolge:
%
%       Pos1: (0,0,0)       (x,y,z)
%       Pos2: (0,180,90)    (x,y,z)
%       Pos3: (90,0,90)     (x,y,z)
%       Pos4: (-90,0,180)   (x,y,z)
%       Pos5: (-90,-90,0)   (x,y,z)
%       Pos6: (180,90,0)    (x,y,z)
%
% Die Parameter muessen bereits gemittelt sein, damit eine Berechnung
% durchgefuehrt werden kann. Dabei wird immer um die Z-Achste gedreht.
%
%[b,Mb]=gyro_error(f_omega,phi)
%
%   f_omega=[f_omega_pose1x;f_omega_pose1y;f_omega_pose1z,f_omega_pose2x;
%   f_omega_pose2y;f_omega_pose2z, ...]
%   phi => Breitengrad  (Rotationstisch)    default=53[°]
%   rot => Rotationsgeschwindigkeit [°/s] , default=40°/s
switch nargin    
    case 4
        d_theta=0;
    case 3
        d_theta=0;
        d_phi=0;
    case 2
        d_phi=0;
        d_theta=0;
        rot=40;
    case 1
        rot=40;
        phi=53;
        d_phi=0;
        d_theta=0;
end
omega_earth=7.292115e-5;

omega=[cosd(phi);0;-sind(phi)]*omega_earth;
omega_ie=omega*180/pi;
omega_turn_table=[0;0;rot];

%Rotationsvektoren
pose.pose1=rotxyz(0+d_phi,0+d_theta,0);
pose.pose2=rotxyz(0+d_phi,180+d_theta,90);
pose.pose3=rotxyz(90+d_phi,0+d_theta,90);
pose.pose4=rotxyz(-90+d_phi,0+d_theta,180);
pose.pose5=rotxyz(-90+d_phi,-90+d_theta,0);
pose.pose6=rotxyz(180+d_phi,90+d_theta,0);
%Rotationsgravitationsvektoren

Cpose1=pose.pose1*(omega_ie+omega_turn_table);
Cpose2=pose.pose2*(omega_ie+omega_turn_table);
Cpose3=pose.pose3*(omega_ie+omega_turn_table);
Cpose4=pose.pose4*(omega_ie+omega_turn_table);
Cpose5=pose.pose5*(omega_ie+omega_turn_table);
Cpose6=pose.pose6*(omega_ie+omega_turn_table);


omega_ib1=[f_omega(1,1);f_omega(2,1);f_omega(3,1)];
omega_ib2=[f_omega(1,2);f_omega(2,2);f_omega(3,2)];
omega_ib3=[f_omega(1,3);f_omega(2,3);f_omega(3,3)];
omega_ib4=[f_omega(1,4);f_omega(2,4);f_omega(3,4)];
omega_ib5=[f_omega(1,5);f_omega(2,5);f_omega(3,5)];
omega_ib6=[f_omega(1,6);f_omega(2,6);f_omega(3,6)];

%Fehlerberechnung
A1=[eye(3,3),kron(omega_ib1',eye(3,3))];
A2=[eye(3,3),kron(omega_ib2',eye(3,3))];
A3=[eye(3,3),kron(omega_ib3',eye(3,3))];
A4=[eye(3,3),kron(omega_ib4',eye(3,3))];
A5=[eye(3,3),kron(omega_ib5',eye(3,3))];
A6=[eye(3,3),kron(omega_ib6',eye(3,3))];
Y1=Cpose1-omega_ib1;
Y2=Cpose2-omega_ib2;
Y3=Cpose3-omega_ib3;
Y4=Cpose4-omega_ib4;
Y5=Cpose5-omega_ib5;
Y6=Cpose6-omega_ib6;
A=[A1;A2;A3;A4;A5;A6];
Y=[Y1;Y2;Y3;Y4;Y5;Y6];
eta1=inv(A'*A)*A'*Y;            %Eta berechnen A*eta1=Y
b=[eta1(1);eta1(2);eta1(3)]*3;    %Biaswerte extrahieren
Mb=(eta1(4:12));                %Mb Matrix extrahieren
Mb=reshape(Mb,[3 3]);           %Mb Matrix zurücktransponieren
end


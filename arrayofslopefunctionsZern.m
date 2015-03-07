function [dx,dy] = arrayofslopefunctionsZern(pupildiam, x, y, maxterm)
% function computes slopes in x and y directions based on the Zernike
% polynomials
%
% copyright by:
% Steffen Mauch, (c) 03/2015
% email: steffen.mauch (at) gmail.com
%
% You can redistribute it and/or modify it under the terms of the GNU 
% General Public License as published by the 
% Free Software Foundation, version 2.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

if( maxterm > 28 )
   error( 'maxterm must be smaller or equal 28!') 
end

x = x / (pupildiam/2);
y = y / (pupildiam/2);

dx = zeros(length(x), maxterm);
dy = zeros(length(y), maxterm);

zero = zeros(size(x));
one  = ones(size(x));

for whichone = 1:maxterm
    switch whichone
        case 1
            dx(:,whichone)=zero;
            dy(:,whichone)= zero;
        case 2
            dx(:,whichone)=zero; 
            dy(:,whichone)= one;
        case 3
            dx(:,whichone)= one;
            dy(:,whichone)= zero;
        case 4
            dx(:,whichone)= 2.*y;
            dy(:,whichone)= 2.*x;
        case 5
            dx(:,whichone)= 4.*x; 
            dy(:,whichone)= 4.*y;
        case 6
            dx(:,whichone)= 2.*x; 
            dy(:,whichone)= -2.*y;
        case 7
            dx(:,whichone)= 6.*x.*y; 
            dy(:,whichone)= 3.*realpow(x,2) -3.*realpow(y,2);
        case 8
            dx(:,whichone)= 6.*x.*y; 
            dy(:,whichone)= 3.*realpow(x,2) +9.*realpow(y,2) -2;
        case 9
            dx(:,whichone)= 9.*realpow(x,2) + 3.*realpow(y,2) -2; 
            dy(:,whichone)= 6.*x.*y;
        case 10
            dx(:,whichone)= 3.*realpow(x,2) -3.*realpow(y,2); 
            dy(:,whichone)= -6.*x.*y;
        case 11
            dx(:,whichone)= 12.*realpow(x,2).*y -4.*realpow(y,3); 
            dy(:,whichone)= 4.*realpow(x,3) -12.*x.*realpow(y,2);
        case 12
            dx(:,whichone)= 24.*realpow(x,2).*y +8.*realpow(y,3) -6.*y; 
            dy(:,whichone)= 8.*realpow(x,3) +24.*x.*realpow(y,2) -6.*x;
        case 13
            dx(:,whichone)= 24.*realpow(x,3) +24.*x.*realpow(y,2) -12.*x; 
            dy(:,whichone)= 24.*realpow(x,2).*y +24.*realpow(y,3) -12.*y;
        case 14
            dx(:,whichone)= 16.*realpow(x,3) +8.*x.*realpow(y,2) -6.*x -8.*x.*realpow(y,2); 
            dy(:,whichone)= 8.*realpow(x,2).*y -8.*realpow(x,2).*y -16.*realpow(y,3) +6.*y;
        case 15
            dx(:,whichone)= 4.*realpow(x,3) -12.*x.*realpow(y,2); 
            dy(:,whichone)= -12.*realpow(x,2).*y +4.*realpow(y,3);
        case 16
            dx(:,whichone)= 20.*realpow(x,3).*y -20.*x.*realpow(y,3); 
            dy(:,whichone)= 5.*realpow(x,4) -30.*realpow(x,2).*realpow(y,2) +5.*realpow(y,4);
        case 17
            dx(:,whichone)= 60.*realpow(x,3).*y +30.*x.*realpow(y,3) -24.*x.*y -10.*x.*realpow(y,3); 
            dy(:,whichone)= 15.*realpow(x,4) +45.*realpow(x,2).*realpow(y,2) -12.*realpow(x,2) -15.*realpow(x,2).*realpow(y,2) -25.*realpow(y,4) +12.*realpow(y,2);
        case 18
            dx(:,whichone)= 40.*realpow(x,3).*y +40.*x.*realpow(y,3) -24.*x.*y; 
            dy(:,whichone)= 10.*realpow(x,4) +60.*realpow(x,2).*realpow(y,2) +50.*realpow(y,4) -12.*realpow(x,2) -36.*realpow(y,2) +3;
        case 19
            dx(:,whichone)= 50.*realpow(x,4) +60.*realpow(x,2).*realpow(y,2) +10.*realpow(y,4) -36.*realpow(x,2) -12.*realpow(y,2) +3; 
            dy(:,whichone)= 40.*realpow(x,3).*y +40.*x.*realpow(y,3) -24.*x.*y;
        case 20
            dx(:,whichone)= 25.*realpow(x,4) +15.*realpow(x,2).*realpow(y,2) -12.*realpow(x,2) -45.*realpow(x,2).*realpow(y,2) -15.*realpow(y,4) +12.*realpow(y,2); 
            dy(:,whichone)= 10.*realpow(x,3).*y -30.*realpow(x,3).*y -60.*x.*realpow(y,3) +24.*x.*y;
        case 21
            dx(:,whichone)= 5.*realpow(x,4) -30.*realpow(x,2).*realpow(y,2) +5.*realpow(y,4); 
            dy(:,whichone)= -20.*realpow(x,3).*y +20.*x.*realpow(y,3);
        case 22
            dx(:,whichone)= 30.*realpow(x,4).*y -60.*realpow(x,2).*realpow(y,3) +6.*realpow(y,5); 
            dy(:,whichone)= 6.*realpow(x,5) -60.*realpow(x,3).*realpow(y,2) +30.*x.*realpow(y,4);
        case 23
            dx(:,whichone)= 120.*realpow(x,4).*y +72.*realpow(x,2).*realpow(y,3) -60.*realpow(x,2).*y -72.*realpow(x,2).*realpow(y,3) -24.*realpow(y,5) +20.*realpow(y,3); 
            dy(:,whichone)= 24.*realpow(x,5) +72.*realpow(x,3).*realpow(y,2) -20.*realpow(x,3) -72.*realpow(x,3).*realpow(y,2) -120.*x.*realpow(y,4) +60.*x.*realpow(y,2);
        case 24
            dx(:,whichone)= 150.*realpow(x,4).*y +180.*realpow(x,2).*realpow(y,3) +30.*realpow(y,5) -120.*realpow(x,2).*y -40.*realpow(y,3) +12.*y; 
            dy(:,whichone)= 30.*realpow(x,5) +180.*realpow(x,3).*realpow(y,2) +150.*x.*realpow(y,4) -40.*realpow(x,3) -120.*x.*realpow(y,2) +12.*x;
        case 25
            dx(:,whichone)= 120.*realpow(x,5) +240.*realpow(x,3).*realpow(y,2) +120.*x.*realpow(y,4) -120.*realpow(x,3) -120.*x.*realpow(y,2) +24.*x; 
            dy(:,whichone)= 120.*realpow(x,4).*y +240.*realpow(x,2).*realpow(y,3) +120.*realpow(y,5) -120.*realpow(x,2).*y -120.*realpow(y,3) +24.*y;
        case 26
            dx(:,whichone)= 90.*realpow(x,5) +120.*realpow(x,3).*realpow(y,2) +30.*x.*realpow(y,4) -80.*realpow(x,3) -40.*x.*realpow(y,2) +12.*x -60.*realpow(x,3).*realpow(y,2) -60.*x.*realpow(y,4) +40.*x.*realpow(y,2); 
            dy(:,whichone)= 60.*realpow(x,4).*y +60.*realpow(x,2).*realpow(y,3) -40.*realpow(x,2).*y -30.*realpow(x,4).*y -120.*realpow(x,2).*realpow(y,3) -90.*realpow(y,5) +40.*realpow(x,2).*y +80.*realpow(y,3) -12.*y;
        case 27
            dx(:,whichone)= 36.*realpow(x,5) +24.*realpow(x,3).*realpow(y,2) -20.*realpow(x,3) -144.*realpow(x,3).*realpow(y,2) -72.*x.*realpow(y,4) +60.*x.*realpow(y,2) +12.*x.*realpow(y,4); 
            dy(:,whichone)= 12.*realpow(x,4).*y -72.*realpow(x,4).*y -144.*realpow(x,2).*realpow(y,3) +60.*realpow(x,2).*y +24.*realpow(x,2).*realpow(y,3) +36.*realpow(y,5) -20.*realpow(y,3);
        case 28
            dx(:,whichone)= 6.*realpow(x,5) -60.*realpow(x,3).*realpow(y,2) +30.*x.*realpow(y,4); 
            dy(:,whichone)= -30.*realpow(x,4).*y +60.*realpow(x,2).*realpow(y,3) -6.*realpow(y,5);
    end
end
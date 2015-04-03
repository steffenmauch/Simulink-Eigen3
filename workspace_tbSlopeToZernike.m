% script to compare Eigen3 Simulink Model with Matlab calculation of the
% corresponding Zernike coefficients
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

slopeX = randn(14*14,1);
slopeY = randn(14*14,1);
orderZ = 28;

% slopeX, slopeY and orderZ are used inside the
% tbSlopesToZernike.slx Simulink model
sim('tbSlopesToZernike.slx')

tempX = slopeX;
tempY = slopeY;

% starting fitting zernike polynoms
[testCoordX,testCoordY] = meshgrid(-1:2/(14-1):1);

zernikeMatlab = fitslopesZern(tempX(:),tempY(:),...
    1, testCoordX(:), testCoordY(:), orderZ);

error = zernikeMatlab - zernike.Data(1,:)';

maxError = max( abs(error));
fprintf( 'maxError is %d\n',maxError );
function [zernikes] = fitslopesZern(xslopes, yslopes, pupildiam, xs, ys, maxterm)
% function takes slopes of e.g. an SHWFS and calculates the zernike
% coefficients till order=maxterm
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

[x_basis,y_basis] = arrayofslopefunctionsZern(pupildiam,xs,ys,maxterm);
A = [x_basis(:,:); y_basis(:,:)];
[u,w,v] = svd(A,0);
w = pinv(w);
zernikes = pupildiam/2*v*w*u'*[xslopes; yslopes];

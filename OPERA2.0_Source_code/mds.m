function [positions, error, ev] = mds(d, order)
% function [positions, error] = mds(d, order)
% Compute the multidimensional scaling that produces the best positions
% given the two dimensional (symmetric) distance matrix supplied as an input.
% Return the order-dimension position data.

% Malcolm Slaney and Michele Covell, "Matlab Multidimensional Scaling Tools,"
% Interval Technical Report #2000-025, 2000 (also available at
% http://web.interval.com/papers/2000-025/).

% This routine written by Malcolm Slaney - Interval Research Corporation - 
% May 1998. (c) Copyright Interval Research, May 1998.

% This is experimental software and is being provided to Licensee
% 'AS IS.'  Although the software has been tested on Macintosh, SGI, 
% Linux, and Windows machines, Interval makes no warranties relating
% to the software's performance on these or any other platforms.
%
% Disclaimer
% THIS SOFTWARE IS BEING PROVIDED TO YOU 'AS IS.'  INTERVAL MAKES
% NO EXPRESS, IMPLIED OR STATUTORY WARRANTY OF ANY KIND FOR THE
% SOFTWARE INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY OF
% PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
% IN NO EVENT WILL INTERVAL BE LIABLE TO LICENSEE OR ANY THIRD
% PARTY FOR ANY DAMAGES, INCLUDING LOST PROFITS OR OTHER INCIDENTAL
% OR CONSEQUENTIAL DAMAGES, EVEN IF INTERVAL HAS BEEN ADVISED OF
% THE POSSIBLITY THEREOF.
%
%   This software program is owned by Interval Research
% Corporation, but may be used, reproduced, modified and
% distributed by Licensee.  Licensee agrees that any copies of the
% software program will contain the same proprietary notices and
% warranty disclaimers which appear in this software program.
%

% If no input distance.... generate a set of test points (x) and then compute
% the distances between points.  Use this as input data.

% By Malcolm, straight from the book.

if nargin < 1
	x=[1 2 1 2;
	   1 1 3 4];
	n = size(x,2);
	x(1,:) = x(1,:)-sum(x(1,:))/n;
	x(2,:) = x(2,:)-sum(x(2,:))/n;

	d = zeros(n,n);
	for i=1:n
		for j=1:n
			d(i,j) = sqrt(sum((x(:,i)-x(:,j)).^2));
		end
	end
end

if nargin < 2
	order = 2;
end

n = size(d,1);
A = -1/2*d.^2;
ardot = 1/n * ones(n,1) * sum(A);
asdot = 1/n * sum(A')' * ones(1,n);
adotdot = 1/n/n * sum(sum(A)) * ones(n,n);

B = A - ardot - asdot + adotdot;
[u,s,v] = svd(B);

positions = (v*s.^.5)';

figure;

if size(positions,1) > 1
	plot(positions(1,:),positions(2,:),'bx');
	title('2D MDS Solutions');
	if nargin < 1
		hold on
		plot(x(1,:),x(2,:),'r+')
		title('2D MDS Solutions (red are test input locations)');
		hold off
	end
	drawnow;
end
ev = diag(s);

error = [];
if nargout > 1
	for order=1:min(10,size(positions,1))
		error(order) = 0;
		for i=1:n
			for j=1:n
				this_d = sqrt(sum((positions(1:order,i) - ...
						positions(1:order,j)).^2));
				error(order) = error(order) + ...
					(this_d - d(i,j)).^2;
			end
		end
	end
end
positions = positions(1:order,:);

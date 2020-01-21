function [yc,yc_weighted,w] = nnrcalcy2(neighbors_y,neighbors_distance,K)


% version 2.0 - June 2016
% Kamel Mansouri
% NCCT US EPA
% mansourikamel@gmail.com


% calc response on the basis of K neighbors
yc = mean(neighbors_y);

% old weighted version
% if sum(neighbors_distance) > 0
%     w = ones(1,length(neighbors_y))./(0.1 + neighbors_distance);
%     w = w./sum(w);
% else
%     w = ones(1,length(neighbors_y))./length(neighbors_y);
% end

% new weighted version

w=zeros(1,K);
% 
if any(neighbors_distance<1e-5)
    f=1e-10;
else
      f=0.05;
end

%  f=0.05;

if length(find(isnan(neighbors_distance))) < (K-1)
    w(~isnan(neighbors_distance)) = ones(1,length(find(~isnan(neighbors_distance))))./(f + neighbors_distance(~isnan(neighbors_distance)));
    w = w./sum(w);

% elseif any(neighbors_distance==0) && length(find(isnan(neighbors_distance))) < (K-1)
%     w(~isnan(neighbors_distance)) = ones(1,length(find(~isnan(neighbors_distance))))./(0.05 + neighbors_distance(~isnan(neighbors_distance)));
%     w = w./sum(w);
else
%     w(find(neighbors_distance==0)) = ones(1,length(find(neighbors_distance==0)))./length(find(neighbors_distance==0));
%     w(find(neighbors_distance~=0)) = zeros(1,length(find(neighbors_distance~=0)));
    w = ones(1,K)./K;
    
end

yc_weighted = neighbors_y'*w';

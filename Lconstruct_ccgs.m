% VX,VY record the edge connecting two nodes, X is the row info, Y is the
% column info
% x, y record the coordination of the starting node
% this is a speed up version
function [VX,VY,x,y,edges,params] = Lconstruct_ccgs(bounds,alpha,r)

if nargin < 2,
    alpha = 0.5;
end
if nargin < 3
    %r = rand(1);
    r = 0.2;
end
disp('Constructing the CCG...');
tic

X = [bounds(:).centroid_r; bounds(:).centroid_c];

%distance matrix
D = pdist(X','euclidean');

% probability matrix
P = D.^-alpha;

%define edges
edges = triu(true(length(X)), 1) & squareform(r < P);

% get edge locations
[xx, yy] = find(edges);
VX = [bounds.centroid_r(xx); bounds.centroid_r(yy)]';
VY = [bounds.centroid_c(xx); bounds.centroid_c(yy)]';

% VX = [bounds(xx).centroid_r; bounds(yy).centroid_r]';
% VY = [bounds(xx).centroid_c; bounds(yy).centroid_c]';

% get node locations
idx = unique([xx, yy],'rows', 'first');
x = bounds.centroid_r(idx);
y = bounds.centroid_c(idx);

params.r = r;
params.alpha = alpha;
toc
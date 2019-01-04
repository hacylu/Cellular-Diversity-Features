function feats = haralick_no_img_correct(SGLD)
feats.names = {'contrast_energy','contrast_inverse_moment', ...
    'contrast_ave','contrast_var','contrast_entropy', ...
    'intensity_ave','intensity_variance','intensity_entropy', ...
    'entropy,','energy', ...
    'correlation', ...
    'information_measure1','information_measure2'};

% Initialize output
%[nrows,ncols,nsl] = size(img);
%feats.img3 = zeros(length(feats.names),nrows,ncols);

%%% Calculate Statistics %%%%
[pi,pj,p] = find(SGLD);
if length(p) <= 1
    return;
end
p = p/sum(p); pi = pi-1; pj = pj-1;

%marginal of x
px_all = sum(SGLD,2);
[pxi,junk,px] = find(px_all);
px = px/sum(px); pxi = pxi-1;
% 	            if length(px) <=3
% 	                continue;
% 	            end

%marginal of y
py_all = sum(SGLD,1)';
[pyi,junk,py] = find(py_all);
py = py/sum(py); pyi = pyi-1;


%%% Calculate Contrast Features %%%%
all_contrast = abs(pi-pj);
[sorted_contrast,sind] = sort(all_contrast);
ind = [find(diff(sorted_contrast)); length(all_contrast)];
contrast = sorted_contrast(ind);
pcontrast = cumsum(p(sind));
pcontrast = diff([0; pcontrast(ind)]);

contrast_energy = sum( contrast.^2 .* pcontrast );
contrast_inverse_moment = sum( (1./(1+contrast.^2)) .* pcontrast );
contrast_ave = sum( contrast .* pcontrast );
contrast_var = sum( (contrast - contrast_ave).^2 .* pcontrast );
contrast_entropy = -sum( pcontrast.*log(pcontrast) );

%%% Calculate Intensity Features %%%%
all_intensity = (pi+pj)/2;
[sorted_intensity,sind] = sort(all_intensity);
ind = [find(diff(sorted_intensity)); length(all_intensity)];
intensity = sorted_intensity(ind);
pintensity = cumsum(p(sind));
pintensity = diff([0; pintensity(ind)]);

intensity_ave = sum( intensity .* pintensity );
intensity_variance = sum( (intensity-intensity_ave).^2 .* pintensity );
intensity_entropy = -sum( pintensity.*log(pintensity) );

%%% Calculate Probability Features %%%%
entropy = -sum( p.*log(p) );
energy = sum( p.*p );

%%% Calculate Correlation Features %%%%
mu_x = sum(pxi.*px);
sigma_x = sqrt(sum( (pxi-mu_x).^2 .* px ));
mu_y = sum(pyi.*py);
sigma_y = sqrt(sum( (pyi-mu_y).^2 .* py));

if sigma_x==0 || sigma_y==0
    warning('Zero standard deviation.');
else
    correlation = sum( (pi-mu_x).*(pj-mu_y).* p ) / (sigma_x*sigma_y);
end

%%% Calculate Information Features %%%%
[px_grid,py_grid] = meshgrid(px,py);
[log_px_grid,log_py_grid] = meshgrid(log(px),log(py));
h1 = -sum( p .* log(px_all(pj+1).*py_all(pi+1)) );
h2 = -sum( px_grid(:).*py_grid(:).*(log_px_grid(:)+log_py_grid(:)) );
hx = -sum( px.*log(px) );
hy = -sum( py.*log(py) );

information_measure1 = (entropy-h1)/max(hx,hy);
information_measure2 = sqrt(1-exp(-2*(h2-entropy)));

for k = 1:length(feats.names)
    feats.val(k) = eval(feats.names{k});
end
addpath(genpath('./veta_watershed'));
addpath(genpath('./staining_normalization'));

%% step 1: load image
img = imread('test.png');
curIM=img(1:2000,1:2000,:);%imshow(img)
%% step 2: segment nuclei and save boundaries
close all

curIMsize=size(curIM);
[curIM_norm] = normalizeStaining(curIM);
curIM_normRed=curIM_norm(:,:,1);

p.scales=[6:4:16];
disp('begin nuclei segmentation using watershed');
[nuclei, properties] = nucleiSegmentationV2(curIM_normRed,p);

figure;imshow(curIM);hold on;
for k = 1:length(nuclei)
    plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
end
hold off;

%%% or you can load precompuated nuclei boundaries and properties
% load('nuclei_properties.mat', 'nuclei', 'properties');


ctemp=[properties.Centroid];
bounds.centroid_c=ctemp(1:2:end);
bounds.centroid_r=ctemp(2:2:end);
%% step 3: extract features
% cellular diversity features (1*20*13*6=1560)
para.CGalpha_min=0.44; para.CGalpha_max=0.44;% larger the alpha, sparse the local cell graph
para.alpha_res=0.02;
para.radius=0.2;

CGinfo=[]; % cell graph information
feats=[];
feats_description=[];
set_alpha=[para.CGalpha_min:para.alpha_res:para.CGalpha_max];
for f=1:length(set_alpha)
    curpara.alpha=set_alpha(f); curpara.radius=para.radius;
    [feat,feat_description,CGinfo{f}] = Lextract_CellularDiversity_features(bounds,properties,...
        curpara.alpha,curpara.radius);
    %%% get feature description
    temp= feat_description;
    for i=1:length(temp)
        cur=temp{i};
        str=sprintf('-a=%.2f',curpara.alpha);
        cur(end+1:end+length(str))=str;
        temp{i}=cur;
    end

    feats=cat(2,feats,feat); 
    feats_description=cat(2,feats_description,temp);
end
disp('feature extraction done. please check the variable feats and feats_description for cellular diversity features and their feature names');
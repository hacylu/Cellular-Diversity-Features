%
% Input
%       -bounds: structure contains the feilds 'centroid_c' and 'centroid_r'  for the coordinates of pre-segmented nuclei? check nuclei_properties.mat for example
%       -properties: n-by-1 struct contains the nuclear morphology
%                   information. check nuclei_properties.mat for example
%       -a?parameter for the local cell cluster graph (CCG), larger value =
%           sparser the graph. Normally set to 0.44 for 40x magnification
%           histology images.
%       -r?parameter for the local cell cluster graph (CCG), set to 0.2
%           normally.
%     
% Program written by Cheng LU
% Case Western Reserve University, email:hacylu@gmail.com
% 2018 March 24th

% check runme.m as a demo for how to use this function for cellular
% diversity feature extraction.

% Reference: 
% This function implement a version of the method proposed in manuscripts
% 1-Cheng Lu, James S. Lewis Jr.*, William D. Dupont, W. Dale Plummer Jr., Andrew Janowczyk, Anant Madabhushi*?
% “An Oral Cavity Squamous Cell Carcinoma Quantitative Histomorphometric-Based Image Classifier (OHbIC) of
%  Nuclear Morphology Can Risk Stratify Patients for Disease Specific Survival”, Modern Pathology. 2017 Aug 4. 

function [allCellDiv,CellDiv_description,CCGinfo] = Lextract_CellularDiversity_features(bounds,properties,a,r)
%% buid local cell cluster graph (CCG)
% build graph
alpha = a;

% build graph
[VX,VY,x,y,edges] = Lconstruct_ccgs(bounds,alpha, r);

CCGinfo.VX=VX;
CCGinfo.VY=VY;

% check the CCG configuration here
% nodes=bounds;
% figure(1);
% hold on;
% for x = 1: numel(nodes)
%     plot(nodes(x).centroid_c, nodes(x).centroid_r,'yo', 'MarkerSize',2, 'MarkerFaceColor', 'g');
% end
% plot(VY', VX', 'g-', 'LineWidth', 1);
% hold off;

%% construct the edge matrix
for j = 1:length(bounds.centroid_c)-1
    for k = j+1:length(bounds.centroid_c)
        edges(k,j) = edges(j,k);
    end
end

% find gland networks
[numcomp,group] = graphconncomp(sparse(edges));

% remove single gland networks
temp_network = hist(group,numcomp);
[a,group_ind] = find(temp_network > 1);

%% parameters settings
allCellDiv=[];
feats=[];
% feature_names'; can expand to include more nuclear morphology
feature_names={'Area','MajorAxisLength','MinorAxisLength','Eccentricity','Orientation','EquivDiameter','Solidity','Perimeter','Circularity','EllipticalDeviation',...
    'MeanIntensity','IntensityDeviation','IntensityRange','MeanInsideBoundaryIntensity','InsideBoundaryIntensityDeviation','InsideBoundaryIntensityRange',...
    'MeanOutsideBoundaryIntensity','OutsideBoundaryIntensityDeviation','OutsideBoundaryIntensityRange','BoundarySaliency'};
feature_max=[2000, 80, 30,   1,    90, 60,    1, 180,    1, 0.6, 180,  60, 180, 180, 60, 180, 180,60, 180, 90];
feature_min=[ 200, 20,  8, 0.4,   -90, 20, 0.86,  50, 0.85,   0,  40,  10,  50,  50, 10,  40, 50, 10,  40, 20];

for i=1:length(feature_names)
    cur_f=eval(['[properties.' feature_names{i} ']']);    
    info.feature_max(i)=feature_max(i);
    info.feature_min(i)=feature_min(i);
    info.num_level(i)=6; % quantification level
    info.bin_size(i)=(info.feature_max(i)-info.feature_min(i))/info.num_level(i);% width of bin
    w = info.bin_size(i);
    % let the out of range value to be the boundary value
    cur_f(cur_f>info.feature_max(i))=info.feature_max(i);
    cur_f(cur_f<info.feature_min(i))=info.feature_min(i);
    % discretization
    feature_discret(:,i) = floor((cur_f-info.feature_min(i))/w)+1; %+1 to make them start from 1 instead of 0
    
    % initialize co-occurence
    bin(i,:) = [info.feature_min(i):w:info.feature_max(i)]; % discretizations!!!
    bin(i,:) = [0:info.num_level(i)]+1; % discretizations!!!
    
    idx_c=1;
    c=cell(1);
    feats=cell(1);
    for ii = 1:length(group_ind) % define a neighborhood for each gland network (number of neighborhoods = number of networks)
        p = zeros(size(bin(i,:),2), size(bin(i,:),2));
        neighborhood_angles = feature_discret(group == group_ind(ii),i); % aggregate angles in cluster cell network
        cur_g=unique(neighborhood_angles);
        cur_g=cur_g(~isnan(cur_g));
        if length(cur_g)>2
            C= nchoosek(cur_g,2);
            for jj=1:size(C,1)
                p(C(jj,1),C(jj,2))=sum(neighborhood_angles==C(jj,2))*sum(neighborhood_angles==C(jj,1));
            end
            p=p+p';
            % account for the identical bins
            for jj=1:length(cur_g)
                 p(jj,jj)=sum(neighborhood_angles==jj);
            end
            
            c(idx_c) = {p./sum(p(:))}; % normalize co-occurence matrix % already check : orignal matrix and the normalized matrix only affect the information measurement 1 feature
            feats{idx_c} = haralick_no_img_correct(c{idx_c}); %haralick_no_img_correct(p);
            idx_c=idx_c+1;
        end
    end
    %% statistics: mean, median, standard deviation, range, kurtosis, skewness, across bounds for each haralick feature
    % check if the current nuclear morpholgy has no haralick features
    if ~isempty(feats) && length(feats)>3
        all_feat=[];% aggregating all haralick features for current nuclear morpholgy for all CCG in an image
        for fi=1:length(feats)
            cur_feat=feats{fi}.val;
            all_feat=[all_feat;cur_feat];
        end
        
        curCellDiv=[];
        curCellDiv=[mean(all_feat) median(all_feat) std(all_feat)  range(all_feat) kurtosis(all_feat) skewness(all_feat)];
        allCellDiv=[allCellDiv curCellDiv];
    else
        allCellDiv=[allCellDiv zeros(1,78)];
    end
end

%% feature names organization
count = 1;
modifier = [{'mean'} {'median'} {'std'} {'range'} {'kurtosis'} {'skewness'} ];
%the format is like this: morphofeaturenames_haralickfeatuename_statistcname
haralick_feature_names={'contrast-energy','contrast-inverse-moment','contrast-ave','contrast-var',...
    'contrast-ent','intensity-ave','intensity-var','intensity-ent','entropy','energy',...
    'correlation','info-measure1','info-measure2'};
for m=1:length(feature_names);
    for mi = 1:numel(modifier)
        for j = 1:numel(haralick_feature_names)
            CellDiv_description{count} = ['CellDiv-' feature_names{m} ':' modifier{mi} '(' haralick_feature_names{j} ')'  ];
            count = count + 1;
        end
    end
end
end
classdef SegmenterModel
    properties
        img_idx
        init_img
        label
        liver_mask
        mask
        adj_img
        phi
        g
        curv
        S
        tumor_mask
        label_mask
        region_stats
        eval_metrics
        p
        results_path
    end
    methods
        %% Constructor
        function obj = SegmenterModel(config_path,results_path)
            fprintf("Initializing the segmenter model from %s: ",config_path);
            if nargin == 2
                % Read config file
                p = fileread(config_path);
                p = jsondecode(p);
                p.num_tiles = reshape(p.num_tiles,1,[]);
                obj.p = p;
            end
            obj.results_path = results_path;
            addpath('helpers\');
            fprintf("Done!\n");
        end

        %% Main methods
        
        % Set the initial image and its label
        function obj = setImage(obj,idx)
            fprintf("Importing image %d and its own label: ",idx);
            obj.img_idx = idx;
            obj.init_img = obj.readImageFromDataset(idx);
            obj.label = obj.readLabelFromDataset(idx);
            obj = obj.createLabelMask();
            obj.liver_mask = [];
            obj.mask = [];
            obj.adj_img = [];
            obj.phi = [];
            obj.g = [];
            obj.curv = [];
            obj.S = [];
            obj.tumor_mask = [];
            obj.region_stats = [];
            obj.eval_metrics = [];
            fprintf("Done!\n");
        end

        % Perform the preprocessing
        function obj = preprocessImage(obj)
            fprintf("Extracting liver mask: ");
            obj.liver_mask = obj.init_img > obj.p.binary_thresh;
            fprintf("Done!\nApplying CLAHE: ");
            level = prctile(obj.init_img(:),obj.p.prctile);
            obj.adj_img = obj.init_img;
            obj.adj_img(obj.init_img > level) = level;
            obj.adj_img = adapthisteq(obj.adj_img, ...
                "NumTiles",obj.p.num_tiles, ...
                "ClipLimit",obj.p.clip_limit);
            level = graythresh(obj.adj_img);
            fprintf("Done!\nApplying thresholding: ");
            obj.mask = imbinarize(obj.adj_img,level);
            obj.mask(~obj.liver_mask) = 1;
            obj.mask = imclose(obj.mask,strel('disk',obj.p.close_r));
            obj.mask = ~bwareaopen(~obj.mask,obj.p.noise_area,obj.p.conn);
            obj.mask = bwareaopen(obj.mask,obj.p.noise_area,obj.p.conn);
            fprintf("Done!\n");
        end
        
        function obj = levelSetEvolution(obj)
            fprintf("Level Set Evolution:\n");
            obj.phi = cell([obj.p.n_iters+1,1]);
            obj.S = cell([obj.p.n_iters,1]);
            obj.phi{1} = bwdist(~obj.mask) - bwdist(obj.mask);
            obj.g = obj.computeG();
            [g_x,g_y] = gradCentral(obj.g);
            for k = 1:obj.p.n_iters
                if mod(k,50)==0
                    fprintf("\tIteration %d\n",k);
                end
                [phi_x,phi_y] = gradCentral(obj.phi{k});
                obj.curv = computeCurvature(phi_x,phi_y);
                F_d = obj.computeDataForce(k);
                F_geo = obj.computeGeoForce(k);
                F_edge = obj.p.beta*(g_x.*phi_x + g_y.*phi_y);
                F_ldtp = obj.p.mu_ldtp * laplacian(obj.phi{k});
                obj.S{k} = F_d + F_geo - F_edge + F_ldtp;
                obj.S{k}(~obj.liver_mask) = 0;
                inside_mask = (obj.phi{k} < 0);
                S_pos = max(obj.S{k}, 0);
                S_neg = min(obj.S{k}, 0);
                obj.S{k}(inside_mask) = S_neg(inside_mask) + ...
                    obj.p.eta * S_pos(inside_mask);
                obj.phi{k+1} = obj.phi{k} + obj.p.dt .* obj.S{k};
                if obj.p.reinit_int > 0 && mod(k, obj.p.reinit_int) == 1
                    obj.phi{k+1} = sussmanReinit(obj.phi{k+1}, ...
                        obj.p.reinit_iters,obj.p.reinit_dt);
                end
            end
            obj.mask = obj.phi{k+1} < 0;
            obj.mask = bwareaopen(obj.mask,obj.p.noise_area);
            obj.mask = imfill(obj.mask,'holes');
            fprintf("Done!\n")
        end

        function obj = extractRegionProps(obj)
            fprintf("Computing the statistics of eligible regions: ");
            CC = bwconncomp(obj.mask);
            N = CC.NumObjects;
            obj.tumor_mask = false(size(obj.mask));
            if N == 0
                obj.region_stats = table([],[],[],[],[],[],'VariableNames', ...
                    {'Area','Circularity','AverageIntensity','PixelIdxList'});
                warning('No connected regions found in regions_mask.');
                return;
            end
            props = regionprops(CC,'Area','Perimeter','Centroid','PixelIdxList');
            pixel_idx_list = {props.PixelIdxList}';
            areas = [props.Area]';
            perimeters = [props.Perimeter]';
            circularities = zeros(N,1);
            for i = 1:N
                circularities(i) = 4*pi*areas(i) / (perimeters(i)^2);
            end
            averages = zeros([N,1]);
            for i = 1:N
                averages(i) = mean(obj.adj_img(pixel_idx_list{i}));
            end
            obj.region_stats = table(areas,circularities, ...
                averages,pixel_idx_list,'VariableNames',{'Area', ...
                'Circularity','AverageIntensity','PixelIdxList'});
            obj.region_stats = sortrows(obj.region_stats,'Area','descend');
            fprintf("Done!\n");
        end

        function obj = tumorIdentification(obj)
            fprintf("Choosing the tumor region: ");
            tumor_idx = [];
            [~, order] = sort(obj.region_stats.Area,'descend');
            for k = 1:length(order)
                idx = order(k);
                if obj.region_stats.AverageIntensity(idx) <= obj.p.max_intensity && ...
                        obj.region_stats.Circularity(idx) >= obj.p.min_circ
                    tumor_idx = idx;
                    obj.tumor_mask(obj.region_stats.PixelIdxList{idx}) = true;
                    break;
                end
            end
            se_expand = strel('disk',obj.p.expand_radius);
            obj.tumor_mask = imdilate(obj.tumor_mask,se_expand);
            if isempty(tumor_idx)
                warning("No tumor region detected.\n");
                return
            else
                fprintf("Done!\n");
            end
        end

        function obj = calculateEvalMetrics(obj)
            fprintf("Calculating evaluation metrics: ");
            [rows,cols] = size(obj.init_img);
            tumor_area = nnz(obj.tumor_mask);
            label_area = nnz(obj.label_mask);
            inter_area = nnz(obj.tumor_mask & obj.label_mask);
            union_area = nnz(obj.tumor_mask | obj.label_mask);
            precision = inter_area/max(tumor_area,eps);
            recall = inter_area/max(label_area,eps);
            fpr = (nnz(obj.tumor_mask & ~obj.label_mask)) / ...
                max((rows*cols - label_area),eps);
            rae = (tumor_area - label_area) / max(label_area,eps);
            if union_area == 0
                dice = NaN;
                jaccard = NaN;
            else
                dice = 2*inter_area/(tumor_area+label_area);
                jaccard = inter_area/union_area;
            end
            try
                seg_rp = regionprops(obj.tumor_mask,'Perimeter','Centroid');
                seg_perim = sum([seg_rp.Perimeter]);
                seg_centroid = seg_rp.Centroid;
            catch
                seg_perim = NaN; seg_centroid = [NaN NaN];
            end
            try
                label_rp = regionprops(obj.label_mask,'Perimeter','Centroid');
                label_perim = sum([label_rp.Perimeter]);
                label_centroid = label_rp.Centroid;
            catch
                label_perim = NaN; label_centroid = [NaN NaN];
            end
            if all(~isnan(seg_centroid)) && all(~isnan(label_centroid))
                centr_dist = sqrt((seg_centroid(1)-label_centroid(1))^2 + ...
                    (seg_centroid(2)-label_centroid(2))^2) * ...
                    obj.p.pixel_spacing;
            else
                centr_dist = NaN;
            end
            obj.eval_metrics = struct();
            obj.eval_metrics.AreaSeg = tumor_area;
            obj.eval_metrics.AreaGT = label_area;
            obj.eval_metrics.PerimeterSeg = seg_perim;
            obj.eval_metrics.PerimeterGT = label_perim;
            obj.eval_metrics.Precision = precision;
            obj.eval_metrics.Recall = recall;
            obj.eval_metrics.FalsePositiveRate = fpr;
            obj.eval_metrics.RelAreaError = rae;
            obj.eval_metrics.Dice = dice;
            obj.eval_metrics.Jaccard = jaccard;
            obj.eval_metrics.CentroidEuclDistance = centr_dist;
            fprintf("Done!\n");
            obj.storeResults();
        end

        %% Helpers
        function img = readImageFromDataset(obj,idx)
            img = imread(strcat(obj.p.img_folder,num2str(idx),'.jpg'));
            if size(img,3) == 3
                img = rgb2gray(img);
            end
            img = im2double(img);
        end

        function label = readLabelFromDataset(obj,idx)
            label = fileread(strcat(obj.p.label_folder,num2str(idx),'.json'));
            label = jsondecode(label);
        end

        function obj = createLabelMask(obj)
            label_x = obj.label(:,1);
            label_y = obj.label(:,2);
            [rows,cols] = size(obj.init_img);
            obj.label_mask = poly2mask(label_x,label_y,rows,cols);
            smt_se = strel('disk', obj.p.smt_disk);
            obj.label_mask = imopen(obj.label_mask, smt_se);
            obj.label_mask = imclose(obj.label_mask, smt_se);
        end

        function g = computeG(obj)
            smt_img = imgaussfilt(obj.adj_img,obj.p.sigma);
            [Ig_x,Ig_y] = gradCentral(smt_img);
            [Ig_xy,~]  = gradCentral(Ig_y);
            [~,Ig_yx] = gradCentral(Ig_x);
            Ig_x_p = Ig_x + (Ig_xy - Ig_yx) * cos(pi/4);
            Ig_y_p = Ig_y + (Ig_xy - Ig_yx) * sin(pi/4);
            mag = sqrt(Ig_x_p.^2 + Ig_y_p.^2);
            g = 1 ./ (1 + mag);
        end
        
        function F_d = computeDataForce(obj,k)
            H = 0.5 * (1 + (2/pi) * atan(obj.phi{k} ./ obj.p.eps_h));
            mask_inside = H.*obj.mask;
            mask_outside = (1-H).*obj.mask;
            inside_area  = sum(mask_inside(:));
            outside_area = sum(mask_outside(:));
            inside_mean = sum(obj.adj_img(:) .* mask_inside(:)) / ...
                (inside_area + eps);          
            outside_mean = sum(obj.adj_img(:) .* mask_outside(:)) / ...
                (outside_area + eps);
            F_d = obj.p.lambda_1*(obj.adj_img-inside_mean).^2 - ...
                obj.p.lambda_2*(obj.adj_img-outside_mean).^2;
        end

        function F_geo = computeGeoForce(obj,k)
            near_iface = abs(obj.phi{k}) <= obj.p.band_width;
            F_geo = obj.p.alpha .* (obj.g .* obj.curv);
            F_geo(~near_iface) = 0;
        end

        %% Data Visualization
        function [] = showPreprocessing(obj)
            if ~isempty(obj.adj_img)
                fig_title = "Preprocessing: binary thresholding, CLAHE, ";
                fig_title = strcat(fig_title,"Otsu thresholding");
                figure('Name',fig_title,'NumberTitle','off')
                subplot(1,2,1);
                imshow(obj.init_img);
                title("Initial image");
                subplot(1,2,2);
                imshow(obj.adj_img);
                title("Adjusted image");
            end
        end

        function [] = showEvolution(obj,interval,time)
            figure_title = strcat("Image ",string(obj.img_idx), ...
                ' - Level Set Evolution');
            f1 = figure('Name',figure_title,'NumberTitle','off');
            for k = 1:obj.p.n_iters
                if k==1 || k==obj.p.n_iters || mod(k,interval)==0
                    clf(f1);
                    imshow(obj.adj_img);
                    hold on;
                    contour(obj.phi{k},[0 0],'r','LineWidth',1.2);
                    title(sprintf('Level Set: step %d',k));
                    hold off;
                    drawnow;
                    pause(time);
                end
            end
        end

        function [] = showSegmentation(obj)
            figure_title = strcat('Image ',string(obj.img_idx), ...
                ': Segmentation vs Ground Truth');
            figure('Name',figure_title,'NumberTitle','off');
            imshow(obj.init_img,[]);
            hold on;
            contour(obj.label_mask,'g-','LineWidth',1.5);
            contour(obj.tumor_mask,'r--','LineWidth',1.5);
        end

        function [] = storeResults(obj)
            data = readstruct(obj.results_path);
            img_name = sprintf('Img_%d',obj.img_idx);
            data.(img_name) = obj.eval_metrics;
            content = jsonencode(data,"PrettyPrint",true);
            fid = fopen(obj.results_path,'w');
            fprintf(fid,"%s\n",content);
            fclose(fid);
        end
    end
end
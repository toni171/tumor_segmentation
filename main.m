close all;
clear;
clc;
i = 1;
while i<=1
m = SegmenterModel('config.json','results.json');
m = m.setImage(i);
m = m.preprocessImage();
m = m.levelSetEvolution();
m = m.extractRegionProps();
m = m.tumorIdentification();
m = m.calculateEvalMetrics();
m.showEvolution(10,0.2);
%m.showSegmentation();
fprintf("Accuracy %f\n",m.eval_metrics.Precision);
i = i+1;
end
results_path = "results.json";
results = readstruct(results_path);
img_names = fieldnames(results);
metric_names = fieldnames(results.(img_names{1}));
metrics = struct();
average = struct();
for j=1:numel(img_names)
    for i=1:numel(metric_names)
        if j == 1
            metrics.(metric_names{i}) = [];
        end
        metrics.(metric_names{i}) = ...
            [metrics.(metric_names{i}) (results.(img_names{j}).(metric_names{i}))];
        if j == numel(img_names)
            average.(metric_names{i}) = mean(metrics.(metric_names{i}));
        end
    end
end
results.overall_results = average;
content = jsonencode(results,"PrettyPrint",true);
fid = fopen(results_path,'w');
fprintf(fid,"%s\n",content);
fclose(fid);
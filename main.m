close all;
clear;
clc;
i = 1;
while i<=11
m = SegmenterModel('config.json','results.json');
m = m.setImage(i);
m = m.preprocessImage();
m = m.levelSetEvolution();
m = m.extractRegionProps();
m = m.tumorIdentification();
m = m.calculateEvalMetrics();
%m.showEvolution(10,0.2);
%m.showSegmentation();
fprintf("Accuracy %f\n",m.eval_metrics.Precision);
i = i+1;
end

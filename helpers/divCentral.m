function d = divCentral(p_x,p_y)
[dpxdx,~] = gradCentral(p_x);
[~,dpydy] = gradCentral(p_y);
d = dpxdx + dpydy;
end
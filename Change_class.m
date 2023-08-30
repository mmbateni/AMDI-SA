y_edges2_p=[0 0.05 0.2 0.3 1];
%[0.000	0.023	0.067	0.159	0.841	0.933	0.977	1.000]%for 7 classes
% for 4 classes[0.000 0.05 0.2 0.3 1]which woud be  -Inf   -1.65   -0.84   -0.52       Inf] in normal space
y_edges2=norminv(y_edges2_p,0,1);%for 7 classes
index_cat_mc=discretize(best_cop_indx,y_edges2,'IncludedEdge','right');
clear index_cat
index_cat=index_cat_mc;
muu=length(y_edges2)/2;
index_cat=((-1)*(index_cat-muu))+muu;%Now, 1 means wettest and 4(end) means driest
save('F:\Program Files\MATLAB\MATLAB Production Server\R2015a\bin\work_thesis\Res_MAT\coupla_derived.mat','index_cat','-append');
M_head=readtable('/fs/project/PAS0272/Tara/DNase_SOM/Brain/training_csv/0.5_0_shifted/1.csv');
sz = size(M_head);
M_mat = table2array(M_head(:,4:sz(2)));
Xnorm = std(M_mat, 0, 2);
minimum = 0.0001;    
signal = M_mat(find(Xnorm > minimum),:);
intervalData.chr = M_head(find(Xnorm > minimum),1);
intervalData.start = M_head(find(Xnorm > minimum),2);
intervalData.stop = M_head(find(Xnorm > minimum),3);
sz_rem = size(find(Xnorm > minimum));
C    = cell(1, sz_rem(1));
C(:) = {'+'};
intervalData.strand = transpose(C);
intervalData.Properties.VarNames = int2str([1:sz_rem(1)])
save 'cagt_1.mat' signal intervalData;
addpath("/fs/project/PAS0272/Tara/DNase_SOM/Brain/cagt/matlab/src")
cagt('cagt_1.mat', 'od', '../data/test', 'op', 'nucleo_around_ctcf_', 'tt', 'CTCF', 'st', 'NUCLEOSOME', ...
  'merge', true, 'bed', true, 'txt', true, 'maxiter', 100, 'replicates', 2, 'start', 'plus', ...
  'overwrite', true, 'flip', true, 'mergeDist', 0.8, 'mergeK', 1, 'distance', 'xcorr', 'maxlag', 20, ...
  'emptyaction', 'drop')
res = load("../data/test/nucleo_around_ctcf_results.mat")
csvwrite("cagt_shapes_1.csv",res.kmeansResults.centroids)
M_head=readtable('/fs/project/PAS0272/Tara/DNase_SOM/A549/training_shifted/chrom22');
sz = size(M_head);
M_mat = table2array(M_head(:,4:sz(2)));

Xnorm = std(M_mat, 0, 2);    
    
signal = M_mat(find(Xnorm),:);
intervalData.chr = M_head(find(Xnorm),1);
intervalData.start = M_head(find(Xnorm),2);
intervalData.stop = M_head(find(Xnorm),2);
sz_rem = size(find(Xnorm));
C    = cell(1, sz_rem(1));
C(:) = {'+'};
intervalData.strand = transpose(C);
intervalData.Properties.VarNames = int2str([1:sz_rem(1)])
save 'cgct_22.mat' signal intervalData;

cagt('cgct_22.mat', 'od', '../data/test', 'op', 'nucleo_around_ctcf_', 'tt', 'CTCF', 'st', 'NUCLEOSOME', ...
  'merge', true, 'bed', true, 'txt', true, 'maxiter', 100, 'replicates', 2, 'start', 'plus', ...
  'overwrite', true, 'flip', true, 'mergeDist', 0.8, 'mergeK', 1, 'distance', 'correlation')
res = load("../data/test/nucleo_around_ctcf_results.mat")
csvwrite("ctcf_centroids_22.csv",res.kmeansResults.centroids)
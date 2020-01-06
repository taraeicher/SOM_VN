function X=run_cagt(in_file, mat_loc, chrom, csv_file, cagt_path, merge_dist, iterations, k_val)
    
    % Read the table of input training regions. Retain only the signals, not the labels.
    M_head=readtable(in_file);
    sz = size(M_head);
    M_mat = table2array(M_head(:,4:sz(2)));
    
    % Filter out all low-variance regions and build the format
    % needed for CAGT. Note: This is the criterion CAGT tests for
    % before the learning procedure.
    signal = M_mat;
    chr = M_head(:,1);
    start = M_head(:,2);
    stop = M_head(:,3);
    sz_rem = size(M_mat);
    C    = cell(sz_rem(1), 1);
    C(:) = {'+'};
    intervalData = dataset(chr, start, stop, C);
    intervalData.Properties.VarNames = ["chr", "start", "stop", "strand"];
    
    % Read the CAGT file, run CAGT, and save the shapes found in a CSV file.
    mat_loc_file = strcat(mat_loc, "/", chrom, ".mat")
    save(mat_loc_file,'intervalData', 'signal')
    addpath(cagt_path)
    cagt(char(mat_loc_file), 'od', strcat(mat_loc, '/data/test'), 'op', '/nucleo_around_ctcf_', 'tt', 'CTCF', 'st', 'NUCLEOSOME', ...
      'merge', true, 'bed', true, 'txt', true, 'maxiter', str2num(iterations), 'replicates', 1, 'k', str2num(k_val), 'start', 'plus', ...
      'overwrite', true, 'flip', true, 'mergeDist', str2num(merge_dist), 'mergeK', 1, 'lowvarcut', 0.01, 'distance', 'correlation', 'maxlag', 20, ...
      'emptyaction', 'drop', 'avgFun', 'median')
    res = load(strcat(mat_loc, "/data/test/nucleo_around_ctcf_results.mat"))
    csvwrite(csv_file,res.kmeansResults.centroids)
    quit()
end
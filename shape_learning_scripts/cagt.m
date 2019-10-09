function cagt(infile, varargin)
%CAGT(INFILE, VARARGIN) Runs CAGT clustering on a set of signals.
%
%   CAGT(INFILE) runs CAGT on data read from INFILE and writes the results
%   in the current directory. INFILE should be a mat file containing the
%   following variables:
%
%      intervalData: dataset with fields
%        chr[nominal] : chromosome names
%        start[double]: start positions (1-based)
%        stop[double]: stop positions (1-based)
%        strand[nominal]: +/-
%      signal: N-by-P matrix where N is the number of intervals. Clustering will be applied on the rows of signal. 
%
%   CAGT writes the following outputs:
%     cagt_params.mat: Calling parameters
%     cagt_results.mat: Results of clusterSignal.
%     cagt_clusters.bed: BED file with the assignments of elements to
%     clusters (See wiki for details).
%     cagt_clusterData.txt: Cluster signal (See wiki for details).
%
%   CAGT(INFILE, ...) can take additional optional parameter/value pairs:
%
%   Output parameters:
%
%     od: output directory. Default: working directory. 
%     op: prefix of output files. Default: 'cagt_'. Output files will be
%     <op>params.mat, <op>results.mat, [<op>clusterData.txt,
%     <op>clusters.bed].
%     tt: target type. Used for labeling the results in the clusterData
%       output file. Default: 'TARGET'.
%     st: signal type. Used for labeling the results in the clusterData output
%       file. Default: 'SIGNAL'.
%     bed: Write BED output. Default: true.
%     txt: Write clusterData.txt output. Default: true.
%     overwrite: If set to false, then CAGT won't overwrite any output
%       files (result, BED, or clusterData), if files with the same name
%       already exist in the output directory. Default false.
%     display: Determines how much output will be printed. Choices are:
%          0 - No output.
%          1 - Output only after the last iteration (default).
%          2 - Output at each iteration.
%
%   Signal filtering parameters:
%
%     maxNan: rows with more than that number of NaN's will be removed.
%       Set to 0 if you want to keep all rows. Default: ceil(.5 * size(signal, 2)).
%     lowSignalCut/lowSignalPrc: Remove rows whose lowSignalPrc-th
%       percentile is less than lowSignalCut. Defaults are 0 for
%       lowSingalCut and 90 for lowSignalPrc.
%     lowVarCut: rows with variance less than lowVarCut will be removed.
%       Useful when the distance metric used is correlation-based. Default: 0.
%     nanTreat: how to treat NaN's in the signal. Possible values:
%         'zero': replace with zeros.
%         'interpolate': use linear interpolation to interpolate missing
%         values (default).
%
%   Distance metric parameters (apply to both clustering stages):
%
%     distance: Distance function, that K-means should minimize with
%     respect to.  Choices are:
%          'sqeuclidean'  - Squared Euclidean distance (default)
%          'correlation'  - One minus the sample correlation between points
%                           (treated as sequences of values)
%          'xcorr'        - One minus the maximum sample cross-correlation
%                           over all possible lags (see the option
%                           'maxlag')
%     avgFun: Method for computing cluster centroids. Choices are 'mean',
%       'median'. Default: mean.
%     maxlag: Maximum lag when xcorr is used as the distance measure.
%       Otherwise maxlag will be set to 0 and the value passed will be ignored.
%
%   K-means/medians parameters:
%
%     k: number of clusters for k-means/medians. Default: 40.
%     start: Method used to choose initial cluster centroid positions, 
%     sometimes known as "seeds".  Choices are:
%         'plus'    -  k-means++ initialization (default). 
%         'sample'  -  Select k observations from X at random.
%          matrix   -  A k-by-P matrix of starting locations.  In this case,
%                      the k parameter can be omitted, and K will be
%                      inferred from the first dimension of the matrix.
%                      You can also supply a 3D array, implying a value for
%                      replicates from the array's third dimension.
%     replicates: Number of times to repeat the clustering, each with a
%       new set of initial centroids.  Default is 1.
%     emptyaction: Action to take if a cluster loses all of its member
%     observations.  Choices are:
%          'error'     - Treat an empty cluster as an error (default)
%          'drop'      - Remove any clusters that become empty, and set
%                        the corresponding values in C and D to NaN.
%          'singleton' - Create a new cluster consisting of the one
%                        observation furthest from its centroid.
%     maxiter: Maximum number of iterations. Default: 100.
%     online: Flag indicating whether an "on-line update phase should be
%       performed in addition to a "batch update" phase.  The on-line phase
%       can be time consuming for large data sets, but guarantees a solution
%       that is a local minimum of the distance criterion, i.e., a partition
%       of the data where moving any single point to a different cluster
%       increases the total sum of distances.  NOT SUPPORTED YET. Default: false.
%
%   Hierarchical clustering parameters:
%
%     merge: Set to false if you do not want to apply hierarchical
%       agglomerative clustering. Default:true.
%     mergeK: final number of clusters after hierarchical agglomerative clustering. 
%       Default is 1 so merging will continue until there's only 1 cluster.
%       k cannot be greater than m, the number of input clusters.
%     mergeDist: maximum distance for merging. Merging will stop if the
%       distance between the two closest centroids is greater than maxDist,
%       or if the number of clusters reaches k. Default Inf.
%     flip: Whether flipping will be allowed or not. If flipping is allowed 
%       both an observation and its reversed (i.e. taking the P features
%       from end to start) will be considered when looking for the closest
%       cluster. Default: true (flipping allowed).
%
%   Author: Sofia Kyriazopoulou (sofiakp@stanford.edu)

% Parse arguments
parser = inputParser;
parser.addRequired('infile', @(x) validateattributes(x, {'char'}, {'nonempty'}, 'cagt', 'infile'));
parser.addParamValue('od', '.', @(x) validateattributes(x, {'char'}, {}, 'cagt', 'od'));
parser.addParamValue('op', 'cagt_', @(x) validateattributes(x, {'char'}, {'nonempty'}, 'cagt', 'op'));
parser.addParamValue('st', 'SIGNAL', @(x) validateattributes(x, {'char'}, {'nonempty'}, 'cagt', 'st'));
parser.addParamValue('tt', 'TARGET', @(x) validateattributes(x, {'char'}, {'nonempty'}, 'cagt', 'tt'));
parser.addParamValue('bed', true, @(x) validateattributes(x, {'logical'}, {}, 'cagt', 'bed'));
parser.addParamValue('txt', true, @(x) validateattributes(x, {'logical'}, {}, 'cagt', 'txt'));
parser.addParamValue('overwrite', false, @(x) validateattributes(x, {'logical'}, {}, 'cagt', 'overwrite'));
% pre-processing parameters
parser.addParamValue('maxNan', NaN);
parser.addParamValue('lowSignalCut', 0);
parser.addParamValue('lowSignalPrc', 90);
parser.addParamValue('lowVarCut', 0);
parser.addParamValue('nanTreat', 'interpolate');
% distance function parameters
parser.addParamValue('distance', 'sqeuclidean');
parser.addParamValue('avgFun', 'mean');
parser.addParamValue('maxlag', 0);
% k-means parameters
parser.addParamValue('k', 40);
parser.addParamValue('start', 'plus');
parser.addParamValue('replicates', 10);
parser.addParamValue('emptyaction', 'error');
parser.addParamValue('maxiter', 100);
parser.addParamValue('online', false);
parser.addParamValue('display', 1);
% merging parameters
parser.addParamValue('merge', true);
parser.addParamValue('flip', true);
parser.addParamValue('mergeK', 1);
parser.addParamValue('mergeDist', Inf);

parser.parse(infile, varargin{:});

%assert(exist(infile, 'file') > 0, 'Input signal file does not exist.');
load(infile, 'signal', 'intervalData');
display = parser.Results.display;
od = parser.Results.od;
op = parser.Results.op;
overwrite = parser.Results.overwrite;

filterParams = struct();
filterParams.maxNan = parser.Results.maxNan;
filterParams.lowSignalCut = parser.Results.lowSignalCut;
filterParams.lowSignalPrc = parser.Results.lowSignalPrc;
filterParams.lowVarCut = parser.Results.lowVarCut;
filterParams.nanTreat = parser.Results.nanTreat;
filterParams = validateFilterParams(signal, filterParams);

distParams = struct();
distParams.avgFun = parser.Results.avgFun;
distParams.distance = parser.Results.distance;
distParams.maxlag = parser.Results.maxlag;
distParams = validateDistParams(distParams);

kmeansParams = struct();
kmeansParams.k = parser.Results.k;
kmeansParams.start = parser.Results.start;
kmeansParams.replicates = parser.Results.replicates;
kmeansParams.emptyaction = parser.Results.emptyaction;
kmeansParams.maxiter = parser.Results.maxiter;
kmeansParams.online = parser.Results.online;
kmeansParams.display = display;
kmeansParams = validateKmeansParams(signal, kmeansParams, distParams);

hcParams = struct();
hcParams.merge = parser.Results.merge;
hcParams.flip = parser.Results.flip;
hcParams.k = parser.Results.mergeK;
hcParams.maxDist = parser.Results.mergeDist;
validateHcParams(hcParams);

if display
  printParams(infile, parser.Results, filterParams, distParams, kmeansParams, hcParams);
end

if ~isdir(od)
  mkdir(od);
end

if ~exist(fullfile(od, [op, 'results.mat']), 'file') || overwrite
  params = filterParams();
  params.kmeansParams = kmeansParams;
  if hcParams.merge
    params.hcParams = hcParams;
  end
  params.distParams = distParams;
  
  results = clusterSignal(signal, params);
  save(fullfile(od, [op, 'params.mat']), 'params');
  save(fullfile(od, [op, 'results.mat']), '-struct', 'results');
  %saveResults(od, params, infile, signal, results);
else
  paramStruct = load(fullfile(od, [op, 'params.mat']));
  params = paramStruct.params;
  results = load(fullfile(od, [op, 'results.mat']));
end

if parser.Results.bed
  bedfile = fullfile(od, [op, 'clusters.bed']);
  if ~exist(bedfile, 'file') || overwrite
    if any(strcmp(intervalData.Properties.VarNames, 'name'))
      names = intervalData.name;
    else
      names = strcat(op, arrayfun(@(x) {num2str(x)}, 1:size(signal, 1)));
    end
    writeTextResults(bedfile, results, intervalData, names, 0);
  end
end

if parser.Results.txt
  outparams = struct();
  outparams.signalType = parser.Results.st;
  outparams.targetType = parser.Results.tt;
  txtfile = fullfile(od, [op, 'clusterData.txt']);
  if ~exist(txtfile, 'file') || overwrite
    outparams.merged = true;
    makeSignalTable(txtfile, results, signal, outparams);
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
function newParams = validateFilterParams(signal, filterParams)
newParams = filterParams;
if isnan(filterParams.maxNan) % if it was left unset
    newParams.maxNan = ceil(.5 * size(signal, 2));
else
  validateattributes(filterParams.maxNan, {'numeric'}, {'scalar', 'integer', 'nonnegative', '<=', size(signal, 2)}, 'cagt', 'maxNan');
end
validateattributes(filterParams.lowSignalCut, {'numeric'}, {'scalar'}, 'cagt', 'lowSignalCut');
validateattributes(filterParams.lowSignalPrc, {'numeric'}, {'scalar', 'nonnegative', '<=', 100}, 'cagt', 'lowSignalPrc');
validateattributes(filterParams.lowVarCut, {'numeric'}, {'scalar', 'nonnegative'}, 'cagt', 'lowVarCut');
validatestring(filterParams.nanTreat, {'zero', 'interpolate'}, 'cagt', 'nanTreat');
end

function newParams = validateDistParams(distParams)
newParams = distParams;
validatestring(distParams.distance, {'sqeuclidean', 'correlation', 'xcorr'}, 'validateDistParams', 'distance');
validatestring(distParams.avgFun, {'mean', 'median'}, 'validateDistParams', 'avgFun');
validateattributes(distParams.maxlag, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, 'validateDistParams', 'maxlag');  
if distParams.maxlag > 0 && ~strcmp(distParams.distance, 'xcorr')
    warning('validateDistParams:UnusedArgument', '''maxlag'' only applies to xcorr. Input value will be replaced by 0.');
    newParams.maxlag = 0;
end
end

function newParams = validateKmeansParams(signal, params, distParams)

newParams = params;
validatestring(params.emptyaction, {'error', 'drop', 'singleton'}, 'validateKmeansParams', 'emptyaction');
validateattributes(params.maxiter, {'numeric'}, {'scalar', 'integer', '>=', 1}, 'validateKmeansParams', 'maxiter');
validateattributes(params.online, {'logical'}, {}, 'validateKmeansParams', 'online');
validateattributes(params.display, {'numeric'}, {'integer', 'scalar', 'nonnegative', '<=', 2}, 'validateKmeansParams', 'display');

if params.online
  warning('cagt:InvalidArgument', 'Online updates are not supported yet. Setting online to false.');
  newParams.online = false;
end

% n points in p dimensional space
[n, p] = size(signal);

validateattributes(params.k, {'numeric'}, {'scalar', 'integer', '>=', 1}, 'validateKmeansParams', 'k');
assert(params.k <= n, 'The number of clusters must be at most equal to the number of input elements');

start = params.start;
if ischar(start)
    validatestring(start, {'sample', 'plus'}, 'validateKmeansParams', 'start');
    validateattributes(params.replicates, {'numeric'}, {'integer', 'scalar', '>=', 1}, 'validateKmeansParams', 'replicates');
elseif isnumeric(start)
    validateattributes(params.start, {'numeric'}, {'nrows', params.k, 'ncols', p}, 'validateKmeansParams', 'start');
    if size(start, 3) ~= params.replicates
      warning('validateKmeansParams:IncompatibleArguments', 'Setting the number of replicates to the size of the third dimension of ''start''.');
    end
    newParams.replicates = size(start, 3);
else
    error('validateKmeansParams:InvalidArgument', ...
          'The ''start'' parameter value must be a string or a numeric matrix or array.');
end
end

function validateHcParams(hcParams)
validateattributes(hcParams.merge, {'logical'}, {}, 'validateHcParams', 'merge');
if hcParams.merge
  validateattributes(hcParams.flip, {'logical'}, {}, 'validateHcParams', 'flip');
  validateattributes(hcParams.k, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, 'validateHcParams', 'k');
  validateattributes(hcParams.maxDist, {'numeric'}, {'scalar'}, 'validateHcParams', 'maxDist');
end
end

function printParams(infile, p, filterParams, distParams, kmeansParams, hcParams)
  fprintf('Input file:\t%s\n', infile);
  fprintf('Output files will be:\t%s.*\n', fullfile(p.od, p.op));
  if ~p.bed
    fprintf('BED file will NOT be written.\n');
  end
  if ~p.txt
    fprintf('clusterData.txt file will NOT be written.\n');
  else
    fprintf('Target type:\t%s\n', p.tt);
    fprintf('Signal type:\t%s\n', p.st);
  end
  fprintf('\nFiltering parameters:\n');
  disp(filterParams);
  printDistParams(distParams);
  printKmeansParams(kmeansParams);
  if hcParams.merge
    printHcParams(hcParams);
  else
    fprintf('Hierarchical clustering will be skipped\n');
  end
end

function printDistParams(distParams)
  fprintf('Distance metric parameters:\n');
  disp(distParams);
end
function printKmeansParams(kmeansParams)
  fprintf('k-means/medians parameters:\n');
  disp(kmeansParams); 
end
function printHcParams(hcParams)
  fprintf('Hierarchical clustering parameters:\n');
  disp(hcParams); 
end
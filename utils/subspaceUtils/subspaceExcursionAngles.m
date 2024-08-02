function angles = subspaceExcursionAngles(X, nAngles, nPCADims)
% angles = subspaceExcursionAngles(X, nAngles [, nPCADims])
% 
% This function is designed to characterize dimensionality in a new way.
% For a collection of subspaces in an embedding space, it will give you a
% series of angles telling you how far your subspaces swing through
% additional dimensions.
% 
% More specifically, it asks: what is the largest subspace angle from one
% of my subspaces to the farthest subspace? That is, if for every pair of
% subspaces I find all of the principal angles, what's the biggest one?
% Once we have this, we can ask: what is the largest principal angle from
% the subspace defined by augmenting my initial subspace with this farthest
% vector? This question is iterated to produce a series of angles.
% 
% To improve noise robustness, angles within a tolerance of one another are
% considered ties, and the next angle is considered for all options that
% tie. If this breaks the tie, we take the largest one; if not, we recurse
% further. The tolerance is chosen somewhat arbitrarily to be 0.5 degrees.
% If the max angle at any step is smaller than this tolerance, the
% tolerance is chosen to be maxAngle / 10.
% 
% If multiple possible vector sequences end up tying, the largest sum of
% angles is returned.
% 
% INPUTS:
% 
% X        -- Subspaces to test. nDims x nVectorsPerSubspace x nSubspaces
% nAngles  -- The number of angles you want returned
% nPCADims -- Optional. If supplied with a natural number, will PCA X to
%             this dimensionality first. Default is no PCA.
% 
% OUTPUTS:
% 
% angles -- the sequence of angles in radians.
% 
% 
% Implemented with a fast semi-greedy algorithm (where ties garner
% recursion) using dynamic programming. Requires subspacea.m by Andrew
% Knyazev and Rico Argentati: 
% https://www.mathworks.com/matlabcentral/fileexchange/55-subspacea-m
% 
% Matt Kaufman 2020


maxTol = 0; % 0.5 degree tolerance for which is max angle

%% Error checking

% Trivial case
if nAngles < 1 || isnan(nAngles)
  angles = [];
  return;
end

if ndims(X) ~= 3
  error('subspaceExcursionAngle:tooFewModes', 'X should be nDims x nVectorsPerSubspace x nSubspaces');
end

if size(X, 3) == 1
  error('subspaceExcursionAngle:tooFewSubspaces', 'Supply multiple subspaces to test (size(X, 3) > 1)');
end

if any(isnan(X(:)))
  error('subspaceExcursionAngle:NaNs', 'No elements of X may be NaN');
end

if nAngles > size(X, 2) * size(X, 3) - 1
  warning('subspaceExcursionAngle:tooManyAngles', 'Requesting more angles than makes sense; reducing nAngles');
  nAngles = size(X, 2) * size(X, 3) - 1;
end


%% PCA preprocessing if requested

[nDims, nVecs, nSubs] = size(X);

if exist('nPCADims', 'var') && ~isempty(nPCADims) && ~isnan(nPCADims) && nPCADims > 0
  XRe = reshape(X, [nDims, nVecs * nSubs]);
  [~, X] = pca(XRe', 'NumComponents', nPCADims);
  X = X';
  X = reshape(X, [nPCADims, nVecs, nSubs]);
  nDims = nPCADims;
end


%% Initialize: compute subspace angles between all pairs

firstAngs = zeros(nSubs);
for s1 = 1 : nSubs - 1
  for s2 = s1 + 1 : nSubs
    theseAngs = subspacea(X(:, :, s1), X(:, :, s2));
    firstAngs(s1, s2) = theseAngs(end);
    firstAngs(s2, s1) = theseAngs(end);
  end
end


%% Find candidate "seed" subspaces: those that make angles as large as possible to another subspace
% We do this as a special case because (1) doing it like the general case
% below would be wasteful of memory; (2) it would take twice as long (since
% it wouldn't know about the symmetry that we have only for the first
% element of the sequence); (3) it lets us do some extra error checking /
% warning

% Data structures for keeping track of things:
% 
% candSeqs: nCandidateSeqs x seqDepth, this are the candidate sequences
%           we're still considering
% candAngs: nCandidateSeqs x seqDepth-1, these are the angles for each
%           sequence
% candSeqVecs: nDims x seqDepth x nCandidateSeqs, these are the vectors
%           used in each sequence

maxAng = max(firstAngs(:));

% Handle perfectly aligned subspaces
% There are SVD's inside subspacea, so we'll call anything within 1000x of
% numerical precision "perfectly aligned"
if maxAng < 1e3 * eps
  angles = zeros(1, nAngles);
  warning('subspaceExcursionAngles:alignedSubspaces', 'No subspaces have a non-trivial angle');
  return;
end

% Handle having too small a tolerance relative to largest subspace angle
if maxTol > maxAng / 10
  warning('subspaceExcursionAngles:smallAngles', 'Tolerance for angles is small relative to subspace angles; reducing tolerance to 10%% of max angle');
  maxTol = maxAng / 10;
end

% Identify candidate seed subspace sequences. These are length 2, because
% we know both the seed and who it makes the biggest angle with
candSeeds = (firstAngs >= maxAng - maxTol);
nSeqs = sum(candSeeds(:));
[candSeqs(:, 1), candSeqs(:, 2)] = find(candSeeds);

% Angles for each candidate sequence. Only length 1 here, because we're
% only looking at the first two subspaces, will grow as sequence does
candAngs = firstAngs(candSeeds);

% Vectors for each candidate sequence. Length 2 because 2 subspaces, will
% grow as sequence does.
candSeqVecs = NaN(nDims, nVecs + 1, nSeqs);
for seq = 1:nSeqs
  [~, U, V] = subspacea(X(:, :, candSeqs(seq, 1)), X(:, :, candSeqs(seq, 2)));
  candSeqVecs(:, 1:nVecs, seq) = U;
  candSeqVecs(:, end, seq) = V(:, end);
end


%% Candidate seeds and second subspace are set, now main loop to fill out sequence

for a = 2:nAngles
  
  %% Evaluate one level of sequences deeper
  
  % Pre-allocate. Use cell arrays because this lets us do much less
  % resizing
  nSeqs = size(candSeqs, 1);
  nextSeqs = cell(1, nSeqs);
  nextAngs = cell(1, nSeqs);
  prevSeqID = cell(1, nSeqs);
  
  % Loop through candidate sequences
  for seq = 1:nSeqs
    % Try this sequence against every possible partner
    nextAngs{seq} = NaN(1, nSubs);
    for nextInSeq = 1:nSubs
      theseAngs = subspacea(candSeqVecs(:, :, seq), X(:, :, nextInSeq));
      nextAngs{seq}(nextInSeq) = theseAngs(end);
    end
    
    % Find largest angle
    maxAng = max(nextAngs{seq});
    
    % Handle having too small a tolerance relative to largest subspace angle
    angTol = min([maxTol, maxAng / 10]);
    
    % Generate the candidate sequences, keep their angles, and keep which
    % sequence they came from (so we can later recover the vectors in the
    % sequence)
    nextCands = find(nextAngs{seq} >= maxAng - angTol);
    if ~isempty(nextCands)
      nextSeqs{seq} = [repmat(candSeqs(seq, :), [length(nextCands), 1]), nextCands'];
      nextAngs{seq} = [repmat(candAngs(seq, :), [length(nextCands), 1]), nextAngs{seq}(nextCands)'];
      prevSeqID{seq} = repmat(seq, [length(nextCands), 1]);
    end
  end
  
  % Concatenate the results for every sequence into single lists
  candSeqs = vertcat(nextSeqs{:});
  candAngs = vertcat(nextAngs{:});
  prevSeqID = vertcat(prevSeqID{:});
  
  
  %% Identify the sequences to pursue further
  
  % Find the biggest last-angles over all sequences (by construction, all
  % but the last angle are within tolerance of one another, so only need to
  % check the last one)
  maxAng = max(candAngs(:, end));
  
  % Filter down to the sequences with big-enough angles
  angTol = min([maxTol, maxAng / 10]);
  bigEnough = candAngs(:, end) >= maxAng - angTol;
  candSeqs = candSeqs(bigEnough, :);
  candAngs = candAngs(bigEnough, :);
  prevSeqID = prevSeqID(bigEnough);
  
  % Filter for ordering degeneracy
  % Note that the very first subspace is special. So, [2 1 5] and [2 5 1]
  % are the same, but [1 2 5] and [2 1 5] are not. We therefore keep the
  % first element of the sequence fixed.
  sortedSeqs = [candSeqs(:, 1), sort(candSeqs(:, 2:end), 2)];
  [~, uSeqs] = unique(sortedSeqs, 'rows');
  candSeqs = candSeqs(uSeqs, :);
  candAngs = candAngs(uSeqs, :);
  prevSeqID = prevSeqID(uSeqs);
  
  
  %% Add new vector for each candidate sequence
  
  nSeqs = size(candSeqs, 1);
  newCandSeqVecs = NaN(nDims, a+nVecs, nSeqs);
  for seq = 1:nSeqs
    [~, ~, V] = subspacea(candSeqVecs(:, :, prevSeqID(seq)), X(:, :, candSeqs(seq, end)));
    newCandSeqVecs(:, 1:end-1, seq) = candSeqVecs(:, :, prevSeqID(seq));
    newCandSeqVecs(:, end, seq) = V(:, end);
  end
  candSeqVecs = newCandSeqVecs;
end


%% Finally, choose an overall winning sequence if we have a tie

if size(candAngs, 1) == 1
  angles = candAngs;
else
  sums = sum(candAngs, 2);
  [~, winner] = max(sums);
  angles = candAngs(winner, :);
end
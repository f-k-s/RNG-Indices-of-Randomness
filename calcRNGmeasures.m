function out = calcRNGmeasures(respSeqsIn, respSet)
%calcRNGmeasures:
%   out = calcRNGmeasures(respSeqsIn, respSet)
%
%
%--------------------------------------------------------------------------
% REFERENCES:
% RNG measure were taken from:
%
%   Towse, J. N., & Neil, D. (1998). Analyzing human random generation 
%       behavior: A review of methods used and a computer program for 
%       describing performance. 
%       Behavior Research Methods, Instruments, & Computers, 30(4), 
%       583-591. doi:10.3758/BF03209475
%   Abbreviated as T&N.98
%
%   Jahanshahi, M., Saleem, T., Ho, A. K., Dirnberger, G., & Fuller, R. 
%       (2006). Random number generation as an index of controlled 
%       processing. Neuropsychology, 20(4), 391-9.
%       doi:10.1037/0894-4105.20.4.391
%   Abbreviated as J&al.06
%
%   Sheskin, David J. Handbook of Parametric and Nonparametric Statistical
%       Procedures. 3rd ed. CRC press 2003.
%   Abbreviated as DJS.03
%
%--------------------------------------------------------------------------
% INPUTS:
%
%   respSeqsIn  A vector or a matrix containing the response sequence(s) in
%               columns.
%
%   respSet     A vector containing the response alternatives.
%
%
%--------------------------------------------------------------------------
% OUTPUTS - A 1xK structure (where K is the number of response sequences)
%           with the following fields: 
%
%   MISS        Missings. MISS is a structure with fields:
%               -.i: Indices of values in the response sequence that are
%               not memebers of the response set.
%               -.n: Number of missing values (numel(MISS.i)).
%               NOTE: all values classified as missing are removed from the
%               response sequence before calculating the output.
%
%   HIST        Histogram of responses. Calculated as 
%                   HIST = histc(responseSequence,responseSet)
%               HIST is a 1xN vector were N = numel(responseSet) and
%               HIST(i) contains the frequency of responseSet(i) in
%               responseSequence.
%
%   CNT         Count scores. CNT(1): Counting Score 1 measures the
%   (J&al.06)   tendency to count in steps of 1 or -1. CNT(2): Counting
%               Score 2 measures the tendency to count in steps of 2 or -2.
%               The scores are calculated as the sum of all squared
%               counting sequence lengths found in responseSequence.
%               Examples:
%               The sequence 3-4-5-6 would increase CNT(1) by 3^2 = 9
%               The sequence 7-5-3 would increase CNT(2) by 2^2 = 4.
%   
%   TRI         Number of unique triplets. There are M-2 triplet in a
%   (J&al.06)   series of M responses. The fewer the number of unique
%               triplets, the greater the tendency to repeat certain
%               steretyped second-order runs.
%
%   DIG         Digrams. DIG is a structure with the following fields:
%   (J&al.06)   -.freqMat: Frequency matrix of digrams. freqMat is a NxN
%                 matrix containing the frequency of all possible pairs.
%                 This includes the 'circular' pair, that is
%                 responseSequence(end)-responseSequence(1).
%                 freqMat(i,j) is the number of occurences of the pair
%                 responseSet(i)-responseSet(j) in responseSequence.
%               -.achieved: Number of achieved digrams, that is the number
%                 of nonzero cells in the .freqMat-Matrix. (Also includes
%                 the circular pair!)
%               -.repIndex: Digram repetition index. The digram repetition
%                 index is the sum of all nonzero cells in the digram
%                 frequency matrix. (Also includes the circular pair!)
%                 repIndex = sum(freqMat(freqMat~=0)-1)
%               -.redIndex: Redundancy index. The redundancy index is the
%                 number of digrams achieved divided by the number of
%                 possible digrams. redIndex = achieved/numel(freqMat)
%
%   CHI         Chi-square goodnes of fit, calculated with the chi2gof
%   (J&al.06)   function:
%                   [~, p, chiStat] =  chi2gof(respSeq,'ctrs',respSet)
%               A structure with the following fields:
%               -.stat = chiStat.chi2stat       s[chi² statistic]
%               -.p = p                         [p value]
%               -.df = chiStat.df               [degrees of freedom]
%
%   COUPON      Coupon measures the mean number of responses produced
%   (T&N.98)    before all the response alternatives are given. One set is
%               defined here as a sequence of responses containing each
%               response alternative at least once. The length of all
%               non-overlapping sets is averaged to calculate the coupon
%               score. Also see: Ginsburg & Karpiuk (1994). Random
%               generation: Analysis of the responses. Perceptual and Motor
%               Skills, 79(1978), 1059-1067. doi:10.2466/pms.1994.79.3.1059
%
%   NSQ         Guttmann's Null-Score Quotient. The null-score (NS) is the
%   (T&N.98)    total number of digram permutations that do not appear in
%               the response sequence (also includes the circular pair!)
%               and takes values between 0 and N^2-1 with N =
%               numel(responseSet). The null-score quotient is 
%                   NSQ = 100 * NS / (N^2-1)
%
%   RNG         Random Number Generation index, introduced by Evans (1978).
%   (T&N.98)    The RNG index describes the distribution of digrams
%               (including the circular digram). It is defined as
%                   RNG = sum(Fij*log2(Fij)) / sum(Fi*log2(Fi))
%               where Fij is the frequency of the digram
%               responseSet(i)-responseSet(j) and Fi is the frequency of
%               the response alternative responseSet(i).
%
%   TPI         Turning Point Index. Turning points mark the change between
%   (T&N.98)    ascending and descending sequences. The turning point index
%               is the number of observed turning points divided by the
%               nubmer of expected turning points, expressed as a
%               percentage. The number of expected turning points is:
%                   TPexpected = 2/3 * (M-2)
%               where M is the number of responses.
%
%   PL          Phase length. The phase length describes the distribution
%   (T&N.98)    of intervals between turning points. PL is a structure with
%               fields for expected and observed phases, containing a count
%               value for each observed phase length. That means
%               PL.observed(d) is the number of phases of length d in the
%               response sequence. See T&N.98 for the equation of expected
%               PL-frequencies.
%
%   REDUND      Redundancy. Redundancy is a percentage score describing the
%   (T&N.98)    frequency distribution of response alternatives. An R score
%               of 0% indicates no redundancy (perfect equality of response
%               alternative frequencies), and an R score of 100% indicates
%               complete redundancy (the same response choice is used
%               throughout).
%
%   FOD         First Order Difference. FOD is a structure with a field
%   (T&N.98)    -.diffs containing the numbers from rSet(1)-rSet(end) to
%               rSet(end)-rSet(1) and a field
%               -.freqs containing the corresponding frequencies of
%               differences between every two consecutive responses, that
%               is diff(responseSequence) (including the circular pair!).
%
%   REPdist     Repetition distance. Repitition distance is the frequency
%   (T&N.98)    of the gap between two identical responses. REPdist is a
%               structure with vectors of length M-1, as the largest gap is
%               M-1 (where M is the number of responses). It contains the
%               following fields:
%               -.length: the length of the gaps. Basically the numbers
%               from 1 to M-1
%               -.freqs: the observed frequencies of each gap.
%               -.expected: the expected frequencies of each gap. This is
%               approximated by a geometric distribution. See T&N.98 for
%               details.
%
%   REPgap      Repetition gap. REPgap is a structure containing the
%   (T&N.98)    following fields:
%               -.median: The median of the repetition distances.
%               -.mean: The mean of the repetition distances.
%               -.mode: The mode(s) of the repetition distances.
%
%   ADJ         Adjacency. Adjacency measure the frequency of digrams
%   (T&N.98)    containing adjacent items (e.g. 1-2 or 5-4). It is
%               expressed seperately for ascending and descending digrams
%               and calculated as a percentage score:
%                   A = 100 * numel(adjacentPairs)/numel(responsePairs)
%               ADJ is a structure with the fields:
%               -.asc for ascending pairs.
%               -.desc for descending pairs.
%               -.comb for the combination of both.
%                   
%   MSSD        Mean squares of the successive differences. MSSD is a
%   (DJS.03)    structure with the following fields:
%               -.mssd: mean square successive differences, calculated as
%                 sum(deriv1.^2) ./ (2*(N-1)) where deriv1 is the first
%                 derivative of the response sequence and N is the number
%                 of responses.
%               -.C: The C-value, calculated as 1-mssd/var(respSeq)
%               -.p: The probability of observing C under the null
%                 hypthesis that the response sequence is random.
%                 The probability distribution was taken from L. C. Young,
%                 'On Randomness in Ordered Sequences,' Ann. Math. Stat.,
%                 vol. 12, no. 3, pp. 293-300, 1941. 
%               -.z: The z-value, calculated as C/sqrt((N-2)/(N^2-1))
%               
%
% V1.2, 2016-10-19, Konrad Schumacher
%

outVars = {'MISS' 'HIST' 'CNT' 'TRI' 'DIG' 'CHI' 'COUPON' 'NSQ' 'RNG' 'TPI' ...
            'PL' 'REDUND' 'FOD' 'REPdist' 'REPgap' 'ADJ' 'MSSD'};

% check input
if numel(size(respSeqsIn))>2, error('Response sequence has to be a vector or a matrix!'); end
singDim = find(size(respSeqsIn)==1);
if ~isempty(singDim) && singDim==1
    respSeqsIn = respSeqsIn';
end
if numel(size(respSet))>2 || ~any(size(respSet)==1)
    error('Response set has to be a vector!');
end
respSet = respSet(:)'; % make sure its a row

% sort reponse set; just in case ...
respSet = sort(respSet); %#ok

setMin = min(respSet);
setMax = max(respSet);
nSet = numel(respSet);


% loop input response sequences

for iRespSeq = 1:size(respSeqsIn,2)
    respSeq = respSeqsIn(:,iRespSeq);

    membr = ismember(respSeq,respSet);
    MISS.i = find(~membr);
    MISS.n = numel(MISS.i);

    % exclude missings from response sequence
    respSeq = respSeq(membr);
%     respSeq(~membr) = NaN;
    
    nResp = numel(respSeq);
    
    % derivatives
    deriv1 = [diff(respSeq,1,1); respSeq(1)-respSeq(end)];
    deriv2 = [diff(deriv1,1,1); deriv1(1)-deriv1(end)];

    %% repeated pairs (REP)
    nREP = numel(find(deriv1==0));

    %% first order difference FOD
    FOD.diffs = setMin-setMax:setMax-setMin;
    FOD.freqs = histc(deriv1,FOD.diffs);

    %% adjacency
    ADJ.asc = 100 * numel(find(deriv1==1))/numel(deriv1);
    ADJ.desc = 100 * numel(find(deriv1==-1))/numel(deriv1);
    ADJ.comb = 100 * numel(find(abs(deriv1)==1))/numel(deriv1);

    %% histogram
    respFreq = histc(respSeq,respSet);
    HIST = respFreq;

    %% median gap score GAP
    gaps = arrayfun(@(x)diff(find(x==respSeq)), respSet', 'UniformOutput',0);
    gaps = cell2mat(gaps);
    REPgap.median = median(gaps);
    REPgap.mean = mean(gaps);
    REPgap.mode = mode(gaps);

    %% repitions distance
    REPdist.length = 1:nResp-1;
    REPdist.freqs = histc(gaps,1:nResp-1)';
    repFrqFun = @(s)(1-1/nSet)^(s-1)*1/nSet*(nResp-1);
    REPdist.expected = arrayfun(repFrqFun,1:nResp-1);

    %% first order information (Towse&Neil 1998)
    H.singl = log2(nResp) - 1/nResp*sum(respFreq.*log2(respFreq));
    H.max = log2(nSet);
    REDUND = 100*(1-H.singl/H.max);

    %% Turning Point Index
    nTPexpect = 2/3*(nResp-2);
    % rm plateaus
    respNoPlat = respSeq(deriv1~=0);
    respNoPlatD1sig = sign(diff(respNoPlat,1,1));
    noPlatTurnPts = find(abs(diff(respNoPlatD1sig))==2)+1;
    nTPobserv = numel(noPlatTurnPts);
    TPI = 100 * nTPobserv/nTPexpect;

    %% Phase length PL
    PLexpFun = @(n,d)2*(n-d-2)*(d^2+3*d+1) / factorial(d+3);
    %reconstruct original turning points by inserting plateau points again
    turnPts = noPlatTurnPts;
    for i = find(deriv1==0)'
        incPts = i<=turnPts; % on plateaus use last item according to RGCalc by Neil&Towse
        turnPts(incPts) = turnPts(incPts)+1;
    end
    tpD1 = diff(turnPts);
    phases = unique(tpD1)';
    nPhases = numel(phases);
    PL.expected = arrayfun(PLexpFun, repmat(nResp,1,nPhases), phases);
    PL.observed = arrayfun(@(p)numel(find(tpD1==p)),phases);

    %% Pair-frequency matrix
    % preferred pairs
    pairs = [respSeq circshift(respSeq,-1)];
    pairs = pairs(1:end-1,:); % remove circular pair
    [uniPairs, iPairs, iUPairs] = unique(pairs,'rows');%,'stable'
    pairFreq = histc(iUPairs,1:numel(iPairs));

    % preferred triplets
    triplets = [respSeq circshift(respSeq,-1) circshift(respSeq,-2)];
    triplets = triplets(1:end-2,:);% remove circular triplets
    [uniTrip, iTrip, iUTrip] = unique(triplets,'rows');%,'stable'
    triplFreq = histc(iUTrip,1:numel(iTrip));

    % pair-frequency matrix
    pairFreqMat = zeros(nSet);
    for i = 1:numel(pairFreq)
        pairFreqMat(uniPairs(i,1),uniPairs(i,2)) = pairFreq(i);
    end
    % should be wrap-around (Towse&Neil 1998), so increment circular pair (last-first response):
    pairFreqMat(respSeq(end),respSeq(1)) = pairFreqMat(respSeq(end),respSeq(1)) + 1;

    %% RNG score Evans 1978
    matI = pairFreqMat > 1;
    respI = respFreq > 1;
    RNG = sum(pairFreqMat(matI).*log2(pairFreqMat(matI))) / sum(respFreq(respI).*log2(respFreq(respI)));

    %% Null-Score Quotient NSQ
    NSQ = 100 * numel(find(pairFreqMat==0))/(numel(pairFreqMat)-1);

    %% Coupon
    iStart = 1;
    iEnd = nSet-1;
    seqLenSum = 0; nSequ = 0;
    while iStart+iEnd < nResp
        occResp = ismember(respSet,respSeq(iStart:iStart+iEnd));
        while ~all(occResp) && iStart+iEnd < nResp
            iEnd = iEnd + 1;
            occResp = ismember(respSet,respSeq(iStart:iStart+iEnd));
        end
        if all(occResp)
            seqLenSum = seqLenSum + iEnd + 1;
            nSequ = nSequ + 1;
        end
        iStart = iStart+iEnd+1;
        iEnd = nSet-1;
    end
    COUPON = seqLenSum/nSequ;


    %% Jahanshahi 2006: ...................................................
    
    %% Digrams achieved DIG
    DIG.achieved = numel(find(pairFreqMat));
    DIG.freqMat = pairFreqMat;
    
    %% Digram repetition index DRI
    DIG.repIndex = sum(pairFreqMat(pairFreqMat~=0)-1);
    
    %% Redundancy Index RI
    DIG.redIndex = DIG.achieved/numel(pairFreqMat);
    
    %% Unique triplets TRI
    TRI = size(uniTrip,1);

    %% Counting score
    CNT = [0 0];
    for cntDif = [-1 1 -2 2]
        cntPos = 0;
        while cntPos < nResp
            cntPos = find(deriv1(cntPos+1:end)==cntDif, 1, 'first') + cntPos;
            cntLen = 0;
            while ~isempty(cntPos) && cntPos<=nResp && deriv1(cntPos)==cntDif
                cntLen = cntLen + 1;
                cntPos = cntPos + 1;
            end
            CNT(abs(cntDif)) = CNT(abs(cntDif)) + cntLen^2;
        end
    end

    %% chi square goodnes of fit
    [~,CHI.p,chiStat] =  chi2gof(respSeq,'ctrs',respSet);
    CHI.stat = chiStat.chi2stat;
    CHI.df = chiStat.chi2stat;
    
    
    %% Mean square successive differences
    % Sheskin, D.J. 2003, Handbook of Parametric and Nonparametric
    % Statistical Procedures: Third Edition, page 365
    % https://books.google.de/books?id=bmwhcJqq01cC
    % Mean of the Squares of the Successive Differences (excl. circular
    % derivative)
    MSSD.mssd = sum(deriv1(1:end-1).^2) ./ (2*(nResp-1));
    MSSD.C = 1 - (MSSD.mssd/var(respSeq));
    % SHESHKIN '03:
    % "In order to rejct the null hypothesis and conclude that the
    % distribution is nonrandom, the abolute value of C must be equal to or
    % greater than the tabled critical value of the C statistic at the
    % prespecified level of significance. A large absolute C value
    % indicates there is a larger discrepancy between the values s^2
    % (variance) and s^2ms (MSSD). A table of critical values can be found
    % in Zar (1999), who developed his table based on the work of Young
    % (1941). Another table of the sampling districution for this test
    % (although not C values) was derived by Hart (1942), and can be found
    % in Bennett and Frankling (1954)."
    MSSD.p = mssd_pval(MSSD.C, nResp);
    % "In lieu of the exact table for the sampling distribution of C, we
    % will employ a large sample normal approximation of the C statistic
    % which is computed with eq. 10.11 (z = ...)."
    MSSD.z = MSSD.C / sqrt((nResp-2)/(nResp.^2-1));
    
    
    %% ##### produce output var: #####

    for i = 1:numel(outVars)
        out(iRespSeq).(outVars{i}) = eval(outVars{i});
    end
end

end % end main function



%% :::::::::::::::::::::::::::::  SUBFUNCTIONS ::::::::::::::::::::::::::::

function p = mssd_pval(C,n)
% Returns the p-value for given C & given number of items.
%
% C - The C-value: 1-mssd/var(items)
% n - The number of items: numel(items)
%
% taken from L. C. Young, On Randomness in Ordered Sequences,
% Ann. Math. Stat., vol. 12, no. 3, pp. 293-300, Sep. 1941.
%
% K. Schumacher, 2015

a  = sqrt( (n.^2+2.*n-12) .* (n-2) ...
      ./ (n.^3 - 13.*n + 24) );
     
m  = (n.^4 - n.^3 - 13.*n.^2 + 37.*n - 60) ...
      ./ (2.*(n.^3 - 13.*n + 24));

% this gets way to large for n>170, because gamma(172)=Inf ...
    % y0 = gamma(2.*m+2) ...
    %       ./ (a.*2.^(2.*m+1) .* gamma(m+1).^2);
%... so we take the logarithm of this equation (gammaln gives results for a
%    wider range of input values): 
y0log = gammaln(2.*m+2) - (log(a.*2.^(2.*m+1)) + 2.*gammaln(m+1));
%... and get back to y0:
y0 = exp(1) .^ y0log;

% the probability density function (somehow very tiny imaginary parts are
% introduced, so just take the real part): 
y = @(c) real( y0 .* (1-c.^2./(a.^2)).^m );

% compute the cumulative density function of y and return the
% probability of the absolute value of C:
p = (1 - integral(y, -1, abs(C)) );

end
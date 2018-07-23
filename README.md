# RNG-Indices-of-Randomness

Implementation of various indices of randomness as measures of RandomNumberGeneration performance; a task assessing executive functioning.

Code has not been externally reviewed; no warranty!


## REFERENCES:
RNG indices were adapted from:

  Towse, J. N., & Neil, D. (1998). Analyzing human random generation 
      behavior: A review of methods used and a computer program for 
      describing performance. 
      Behavior Research Methods, Instruments, & Computers, 30(4), 
      583-591. doi:10.3758/BF03209475
  Abbreviated as T&N.98

  Jahanshahi, M., Saleem, T., Ho, A. K., Dirnberger, G., & Fuller, R. 
      (2006). Random number generation as an index of controlled 
      processing. Neuropsychology, 20(4), 391-9.
      doi:10.1037/0894-4105.20.4.391
  Abbreviated as J&al.06

  Sheskin, David J. Handbook of Parametric and Nonparametric Statistical
      Procedures. 3rd ed. CRC press 2003.
  Abbreviated as DJS.03


## INDICES:

###  MISS        Missings. 
    MISS is a structure with fields:
              -.i: Indices of values in the response sequence that are
              not memebers of the response set.
              -.n: Number of missing values (numel(MISS.i)).
              NOTE: all values classified as missing are removed from the
              response sequence before calculating the output.

###  HIST        Histogram of responses. 
    Calculated as: 
              HIST = histc(responseSequence,responseSet)
              HIST is a 1xN vector were N = numel(responseSet) and
              HIST(i) contains the frequency of responseSet(i) in
              responseSequence.

###  CNT (J&al.06)        Count scores. 
    CNT(1): Counting Score 1 measures the tendency to count in steps of 1 or -1.
    CNT(2): Counting Score 2 measures the tendency to count in steps of 2 or -2.
              The scores are calculated as the sum of all squared
              counting sequence lengths found in responseSequence.
              Examples:
              The sequence 3-4-5-6 would increase CNT(1) by 3^2 = 9
              The sequence 7-5-3 would increase CNT(2) by 2^2 = 4.
  
###  TRI (J&al.06)        Number of unique triplets. 
    There are M-2 triplet in a series of M responses. The fewer the number of unique
              triplets, the greater the tendency to repeat certain
              steretyped second-order runs.

###  DIG (J&al.06)       Digrams.
    DIG is a structure with the following fields:
     -.freqMat: Frequency matrix of digrams. freqMat is a NxN
                matrix containing the frequency of all possible pairs.
                This includes the 'circular' pair, that is
                responseSequence(end)-responseSequence(1).
                freqMat(i,j) is the number of occurences of the pair
                responseSet(i)-responseSet(j) in responseSequence.
              -.achieved: Number of achieved digrams, that is the number
                of nonzero cells in the .freqMat-Matrix. (Also includes
                the circular pair!)
              -.repIndex: Digram repetition index. The digram repetition
                index is the sum of all nonzero cells in the digram
                frequency matrix. (Also includes the circular pair!)
                repIndex = sum(freqMat(freqMat~=0)-1)
              -.redIndex: Redundancy index. The redundancy index is the
                number of digrams achieved divided by the number of
                possible digrams. redIndex = achieved/numel(freqMat)

###  CHI (J&al.06)        Chi-square goodnes of fit
    Calculated with the chi2gof function:
              [~, p, chiStat] =  chi2gof(respSeq,'ctrs',respSet)
              Returns a structure with the following fields:
              -.stat = chiStat.chi2stat       [chi² statistic]
              -.p = p                         [p value]
              -.df = chiStat.df               [degrees of freedom]

###  COUPON (T&N.98)     
    Coupon measures the mean number of responses produced
              before all the response alternatives are given. One set is
              defined here as a sequence of responses containing each
              response alternative at least once. The length of all
              non-overlapping sets is averaged to calculate the coupon
              score. Also see: Ginsburg & Karpiuk (1994). Random
              generation: Analysis of the responses. Perceptual and Motor
              Skills, 79(1978), 1059-1067. doi:10.2466/pms.1994.79.3.1059

###  NSQ (T&N.98)        Guttmann's Null-Score Quotient. 
    The null-score (NS) is the total number of digram permutations that do not appear in
              the response sequence (also includes the circular pair!)
              and takes values between 0 and N^2-1 with N =
              numel(responseSet). The null-score quotient is 
                  NSQ = 100 * NS / (N^2-1)

###  RNG (T&N.98)        Random Number Generation index, introduced by Evans (1978).
      The RNG index describes the distribution of digrams
              (including the circular digram). It is defined as
                  RNG = sum(Fij*log2(Fij)) / sum(Fi*log2(Fi))
              where Fij is the frequency of the digram
              responseSet(i)-responseSet(j) and Fi is the frequency of
              the response alternative responseSet(i).

###  TPI (T&N.98)        Turning Point Index. 
    Turning points mark the change between ascending and descending sequences. 
              The turning point index is the number of observed turning points 
              divided by the nubmer of expected turning points, expressed as a
              percentage. The number of expected turning points is:
                  TPexpected = 2/3 * (M-2)
              where M is the number of responses.

###  PL (T&N.98)         Phase length. 
    The phase length describes the distribution of intervals between turning points.
              PL is a structure withfields for expected and observed phases, 
              containing a count value for each observed phase length. That means
              PL.observed(d) is the number of phases of length d in the
              response sequence. See T&N.98 for the equation of expected
              PL-frequencies.

###  REDUND (T&N.98)     Redundancy.
    Redundancy is a percentage score describing the frequency distribution
              of response alternatives. An R score of 0% indicates no redundancy 
              (perfect equality of response alternative frequencies), 
              and an R score of 100% indicates complete redundancy 
              (the same response choice is used throughout).

###  FOD (T&N.98)        First Order Difference.
    FOD is a structure with a field
      -.diffs containing the numbers from rSet(1)-rSet(end) to
              rSet(end)-rSet(1) and a field
      -.freqs containing the corresponding frequencies of
              differences between every two consecutive responses, that
              is diff(responseSequence) (including the circular pair!).

###  REPdist (T&N.98)    Repetition distance.
    Repitition distance is the frequency of the gap between two identical responses. 
              REPdist is a structure with vectors of length M-1, as the largest gap is
              M-1 (where M is the number of responses). It contains the
              following fields:
              -.length: the length of the gaps. Basically the numbers from 1 to M-1
              -.freqs: the observed frequencies of each gap.
              -.expected: the expected frequencies of each gap. This is
              approximated by a geometric distribution. See T&N.98 for
              details.

###  REPgap (T&N.98)     Repetition gap. 
    REPgap is a structure containing the following fields:
              -.median: The median of the repetition distances.
              -.mean: The mean of the repetition distances.
              -.mode: The mode(s) of the repetition distances.

###  ADJ (T&N.98)        Adjacency.
    Adjacency measure the frequency of digrams containing adjacent items (e.g. 1-2 or 5-4).
          It is expressed seperately for ascending and descending digrams
          and calculated as a percentage score:
                  A = 100 * numel(adjacentPairs)/numel(responsePairs)
          ADJ is a structure with the fields:
              -.asc for ascending pairs.
              -.desc for descending pairs.
              -.comb for the combination of both.
                  
###  MSSD (DJS.03)       Mean squares of the successive differences. 
    MSSD is a structure with the following fields:
              -.mssd: mean square successive differences, calculated as
                sum(deriv1.^2) ./ (2*(N-1)) where deriv1 is the first
                derivative of the response sequence and N is the number
                of responses.
              -.C: The C-value, calculated as 1-mssd/var(respSeq)
              -.p: The probability of observing C under the null
                hypthesis that the response sequence is random.
                The probability distribution was taken from L. C. Young,
                'On Randomness in Ordered Sequences,' Ann. Math. Stat.,
                vol. 12, no. 3, pp. 293-300, 1941. 
              -.z: The z-value, calculated as C/sqrt((N-2)/(N^2-1))

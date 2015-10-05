# original asm settings

assemble = 0

utgErrorRate = 0.25
utgErrorLimit = 6.5

cnsErrorRate = 0.25
cgwErrorRate = 0.25
ovlErrorRate = 0.25

merSize=14

merylMemory = 128000
merylThreads = 56
merOverlapperThreads=1
merOverlapperExtendConcurrency=56
merOverlapperSeedConcurrency=56

ovlStoreMemory = 8192

# grid info
useGrid = 0
scriptOnGrid = 0
frgCorrOnGrid = 0
ovlCorrOnGrid = 0

ovlHashBits = 24
maxGap=1500
ovlHashBlockLength = 1000000000
ovlRefBlockLength = 1000000000
ovlRefBlockSize = 0
blasr=-noRefineAlign -advanceHalf -noSplitSubreads -minMatch 10 -minPctIdentity 70 -bestn 24 -nCandidates 24

frgCorrThreads = 1
frgCorrConcurrency = 56
frgCorrBatchSize = 100000

ovlCorrBatchSize = 100000

ovlThreads = 1
ovlConcurrency = 56
cnsConcurrency = 56

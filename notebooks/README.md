This repository of notebooks is my exploratory analysis of the NINDS paRemap data. The files can be viewed in the following order:

1. wordmatch Distances Analysis: This notebook looked at the different word match pairs in the paRemap task. I computed various distance metrics/features from each of these events to somewhat characterize a quantitative distance between brain states when each of these word pairs were encoded.
2. Classificaion Analysis: This notebook performed different classification testing with different algorithms. Most performed at chance given the features we input as the distances computed from the previous notebook. It would be good to retest when our features become more fine tuned.
3. exploratory_analysis: This notebook looked at the differences in response times for each of the probe words. Also planning on looking at each of the word pairs and looking at the different response times. 

Questions to Answer:
1. Are different word pairs encoded differently?
2. What happens when you filter out events by response times (e.g. fast and slow)?

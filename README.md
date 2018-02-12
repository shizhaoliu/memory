# memory
Demo code for paper "Methods for Assessment of Memory Reactivation"



This MATLAB software includes custom functions written for detecting significant memory reactivation events
(assuming the sorted population spike activity during the candidate event is given).
We include a demon code with examples.

The data examples contains two candidate events.

In the demo, we show both results based on the supervsied decoding method (both Rw and Rwd are calculated)
and the unsupervised decoding method (only Rw are calculated). For both methods, supervised Rw
methods returns negative results while supervised and unsupervised returns positive results.
Detailed documentation can be found in the functions.


Citation: 
  Shizhao Liu, Andres D. Grosmark, Zhe Chen,  "Methods for Assessment of Memory Reactivation",  Neural Computation, 2018.

Notation: 
  Rw (linear weighted correlation), Rwd (weighted distance correlation), Zw, Zwd (their corresponding Z-score) 

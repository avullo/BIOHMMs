# BIOHMMs
An implementation of non-causal input-output HMMs for probabilistic sequence modeling.

## Introduction

A bidirectional IOHMM is a non-causal model of a stochastic translation defined on a space
of finite sequences.

The motivation for this architecture is the ability to model non temporal sequences which appear in some
important domains like computational biology. Causal models (i.e. RNNs, IOHMMs) are *causal* in the sense
that the output at time t does not depend on future inputs. Causality is easy to justify in dynamics that
attempt to model the behavior of physical systems, or that need to operate in real time. Clearly, in these
cases the response at time t cannot depend on stimulae that the system has not yet encountered. But biological
sequences are not temporal: the conformation and the function of a region in a sequence may strongly depend
on events located both upstream and downstream. 

A bidirectional IOHMMs extends IOHMMs by introducing non-causal bidirectional dynamics to capture both upstream
and downstream information. Like IOHMMs, the model describes the conditional probability distribution
P(Y|U), where U = U1, U2, ..., UT is a generic input sequence and Y = Y1, Y2, ..., YT the corresponding output
sequence.

The model is based on two Markov sequences of hidden state variables, denoted by F and B, respectively.
For each time step, F(t) and B(t) are discrete variables with realizations (states) in {f1 , ..., fn }
and {b1 , ..., bm }, respectively. As in HMMs, F(t) is assumed to have a causal impact on the next state
F(t+1). Hence, F(t) stores contextual information contained on the left of t (propagated in the forward
direction). Symmetrically, B(t) is assumed to have a causal impact on the state B(t−1) , thus summarizing
information contained on the right of t (propagated in the backward direction).

## References

P. Baldi, S. Brunak, P. Frasconi, G. Soda and G. Pollastri. "Exploiting the Past and the Future in Protein Secondary Structure Prediction", *Bioinformatics*, 15, 937-946, 1999.
Y. Bengio and P. Frasconi. Input-output HMM’s for sequence processing. *IEEE Trans. on Neural Networks*, 7:1231–1249, 1996.
P. Smyth, D. Heckerman, and M. I. Jordan. Probabilistic independence networks for hidden markov probability models. *Neural Computation*, 9(2):227–269, 1997.



# Morphogenesis dynamics (Shh, Gli & Ptc)
Masters degree thesis about the sonic Hedgehog signaling system.

During human development, cells are exposed to a complex network of regulatory signals, which must  be interpreted correctly in order to success doing functions necessary for the organism. Therefore, the transition of signals and cascades of genetic regulation can be understood as mechanisms of information processing that translate extracellular information into intracellular decisions.
	
The present work aims to show the differences, in terms of a qualitative behavior, that can be found through these complex systems of biological regulation process through different theoretical approaches.
	
We were particularly interested in the Shh signaling system, among its many roles during development, patterns spinal cord and limb bud tissue differentiation and controls midbrain and ventral forebrain neuronal differentiation. 
	
There are a few models that has been described in order to understand its behaviour. The most important was developed in \cite{schaffer}. In particular, these mathematicians use the thermodinamic approach to Shh's gene expression mechanish. This approach, in rough outlines, aim to model the gene transcription systems enumerating all the possible states of the promoter and enchancers of gene transcription activation, and then, relating all of them with their theoretically calculate transcriptional activation level. While these steps can be done in multiple ways, they focused their work on the stimulated approach, which links the transcriptional activity to the transcriptional factors (uniquely).
	
However, some discrepancies has been observed in biological experiments, mainly focused in the existence of an unique stable state, casting aside the biestable swicht behavior shown in the classic model.
	
As a result of these state of the art, we thought it would be necessary to offer a new model that update the classic one, but at the same time, we still want to use the BEWARE approach.
	
Our main goal has been to develop and study a new model based on the same BEWARE strategy but approaching it by the recruitment perspective, that is, taking into account the fundamental part of the RNAp in this whole biological process and its hability of produce the transcriptional activation. 
	 
This work presents a deduction of our new model based on the approach made in \cite{cambon1} and a qualitative study of both models (old and new) through numerical simulations (numerical integration, parameter tunning, bifurcation diagrams, steady states, etc), our discoveries during these and our conclusions about the new and old models. 
	 
Specifically, we linked the biestable behaviour of the old model not only to the Shh/$K_{Shh}$ ratio in our cell but with the basal rate of Gli3 (a common trancription factor in this process). We haven't seen that last link in any paper.
	 
 Furthermore, we found that our new model offers a surprising searched conclusion. Even though the main thermodinamic operator exhibit a highly similar behaviour, our global simulations show a single steady state, no matter how far we alter our parameters. 
	 
We hope that these experiments should motivate a deep analytic research of our new model, because our results suggest it could fit the actual biological paradigm, helping us to understand quite a lot about this important topic.

## Code

![Python scripts](https://github.com/thebooort/shh-signal-model/tree/master/scripts/python_codes)


![Octave scripts](https://github.com/thebooort/shh-signal-model/tree/master/scripts/octave_codes)


![XPPAUT/AUTO07-P scripts](https://github.com/thebooort/shh-signal-model/tree/master/scripts/xppaut_codes)

![Jupyter notebook for symbolic calculations scripts](https://github.com/thebooort/shh-signal-model/tree/master/scripts/Jupyter_Notebook_For_Symbolic_Calculations)

## Document (not finished, errors and similar expected)

![Pdf](https://github.com/thebooort/shh-signal-model/blob/master/thesis.pdf)

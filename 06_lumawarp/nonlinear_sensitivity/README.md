# 31_nonmonotone_54node: the paper's generator with curved mechanisms (experiment 1 of the rework note)

The code change lives in the pipeline repository, branch `nonmonotone-variant`
(`../../causal-dag-eval-pipeline`). This folder holds the pre-registration and, once run, the scores.

    cd ../../causal-dag-eval-pipeline && git checkout nonmonotone-variant
    make recalibrate && make data && make shape      # needs renv restored (dagitty, rCausalMGM)

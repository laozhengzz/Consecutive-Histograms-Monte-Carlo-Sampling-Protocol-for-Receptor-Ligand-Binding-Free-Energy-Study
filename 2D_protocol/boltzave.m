function ave = boltzave(vec)
expvec = exp(vec/-0.5918);
ave = sum(vec.*expvec)/sum(expvec);
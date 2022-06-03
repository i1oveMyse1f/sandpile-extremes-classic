# Prediction of Extremes in BTW and Manna Sandpiles

## Abstract

The state-of-the-art in the theory of self-organized criticality exposes that a certain quiescence precedes events that are located on the tail of the probability distribution of events with respect to their sizes. The existence of the quiescence allows us to predict the occurrence of these events in advance. In this work, we estimate the predictability of the Bak-Tang-Wiesenfeld (BTW) and Manna models on the square lattice in the thermodynamic limit defined by the tendency of the system volume to infinity. For both models, we define an algorithm that forecasts the occurrence of large events after a fall in activity. The collapse of the algorithm efficiency computed with various lattices is found if the size of events is normalized by a power-law function of the lattice length. The power-law exponents are $2.75$ and $3$ for the Manna and BTW models respectively. This yields that the prediction in thermodynamic limit does not exist in the BTW but not in the Manna model, at least based on the quiescence.

**Keywords:** Self-organized criticality, Manna model, BTW model, Scalability.

## Results

All results provided in [results](./results/) directory.

## How to run?

To initialize repository, run the following script in the shell:

```bash
python ./init.py
```

To generate data, run the following script in the shell:

```bash
./exec_files/generation L model
```

where `L` is a lattice size., for example, 64, and `model` is a one of `rand` and `determ` - a type of Sandpile model (`rand` - Manna Sandpile, `determ` - BTW Sandpile)

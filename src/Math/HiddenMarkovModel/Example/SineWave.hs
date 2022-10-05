{- |
Example of an HMM with continuous emissions.
We train a model to accept sine waves of a certain frequency.

There are four hidden states: 'Rising', 'High', 'Falling', 'Low'.
-}
module Math.HiddenMarkovModel.Example.SineWave
{-# WARNING "do not import that module, it is only intended for demonstration" #-}
   (module Math.HiddenMarkovModel.Example.SineWavePrivate)
   where

import Math.HiddenMarkovModel.Example.SineWavePrivate

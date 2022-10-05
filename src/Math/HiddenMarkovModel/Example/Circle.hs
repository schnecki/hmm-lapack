{- |
Example of an HMM with continuous emissions with two-dimensional observations.
We train a model to accept a parametric curve of a circle with a certain speed.
This is like "Math.HiddenMarkovModel.Example.SineWave" but in two dimensions.

The four hidden states correspond to the four quadrants.
-}
module Math.HiddenMarkovModel.Example.Circle
{-# WARNING "do not import that module, it is only intended for demonstration" #-}
   (module Math.HiddenMarkovModel.Example.CirclePrivate) where

import Math.HiddenMarkovModel.Example.CirclePrivate

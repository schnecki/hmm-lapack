{- |
This is an example of an HMM with discrete emissions.
We model a traffic light consisting of the colors red, yellow, green,
where only one lamp can be switched on at every point in time.
This way, when it is yellow you cannot tell immediately
whether it will switch to green or red.
We can only infer this from the light seen before.

There are four hidden states:
'StateRed' emits red, 'StateYellowRG' emits yellow between red and green,
'StateGreen' emits green, 'StateYellowGR' emits yellow between green and red.

We quantise time in time steps.
The transition matrix of the model 'hmm' encodes
the expected duration of every state counted in time steps
and what states follow after each other.
E.g. transition probability of 0.8 of a state to itself means
that the expected duration of the state is 5 time steps (1/(1-0.8)).
However, it is a geometric distribution,
that is, shorter durations are always more probable.

The distribution of 'hmm' encodes which lights a state activates.
In our case everything is deterministic:
Every state can switch exactly one light on.

Given a sequence of observed lights
the function 'HMM.reveal' tells us the most likely sequence of states.
We test this with the light sequences in 'stateSequences'
where we already know the hidden states
as they are stored in 'labeledSequences'.
'verifyRevelation' compares the computed state sequence with the given one.

We also try some trainings in 'hmmTrainedSupervised' et.al.
-}
module Math.HiddenMarkovModel.Example.TrafficLight
{-# WARNING "do not import that module, it is only intended for demonstration" #-}
   (
   HMM,
   Color(..),
   hmm,
   hmmDisturbed,
   red, yellowRG, green, yellowGR,
   labeledSequences,
   hmmTrainedSupervised,
   stateSequences,
   hmmTrainedUnsupervised,
   hmmIterativelyTrained,
   ) where

import Math.HiddenMarkovModel.Example.TrafficLightPrivate

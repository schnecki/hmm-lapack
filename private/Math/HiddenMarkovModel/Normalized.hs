{-# LANGUAGE TypeFamilies #-}
{- |
Counterparts to functions in "Math.HiddenMarkovModel.Private"
that normalize interim results.
We need to do this in order to prevent
to round very small probabilities to zero.
-}
module Math.HiddenMarkovModel.Normalized where

import qualified Math.HiddenMarkovModel.Public.Distribution as Distr
import Math.HiddenMarkovModel.Private
          (T(..), Trained(..), emission,
           biscaleTransition, revealGen, sumTransitions)
import Math.HiddenMarkovModel.Utility (normalizeFactor, normalizeProb)

import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix ((-*#), (#*|))
import Numeric.LAPACK.Vector (Vector)

import qualified Numeric.Netlib.Class as Class

import qualified Control.Functor.HT as Functor

import qualified Data.Array.Comfort.Storable as StorableArray
import qualified Data.Array.Comfort.Shape as Shape

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Foldable as Fold
import Data.Traversable (Traversable)


{- $setup
>>> import qualified Data.NonEmpty as NonEmpty
-}


{- |
Logarithm of the likelihood to observe the given sequence.
We return the logarithm because the likelihood can be so small
that it may be rounded to zero in the choosen number type.
-}
logLikelihood ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh, Floating prob,
    Class.Real prob, Distr.Emission typ prob ~ emission,
    Traversable f) =>
   T typ sh prob -> NonEmpty.T f emission -> prob
logLikelihood hmm = Fold.sum . fmap (log . fst) . alpha hmm

alpha ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh,
    Class.Real prob, Distr.Emission typ prob ~ emission,
    Traversable f) =>
   T typ sh prob ->
   NonEmpty.T f emission -> NonEmpty.T f (prob, Vector sh prob)
alpha hmm (NonEmpty.Cons x xs) =
   let normMulEmiss y = normalizeFactor . Vector.mul (emission hmm y)
   in  NonEmpty.scanl
          (\(_,alphai) xi -> normMulEmiss xi (transition hmm #*| alphai))
          (normMulEmiss x (initial hmm))
          xs

beta ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh,
    Class.Real prob, Distr.Emission typ prob ~ emission,
    Traversable f, NonEmptyC.Reverse f) =>
   T typ sh prob ->
   f (prob, emission) -> NonEmpty.T f (Vector sh prob)
beta hmm =
   nonEmptyScanr
      (\(ci,xi) betai ->
         Vector.scale (recip ci) $
         Vector.mul (emission hmm xi) betai -*# transition hmm)
      (Vector.one $ StorableArray.shape $ initial hmm)

alphaBeta ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh,
    Class.Real prob, Distr.Emission typ prob ~ emission,
    Traversable f, NonEmptyC.Zip f, NonEmptyC.Reverse f) =>
   T typ sh prob ->
   NonEmpty.T f emission ->
   (NonEmpty.T f (prob, Vector sh prob), NonEmpty.T f (Vector sh prob))
alphaBeta hmm xs =
   let calphas = alpha hmm xs
   in  (calphas,
        beta hmm $ NonEmpty.tail $ NonEmptyC.zip (fmap fst calphas) xs)


xiFromAlphaBeta ::
   (Distr.EmissionProb typ, Shape.C sh, Eq sh,
    Class.Real prob, Distr.Emission typ prob ~ emission,
    Traversable f, NonEmptyC.Zip f) =>
   T typ sh prob ->
   NonEmpty.T f emission ->
   NonEmpty.T f (prob, Vector sh prob) ->
   NonEmpty.T f (Vector sh prob) ->
   f (Matrix.Square sh prob)
xiFromAlphaBeta hmm xs calphas betas =
   let (cs,alphas) = Functor.unzip calphas
   in  NonEmptyC.zipWith4
          (\x alpha0 c1 beta1 ->
             Matrix.scale (recip c1) $ biscaleTransition hmm x alpha0 beta1)
          (NonEmpty.tail xs)
          (NonEmpty.init alphas)
          (NonEmpty.tail cs)
          (NonEmpty.tail betas)

zetaFromAlphaBeta ::
   (Shape.C sh, Eq sh, Class.Real prob, NonEmptyC.Zip f) =>
   NonEmpty.T f (prob, Vector sh prob) ->
   NonEmpty.T f (Vector sh prob) ->
   NonEmpty.T f (Vector sh prob)
zetaFromAlphaBeta calphas betas =
   NonEmptyC.zipWith (Vector.mul . snd) calphas betas


{- |
Reveal the state sequence
that led most likely to the observed sequence of emissions.
It is found using the Viterbi algorithm.
-}
reveal ::
   (Distr.EmissionProb typ, Shape.InvIndexed sh, Eq sh, Shape.Index sh ~ state,
    Distr.Emission typ prob ~ emission, Class.Real prob, Traversable f) =>
   T typ sh prob -> NonEmpty.T f emission -> NonEmpty.T f state
reveal = revealGen normalizeProb


{- |
Variant of NonEmpty.scanr with less stack consumption.

prop> \x xs -> nonEmptyScanr (-) x xs == NonEmpty.scanr (-) x (xs::[Int])
-}
nonEmptyScanr ::
   (Traversable f, NonEmptyC.Reverse f) =>
   (a -> b -> b) -> b -> f a -> NonEmpty.T f b
nonEmptyScanr f x =
   NonEmptyC.reverse . NonEmpty.scanl (flip f) x . NonEmptyC.reverse


{- |
Consider a superposition of all possible state sequences
weighted by the likelihood to produce the observed emission sequence.
Now train the model with respect to all of these sequences
with respect to the weights.
This is done by the Baum-Welch algorithm.
-}
trainUnsupervised ::
   (Distr.Estimate typ, Shape.C sh, Eq sh,
    Class.Real prob, Distr.Emission typ prob ~ emission) =>
   T typ sh prob -> NonEmpty.T [] emission -> Trained typ sh prob
trainUnsupervised hmm xs =
   let (alphas, betas) = alphaBeta hmm xs
       zetas = zetaFromAlphaBeta alphas betas
       zeta0 = NonEmpty.head zetas

   in  Trained {
          trainedInitial = zeta0,
          trainedTransition =
             sumTransitions hmm $ xiFromAlphaBeta hmm xs alphas betas,
          trainedDistribution =
             Distr.accumulateEmissionVectors $ NonEmptyC.zip xs zetas
       }

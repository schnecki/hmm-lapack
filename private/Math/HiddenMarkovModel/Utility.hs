module Math.HiddenMarkovModel.Utility where

import qualified Numeric.LAPACK.Matrix.Hermitian as Hermitian
import qualified Numeric.LAPACK.Matrix.Layout as Layout
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix.Array as ArrMatrix
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Array (ArrayMatrix)
import Numeric.LAPACK.Vector (Vector, (.*|))

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as StorableArray
import qualified Data.Array.Comfort.Boxed as Array
import qualified Data.Array.Comfort.Shape as Shape

import Foreign.Storable (Storable)

import qualified System.Random as Rnd

import qualified Control.Monad.Trans.State as MS


normalizeProb :: (Shape.C sh, Class.Real a) => Vector sh a -> Vector sh a
normalizeProb = snd . normalizeFactor

normalizeFactor :: (Shape.C sh, Class.Real a) => Vector sh a -> (a, Vector sh a)
normalizeFactor xs =
   let c = Vector.sum xs
   in  (c, recip c .*| xs)

-- see htam:Stochastic
randomItemProp ::
   (Rnd.RandomGen g, Rnd.Random b, Num b, Ord b) =>
   [(a,b)] -> MS.State g a
randomItemProp props =
   let (keys,ps) = unzip props
   in  do p <- MS.state (Rnd.randomR (0, sum ps))
          return $
             fst $ head $ dropWhile ((0<=) . snd) $
             zip keys $ tail $ scanl (-) p ps

attachOnes :: (Num b) => [a] -> [(a,b)]
attachOnes = map (flip (,) 1)


vectorDim :: Shape.C sh => Vector sh a -> Int
vectorDim = Shape.size . StorableArray.shape


hermitianFromList ::
   (Shape.C sh, Class.Floating a) => sh -> [a] -> Hermitian.Hermitian sh a
hermitianFromList = Hermitian.fromList Layout.RowMajor


squareConstant ::
   (Shape.C sh, Class.Real a) => sh -> a -> Matrix.Square sh a
squareConstant =
   (ArrMatrix.fromVector .) .
   Vector.constant . Layout.square Layout.RowMajor

squareFromLists ::
   (Shape.C sh, Eq sh, Storable a) => sh -> [Vector sh a] -> Matrix.Square sh a
squareFromLists sh =
   Square.fromFull . Matrix.fromRowArray sh . Array.fromList sh

diagonal :: (Shape.C sh, Class.Real a) => Vector sh a -> Matrix.Diagonal sh a
diagonal = Matrix.diagonal Layout.RowMajor


newtype Distance f a = Distance {getDistance :: f a -> f a -> a}

distance ::
   (Shape.C sh, Eq sh, Class.Real a) =>
   Vector sh a -> Vector sh a -> a
distance =
   getDistance $
   Class.switchReal
      (Distance $ (Vector.normInf .) . Vector.sub)
      (Distance $ (Vector.normInf .) . Vector.sub)

matrixDistance ::
   (Extent.Measure meas, Extent.C vert, Extent.C horiz) =>
   (Shape.C height, Shape.C width, Eq height, Eq width, Class.Real a) =>
   ArrayMatrix pack prop lower upper meas vert horiz height width a ->
   ArrayMatrix pack prop lower upper meas vert horiz height width a ->
   a
matrixDistance a b = distance (ArrMatrix.unwrap a) (ArrMatrix.unwrap b)

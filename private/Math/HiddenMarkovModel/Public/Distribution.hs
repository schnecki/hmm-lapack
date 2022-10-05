{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE EmptyDataDecls #-}
module Math.HiddenMarkovModel.Public.Distribution (
   T(..), Trained(..), Emission,
   Show(..), NFData(..), Format(..),
   Info(..), Generate(..), EmissionProb(..),
   Estimate(..), accumulateEmissionVectors,

   Discrete, discreteFromList,
   Gaussian, gaussian, gaussianTrained,

   ToCSV(..), FromCSV(..), HMMCSV.CSVParser, CSVSymbol(..),
   ) where

import qualified Math.HiddenMarkovModel.CSV as HMMCSV
import Math.HiddenMarkovModel.Utility (randomItemProp, vectorDim)

import qualified Numeric.LAPACK.Matrix.HermitianPositiveDefinite as HermitianPD
import qualified Numeric.LAPACK.Matrix.Hermitian as Hermitian
import qualified Numeric.LAPACK.Matrix.Triangular as Triangular
import qualified Numeric.LAPACK.Matrix.Layout as Layout
import qualified Numeric.LAPACK.Matrix.Array as ArrMatrix
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Format as Format
import qualified Numeric.LAPACK.Output as Output
import Numeric.LAPACK.Matrix ((-*#), (-/#), (#/\), (|*-), (#!))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Format (FormatArray)
import Numeric.LAPACK.Output (Output)

import qualified Type.Data.Bool as TBool

import qualified Numeric.Netlib.Class as Class
import Foreign.Storable (Storable)

import qualified Data.Array.Comfort.Storable as StorableArray
import qualified Data.Array.Comfort.Shape as Shape
import qualified Data.Array.Comfort.Boxed as Array
import Data.Array.Comfort.Boxed (Array, (!))
import Data.Array.Comfort.Shape ((::+)((::+)))

import qualified System.Random as Rnd

import qualified Text.CSV.Lazy.String as CSV
import Text.Read.HT (maybeRead)
import Text.Printf (printf)

import qualified Control.Monad.Exception.Synchronous as ME
import qualified Control.Monad.Trans.Class as MT
import qualified Control.Monad.Trans.State as MS
import qualified Control.DeepSeq as DeepSeq
import Control.Monad (liftM2)
import Control.Applicative (liftA2)

import qualified Data.NonEmpty.Map as NonEmptyMap
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Semigroup as Sg
import qualified Data.Map as Map
import qualified Data.List.HT as ListHT
import qualified Data.List as List
import Data.Functor.Identity (Identity(Identity), runIdentity)
import Data.Tuple.HT (snd3)
import Data.Set (Set)
import Data.Maybe (listToMaybe)

import qualified Prelude as P
import Prelude2010 hiding (Show, showsPrec)



data family T typ sh prob
data family Trained typ sh prob

type family Emission typ prob


class Show typ where
   showsPrec ::
      (Shape.C sh, P.Show sh, P.Show prob, Storable prob) =>
      Int -> T typ sh prob -> ShowS
   showsPrecTrained ::
      (Shape.C sh, P.Show sh, P.Show prob, Storable prob) =>
      Int -> Trained typ sh prob -> ShowS

instance
   (Show typ, Shape.C sh, P.Show sh, P.Show prob, Storable prob) =>
      P.Show (T typ sh prob) where
   showsPrec = showsPrec

instance
   (Show typ, Shape.C sh, P.Show sh, P.Show prob, Storable prob) =>
      P.Show (Trained typ sh prob) where
   showsPrec = showsPrecTrained


class NFData typ where
   rnf ::
      (DeepSeq.NFData sh, DeepSeq.NFData prob, Shape.C sh) =>
      T typ sh prob -> ()
   rnfTrained ::
      (DeepSeq.NFData sh, DeepSeq.NFData prob, Shape.C sh) =>
      Trained typ sh prob -> ()

instance
   (NFData typ, DeepSeq.NFData sh, DeepSeq.NFData prob, Shape.C sh) =>
      DeepSeq.NFData (T typ sh prob) where
   rnf = rnf

instance
   (NFData typ, DeepSeq.NFData sh, DeepSeq.NFData prob, Shape.C sh) =>
      DeepSeq.NFData (Trained typ sh prob) where
   rnf = rnfTrained


class Format typ where
   format ::
      (Shape.C sh, Output out, Class.Real prob) =>
      Format.Config -> T typ sh prob -> out

instance
   (Format typ, Shape.C sh, Class.Real prob) =>
      Format.Format (T typ sh prob) where
   format = format



class Info typ where
   statesShape :: (Shape.C sh) => T typ sh prob -> sh
   statesShapeTrained :: (Shape.C sh) => Trained typ sh prob -> sh

class Generate typ where
   generate ::
      (Shape.Indexed sh, Class.Real prob, Rnd.Random prob, Rnd.RandomGen g) =>
      T typ sh prob -> Shape.Index sh -> MS.State g (Emission typ prob)

class EmissionProb typ where
   mapStatesShape ::
      (Shape.C sh0, Shape.C sh1) =>
      (sh0 -> sh1) -> T typ sh0 prob -> T typ sh1 prob
   {-
   This function could be implemented generically in terms of emissionStateProb
   but that would require an Info constraint.
   -}
   emissionProb ::
      (Shape.C sh, Class.Real prob) =>
      T typ sh prob -> Emission typ prob -> Vector sh prob
   emissionStateProb ::
      (Shape.Indexed sh, Class.Real prob) =>
      T typ sh prob -> Emission typ prob -> Shape.Index sh -> prob
   emissionStateProb distr e s = emissionProb distr e StorableArray.! s

class (EmissionProb typ) => Estimate typ where
   accumulateEmissions ::
      (Shape.Indexed sh, Class.Real prob, Shape.Index sh ~ state) =>
      sh -> NonEmpty.T [] (state, Emission typ prob) -> Trained typ sh prob
   trainVector ::
      (Shape.C sh, Eq sh, Class.Real prob) =>
      Emission typ prob -> Vector sh prob -> Trained typ sh prob
   combine ::
      (Shape.C sh, Eq sh, Class.Real prob) =>
      Trained typ sh prob -> Trained typ sh prob -> Trained typ sh prob
   normalize ::
      (Shape.C sh, Eq sh, Class.Real prob) =>
      Trained typ sh prob -> T typ sh prob

accumulateEmissionVectors ::
   (Estimate typ, Shape.C sh, Eq sh, Class.Real prob) =>
   NonEmpty.T [] (Emission typ prob, Vector sh prob) -> Trained typ sh prob
accumulateEmissionVectors = NonEmpty.foldl1Map combine (uncurry trainVector)

instance
   (Estimate typ, Shape.C sh, Eq sh, Class.Real prob) =>
      Sg.Semigroup (Trained typ sh prob) where
   (<>) = combine


data Discrete symbol

newtype instance T (Discrete symbol) sh prob =
      Discrete (Matrix.General (Set symbol) sh prob)

newtype instance Trained (Discrete symbol) sh prob =
      DiscreteTrained (NonEmptyMap.T symbol (Vector sh prob))

type instance Emission (Discrete symbol) prob = symbol


instance (P.Show symbol, Ord symbol) => Show (Discrete symbol) where
   showsPrec prec (Discrete m) = P.showsPrec prec m
   showsPrecTrained prec (DiscreteTrained m) = P.showsPrec prec m

instance (DeepSeq.NFData symbol) => NFData (Discrete symbol) where
   rnf (Discrete m) = DeepSeq.rnf m
   rnfTrained (DiscreteTrained m) = DeepSeq.rnf m

instance (P.Show symbol, Ord symbol) => Format (Discrete symbol) where
   format fmt (Discrete m) =
      Output.formatAligned $
      map (\(sym,v) ->
            map (Identity . Output.text) $
            (show sym ++ ":") :
               map (printFmt $ Format.configFormat fmt) (Vector.toList v)) $
      Array.toAssociations $ Matrix.toRowArray m

-- cf. Data.Bifunctor.Flip
newtype Flip f b a = Flip {getFlip :: f a b}

printFmt :: (Class.Real a) => String -> a -> String
printFmt fmt =
   getFlip $ Class.switchReal (Flip $ printf fmt) (Flip $ printf fmt)

instance (Ord symbol) => Info (Discrete symbol) where
   statesShape (Discrete m) = Matrix.width m
   statesShapeTrained (DiscreteTrained m) = discreteStateShape m

instance (Ord symbol) => Generate (Discrete symbol) where
   generate (Discrete m) =
      randomItemProp . StorableArray.toAssociations . Matrix.takeColumn m

instance (Ord symbol) => EmissionProb (Discrete symbol) where
   mapStatesShape f (Discrete m) = Discrete $ Matrix.mapWidth f m
   emissionProb (Discrete m) = Matrix.takeRow m
   emissionStateProb (Discrete m) x s = m #! (x,s)

instance (Ord symbol) => Estimate (Discrete symbol) where
   accumulateEmissions sh =
      DiscreteTrained .
      NonEmptyMap.map
         (StorableArray.reshape sh .
          StorableArray.fromAssociations 0 (Shape.Deferred sh) .
          Map.toList) .
      NonEmptyMap.fromListWith (Map.unionWith (+)) .
      fmap (\(state,sym) -> (sym, Map.singleton (Shape.deferIndex sh state) 1))
   trainVector sym = DiscreteTrained . NonEmptyMap.singleton sym
   combine (DiscreteTrained distr0) (DiscreteTrained distr1) =
      DiscreteTrained $ NonEmptyMap.unionWith Vector.add distr0 distr1
   normalize (DiscreteTrained distr) =
      Discrete $ normalizeProbColumns $ discreteFromMap distr

normalizeProbColumns ::
   (Shape.C height, Shape.C width, Eq width, Class.Real a) =>
   Matrix.General height width a -> Matrix.General height width a
normalizeProbColumns m = m #/\ Matrix.columnSums m

discreteStateShape ::
   (Shape.C sh) => NonEmptyMap.T symbol (Vector sh prob) -> sh
discreteStateShape =
   StorableArray.shape . snd . fst . NonEmptyMap.minViewWithKey

discreteFromMap ::
   (Ord symbol, Shape.C sh, Eq sh, Class.Real prob) =>
   NonEmptyMap.T symbol (Vector sh prob) -> Matrix.General (Set symbol) sh prob
discreteFromMap m =
   Matrix.fromRowArray (discreteStateShape m) $
   Array.fromMap $ NonEmptyMap.flatten m

discreteFromList ::
   (Ord symbol, Shape.C sh, Eq sh, Class.Real prob) =>
   NonEmpty.T [] (symbol, Vector sh prob) -> T (Discrete symbol) sh prob
discreteFromList = Discrete . discreteFromMap . NonEmptyMap.fromList



data Gaussian emiSh

newtype instance T (Gaussian emiSh) stateSh a =
   Gaussian (Array stateSh (a, Vector emiSh a, Triangular.Upper emiSh a))

newtype instance Trained (Gaussian emiSh) stateSh a =
   GaussianTrained
      (StorableArray.Array (stateSh, Layout.Hermitian (()::+emiSh)) a)

type instance Emission (Gaussian emiSh) a = Vector emiSh a


instance (Shape.C emiSh, P.Show emiSh) => Show (Gaussian emiSh) where
   showsPrec prec (Gaussian m) = P.showsPrec prec m
   showsPrecTrained prec (GaussianTrained m) = P.showsPrec prec m

instance (DeepSeq.NFData emiSh) => NFData (Gaussian emiSh) where
   rnf (Gaussian params) = DeepSeq.rnf params
   rnfTrained (GaussianTrained params) = DeepSeq.rnf params


instance (FormatArray emiSh) => Format (Gaussian emiSh) where
   format = runFormatGaussian $ Class.switchReal formatGaussian formatGaussian

newtype FormatGaussian out emiSh stateSh a =
   FormatGaussian {
      runFormatGaussian ::
         Format.Config -> T (Gaussian emiSh) stateSh a -> out
   }

formatGaussian ::
   (FormatArray emiSh, Shape.C stateSh,
    Class.Real a, Format.Format a, Output out) =>
   FormatGaussian out emiSh stateSh a
formatGaussian =
   FormatGaussian $ \fmt (Gaussian params) ->
      Format.format fmt $ Array.toList params


instance Info (Gaussian emiSh) where
   statesShape (Gaussian params) = Array.shape params
   statesShapeTrained (GaussianTrained params) =
      fst $ StorableArray.shape params

instance (Shape.C emiSh, Eq emiSh) => Generate (Gaussian emiSh) where
   generate (Gaussian allParams) state = do
      let (_c, center, covarianceChol) = allParams ! state
      seed <- MS.state Rnd.random
      return $
         Vector.add center $
         Vector.random Vector.Normal (StorableArray.shape center) seed
            -*# covarianceChol

instance (Shape.C emiSh, Eq emiSh) => EmissionProb (Gaussian emiSh) where
   mapStatesShape f (Gaussian m) = Gaussian $ Array.mapShape f m
   emissionProb (Gaussian allParams) x =
      StorableArray.fromBoxed $ fmap (gaussianEmissionProb x) allParams
   emissionStateProb (Gaussian allParams) x s =
      gaussianEmissionProb x $ allParams ! s

gaussianEmissionProb ::
   (Shape.C emiSh, Eq emiSh, Class.Real a) =>
   Vector emiSh a -> (a, Vector emiSh a, Triangular.Upper emiSh a) -> a
gaussianEmissionProb x (c, center, covarianceChol) =
   c * expSquared (Vector.sub x center -/# covarianceChol)

expSquared :: (Shape.C sh, Class.Real a) => Vector sh a -> a
expSquared =
   getNorm $ Class.switchReal (Norm expSquaredAux) (Norm expSquaredAux)

newtype Norm f a = Norm {getNorm :: f a -> a}

expSquaredAux ::
   (Shape.C sh, Class.Floating a, Vector.RealOf a ~ ar, Class.Real ar) =>
   Vector sh a -> ar
expSquaredAux x = exp ((-1/2) * Vector.norm2Squared x)


instance (Shape.C emiSh, Eq emiSh) => Estimate (Gaussian emiSh) where
   accumulateEmissions sh xs =
      let emiSh = StorableArray.shape $ snd $ NonEmpty.head xs
          hermSh = Layout.hermitian Layout.RowMajor (()::+emiSh)
      in GaussianTrained $
         Matrix.toRowMajor . Matrix.fromRowArray hermSh . Array.reshape sh .
         Array.accumulate Vector.add
            (Array.replicate (Shape.Deferred sh) (Vector.zero hermSh)) .
         map (\(state,v) -> (Shape.deferIndex sh state, extendedHermitian v)) .
         NonEmpty.flatten
            $ xs
   trainVector xs probs =
      GaussianTrained $ Matrix.toRowMajor $ probs |*- extendedHermitian xs
   combine (GaussianTrained m0) (GaussianTrained m1) =
      GaussianTrained $ Vector.add m0 m1
   {-
     Sum_i (xi-m) * (xi-m)^T
   = Sum_i xi*xi^T + Sum_i m*m^T - Sum_i xi*m^T - Sum_i m*xi^T
   = Sum_i xi*xi^T - Sum_i m*m^T
   = Sum_i xi*xi^T - n * m*m^T
   -}
   normalize (GaussianTrained m) =
      let params (weight, centerSum, covarianceSum) =
             let c = recip (weight#!((),()))
                 center = Vector.scale c $ Matrix.flattenRow centerSum
             in  (center,
                  HermitianPD.assurePositiveDefiniteness $
                  Matrix.sub
                     (Matrix.scaleRealReal c covarianceSum)
                     (Hermitian.relaxIndefinite $
                      Hermitian.outer Layout.RowMajor center))
      in Gaussian $
         fmap (gaussianParameters . params .
               Hermitian.split . ArrMatrix.fromVector) $
         Matrix.toRowArray $ Matrix.fromRowMajor m

extendedHermitian ::
   (Shape.C emiSh, Class.Floating a) =>
   StorableArray.Array emiSh a ->
   StorableArray.Array (Layout.Hermitian (()::+emiSh)) a
extendedHermitian =
   ArrMatrix.toVector .
   Hermitian.outer Layout.RowMajor . Vector.append (Vector.one ())

{- |
input array must be non-empty
-}
gaussianTrained ::
   (TBool.C zero, Shape.C emiSh, Eq emiSh, Shape.C stateSh, Class.Real prob) =>
   Array stateSh
      (prob, Vector emiSh prob,
       Matrix.FlexHermitian TBool.False zero TBool.True emiSh prob) ->
   Trained (Gaussian emiSh) stateSh prob
gaussianTrained =
   GaussianTrained . Matrix.toRowMajor .
   matrixFromRowArray "HMM.Distribution.gaussianTrained" .
   fmap
      (\(weight, center, covariance) ->
         ArrMatrix.toVector $
         Hermitian.stack
            (Hermitian.fromList Layout.RowMajor () [weight])
            (Matrix.singleRow Layout.RowMajor center)
            (Hermitian.relaxIndefinite covariance))

matrixFromRowArray ::
   (Shape.C width, Eq width, Shape.C height, Class.Real a) =>
   String ->
   Array height (StorableArray.Array width a) ->
   Matrix.General height width a
matrixFromRowArray name xs =
   case Array.toList xs of
      [] -> error $ name ++ ": empty array"
      x:_ -> Matrix.fromRowArray (StorableArray.shape x) xs

gaussian ::
   (Shape.C emiSh, Shape.C stateSh, Class.Real prob) =>
   Array stateSh (Vector emiSh prob, Matrix.HermitianPosDef emiSh prob) ->
   T (Gaussian emiSh) stateSh prob
gaussian = Gaussian . fmap gaussianParameters

gaussianParameters ::
   (Shape.C emiSh, Class.Real prob) =>
   (Vector emiSh prob, Matrix.HermitianPosDef emiSh prob) ->
   (prob, Vector emiSh prob, Triangular.Upper emiSh prob)
gaussianParameters (center, covariance) =
   gaussianFromCholesky center $ HermitianPD.decompose covariance

gaussianFromCholesky ::
   (Shape.C emiSh, Class.Real prob) =>
   Vector emiSh prob -> Triangular.Upper emiSh prob ->
   (prob, Vector emiSh prob, Triangular.Upper emiSh prob)
gaussianFromCholesky center covarianceChol =
   let covarianceSqrtDet =
         Vector.product $ Triangular.takeDiagonal covarianceChol
   in  (recip (sqrt2pi ^ vectorDim center * covarianceSqrtDet),
        center, covarianceChol)

sqrt2pi :: (Class.Real a) => a
sqrt2pi = runIdentity $ Class.switchReal sqrt2piAux sqrt2piAux

sqrt2piAux :: (Floating a) => Identity a
sqrt2piAux = Identity $ sqrt (2*pi)


class ToCSV typ where
   toCells ::
      (Shape.C sh, Class.Real prob, P.Show prob) =>
      T typ sh prob -> [[String]]

class FromCSV typ where
   parseCells ::
      (Shape.C sh, Eq sh, Class.Real prob, Read prob) =>
      sh -> HMMCSV.CSVParser (T typ sh prob)

class (Ord symbol) => CSVSymbol symbol where
   cellFromSymbol :: symbol -> String
   symbolFromCell :: String -> Maybe symbol

instance CSVSymbol Char where
   cellFromSymbol = (:[])
   symbolFromCell = listToMaybe

instance CSVSymbol Int where
   cellFromSymbol = show
   symbolFromCell = maybeRead


instance (CSVSymbol symbol) => ToCSV (Discrete symbol) where
   toCells (Discrete m) =
      map
         (\(symbol, probs) ->
            cellFromSymbol symbol : HMMCSV.cellsFromVector probs) $
      Array.toAssociations $ Matrix.toRowArray m

instance (CSVSymbol symbol) => FromCSV (Discrete symbol) where
   parseCells n =
      let p = parseSymbolProb n
      in fmap discreteFromList $
         liftA2 NonEmpty.Cons (HMMCSV.getRow >>= p) (HMMCSV.manyRowsUntilEnd p)

parseSymbolProb ::
   (Shape.C sh, Class.Real prob, Read prob, CSVSymbol symbol) =>
   sh -> CSV.CSVRow -> HMMCSV.CSVParser (symbol, Vector sh prob)
parseSymbolProb sh row =
   case row of
      [] -> MT.lift $ ME.throw "missing symbol"
      c:cs ->
         liftM2 (,)
            (let str = CSV.csvFieldContent c
             in  MT.lift $ ME.fromMaybe (printf "unknown symbol %s" str) $
                 symbolFromCell str)
            (do v <- HMMCSV.parseVectorFields cs
                let n = Shape.size sh
                let m = vectorDim v
                HMMCSV.assert (n == m)
                   (printf "number of states (%d) and size of probability vector (%d) mismatch"
                      n m)
                return $ StorableArray.reshape sh v)


instance (Shape.Indexed emiSh) => ToCSV (Gaussian emiSh) where
   toCells (Gaussian params) =
      List.intercalate [[]] $
      map
         (\(_, center, covarianceChol) ->
            HMMCSV.cellsFromVector center :
            HMMCSV.cellsFromSquare (Triangular.toSquare covarianceChol)) $
      Array.toList params

instance (emiSh ~ Matrix.ShapeInt) => FromCSV (Gaussian emiSh) where
   parseCells sh = do
      let n = Shape.size sh
      gs <- HMMCSV.manySepUntilEnd parseSingleGaussian
      HMMCSV.assert (length gs == n) $
         printf "number of states (%d) and number of Gaussians (%d) mismatch"
            n (length gs)
      let sizes = map (vectorDim . snd3) gs
      HMMCSV.assert (ListHT.allEqual sizes) $
         printf "dimensions of emissions mismatch: %s" (show sizes)
      return $ Gaussian $ Array.fromList sh gs

parseSingleGaussian ::
   (emiSh ~ Matrix.ShapeInt, Class.Real prob, Eq prob, Read prob) =>
   HMMCSV.CSVParser (prob, Vector emiSh prob, Triangular.Upper emiSh prob)
parseSingleGaussian = do
   center <- HMMCSV.parseNonEmptyVectorCells
   covarianceCholSquare <-
      HMMCSV.parseSquareMatrixCells $ StorableArray.shape center
   let covarianceChol = Triangular.takeUpper covarianceCholSquare
   HMMCSV.assert
      (isUpperTriang covarianceCholSquare covarianceChol)
      "matrices must be upper triangular"
   return $ gaussianFromCholesky center covarianceChol


{-
Maybe this test is too strict.
It would also be ok, and certainly more intuitive
to use an orthogonal but not normalized matrix.
We could get such a matrix from the eigensystem.
-}
isUpperTriang ::
   (Shape.C sh, Class.Real a, Eq a) =>
   Matrix.Square sh a -> Triangular.Upper sh a -> Bool
isUpperTriang m mt =
   Vector.toList (ArrMatrix.toVector m)
   ==
   Vector.toList (ArrMatrix.toVector (Triangular.toSquare mt))

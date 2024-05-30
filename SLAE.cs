namespace VectorFEM3D;
public abstract class SLAE
{
   protected SparseMatrix matrix = default!;
   protected Vector vector = default!;
   public Vector solution = default!;
   public double time;
   protected double eps;
   protected int maxIters;
   public int lastIter;

   public SLAE() 
   {
      eps = 1e-16;
      maxIters = 2000;
   }

   public SLAE(double eps, int maxIters)
   {
      this.eps = eps;
      this.maxIters = maxIters;
   }

   public void SetSLAE(Vector vector, SparseMatrix matrix, Vector solution)
   {
      this.vector = vector;
      this.matrix = matrix;
      this.solution = solution;
   }

   public abstract Vector Solve();


   protected void LU()
   {
      for (int i = 0; i < matrix.Size; i++)
      {

         for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
         {
            int jCol = matrix.Jg[j];
            int jk = matrix.Ig[jCol];
            int k = matrix.Ig[i];

            int sdvig = matrix.Jg[matrix.Ig[i]] - matrix.Jg[matrix.Ig[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < matrix.Ig[jCol + 1]; k++, jk++)
            {
               sumL += matrix.Ggl[k] * matrix.Ggu[jk];
               sumU += matrix.Ggu[k] * matrix.Ggl[jk];
            }

            matrix.Ggl[j] -= sumL;
            matrix.Ggu[j] -= sumU;
            matrix.Ggu[j] /= matrix.Di[jCol];
         }

         double sumD = 0.0;
         for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
            sumD += matrix.Ggl[j] * matrix.Ggu[j];

         matrix.Di[i] -= sumD;
      }
   }

   protected void ForwardElimination()
   {
      for (int i = 0; i < matrix.Size; i++)
      {
         for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
         {
            solution[i] -= matrix.Ggl[j] * solution[matrix.Jg[j]];
         }

         solution[i] /= matrix.Di[i];
      }
   }

   protected void BackwardSubstitution()
   {
      for (int i = matrix.Size - 1; i >= 0; i--)
      {
         for (int j = matrix.Ig[i + 1] - 1; j >= matrix.Ig[i]; j--)
         {
            solution[matrix.Jg[j]] -= matrix.Ggu[j] * solution[i];
         }
      }
   }

   public void PrintSolution()
   {
      for(int i = 0; i < solution.Length; i++)
      {
         Console.WriteLine(solution[i]);
      }
   }
}

public class LOSLTSolver : SLAE
{
   public LOSLTSolver(double eps, int maxIters) : base(eps, maxIters) { }
   
   public override Vector Solve()
   {
      SparseMatrix ltMatrix = new SparseMatrix(matrix.Size, matrix.Jg.Length, matrix.Symmetric);
      SparseMatrix.Copy(matrix, ltMatrix);
      
      ConvertLT(ltMatrix);
      
      int i;
      
      Vector r = DirElim(vector - matrix * solution, ltMatrix);
      Vector z = BackSub(r, ltMatrix);
      Vector p = DirElim(matrix * z, ltMatrix);
      
      double error = r * r;

      for (i = 1; i <= maxIters && error > eps; i++)
      {
         var alpha = p * r / (p * p);
         solution += alpha * z;
         error = (r * r) - (alpha * alpha * (p * p));
         r -= alpha * p;

         var tmp = DirElim(matrix * BackSub(r, ltMatrix), ltMatrix);

         var beta = -(p * tmp) / (p * p);
         z = BackSub(r, ltMatrix) + (beta * z);
         p = tmp + (beta * p);
         
         //Console.WriteLine($"{i}: {error}");
      }

      Console.WriteLine($"Iters: {i}\tError: {error}");
      return solution;
   }

   private void ConvertLT(SparseMatrix _matrix)
   {
      for (int i = 0; i < _matrix.Size; i++)
      {
         for (int j = _matrix.Ig[i]; j < _matrix.Ig[i + 1]; j++)
         {
            int jCol = _matrix.Jg[j];
            int jk = _matrix.Ig[jCol];
            int k = _matrix.Ig[i];

            double sumL = 0.0;

            for (; k < j && jk < _matrix.Ig[jCol + 1]; )
            {
               if (_matrix.Jg[k] == _matrix.Jg[jk])
               {
                  sumL += _matrix.Gg[k] * _matrix.Gg[jk];
                  k++;
                  jk++;
               }

               else
               {
                  if (_matrix.Jg[k] > _matrix.Jg[jk])
                  {
                     jk++;
                  }

                  else
                  {
                     k++;
                  }
               }
            }

            _matrix.Gg[j] -= sumL;
            _matrix.Gg[j] /= _matrix.Di[jCol];
         }

         double sumD = 0.0;
         for (int j = _matrix.Ig[i]; j < _matrix.Ig[i + 1]; j++)
            sumD += _matrix.Gg[j] * _matrix.Gg[j];

         _matrix.Di[i] -= sumD;
         _matrix.Di[i] = Math.Sqrt(Math.Abs(_matrix.Di[i]));
      } 
   }

   private Vector DirElim(Vector vector, SparseMatrix lltMatrix)
   {
      Vector result = new Vector(vector.Length);
      Vector.Copy(vector, result);

      for (int i = 0; i < vector.Length; i++)
      {
         for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
         {
            result[i] -= lltMatrix.Gg[j] * result[matrix.Jg[j]];
         }

         result[i] /= lltMatrix.Di[i];
      }

      return result;
   }

   private Vector BackSub(Vector vector, SparseMatrix lltMatrix)
   {
      Vector result = new Vector(vector.Length);
      Vector.Copy(vector, result);

      for (int i = matrix.Size - 1; i >= 0; i--)
      {
         result[i] /= lltMatrix.Di[i];
         for (int j = matrix.Ig[i + 1] - 1; j >= matrix.Ig[i]; j--)
         {
            result[matrix.Jg[j]] -= lltMatrix.Gg[j] * result[i];
         }
      }

      return result;
   }
}

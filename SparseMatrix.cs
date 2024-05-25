namespace VectorFEM3D;

public class SparseMatrix
{
    public int[] Ig { get; }
    public int[] Jg { get; }
    public double[] Di { get; }
    public double[] Gg { get; }
    public int Size { get; }
    public double[] Ggl { get; }
    public double[] Ggu { get; }
    public bool Symmetric { get; }

    public SparseMatrix(int size, int sizeOffDiag, bool flag)
    {
        Size = size;
        Ig = new int[size + 1];
        Jg = new int[sizeOffDiag];
        Gg = new double[sizeOffDiag];
        Di = new double[size];
        Ggl = new double[sizeOffDiag];
        Ggu = new double[sizeOffDiag];
        Symmetric = flag;
    }

    public static Vector operator *(SparseMatrix matrix, Vector vector)
    {
        Vector product = new(vector.Length);

        if (matrix.Symmetric)
        {
            for (int i = 0; i < vector.Length; i++)
            {
                product[i] = matrix.Di[i] * vector[i];

                for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
                {
                    product[i] += matrix.Gg[j] * vector[matrix.Jg[j]];
                    product[matrix.Jg[j]] += matrix.Gg[j] * vector[i];
                }
            }
        }

        else
        {
            for (int i = 0; i < vector.Length; i++)
            {
                product[i] += matrix.Di[i] * vector[i];

                for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
                {
                    product[i] += matrix.Ggl[j] * vector[matrix.Jg[j]];
                    product[matrix.Jg[j]] += matrix.Ggu[j] * vector[i];
                }
            }
        }

        return product;
    }

    public void PrintDense(string path)
    {
        double[,] A = new double[Size, Size];

        for (int i = 0; i < Size; i++)
        {
            A[i, i] = Di[i];

            for (int j = Ig[i]; j < Ig[i + 1]; j++)
            {
                A[i, Jg[j]] = Gg[j];
                A[Jg[j], i] = Gg[j];
            }
        }

        using var sw = new StreamWriter(path);
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size; j++)
            {
                sw.Write(A[i, j].ToString("0.00") + "\t");
            }

            sw.WriteLine();
        }
    }

    public void Clear()
    {
        for (int i = 0; i < Size; i++)
        {
            Di[i] = 0.0;

            for (int k = Ig[i]; k < Ig[i + 1]; k++)
            {
                Gg[k] = 0.0;
            }
        }
    }
    
    public static void Copy(SparseMatrix source, SparseMatrix product)
    {
        for (int i = 0; i < product.Size + 1; i++)
        {
            product.Ig[i] = source.Ig[i];
        }

        for (int i = 0; i < product.Size; i++)
        {
            product.Di[i] = source.Di[i];
        }

        for (int i = 0; i < product.Jg.Length; i++)
        {
            product.Jg[i] = source.Jg[i];
            product.Ggl[i] = source.Ggl[i];
            product.Ggu[i] = source.Ggu[i];
            product.Gg[i] = source.Gg[i];
        }
    }
}

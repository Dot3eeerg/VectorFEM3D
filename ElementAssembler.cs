namespace VectorFEM3D;

public record ElementAssembler
{
    public Matrix StiffnessMatrix;
    public Matrix MassMatrix;
    public Vector LocalVector;
    
    public Parallelepiped Element;
    
    private static Integration _integration;
    private static IBasis3D _basis;

    private readonly double[][] _g1 =
    {
        new double[] { 2, 1, -2, -1 },
        new double[] { 1, 2, -1, -2 },
        new double[] { -2, -1, 2, 1 },
        new double[] { -1, -2, 1, 2 }
    };

    private double[][] _g2 =
    {
        new double[] { 2, -2, 1, -1 },
        new double[] { -2, 2, -1, 1 },
        new double[] { 1, -1, 2, -2 },
        new double[] { -1, 1, -2, 2 }
    };

    private double[][] _g3 =
    {
        new double[] { -2, 2, -1, 1 },
        new double[] { -1, 1, -2, 2 },
        new double[] { 2, -2, 1, -1 },
        new double[] { 1, -1, 2, -2 }
    };

    private double[][] _g3T =
    {
        new double[] { -2, -1, 2, 1 },
        new double[] { 2, 1, -2, -1 },
        new double[] { -1, -2, 1, 2 },
        new double[] { 1, 2, -1, -2 }
    };

    private double[][] _m =
    {
        new double[] { 4, 2, 2, 1 },
        new double[] { 2, 4, 1, 2 },
        new double[] { 2, 1, 4, 2 },
        new double[] { 1, 2, 2, 4 }
    };
    
    public ElementAssembler(IBasis3D basis)
    {
        _basis = basis;
        _integration = new Integration(new SegmentGaussOrder9());

        StiffnessMatrix = new Matrix(_basis.Size);
        MassMatrix = new Matrix(_basis.Size);
        LocalVector = new Vector(_basis.Size);

        Element = new Parallelepiped();
    }

    public void AnalyticalAssembly(double time, Func<Point3D, double, int, double, double> f)
    {
        // Mass matrix
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                MassMatrix[i, j] = Element.hx * Element.hy * Element.hz / 36 * _m[i % 4][j % 4];
            }
        }
        
        for (int i = 4; i < 8; i++)
        {
            for (int j = 4; j < 8; j++)
            {
                MassMatrix[i, j] = Element.hx * Element.hy * Element.hz / 36 * _m[i % 4][j % 4];
            }
        }
        
        for (int i = 8; i < 12; i++)
        {
            for (int j = 8; j < 12; j++)
            {
                MassMatrix[i, j] = Element.hx * Element.hy * Element.hz / 36 * _m[i % 4][j % 4];
            }
        }
        
        // Stiffness matrix
        // G1 + G2
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                StiffnessMatrix[i, j] = Element.hx * Element.hy / (6 * Element.hz) * _g1[i % 4][j % 4] +
                                        Element.hx * Element.hz / (6 * Element.hy) * _g2[i % 4][j % 4];
            }
        }

        // G2
        for (int i = 4; i < 8; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                StiffnessMatrix[i, j] = StiffnessMatrix[j, i] = -Element.hz / 6 * _g2[i % 4][j % 4];
            }
        }
        
        // G3
        for (int i = 8; i < 12; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                StiffnessMatrix[i, j] = Element.hy / 6 * _g3T[i % 4][j % 4];
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 8; j < 12; j++)
            {
                StiffnessMatrix[i, j] = Element.hy / 6 * _g3[i % 4][j % 4];
            }
        }
        
        // G1 + G2
        for (int i = 4; i < 8; i++)
        {
            for (int j = 4; j < 8; j++)
            {
                StiffnessMatrix[i, j] = Element.hx * Element.hy / (6 * Element.hz) * _g1[i % 4][j % 4] +
                                        Element.hy * Element.hz / (6 * Element.hx) * _g2[i % 4][j % 4];
            }
        }
        
        // G1
        for (int i = 8; i < 12; i++)
        {
            for (int j = 4; j < 8; j++)
            {
                StiffnessMatrix[i, j] = StiffnessMatrix[j, i] = -Element.hx / 6 * _g1[i % 4][j % 4];
            }
        }
        
        // G1 + G2
        for (int i = 8; i < 12; i++)
        {
            for (int j = 8; j < 12; j++)
            {
                StiffnessMatrix[i, j] = Element.hx * Element.hz / (6 * Element.hy) * _g1[i % 4][j % 4] +
                                        Element.hy * Element.hz / (6 * Element.hx) * _g2[i % 4][j % 4];
            }
        }

        for (int i = 0; i < _basis.Size; i++)
        {
            //Func<Point3D, double> kek;
            //Vector3D psi1 = new(0, 0, 0);
            //Vector3D func = new(0, 0, 0);

            //var ik = i;
            //kek = point =>
            //{
            //    psi1.Copy(_basis.GetPsi(ik, point));
            //    switch (ik / 4)
            //    {
            //        case 0:
            //            func.UpdateVector(f(Element.PointList[ik], time, ik, Element.Sigma), 0, 0);
            //            break;

            //        case 1:
            //            func.UpdateVector(0, f(Element.PointList[ik], time, ik, Element.Sigma), 0);
            //            break;

            //        case 2:
            //            func.UpdateVector(0, 0, f(Element.PointList[ik], time, ik, Element.Sigma));
            //            break;
            //    }

            //    return func * psi1;
            //};
            //LocalVector[i] = Element.hx * Element.hy * Element.hz * _integration.Gauss3D(kek);
            
            LocalVector[i] = f(Element.PointList[i], time, i, Element.Sigma);
        }

        LocalVector = MassMatrix * LocalVector;
    }
}

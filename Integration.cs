namespace VectorFEM3D;

public class Integration
{
    private readonly SegmentGaussOrder9 _quadratures;

    public Integration(SegmentGaussOrder9 quadratures)
    {
        _quadratures = quadratures;
    }

    public double Gauss3D(Func<Point3D, double> psi)
    {
        double result = 0;
        Point3D point = new(0, 0, 0);

        for (int i = 0; i < _quadratures.Size; i++)
        {
            point.X = (_quadratures.GetPoint(i) + 1) / 2.0;

            for (int j = 0; j < _quadratures.Size; j++)
            {
                point.Y = (_quadratures.GetPoint(j) + 1) / 2.0;

                for (int k = 0; k < _quadratures.Size; k++)
                {
                    point.Z = (_quadratures.GetPoint(k) + 1) / 2.0;

                    result += psi(point) * _quadratures.GetWeight(i) * _quadratures.GetWeight(j) *
                              _quadratures.GetWeight(k);
                }
            }
        }

        return result / 8.0;
    }
}
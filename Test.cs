using System.Runtime.CompilerServices;

namespace VectorFEM3D;

public abstract class Test
{
    protected double lambda;
    protected double sigma;
    protected double epsilon;

    public Test(Grid grid)
    {
        lambda = grid.Mu;
        sigma = grid.Sigma;
        epsilon = grid.Epsilon;
    }
    
    public abstract double UValue(Point3D point, double t, int i);
    protected abstract double FValue(Point3D point, double t, int i, double sigma);

    public abstract double Theta(Point3D point, double t, ElementSide elementSide);
    
    public double F(Point3D point, double t, int i, double sigma)
        => i switch
        {
            0 => FValue(point, t, 0, sigma),
            1 => FValue(point, t, 0, sigma),
            2 => FValue(point, t, 1, sigma),
            3 => FValue(point, t, 1, sigma),
            4 => FValue(point, t, 2, sigma),
            5 => FValue(point, t, 2, sigma),
            6 => FValue(point, t, 2, sigma),
            7 => FValue(point, t, 2, sigma),
            8 => FValue(point, t, 0, sigma),
            9 => FValue(point, t, 0, sigma),
            10 => FValue(point, t, 1, sigma),
            11 => FValue(point, t, 1, sigma),
        };
}

public class Test1 : Test
{
    public Test1(Grid grid) : base(grid) { }
    
    public override double UValue(Point3D point, double t, int i)
        => i switch
        {
            0 => t,
            1 => t,
            2 => t,
            _ => throw new Exception("Can't find UValue type")
        };

    protected override double FValue(Point3D point, double t, int i, double sigma)
        => i switch
        {
            0 => sigma * 1,
            1 => sigma * 1,
            2 => sigma * 1,
            _ => throw new Exception("Can't find FValue type")
        };

    public override double Theta(Point3D point, double t, ElementSide elementSide)
    {
        switch (elementSide)
        {
            case ElementSide.Left:
            case ElementSide.Right:
                return 0;
            
            case ElementSide.Bottom:
            case ElementSide.Upper:
                return 2;
            
            case ElementSide.Rear:
            case ElementSide.Front:
                return 1;
            
            default:
                throw new Exception("Can't find function type");
        }
    }
}
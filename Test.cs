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
    protected abstract double FValue(Point3D point, double t, int i);

    public abstract double Theta(Point3D point, double t, ElementSide elementSide);
    
    public double F(Point3D point, double t, int i)
        => i switch
        {
            0 => FValue(point, t, 0),
            1 => FValue(point, t, 0),
            2 => FValue(point, t, 1),
            3 => FValue(point, t, 1),
            4 => FValue(point, t, 2),
            5 => FValue(point, t, 2),
            6 => FValue(point, t, 2),
            7 => FValue(point, t, 2),
            8 => FValue(point, t, 0),
            9 => FValue(point, t, 0),
            10 => FValue(point, t, 1),
            11 => FValue(point, t, 1),
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

    protected override double FValue(Point3D point, double t, int i)
        => i switch
        {
            0 => 1,
            1 => 1,
            2 => 1,
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
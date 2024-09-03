namespace VectorFEM3D;

public class Grid2D
{
    private readonly List<double> _xStart = new List<double>();
    private readonly List<double> _xEnd = new List<double>();
    private readonly List<int> _xSteps = new List<int>();
    private readonly List<double> _xRaz = new List<double>();
    private readonly List<double> _yStart = new List<double>();
    private readonly List<double> _yEnd = new List<double>();
    private readonly List<int> _ySteps = new List<int>();
    private readonly List<double> _yRaz = new List<double>();
    private readonly int[] _boundaries;

    private readonly double[] _xZones;
    private readonly double[] _yZones;
    private readonly int[][] _zones;
    private readonly double[] _sigmaValues;
    private readonly int[][] _mainZones;
    private readonly double[] _mainSigmaValues;
    
    public double GeneratorX;
    public double GeneratorY;
    private double _receiver;
    public double SourceValue;

    private readonly List<double> _xValues = new List<double>();
    private readonly List<double> _yValues = new List<double>();

    public List<int> SumSteps;
    
    public Point2D[] Nodes { get; private set; }
    public HashSet<int> DirichletBoundaries { get; private set; } 
    public int[][] Elements { get; private set; }
    public double Mu { get; set; }
    public HashSet<int> LeftSideElement = new HashSet<int>();

    public int Generator { get; private set; }
    public int Receiver { get; private set; }

    public Grid2D(string path)
    {
        using (var sr = new StreamReader(path))
        {
            string[] data;
            data = sr.ReadLine()!.Split(" ").ToArray();
            int kek = Convert.ToInt32(data[0]);
            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _xStart.Add(Convert.ToDouble(data[0]));
                _xEnd.Add(Convert.ToDouble(data[1]));
                _xSteps.Add(Convert.ToInt32(data[2]));
                _xRaz.Add(Convert.ToDouble(data[3]));
            }

            _xZones = sr.ReadLine()!.Split(" ").Select(Convert.ToDouble).ToArray();
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            kek = Convert.ToInt32(data[0]);
            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _yStart.Add(Convert.ToDouble(data[0]));
                _yEnd.Add(Convert.ToDouble(data[1]));
                _ySteps.Add(Convert.ToInt32(data[2]));
                _yRaz.Add(Math.Sqrt(Convert.ToDouble(data[3])));
            }
            
            _yZones = sr.ReadLine()!.Split(" ").Select(Convert.ToDouble).ToArray();
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            _boundaries = new int[4];
            _boundaries[0] = Convert.ToInt32(data[0]);
            _boundaries[1] = Convert.ToInt32(data[1]);
            _boundaries[2] = Convert.ToInt32(data[2]);
            _boundaries[3] = Convert.ToInt32(data[3]);
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            Mu = Convert.ToDouble(data[0]);

            kek = Convert.ToInt32(sr.ReadLine());
            _zones = new int[kek].Select(_ => new int[4]).ToArray();
            _sigmaValues = new double[kek];

            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _sigmaValues[i] = Convert.ToDouble(data[0]);
                _zones[i][0] = Convert.ToInt32(data[1]);
                _zones[i][1] = Convert.ToInt32(data[2]);
                _zones[i][2] = Convert.ToInt32(data[3]);
                _zones[i][3] = Convert.ToInt32(data[4]);
            }

            data = sr.ReadLine()!.Split(" ").ToArray();
            GeneratorX = Convert.ToDouble(data[0]);
            GeneratorY = Convert.ToDouble(data[1]);
            SourceValue = Convert.ToDouble(data[2]);
            
            data = sr.ReadLine()!.Split(" ").ToArray();
            _receiver = Convert.ToDouble(data[0]);

            kek = Convert.ToInt32(sr.ReadLine());
            _mainZones = new int[kek].Select(_ => new int[4]).ToArray();
            _mainSigmaValues = new double[kek];

            for (int i = 0; i < kek; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                _mainSigmaValues[i] = Convert.ToDouble(data[0]);
                _mainZones[i][0] = Convert.ToInt32(data[1]);
                _mainZones[i][1] = Convert.ToInt32(data[2]);
                _mainZones[i][2] = Convert.ToInt32(data[3]);
                _mainZones[i][3] = Convert.ToInt32(data[4]);
            }
        }
    }

    public void BuildGrid(bool flag)
    {
        SumSteps = new List<int>();
        int count = 1;
        int nodeCount = 1;
        SumSteps.Add(0);
        for (int i = 0; i < _xSteps.Count; i++)
        {
            SumSteps[0] += _xSteps[i];
        }
        count *= SumSteps[0];
        nodeCount *= SumSteps[0] + 1;
        int nodesInRow = SumSteps[0] + 1;
        
        SumSteps.Add(0);
        for (int i = 0; i < _ySteps.Count; i++)
        {
            SumSteps[1] += _ySteps[i];
        }
        count *= SumSteps[1];
        nodeCount *= SumSteps[1] + 1;
        
        Elements = new int[count].Select(_ => new int[4]).ToArray();
        Nodes = new Point2D[nodeCount];

        List<double> sumRazX = new(), sumRazY = new();
        for (int i = 0; i < _xSteps.Count; i++)
        {
            sumRazX.Add(0);
            for (int j = 0; j < _xSteps[i]; j++)
            {
                sumRazX[i] += Math.Pow(_xRaz[i], j);
            }
        }
        
        for (int i = 0; i < _ySteps.Count; i++)
        {
            sumRazY.Add(0);
            for (int j = 0; j < _ySteps[i]; j++)
            {
                sumRazY[i] += Math.Pow(_yRaz[i], j);
            }
        }
        
        DirichletBoundaries = new();

        double x = 0, y = 0;
        double xStep = 0, yStep = 0;

        for (int i = 0; i < _xSteps.Count; i++)
        {
            x = _xStart[i];
            xStep = (_xEnd[i] - _xStart[i]) / sumRazX[i];
            
            for (int j = 0; j < _xSteps[i]; j++)
            {
                _xValues.Add(x);
                x += xStep;
                xStep *= _xRaz[i];
            }
            
        }
        _xValues.Add(_xEnd[^1]);

        if (flag)
        {
            for (int i = 0; i < _xValues.Count - 1; i++)
            {
                if (_xValues[i] == _receiver)
                {
                    break;
                }
                
                if (_xValues[i] < _receiver && _xValues[i + 1] > _receiver)
                {
                    if (Math.Abs(_receiver - _xValues[i]) > Math.Abs(_xValues[i + 1] - _receiver))
                    {
                        _xValues[i + 1] = _receiver;
                    }
                    else
                    {
                        _xValues[i] = _receiver;
                    }
                    
                    break;
                }
            }
        }

        if (flag)
        {
            for (int i = 0; i < _xValues.Count - 1; i++)
            {
                if (_xValues[i] == GeneratorX)
                {
                    break;
                }
                
                if (_xValues[i] < GeneratorX && _xValues[i + 1] > GeneratorX)
                {
                    _xValues[i + 1] = GeneratorX;
                    
                    break;
                }
            }
        }
        
        for (int i = 0; i < _ySteps.Count; i++)
        {
            y = _yStart[i];
            yStep = (_yEnd[i] - _yStart[i]) / sumRazY[i];
            
            for (int j = 0; j < _ySteps[i]; j++)
            {
                _yValues.Add(y);
                y += yStep;
                yStep *= _yRaz[i];
            }
            
        }
        _yValues.Add(_yEnd[^1]);

        if (flag)
        {
            for (int i = 0; i < _yValues.Count - 1; i++)
            {
                if (_yValues[i] == GeneratorY)
                {
                    break;
                }
                
                if (_yValues[i] < GeneratorY && _yValues[i + 1] > GeneratorY)
                {
                    _yValues[i + 1] = GeneratorY;
                    
                    break;
                }
            }
        }
        
        for (int j = 0; j < _yValues.Count; j++)
        {
            for (int k = 0; k < _xValues.Count; k++)
            {
                Nodes[j * nodesInRow + k] = new Point2D(_xValues[k], _yValues[j]);

                if (_xValues[k] == GeneratorX && _yValues[j] == GeneratorY)
                {
                    Generator = j * nodesInRow + k;
                }

                if (_xValues[k] == _receiver && _yValues[j] == GeneratorY)
                {
                    Receiver = j * nodesInRow + k;
                }
            }
        }

        int index = 0;
        
        index = 0;

        for (int i = 0; i < SumSteps[1]; i++)
        {
            for (int j = 0; j < SumSteps[0]; j++)
            {
                Elements[index][0] = i * nodesInRow + j;
                Elements[index][1] = i * nodesInRow + j + 1;
                Elements[index][2] = (i + 1) * nodesInRow + j;
                Elements[index++][3] = (i + 1) * nodesInRow + j + 1;
            }
        }
    }

    public void AmletGrid()
    {
        using (var sw = new StreamWriter("rValues.txt"))
        {
            foreach (var kek in _xValues)
            {
               sw.WriteLine(kek);
            }
        }
        
        using (var sw = new StreamWriter("zValues.txt"))
        {
            foreach (var kek in _yValues)
            {
               sw.WriteLine(kek);
            }
        }

        using (var sw = new StreamWriter("rz"))
        {
            foreach (var yValue in _yValues)
            {
                foreach (var xValue in _xValues)
                {
                    sw.WriteLine($"{xValue} {yValue}");
                }
            }
        }
        
        using (var sw = new StreamWriter("nvkat2d"))
        {
            for (int i = 0; i < Elements.Length; i++)
            {
                sw.WriteLine(Elements[i][0] == Generator ? "2" : "1");
            }
        }
        
        using (var sw = new StreamWriter("results.txt"))
        {
            foreach (var node in Nodes)
            {
                sw.WriteLine($"{node.X} {node.Y}");
            }
        }

        using (var sw = new StreamWriter("nvtr"))
        {
            for (int i = 0; i < Elements.Length; i++)
            {
                sw.WriteLine($"{Elements[i][0] + 1} {Elements[i][1] + 1} {Elements[i][2] + 1} {Elements[i][3] + 1} 0 1");
            }
        }

        using (var sw = new StreamWriter("l1"))
        {
            for (int i = 0; i < Nodes.Length; i++)
            {
                if (DirichletBoundaries.Contains(i))
                {
                    sw.WriteLine($"{i + 1}");
                }
            }
        }
    }

    public double GetSigma(Point2D point)
    {
        for (int i = 0; i < _zones.Length; i++)
        {
            if (point.X <= _xZones[_zones[i][1]] && point.Y <= _yZones[_zones[i][3]] &&
                point.X >= _xZones[_zones[i][0]] && point.Y >= _yZones[_zones[i][2]])
            {
                return _sigmaValues[i];
            }
        }

        throw new Exception("Can't find eligible zone for sigma");
    }

    public double GetMainSigma(Point2D point)
    {
        for (int i = 0; i < _mainZones.Length; i++)
        {
            if (point.X <= _xZones[_mainZones[i][1]] && point.Y <= _yZones[_mainZones[i][3]] &&
                point.X >= _xZones[_mainZones[i][0]] && point.Y >= _yZones[_mainZones[i][2]])
            {
                return _mainSigmaValues[i];
            }
        }

        throw new Exception("Can't find eligible zone for main sigma");
    }

    public double GetDifferenceSigma(Point2D point)
    {
        return -(GetMainSigma(point) - GetSigma(point));
    }

    public bool CheckZone(Point2D point)
    {
        if (GetMainSigma(point) - GetSigma(point) != 0)
        {
            return true;
        }

        return false;
    }

    public void AccountBoundaryConditions()
    {
        for (int ielem = 0; ielem < Elements.Length; ielem++)
        {
            if (ielem < SumSteps[0])
            {
                if (_boundaries[2] == 1) DirichletBoundary(ElementSide.Bottom, ielem);
            }

            if (ielem >= SumSteps[0] * (SumSteps[1] - 1))
            {
                if (_boundaries[3] == 1) DirichletBoundary(ElementSide.Upper, ielem);
            }

            if (ielem % SumSteps[0] == 0)
            {
                if (_boundaries[0] == 1) DirichletBoundary(ElementSide.Left, ielem);
            }

            if ((ielem + 1) % SumSteps[0] == 0)
            {
                if (_boundaries[1] == 1) DirichletBoundary(ElementSide.Right, ielem);
            }
        }
    }

    private void DirichletBoundary(ElementSide elementSide, int ielem)
    {
        switch (elementSide)
        {
            case ElementSide.Bottom:
                DirichletBoundaries.Add(Elements[ielem][0]);
                DirichletBoundaries.Add(Elements[ielem][1]);
                break;
            
            case ElementSide.Upper:
                DirichletBoundaries.Add(Elements[ielem][2]);
                DirichletBoundaries.Add(Elements[ielem][3]);
                break;
            
            case ElementSide.Left:
                DirichletBoundaries.Add(Elements[ielem][0]);
                DirichletBoundaries.Add(Elements[ielem][2]);
                LeftSideElement.Add(ielem);
                break;
                
            case ElementSide.Right:
                DirichletBoundaries.Add(Elements[ielem][1]);
                DirichletBoundaries.Add(Elements[ielem][3]);
                break;
        }
    }
}

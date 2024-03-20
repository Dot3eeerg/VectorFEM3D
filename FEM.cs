namespace VectorFEM3D;

public class FEM
{
    private SparseMatrix? _globalMatrix;
    private Vector? _globalVector;
    private Vector[]? _layers;
    private Vector? _solution;
    private Vector? _localVector;
    private Matrix? _stiffnessMatrix;
    private Matrix? _massMatrix;
    private Grid? _grid;
    private ITimeGrid _timeGrid;
    private Test? _test;
    private IBasis3D? _basis;
    private Integration? _integration;
    private SLAE? _slae;
    private Scheme _scheme;
    private Scheme _activeScheme;
    //private Generator _generator = new Generator(-100, -100, 25, 100, 100, 25);
    private Generator _generator = new Generator(-50, -50, 25, 50, 50, 25);
    
    private const double _mu0 = 1.25653706212 * 10e-6;
    //private const double _mu0 = 4 * Math.PI * 10e-7;

    public FEM(Grid grid, ITimeGrid timeGrid)
    {
        _grid = grid;
        _timeGrid = timeGrid;
        _basis = new TriLinearVectorBasis();
        _integration = new Integration(new SegmentGaussOrder9());
        _stiffnessMatrix = new(_basis.Size);
        _massMatrix = new(_basis.Size);
        _localVector = new(_basis.Size);
    }

    public void SetTest(Test test)
    {
        _test = test;
    }

    public void SetSolver(SLAE slae)
    {
        _slae = slae;
    }

    public void SetScheme(Scheme scheme)
    {
        _scheme = scheme;
        _activeScheme = scheme;
    }

    public void Compute()
    {
        BuildPortrait();
        PrepareLayers();

        int itime = 0;
        
        switch (_scheme)
        {
            case Scheme.Natural:
                itime = 0;

                _activeScheme = Scheme.Two_layer_Implicit;
                itime++;
                
                var genNumber = FindElementNumberForGenerator();
                //HashSet<int> genEdges = new HashSet<int>();
                //HashSet<int> genSkip = new HashSet<int>();
                //foreach (var gen in genNumber)
                //{
                //    for (int iedge = 0; iedge < 4; iedge++)
                //    {
                //        
                //        if ((_grid.Edges[_grid.Elements[gen][iedge]].Point0.X < _generator.xStart ||
                //             _grid.Edges[_grid.Elements[gen][iedge]].Point1.X > _generator.xEnd) &&
                //            (_grid.Edges[_grid.Elements[gen][iedge]].Point0.Y < _generator.yStart ||
                //             _grid.Edges[_grid.Elements[gen][iedge]].Point1.Y > _generator.yEnd))
                //        {
                //            genEdges.Add(_grid.Elements[gen][iedge]);
                //        }
                //        else
                //        {
                //            genSkip.Add(_grid.Elements[gen][iedge]);
                //        }
                //    }
                //}
                
                var hx = 5;
                var hy = 5;
                double delta = 1e-10;
                
                var stepsX = _generator.Length / hx;
                var stepsY = _generator.Width / hy;

                for (int i = 0; i < _grid.Edges.Length; i++)
                {
                    double point;
                    double len1;
                    double len2;
                    
                    switch (_grid.Edges[i].GetAxis())
                    {
                        case 0:
                            point = _generator.xStart + hx / 2.0;
                            
                            for (int j = 0; j < stepsX; j++)
                            {
                                len1 = Math.Sqrt(Math.Pow(point - _grid.Edges[i].Point.X, 2) +
                                                 Math.Pow(_generator.yStart - _grid.Edges[i].Point.Y, 2) +
                                                 Math.Pow(_generator.zStart - _grid.Edges[i].Point.Z, 2));
                                
                                if (len1 < 1e-10)
                                {
                                    len1 = delta;
                                }

                                len2 = Math.Sqrt(Math.Pow(point - _grid.Edges[i].Point.X, 2) +
                                                 Math.Pow(_generator.yEnd - _grid.Edges[i].Point.Y, 2) +
                                                 Math.Pow(_generator.zStart - _grid.Edges[i].Point.Z, 2));
                                
                                if (len2 < 1e-10)
                                {
                                    len2 = delta;
                                }
                                
                                _solution[i] += _mu0 / (4 * Math.PI) * hx / len1;
                                _solution[i] += -_mu0 / (4 * Math.PI) * hx / len2;
                                
                                point += hx;
                            }

                            break;
                        
                        case 1:
                            point = _generator.yStart + hy / 2.0;
                            
                            for (int j = 0; j < stepsY; j++)
                            {
                                len1 = Math.Sqrt(
                                    Math.Pow(_generator.xEnd - _grid.Edges[i].Point.X, 2) +
                                    Math.Pow(point - _grid.Edges[i].Point.Y, 2) +
                                    Math.Pow(_generator.zStart - _grid.Edges[i].Point.Z, 2));

                                if (len1 < 1e-10)
                                {
                                    len1 = delta;
                                }

                                len2 = Math.Sqrt(
                                    Math.Pow(_generator.xStart - _grid.Edges[i].Point.X, 2) +
                                    Math.Pow(point - _grid.Edges[i].Point.Y, 2) +
                                    Math.Pow(_generator.zStart - _grid.Edges[i].Point.Z, 2));

                                if (len2 < 1e-10)
                                {
                                    len2 = delta;
                                }
                                
                                _solution[i] += _mu0 / (4 * Math.PI) * hy / len1;
                                _solution[i] += -_mu0 / (4 * Math.PI) * hy / len2;
                                
                                point += hy;
                            }
                            
                            break;
                        
                        case 2:
                            break;
                    }
                    
                }
                
                //foreach (var iGen in genNumber)
                //{
                //    stepsX = (_grid.Edges[_grid.Elements[iGen][11]].Point1.X -
                //             _grid.Edges[_grid.Elements[iGen][0]].Point0.X) / hx;
                //    stepsY = (_grid.Edges[_grid.Elements[iGen][11]].Point1.Y -
                //             _grid.Edges[_grid.Elements[iGen][0]].Point0.Y) / hy;
                //    

                //    for (int i = 0; i < _grid.Edges.Length; i++)
                //    {
                //        var sum1 = 0;
                //        var sum2 = 0;
                //        double point = 0;
                //        if (i == _grid.Elements[iGen][0] || i == _grid.Elements[iGen][1] || i == _grid.Elements[iGen][2] || i == _grid.Elements[iGen][3])
                //        {
                //            delta = 1e-10;
                //        }
                //        else
                //        {
                //            delta = 0;
                //        }

                //        switch (_grid.Edges[i].GetAxis())
                //        {
                //            case 0:
                //                point = _grid.Edges[_grid.Elements[iGen][0]].Point0.X + hx / 2.0;
                //                
                //                for (int j = 0; j < stepsX; j++)
                //                {
                //                    var tm1 = _mu0 / (4 * Math.PI) * hx / (Math.Sqrt(
                //                        Math.Pow(point - _grid.Edges[i].Point.X, 2) +
                //                        Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][0]].Point0.Y - _grid.Edges[i].Point.Y,
                //                            2) + Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][0]].Point0.Z - _grid.Edges[i].Point.Z,
                //                            2)) + delta);
                //                    _solution[i] += tm1;
                //                    
                //                    var tm2 = -_mu0 / (4 * Math.PI) * hx / (Math.Sqrt(
                //                        Math.Pow(point - _grid.Edges[i].Point.X, 2) +
                //                        Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][1]].Point0.Y - _grid.Edges[i].Point.Y,
                //                            2) + Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][1]].Point0.Z - _grid.Edges[i].Point.Z,
                //                            2)) + delta);
                //                    _solution[i] += tm2;
                //                    //if (i == 13104)
                //                    //{
                //                    //    Console.WriteLine($"{tm1} {tm2}");
                //                    //    //Console.WriteLine($"{} \t {} \t {}");
                //                    //}
                //                    
                //                    
                //                    point += hx;
                //                }

                //                //if (i == 13104)
                //                //{
                //                //    Console.WriteLine($"{tm1} {tm2}");
                //                //}
                //                
                //                break;
                //            
                //            case 1:
                //                point = _grid.Edges[_grid.Elements[iGen][0]].Point0.Y + hy / 2.0;
                //                
                //                for (int j = 0; j < stepsY; j++)
                //                {
                //                    var t1 = _mu0 / (4 * Math.PI) * hy / (Math.Sqrt(
                //                        Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][3]].Point0.X - _grid.Edges[i].Point.X,
                //                            2) + Math.Pow(point - _grid.Edges[i].Point.Y, 2) +
                //                        Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][3]].Point0.Z - _grid.Edges[i].Point.Z,
                //                            2)) + delta);
                //                    
                //                    _solution[i] += t1;
                //                    
                //                    var t2 = -_mu0 / (4 * Math.PI) * hy / (Math.Sqrt(
                //                        Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][2]].Point0.X - _grid.Edges[i].Point.X,
                //                            2) + Math.Pow(point - _grid.Edges[i].Point.Y, 2) +
                //                        Math.Pow(
                //                            _grid.Edges[_grid.Elements[iGen][2]].Point0.Z - _grid.Edges[i].Point.Z,
                //                            2)) + delta);
                //                    
                //                    _solution[i] += t2;
                //                    
                //                    //if (i == 13125)
                //                    //{
                //                    //    Console.WriteLine($"{t1} {t2}");
                //                    //    Console.WriteLine($"{Math.Pow(
                //                    //        _grid.Edges[_grid.Elements[iGen][2]].Point0.X - _grid.Edges[i].Point.X,
                //                    //        2)} \t {Math.Pow(point - _grid.Edges[i].Point.Y, 2)} \t {Math.Pow(
                //                    //        _grid.Edges[_grid.Elements[iGen][2]].Point0.Z - _grid.Edges[i].Point.Z,
                //                    //        2)}");
                //                    //}
                //                    
                //                    point += hy;
                //                }
                //                
                //                break;
                //            
                //            case 2:
                //                break;
                //        }
                //    }
                //    
                //    //_solution[_grid.Elements[iGen][0]] = 1e-6;
                //    //_solution[_grid.Elements[iGen][1]] = -1e-6;
                //    //_solution[_grid.Elements[iGen][2]] = -1e-6;
                //    //_solution[_grid.Elements[iGen][3]] = 1e-6;
                //}

                var A1 = GetValue(new Point3D(-10, 0.0, 30.0));
                var A2 = GetValue(new Point3D(0.0, -10, 30.0));
                var modA1 = Math.Sqrt(Math.Pow(A1.X, 2) + Math.Pow(A1.Y, 2) + Math.Pow(A1.Z, 2));
                var modA2 = Math.Sqrt(Math.Pow(A2.X, 2) + Math.Pow(A2.Y, 2) + Math.Pow(A2.Z, 2));
                
                var A3 = GetValue(new Point3D(8, 0.0, 30.0));
                var A4 = GetValue(new Point3D(0.0, 8, 30.0));
                var modA3 = Math.Sqrt(Math.Pow(A3.X, 2) + Math.Pow(A3.Y, 2) + Math.Pow(A3.Z, 2));
                var modA4 = Math.Sqrt(Math.Pow(A4.X, 2) + Math.Pow(A4.Y, 2) + Math.Pow(A4.Z, 2));
                
                var A5 = GetValue(new Point3D(-48, 0.0, -10.0));
                var A6 = GetValue(new Point3D(0.0, -48, -10.0));
                var modA5 = Math.Sqrt(Math.Pow(A5.X, 2) + Math.Pow(A5.Y, 2) + Math.Pow(A5.Z, 2));
                var modA6 = Math.Sqrt(Math.Pow(A6.X, 2) + Math.Pow(A6.Y, 2) + Math.Pow(A6.Z, 2));
                
                var A7 = GetValue(new Point3D(29, 0.0, 5.0));
                var A8 = GetValue(new Point3D(0.0, 29, 5.0));
                var modA7 = Math.Sqrt(Math.Pow(A7.X, 2) + Math.Pow(A7.Y, 2) + Math.Pow(A7.Z, 2));
                var modA8 = Math.Sqrt(Math.Pow(A8.X, 2) + Math.Pow(A8.Y, 2) + Math.Pow(A8.Z, 2));
                
                Vector.Copy(_solution, _layers[0]);
                
                break;
            
            case Scheme.Two_layer_Implicit:
                itime = 1;
                break;
            
            case Scheme.Three_layer_Implicit:
                itime = 2;
                break;
            
            case Scheme.Four_layer_Implicit:
                itime = 3;
                break;
        }
        

        using (var sw = new StreamWriter("Tests/5.csv"))
        {
            for ( ; itime < _timeGrid.TGrid.Length; itime++)
            {
                AssemblySLAE(itime);
                AccountDirichletBoundaries(itime);

                //Vector kek = new Vector(_globalVector.Length);
                //for (int i = 0; i < kek.Length; i++)
                //{
                //    kek[i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[itime], _grid.Edges[i].GetAxis());
                //}

                //_solution = _globalMatrix * kek;
                //break;
                
                _slae.SetSLAE(_globalVector, _globalMatrix);
                _solution = _slae.Solve();

                //var kek = GetValue(new Point3D(0.0, 0.0, 30));
                //var anime = CalculateEMF(new Point3D(0.25, 0.25, 30), itime);
                var A1 = GetValue(new Point3D(2.5, 0.0, 0.0));
                var A2 = GetValue(new Point3D(0.0, 2.5, 0.0));
                var modA1 = Math.Sqrt(Math.Pow(A1.X, 2) + Math.Pow(A1.Y, 2) + Math.Pow(A1.Z, 2));
                var modA2 = Math.Sqrt(Math.Pow(A2.X, 2) + Math.Pow(A2.Y, 2) + Math.Pow(A2.Z, 2));
                
                
                //Console.WriteLine($"{itime} = {_timeGrid[itime]} EMF = {CalculateEMF(new Point3D(0.1, 0.1, 30), itime)}");
                Console.WriteLine($"{itime} = {_timeGrid[itime]} dBz = {-CalculatedBz(new Point3D(0, 0, 25), itime)}");

                switch (_activeScheme)
                {
                    case Scheme.Two_layer_Implicit:
                        if (_scheme == Scheme.Natural)
                        {
                            _activeScheme = Scheme.Three_layer_Implicit;
                            
                            Vector.Copy(_solution, _layers[1]);
                        }
                        
                        else
                        {
                            Vector.Copy(_solution, _layers[0]);
                        }
                        break;
                    
                    case Scheme.Three_layer_Implicit:
                        if (_scheme == Scheme.Natural)
                        {
                            //_activeScheme = Scheme.Four_layer_Implicit;
                            
                            //Vector.Copy(_solution, _layers[2]);
                            
                            Vector.Copy(_layers[1], _layers[0]);
                            Vector.Copy(_solution, _layers[1]);
                        }

                        else
                        {
                            Vector.Copy(_layers[1], _layers[0]);
                            Vector.Copy(_solution, _layers[1]);
                        }
                        
                        break;
                    
                    case Scheme.Four_layer_Implicit:
                        Vector.Copy(_layers[1], _layers[0]);
                        Vector.Copy(_layers[2], _layers[1]);
                        Vector.Copy(_solution, _layers[2]);
                        break;
                }
                
                //double error = 0;
                //for (int i = 0; i < _grid.Edges.Length; i++)
                //{
                //    error += Math.Pow(
                //        _test.UValue(_grid.Edges[i].Point, _timeGrid[itime], _grid.Edges[i].GetAxis()) - _solution[i], 2);
                //}
                //PrintError(itime);
                //sw.WriteLine($"{_timeGrid[itime]},{Math.Sqrt(error / _grid.Edges.Length)}");

                //for (int i = 0; i < _grid.Edges.Length; i++)
                //{
                //    sw.WriteLine(
                //        $"{i + 1},{_solution[i]},{_test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis())},{_test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis()) - _solution[i]}");
                //    
                //    Console.WriteLine(
                //        $"{i + 1},{_solution[i]},{_test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis())},{_test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis()) - _solution[i]}");
                //}
                //break;
            }
        }
    }

    private void AssemblySLAE(int itime)
    {
        _globalVector.Fill(0);
        _globalMatrix.Clear();

        for (int ielem = 0; ielem < _grid.Elements.Length; ielem++)
        {
            AssemblyLocalElement(ielem, itime);

            if (_activeScheme != Scheme.Natural)
            {
                _stiffnessMatrix += SchemeUsage(ielem, itime, _activeScheme, 0) * _massMatrix;
            }

            for (int i = 0; i < _basis.Size; i++)
            {
                for (int j = 0; j < _basis.Size; j++)
                {
                    AddElement(_grid.Elements[ielem][i], _grid.Elements[ielem][j], _stiffnessMatrix[i, j]);
                }
            }
            
            AssemblyGlobalVector(ielem, itime);
            
            _stiffnessMatrix.Clear();
            _massMatrix.Clear();
            _localVector.Fill(0);
        }
    }

    private void AssemblyGlobalVector(int ielem, int itime)
    {
        double[] qj3 = new double[_basis.Size];
        double[] qj2 = new double[_basis.Size];
        double[] qj1 = new double[_basis.Size];
        
        switch (_activeScheme)
        {
            case Scheme.Natural:
              if (_grid.Edges[_grid.Elements[ielem][0]].Point0.X >= -51 &&
                  _grid.Edges[_grid.Elements[ielem][0]].Point0.X <= 51 &&
                  _grid.Edges[_grid.Elements[ielem][0]].Point0.Y >= -51 &&
                  _grid.Edges[_grid.Elements[ielem][0]].Point0.Y <= 51 &&
                  _grid.Edges[_grid.Elements[ielem][0]].Point0.Z <= 51 &&
                  _grid.Edges[_grid.Elements[ielem][11]].Point1.Z >= 30)
              {
                  for (int i = 0; i < _basis.Size; i++)
                  {
                      _globalVector[_grid.Elements[ielem][i]] = 1;
                  }
              }

              break;
            
            case Scheme.Two_layer_Implicit:
                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj1[i] += _massMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }

                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _localVector[i] += SchemeUsage(ielem, itime, _activeScheme, 1) * qj1[i];

                    _globalVector[_grid.Elements[ielem][i]] += _localVector[i];
                }
                
                break;
            
            case Scheme.Three_layer_Implicit:

                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj2[i] += _massMatrix[i, j] * _layers[1][_grid.Elements[ielem][j]];
                        qj1[i] += _massMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }
                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _localVector[i] += SchemeUsage(ielem, itime, _activeScheme, 1) * qj2[i];
                    _localVector[i] += SchemeUsage(ielem, itime, _activeScheme, 2) * qj1[i];
                    
                    _globalVector[_grid.Elements[ielem][i]] += _localVector[i];
                }
                
                break;
            
            case Scheme.Four_layer_Implicit:

                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj3[i] += _massMatrix[i, j] * _layers[2][_grid.Elements[ielem][j]];
                        qj2[i] += _massMatrix[i, j] * _layers[1][_grid.Elements[ielem][j]];
                        qj1[i] += _massMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }
                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _localVector[i] += SchemeUsage(ielem, itime, _activeScheme, 1) * qj3[i];
                    _localVector[i] += SchemeUsage(ielem, itime, _activeScheme, 2) * qj2[i];
                    _localVector[i] += SchemeUsage(ielem, itime, _activeScheme, 3) * qj1[i];
                    
                    _globalVector[_grid.Elements[ielem][i]] += _localVector[i];
                }
                
                break;
        }
    }

    private void AddElement(int i, int j, double value)
    {
        if (i == j)
        {
            _globalMatrix.Di[i] += value;
            return;
        }

        if (i > j)
        {
            for (int icol = _globalMatrix.Ig[i]; icol < _globalMatrix.Ig[i + 1]; icol++)
            {
                if (_globalMatrix.Jg[icol] == j)
                {
                    _globalMatrix.Ggl[icol] += value;
                    return;
                }
            }
        }

        else
        {
            for (int icol = _globalMatrix.Ig[j]; icol < _globalMatrix.Ig[j + 1]; icol++)
            {
                if (_globalMatrix.Jg[icol] == i)
                {
                    _globalMatrix.Ggu[icol] += value;
                    return;
                }
            }
        }
    }

    private void AssemblyLocalElement(int ielem, int itime)
    {
        double hx = _grid.Edges[_grid.Elements[ielem][0]].Length;
        double hy = _grid.Edges[_grid.Elements[ielem][2]].Length;
        double hz = _grid.Edges[_grid.Elements[ielem][4]].Length;

        for (int i = 0; i < _basis.Size; i++)
        {
            for (int j = 0; j < _basis.Size; j++)
            {
                Func<Point3D, double> kek;
                Vector3D psi1 = new(0, 0, 0);
                Vector3D psi2 = new(0, 0, 0);
                Vector3D dPsi1 = new(0, 0, 0);
                Vector3D dPsi2 = new(0, 0, 0);

                int ik = i;
                int jk = j;
                kek = point =>
                {
                    psi1.Copy(_basis.GetPsi(ik, point));
                    psi2.Copy(_basis.GetPsi(jk, point));

                    return psi1 * psi2;
                };

                _massMatrix[i, j] += hx * hy * hz * _integration.Gauss3D(kek);

                kek = point =>
                {
                    dPsi1.Copy(_basis.GetDPsi(ik, point));
                    dPsi2.Copy(_basis.GetDPsi(jk, point));

                    return Vector3D.DotProductJacob(dPsi1, dPsi2, hx, hy, hz);
                };

                _stiffnessMatrix[i, j] += 1 / _grid.Mu * _integration.Gauss3D(kek);
            }

            _localVector[i] = _test.F(_grid.Edges[_grid.Elements[ielem][i]].Point, _timeGrid[itime], i, _grid.GetSigma(
                new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                    _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                    _grid.Edges[_grid.Elements[ielem][11]].Point.Z)));
        }

        _localVector = _massMatrix * _localVector;
    }
    
    private void AccountDirichletBoundaries(int itime)
    {
        foreach (var edge in _grid.DirichletBoundaries)
        {

            _globalMatrix.Di[edge] = 1;
            _globalVector[edge] = _test.UValue(_grid.Edges[edge].Point, _timeGrid[itime], _grid.Edges[edge].GetAxis());

            for (int i = _globalMatrix.Ig[edge]; i < _globalMatrix.Ig[edge + 1]; i++)
                _globalMatrix.Ggl[i] = 0;

            for (int col = edge + 1; col < _globalMatrix.Size; col++)
                for (int j = _globalMatrix.Ig[col]; j < _globalMatrix.Ig[col + 1]; j++)
                    if (_globalMatrix.Jg[j] == edge)
                    {
                        _globalMatrix.Ggu[j] = 0;
                        break;
                    }
        }
    }

    private double SchemeUsage(int ielem, int itime, Scheme scheme, int i)
    {
       double t01;
       double t02;
       double t12;
       
       switch (scheme)
        {
            case Scheme.Two_layer_Implicit:
                t01 = _timeGrid[itime] - _timeGrid[itime - 1];
                switch (i)
                {
                    // тут только параболическая, эпсилон для гиперболики придётся добавить
                    case 0:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                            _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                            _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) / t01);
                    
                    case 1:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                            _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                            _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) / t01);
                }
                
                break;
            
            case Scheme.Three_layer_Implicit:
                t01 = _timeGrid[itime] - _timeGrid[itime - 1];
                t02 = _timeGrid[itime] - _timeGrid[itime - 2];
                t12 = _timeGrid[itime - 1] - _timeGrid[itime - 2];
                
                switch (i)
                {
                    case 0:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                   _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                   _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * (t01 + t02) + 2 * _grid.Epsilon) /
                               (t01 * t02);
                    
                    case 1:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                            _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                            _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t02 + 2 * _grid.Epsilon) / (t01 * t12);
                    
                    case 2:
                        return -(_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                            _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                            _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t01 + 2 * _grid.Epsilon) / (t02 * t12);
                    
                }

                break;
            
            case Scheme.Four_layer_Implicit:
                t01 = _timeGrid[itime] - _timeGrid[itime - 1];
                t02 = _timeGrid[itime] - _timeGrid[itime - 2];
                t12 = _timeGrid[itime - 1] - _timeGrid[itime - 2];
                double t03 = _timeGrid[itime] - _timeGrid[itime - 3];
                double t13 = _timeGrid[itime - 1] - _timeGrid[itime - 3];
                double t23 = _timeGrid[itime - 2] - _timeGrid[itime - 3];

                switch (i)
                {
                    case 0:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                       _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                       _grid.Edges[_grid.Elements[ielem][11]].Point.Z))
                                   * (t01 * t02 + t01 * t03 + t02 * t03) + 2 * _grid.Epsilon * (t01 + t02 + t03)) /
                               (t01 * t02 * t03);
                    
                    case 1:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                    _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                    _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t02 * t03 +
                                2 * _grid.Epsilon * (t02 + t03)) / (t01 * t12 * t13);
                    
                    case 2:
                        return -(_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                     _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                     _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t01 * t03 +
                                 2 * _grid.Epsilon * (t01 + t03)) / (t02 * t12 * t23);
                    
                    case 3:
                        return (_grid.GetSigma(new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
                                    _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
                                    _grid.Edges[_grid.Elements[ielem][11]].Point.Z)) * t01 * t02 +
                                2 * _grid.Epsilon * (t01 + t02)) / (t03 * t13 * t23);
                }

                return 0;
            
            default:
                throw new Exception("Undefined scheme");
        }

        throw new Exception("SchemeUsage can't return a value");
    }

    private void BuildPortrait()
    {
        HashSet<int>[] list = new HashSet<int>[_grid.Edges.Length].Select(_ => new HashSet<int>()).ToArray();
        foreach (var element in _grid.Elements)
            foreach (var pos in element)
                foreach (var node in element)
                    if (pos > node)
                        list[pos].Add(node);

        list = list.Select(childlist => childlist.Order().ToHashSet()).ToArray();
        int count = list.Sum(childlist => childlist.Count);

        _globalMatrix = new(_grid.Edges.Length, count);
        _globalVector = new(_grid.Edges.Length);
        _solution = new(_grid.Edges.Length);
        _layers = new Vector[3].Select(_ => new Vector(_grid.Edges.Length)).ToArray();

        _globalMatrix.Ig[0] = 0;

        for (int i = 0; i < list.Length; i++)
            _globalMatrix.Ig[i + 1] = _globalMatrix.Ig[i] + list[i].Count;

        int k = 0;

        foreach (var childlist in list)
            foreach (var value in childlist)
                _globalMatrix.Jg[k++] = value;
    }
    
    private void PrepareLayers()
    {
        switch (_scheme)
        {
            case Scheme.Two_layer_Implicit:
                for (int i = 0; i < _grid.Edges.Length; i++)
                {
                    _layers[0][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis());
                }
                break;
            
            case Scheme.Three_layer_Implicit:
                for (int i = 0; i < _grid.Edges.Length; i++)
                {
                    _layers[0][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis());
                    _layers[1][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[1], _grid.Edges[i].GetAxis());
                }
                
                break;
            
            case Scheme.Four_layer_Implicit:
                for (int i = 0; i < _grid.Edges.Length; i++)
                {
                    _layers[0][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[0], _grid.Edges[i].GetAxis());
                    _layers[1][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[1], _grid.Edges[i].GetAxis());
                    _layers[2][i] = _test.UValue(_grid.Edges[i].Point, _timeGrid[2], _grid.Edges[i].GetAxis());
                }
                
                break;
        }
    }

    private void PrintError(int itime)
    {
        double error = 0;
        for (int i = 0; i < _grid.Edges.Length; i++)
        {
            error += Math.Pow(
                _test.UValue(_grid.Edges[i].Point, _timeGrid[itime], _grid.Edges[i].GetAxis()) - _solution[i], 2);
        }
        Console.WriteLine($"Layer error {itime} = {Math.Sqrt(error / _grid.Edges.Length)}");
    }

    public Vector3D GetValue(Point3D point)
    {
        var vector = new Vector3D(0, 0, 0);
        var kek = new Point3D(0, 0, 0);

        foreach (var elem in _grid.Elements)
        {
            if (point.X >= _grid.Edges[elem[0]].Point0.X && point.X < _grid.Edges[elem[11]].Point1.X &&
                point.Y >= _grid.Edges[elem[0]].Point0.Y && point.Y < _grid.Edges[elem[11]].Point1.Y &&
                point.Z >= _grid.Edges[elem[0]].Point0.Z && point.Z < _grid.Edges[elem[11]].Point1.Z)
            {
                kek.X = (point.X - _grid.Edges[elem[0]].Point0.X) / _grid.Edges[elem[0]].Length;
                kek.Y = (point.Y - _grid.Edges[elem[2]].Point0.Y) / _grid.Edges[elem[2]].Length;
                kek.Z = (point.Z - _grid.Edges[elem[4]].Point0.Z) / _grid.Edges[elem[4]].Length;
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    vector += _basis.GetPsi(i, kek) * _solution[elem[i]];
                }
                
                break;
            }
        }

        return vector;
    }

    private double CalculateEMF(Point3D point, int itime)
    {
        double result = 0;
        double height = 30;
        double pointX = 0;
        double pointY = 0;
        
        for (int ielem = 0; ielem < _grid.Elements.Length; ielem++)
        {
            if (point.X >= _grid.Edges[_grid.Elements[ielem][0]].Point0.X &&
                point.X < _grid.Edges[_grid.Elements[ielem][11]].Point1.X &&
                point.Y >= _grid.Edges[_grid.Elements[ielem][0]].Point0.Y &&
                point.Y < _grid.Edges[_grid.Elements[ielem][11]].Point1.Y &&
                point.Z >= _grid.Edges[_grid.Elements[ielem][0]].Point0.Z &&
                point.Z < _grid.Edges[_grid.Elements[ielem][11]].Point1.Z)
            {
                var hx = 1e-4;
                var hy = 1e-4;
                
                var stepsX = (_grid.Edges[_grid.Elements[ielem][11]].Point1.X -
                         _grid.Edges[_grid.Elements[ielem][0]].Point0.X) / hx;
                var stepsY = (_grid.Edges[_grid.Elements[ielem][11]].Point1.Y -
                         _grid.Edges[_grid.Elements[ielem][0]].Point0.Y) / hy;

                pointX = _grid.Edges[_grid.Elements[ielem][0]].Point0.X + hx / 2;
                pointY = _grid.Edges[_grid.Elements[ielem][0]].Point0.Y + hy / 2;
                
                for (int i = 0; i < stepsX; i++)
                {
                    result += hx *
                              GetValueFordAdt(
                                  new Point3D(pointX, _grid.Edges[_grid.Elements[ielem][0]].Point0.Y, point.Z), ielem,
                                  itime).X;
                    
                    result += -hx *
                              GetValueFordAdt(
                                  new Point3D(pointX, _grid.Edges[_grid.Elements[ielem][11]].Point1.Y, point.Z), ielem,
                                  itime).X;
                    
                    pointX += hx;
                }
                
                for (int i = 0; i < stepsY; i++)
                {
                    result += hy *
                              GetValueFordAdt(
                                  new Point3D( _grid.Edges[_grid.Elements[ielem][0]].Point0.X, pointY, point.Z), ielem,
                                  itime).Y;
                    
                    result += -hy *
                              GetValueFordAdt(
                                  new Point3D(_grid.Edges[_grid.Elements[ielem][11]].Point1.X, pointY, point.Z), ielem,
                                  itime).Y;
                    
                    pointY += hy;
                }
                
                break;
            }
        }

        return result;
    }

    private Vector3D GetValueFordAdt(Point3D point, int ielem, int itime)
    {
        var vector1 = new Vector3D(0, 0, 0);
        var vector2 = new Vector3D(0, 0, 0);
        var kek = new Point3D(0, 0, 0);
        var dt = _timeGrid[itime] - _timeGrid[itime - 1];
        
        kek.X = (point.X - _grid.Edges[_grid.Elements[ielem][0]].Point0.X) / _grid.Edges[_grid.Elements[ielem][0]].Length;
        kek.Y = (point.Y - _grid.Edges[_grid.Elements[ielem][2]].Point0.Y) / _grid.Edges[_grid.Elements[ielem][2]].Length;
        kek.Z = (point.Z - _grid.Edges[_grid.Elements[ielem][4]].Point0.Z) / _grid.Edges[_grid.Elements[ielem][4]].Length;
        
        for (int i = 0; i < _basis.Size; i++)
        {
            vector1 += _basis.GetPsi(i, kek) * _solution[_grid.Elements[ielem][i]];
            switch (_activeScheme)
            {
                case Scheme.Two_layer_Implicit:
                    vector2 += _basis.GetPsi(i, kek) * _layers[0][_grid.Elements[ielem][i]];
                    break;
                
                case Scheme.Three_layer_Implicit:
                    vector2 += _basis.GetPsi(i, kek) * _layers[1][_grid.Elements[ielem][i]];
                    break;
                
                case Scheme.Four_layer_Implicit:
                    vector2 += _basis.GetPsi(i, kek) * _layers[2][_grid.Elements[ielem][i]];
                    break;
            }
        }

        vector1 -= vector2;
        vector1 /= dt;
        
        return vector1;
    }
    
    private Vector3D GetValueForRotAdt(Point3D point, int ielem, int itime)
    {
        var vector1 = new Vector3D(0, 0, 0);
        var vector2 = new Vector3D(0, 0, 0);
        var kek = new Point3D(0, 0, 0);
        var dt = _timeGrid[itime] - _timeGrid[itime - 1];
        
        kek.X = (point.X - _grid.Edges[_grid.Elements[ielem][0]].Point0.X) / _grid.Edges[_grid.Elements[ielem][0]].Length;
        kek.Y = (point.Y - _grid.Edges[_grid.Elements[ielem][2]].Point0.Y) / _grid.Edges[_grid.Elements[ielem][2]].Length;
        kek.Z = (point.Z - _grid.Edges[_grid.Elements[ielem][4]].Point0.Z) / _grid.Edges[_grid.Elements[ielem][4]].Length;
        
        for (int i = 0; i < _basis.Size; i++)
        {
            if (i is >= 4 and <= 7)
            {
                continue;
            }
            
            vector1 += _basis.GetDPsi(i, kek) * _solution[_grid.Elements[ielem][i]];
            
            switch (_activeScheme)
            {
                case Scheme.Two_layer_Implicit:
                    vector2 += _basis.GetDPsi(i, kek) * _layers[0][_grid.Elements[ielem][i]];
                    break;
                
                case Scheme.Three_layer_Implicit:
                    vector2 += _basis.GetDPsi(i, kek) * _layers[1][_grid.Elements[ielem][i]];
                    break;
                
                case Scheme.Four_layer_Implicit:
                    vector2 += _basis.GetDPsi(i, kek) * _layers[2][_grid.Elements[ielem][i]];
                    break;
            }
        }

        vector1 -= vector2;
        vector1 /= dt;
        
        return vector1;
    }
    
    private double CalculatedBz(Point3D point, int itime)
    {
        Vector3D vec1 = new(0.0, 0.0, 0.0);
        double Axy = 0;
        double Ayx = 0;

        for (int ielem = 0; ielem < _grid.Elements.Length; ielem++)
        {
            if (point.X >= _grid.Edges[_grid.Elements[ielem][0]].Point0.X &&
                point.X < _grid.Edges[_grid.Elements[ielem][11]].Point1.X &&
                point.Y >= _grid.Edges[_grid.Elements[ielem][0]].Point0.Y &&
                point.Y < _grid.Edges[_grid.Elements[ielem][11]].Point1.Y &&
                point.Z >= _grid.Edges[_grid.Elements[ielem][0]].Point0.Z &&
                point.Z < _grid.Edges[_grid.Elements[ielem][11]].Point1.Z)
                {
                    vec1 = GetValueForRotAdt(point, ielem, itime);
                    
                    break;
                }
        }

        //return Axy - Ayx;
        return vec1.Y - vec1.X;
    }

    private List<int> FindElementNumberForGenerator()
    {
        List<int> genNumbers = new List<int>();
        for (int i = 0; i < _grid.Elements.Length; i++)
        {

            if (!(_grid.Edges[_grid.Elements[i][0]].Point0.X > _generator.xEnd || _generator.xStart > _grid.Edges[_grid.Elements[i][0]].Point1.X ||
                  _grid.Edges[_grid.Elements[i][2]].Point0.Y > _generator.yEnd || _generator.yStart > _grid.Edges[_grid.Elements[i][2]].Point1.Y) &&
                  _grid.Edges[_grid.Elements[i][0]].Point0.Z <= _generator.zStart && _grid.Edges[_grid.Elements[i][11]].Point1.Z >= _generator.zEnd)
            {
                genNumbers.Add(i);
            }
        }

        return genNumbers;
    }
}
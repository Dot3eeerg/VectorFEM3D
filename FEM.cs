﻿namespace VectorFEM3D;

public class FEM
{
    private SparseMatrix? _globalMatrix;
    private Vector? _globalVector;
    private Vector[]? _layers;
    private Vector? _solution;
    private Vector? _solution1;
    private Grid? _grid;
    private ITimeGrid _timeGrid;
    private Test? _test;
    private IBasis3D? _basis;
    private IBasis2D? _basis2D;
    private Integration? _integration;
    private SLAE? _slae;
    private Scheme _scheme;
    private Scheme _activeScheme;
    private double _sum = 0;

    private Vector _2DSolution;

    private ElementAssembler _elementAssembler;
    
    private Point3D[] _elementPointList;

    public FEM(Grid grid, ITimeGrid timeGrid)
    {
        _grid = grid;
        _timeGrid = timeGrid;
        _basis = new TriLinearVectorBasis();
        _basis2D = new BiLinearBasis();
        _integration = new Integration(new SegmentGaussOrder9());
        _elementPointList = new Point3D[_basis.Size];
        _elementAssembler = new ElementAssembler(_basis);
        BuildPortrait();
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

    public void Convert2DSolution(Grid2D grid)
    {
        _2DSolution = new Vector(grid.Nodes.Length);
        
        using (var sr = new StreamReader("Tests/Fifth/2DResults"))
        //using (var sr = new StreamReader("Tests/Fifth/EMFResultsLast"))
        //using (var sr = new StreamReader("Tests/Fifth/EMFResultsMiddle"))
        {
            for (int i = 0; i < grid.Nodes.Length; i++)
            {
                _2DSolution[i] = Convert.ToDouble(sr.ReadLine());
            }
        }
        
        Point2D point = new Point2D(0, 0);
        double result;
        
        for (int i = 0; i < _solution.Length; i++)
        {
            switch (_grid.Edges[i].GetAxis())
            {
                case 0:
                    point.X = Math.Sqrt(_grid.Edges[i].Point.X * _grid.Edges[i].Point.X +
                                        _grid.Edges[i].Point.Y * _grid.Edges[i].Point.Y);
                    point.Y = _grid.Edges[i].Point.Z;

                    result = GetValue(point, grid);

                    _solution[i] = -_grid.Edges[i].Point.Y / point.X * result;
                    break;
                
                case 1:
                    point.X = Math.Sqrt(_grid.Edges[i].Point.X * _grid.Edges[i].Point.X +
                                        _grid.Edges[i].Point.Y * _grid.Edges[i].Point.Y);
                    point.Y = _grid.Edges[i].Point.Z;
                    
                    result = GetValue(point, grid);

                    _solution[i] = _grid.Edges[i].Point.X / point.X * result;
                    
                    break;
                
                case 2:
                    _solution[i] = 0;

                    break;
            }
        }
        
        //Point3D anime = new Point3D(0, 0, 0);
        //Point2D kek = new Point2D(0, 0);
        //for (int j = 1; j < _grid.Yvalues.Count; j++)
        //{
        //    //anime.X = (_grid.Xvalues[j] + _grid.Xvalues[j - 1]) / 2;
        //    anime.X = 0;
        //    anime.Y = (_grid.Yvalues[j] + _grid.Yvalues[j - 1]) / 2;
        //    anime.Z = _grid.GeneratorZ;
        //    
        //    //Console.Write($"${anime.Y:F}$ & ");
        //    
        //    //Console.Write($"${GetValue(anime).GetLength():e2}$ & ");
        //    
        //    double x = (_grid.Xvalues[j] + _grid.Xvalues[j - 1]) / 2;
        //    double y = (_grid.Yvalues[j] + _grid.Yvalues[j - 1]) / 2;
        //    //kek.X = Math.Sqrt(x * x + y * y);
        //    kek.X = Math.Abs(y);
        //    if (kek.X == 0)
        //    {
        //        continue;
        //    }
        //    kek.Y = _grid.GeneratorZ;
        //    //Console.WriteLine($"${GetValue(kek, grid):e2}$ & $\\dots$ \\\\ \\hline");
        //}
        
        using (var sw1 = new StreamWriter("Tests/A1/x"))
        {
            Point3D anime1 = new Point3D(0, 0, 0);
            
            for (int i = 0; i < _grid.Zvalues.Count; i++)
            {
                for (int j = 1; j < _grid.Yvalues.Count; j++)
                {
                    //anime1.X = (_grid.Xvalues[j] + _grid.Xvalues[j - 1]) / 2;
                    anime1.X = 0;
                    anime1.Y = (_grid.Yvalues[j] + _grid.Yvalues[j - 1]) / 2;
                    anime1.Z = _grid.Zvalues[i];
                    sw1.WriteLine(GetValue(anime1).GetLength());
                }
            }
        }

        using (var sw1 = new StreamWriter("Tests/A1/r0"))
        {
            Point2D kek1 = new Point2D(0, 0);

            for (int i = 0; i < _grid.Zvalues.Count; i++)
            {
                for (int j = 1; j < _grid.Xvalues.Count; j++)
                {
                    double x = (_grid.Xvalues[j] + _grid.Xvalues[j - 1]) / 2;
                    double y = (_grid.Yvalues[j] + _grid.Yvalues[j - 1]) / 2;
                    //kek1.X = Math.Sqrt(x * x + y * y);
                    kek1.X = Math.Abs(y);
                    if (kek1.X == 0)
                    {
                        continue;
                    }
                    kek1.Y = _grid.Zvalues[i];
                    sw1.WriteLine(GetValue(kek1, grid));
                }
            }
        }
    }

    public void Compute(Grid2D grid)
    {
        PrepareLayers();

        int itime = 0;
        
        switch (_scheme)
        {
            case Scheme.Natural:
                itime = 0;
                
                _activeScheme = Scheme.Two_layer_Implicit;
                itime++;
                
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
        
        using (var sw = new StreamWriter("Tests/Fifth/EMFResultsFull"))
        {
            for ( ; itime < _timeGrid.TGrid.Length; itime++)
            {
                AssemblySLAE(itime);
                AccountDirichletBoundaries(itime);

                switch (_activeScheme)
                {
                    case Scheme.Two_layer_Implicit:
                        _slae.SetSLAE(_globalVector, _globalMatrix, _layers[0]);
                        
                        break;
                    
                    case Scheme.Three_layer_Implicit:
                        _slae.SetSLAE(_globalVector, _globalMatrix, _layers[1]);
                        //_slae.SetSLAE(_globalVector, _globalMatrix, _solution);
                        
                        break;
                    
                    case Scheme.Four_layer_Implicit:
                        _slae.SetSLAE(_globalVector, _globalMatrix, _layers[2]);
                        
                        break;
                        
                }
                
                //Console.Write($"Layer {itime}: ");
                _solution = _slae.Solve();
                
                PrintError(itime);
                
                //using (var sw1 = new StreamWriter("Tests/A0/x"))
                //{
                //    Point3D anime = new Point3D(0, 0, 0);
                //    
                //    //for (int i = 0; i < _grid.Zvalues.Count; i++)
                //    //{
                //        for (int j = 1; j < _grid.Yvalues.Count; j++)
                //        {
                //            //anime.X = (_grid.Xvalues[j] + _grid.Xvalues[j - 1]) / 2;
                //            anime.X = 0;
                //            anime.Y = (_grid.Yvalues[j] + _grid.Yvalues[j - 1]) / 2;
                //            anime.Z = _grid.GeneratorZ;
                //            sw1.WriteLine(GetValue(anime).GetLength());
                //        }
                //    //}
                //}

                //using (var sw1 = new StreamWriter("Tests/A0/r0"))
                //{
                //    Point2D kek = new Point2D(0, 0);

                //    //for (int i = 0; i < _grid.Zvalues.Count; i++)
                //    //{
                //        for (int j = 1; j < _grid.Xvalues.Count; j++)
                //        {
                //            double x = (_grid.Xvalues[j] + _grid.Xvalues[j - 1]) / 2;
                //            double y = (_grid.Yvalues[j] + _grid.Yvalues[j - 1]) / 2;
                //            //kek.X = Math.Sqrt(x * x + y * y);
                //            kek.X = Math.Abs(y);
                //            if (kek.X == 0)
                //            {
                //                kek.X = 0.01;
                //            }
                //            kek.Y = _grid.GeneratorZ;
                //            sw1.WriteLine(GetValue(kek, grid));
                //        }
                //    //}
                //}
                //
                //
                //return;

                //var anime = GetValue(new Point3D(1, 0, _grid.GeneratorZ));
                ////sw.WriteLine($"{anime.X} {anime.Y} {anime.Z} {anime.GetLength()}");
                //sw.WriteLine($"{anime.GetLength()}");

                switch (_activeScheme)
                {
                    case Scheme.Two_layer_Implicit:
                        if (_scheme == Scheme.Natural)
                        {
                            if (itime == 50)
                            {
                                _activeScheme = Scheme.Three_layer_Implicit;
                                
                                Vector.Copy(_solution, _layers[1]);
                                
                                break;
                            }
                            
                            Vector.Copy(_solution, _layers[0]);
                        }
                        
                        else
                        {
                            Vector.Copy(_solution, _layers[0]);
                        }
                        break;
                    
                    case Scheme.Three_layer_Implicit:
                        if (_scheme == Scheme.Natural)
                        {
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
            }
        }

        Console.WriteLine(_sum / (_timeGrid.TGrid.Length - 2));
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
                _elementAssembler.StiffnessMatrix += SchemeUsage(ielem, itime, _activeScheme, 0) * _elementAssembler.MassMatrix;
            }

            for (int i = 0; i < _basis.Size; i++)
            {
                for (int j = 0; j < _basis.Size; j++)
                {
                    AddElement(_grid.Elements[ielem][i], _grid.Elements[ielem][j], _elementAssembler.StiffnessMatrix[i, j]);
                }
            }
            
            AssemblyGlobalVector(ielem, itime);
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
                if (_grid.Edges[_grid.Elements[ielem][0]].Point0.Z < _grid.GeneratorZ &&
                    _grid.Edges[_grid.Elements[ielem][11]].Point1.Z > _grid.GeneratorZ)
                {
                    var point = new Point3D(0, 0, 0);
                    
                    // 1
                    if (_grid.Edges[_grid.Elements[ielem][0]].Point0.X > _grid.GeneratorX.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][3]].Point1.X < _grid.GeneratorX.Item2 &&
                        _grid.Edges[_grid.Elements[ielem][0]].Point0.Y < _grid.GeneratorY.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][0]].Point0.Y > _grid.GeneratorY.Item2)
                    {
                        point.X = (_grid.Edges[_grid.Elements[ielem][0]].Point.X -
                                   _grid.Edges[_grid.Elements[ielem][0]].Point0.X) /
                                  _grid.Edges[_grid.Elements[ielem][0]].Length;
                        
                        point.Y = (_grid.GeneratorY.Item1 - _grid.Edges[_grid.Elements[ielem][0]].Point0.Y) /
                                  _grid.Elements[_grid.Elements[ielem][0]].Length;
                        
                        point.Z = 0;
                        
                        for (int i = 0; i < _basis.Size; i++)
                        {
                            _globalVector[_grid.Elements[ielem][i]] += _grid.SourceValue * _basis.GetPsi(i, point).X;
                        }
                    }
                    
                    // 2
                    if (_grid.Edges[_grid.Elements[ielem][0]].Point0.X > _grid.GeneratorX.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][3]].Point1.X < _grid.GeneratorX.Item2 &&
                        _grid.Edges[_grid.Elements[ielem][3]].Point1.Y < _grid.GeneratorY.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][3]].Point1.Y > _grid.GeneratorY.Item2)
                    {
                        point.X = (_grid.Edges[_grid.Elements[ielem][0]].Point.X -
                                   _grid.Edges[_grid.Elements[ielem][0]].Point0.X) /
                                  _grid.Edges[_grid.Elements[ielem][0]].Length;
                        
                        point.Y = (_grid.GeneratorY.Item2 - _grid.Edges[_grid.Elements[ielem][3]].Point0.Y) /
                                  _grid.Elements[_grid.Elements[ielem][0]].Length;
                        
                        point.Z = 0;
                        
                        for (int i = 0; i < _basis.Size; i++)
                        {
                            _globalVector[_grid.Elements[ielem][i]] += _grid.SourceValue * _basis.GetPsi(i, point).X;
                        }
                    }
                    
                    // 3
                    if (_grid.Edges[_grid.Elements[ielem][0]].Point0.X > _grid.GeneratorX.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][0]].Point0.X < _grid.GeneratorX.Item2 &&
                        _grid.Edges[_grid.Elements[ielem][4]].Point0.Y < _grid.GeneratorY.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][4]].Point1.Y > _grid.GeneratorY.Item2)
                    {
                        point.X = (_grid.GeneratorX.Item1 - _grid.Edges[_grid.Elements[ielem][0]].Point0.X) /
                                  _grid.Elements[_grid.Elements[ielem][4]].Length;

                        point.Y = (_grid.Edges[_grid.Elements[ielem][4]].Point.Y -
                                   _grid.Edges[_grid.Elements[ielem][4]].Point0.Y) /
                                  _grid.Elements[_grid.Elements[ielem][4]].Length;
                        
                        point.Z = 0;
                        
                        for (int i = 0; i < _basis.Size; i++)
                        {
                            _globalVector[_grid.Elements[ielem][i]] += _grid.SourceValue * _basis.GetPsi(i, point).X;
                        }
                    }
                    
                    // 4
                    if (_grid.Edges[_grid.Elements[ielem][1]].Point0.X > _grid.GeneratorX.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][1]].Point0.X < _grid.GeneratorX.Item2 &&
                        _grid.Edges[_grid.Elements[ielem][5]].Point0.Y < _grid.GeneratorY.Item1 &&
                        _grid.Edges[_grid.Elements[ielem][5]].Point1.Y > _grid.GeneratorY.Item2)
                    {
                        point.X = (_grid.GeneratorX.Item2 - _grid.Edges[_grid.Elements[ielem][5]].Point0.X) /
                                  _grid.Elements[_grid.Elements[ielem][5]].Length;

                        point.Y = (_grid.Edges[_grid.Elements[ielem][5]].Point.Y -
                                   _grid.Edges[_grid.Elements[ielem][5]].Point0.Y) /
                                  _grid.Elements[_grid.Elements[ielem][5]].Length;
                        
                        point.Z = 0;
                        
                        for (int i = 0; i < _basis.Size; i++)
                        {
                            _globalVector[_grid.Elements[ielem][i]] += _grid.SourceValue * _basis.GetPsi(i, point).X;
                        }
                    }
                }
                
                //if (_grid.Edges[_grid.Elements[ielem][0]].Point0.X >= -51 &&
                //    _grid.Edges[_grid.Elements[ielem][0]].Point0.X <= 51 &&
                //    _grid.Edges[_grid.Elements[ielem][0]].Point0.Y >= -51 &&
                //    _grid.Edges[_grid.Elements[ielem][0]].Point0.Y <= 51 &&
                //    _grid.Edges[_grid.Elements[ielem][0]].Point0.Z <= 51 &&
                //    _grid.Edges[_grid.Elements[ielem][11]].Point1.Z >= 30)
                //{
                //    for (int i = 0; i < _basis.Size; i++)
                //    {
                //        _globalVector[_grid.Elements[ielem][i]] = 1;
                //    }
                //}

                break;
            
            case Scheme.Two_layer_Implicit:
                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj1[i] += _elementAssembler.MassMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }

                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _elementAssembler.LocalVector[i] += SchemeUsage(ielem, itime, _activeScheme, 1) * qj1[i];

                    _globalVector[_grid.Elements[ielem][i]] += _elementAssembler.LocalVector[i];
                }
                
                break;
            
            case Scheme.Three_layer_Implicit:

                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj2[i] += _elementAssembler.MassMatrix[i, j] * _layers[1][_grid.Elements[ielem][j]];
                        qj1[i] += _elementAssembler.MassMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }
                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _elementAssembler.LocalVector[i] += SchemeUsage(ielem, itime, _activeScheme, 1) * qj2[i];
                    _elementAssembler.LocalVector[i] += SchemeUsage(ielem, itime, _activeScheme, 2) * qj1[i];
                    
                    _globalVector[_grid.Elements[ielem][i]] += _elementAssembler.LocalVector[i];
                }
                
                break;
            
            case Scheme.Four_layer_Implicit:

                for (int i = 0; i < _basis.Size; i++)
                {
                    for (int j = 0; j < _basis.Size; j++)
                    {
                        qj3[i] += _elementAssembler.MassMatrix[i, j] * _layers[2][_grid.Elements[ielem][j]];
                        qj2[i] += _elementAssembler.MassMatrix[i, j] * _layers[1][_grid.Elements[ielem][j]];
                        qj1[i] += _elementAssembler.MassMatrix[i, j] * _layers[0][_grid.Elements[ielem][j]];
                    }
                }
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    _elementAssembler.LocalVector[i] += SchemeUsage(ielem, itime, _activeScheme, 1) * qj3[i];
                    _elementAssembler.LocalVector[i] += SchemeUsage(ielem, itime, _activeScheme, 2) * qj2[i];
                    _elementAssembler.LocalVector[i] += SchemeUsage(ielem, itime, _activeScheme, 3) * qj1[i];
                    
                    _globalVector[_grid.Elements[ielem][i]] += _elementAssembler.LocalVector[i];
                }
                
                break;
        }
    }

    private void AddElement(int i, int j, double value)
    {
        if (_globalMatrix.Symmetric)
        {
            if (i == j)
            {
                _globalMatrix.Di[i] += value;
                return;
            }

            for (int icol = _globalMatrix.Ig[i]; icol < _globalMatrix.Ig[i + 1]; icol++)
            {
                if (_globalMatrix.Jg[icol] == j)
                {
                    _globalMatrix.Gg[icol] += value;
                    return;
                }
            }
        }
        
        else
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
    }

    private void AssemblyLocalElement(int ielem, int itime)
    {
        for (int i = 0; i < _basis.Size; i++)
        {
            _elementPointList[i] = _grid.Edges[_grid.Elements[ielem][i]].Point;
        }
        
        double hx = _grid.Edges[_grid.Elements[ielem][0]].Length;
        double hy = _grid.Edges[_grid.Elements[ielem][4]].Length;
        double hz = _grid.Edges[_grid.Elements[ielem][8]].Length;

        _elementAssembler.Element.SetElement(_elementPointList, hx, hy, hz,
            _grid.GetSigma(_grid.Edges[_grid.Elements[ielem][3]].Point));

        _elementAssembler.AnalyticalAssembly(_timeGrid[itime], _test.F);
        _elementAssembler.StiffnessMatrix = 1 / _grid.Mu * _elementAssembler.StiffnessMatrix;

        //for (int i = 0; i < _basis.Size; i++)
        //{
        //    for (int j = 0; j < _basis.Size; j++)
        //    {
        //        Func<Point3D, double> kek;
        //        Vector3D psi1 = new(0, 0, 0);
        //        Vector3D psi2 = new(0, 0, 0);
        //        Vector3D dPsi1 = new(0, 0, 0);
        //        Vector3D dPsi2 = new(0, 0, 0);

        //        int ik = i;
        //        int jk = j;
        //        kek = point =>
        //        {
        //            psi1.Copy(_basis.GetPsi(ik, point));
        //            psi2.Copy(_basis.GetPsi(jk, point));

        //            return psi1 * psi2;
        //        };

        //        _massMatrix[i, j] += hx * hy * hz * _integration.Gauss3D(kek);

        //        kek = point =>
        //        {
        //            dPsi1.Copy(_basis.GetDPsi(ik, point));
        //            dPsi2.Copy(_basis.GetDPsi(jk, point));

        //            return Vector3D.DotProductJacob(dPsi1, dPsi2, hx, hy, hz);
        //        };

        //        _stiffnessMatrix[i, j] += 1 / _grid.Mu * _integration.Gauss3D(kek);
        //    }

        //    _localVector[i] = _test.F(_grid.Edges[_grid.Elements[ielem][i]].Point, _timeGrid[itime], i, _grid.GetSigma(
        //        new Point3D(_grid.Edges[_grid.Elements[ielem][0]].Point.X,
        //            _grid.Edges[_grid.Elements[ielem][3]].Point.Y,
        //            _grid.Edges[_grid.Elements[ielem][11]].Point.Z)));
        //}

        //_localVector = _massMatrix * _localVector;
    }
    
    private void AccountDirichletBoundaries(int itime)
    {
        if (_globalMatrix.Symmetric)
        {
            foreach (var edge in _grid.DirichletBoundaries)
            {
                _globalMatrix.Di[edge] = 1;
                var value = _test.UValue(_grid.Edges[edge].Point, _timeGrid[itime], _grid.Edges[edge].GetAxis());
                _globalVector[edge] = value;

                for (int i = _globalMatrix.Ig[edge]; i < _globalMatrix.Ig[edge + 1]; i++)
                {
                    _globalVector[_globalMatrix.Jg[i]] -= value * _globalMatrix.Gg[i];
                    _globalMatrix.Gg[i] = 0;
                }

                for (int i = edge + 1; i < _globalMatrix.Size; i++)
                {
                    for (int j = _globalMatrix.Ig[i]; j < _globalMatrix.Ig[i + 1]; j++)
                    {
                        if (_globalMatrix.Jg[j] == edge)
                        {
                            _globalVector[i] -= value * _globalMatrix.Gg[j];
                            _globalMatrix.Gg[j] = 0;
                        }
                    }
                }
            }
        }

        else
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

        _globalMatrix = new(_grid.Edges.Length, count, true);
        _globalVector = new(_grid.Edges.Length);
        _solution = new(_grid.Edges.Length);
        _solution1 = new(_grid.Edges.Length);
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
        double errorX = 0;
        double errorY = 0;
        double errorZ = 0;
        //for (int i = 0; i < _grid.Edges.Length; i++)
        //{
        //    error += Math.Pow(
        //        _test.UValue(_grid.Edges[i].Point, _timeGrid[itime], _grid.Edges[i].GetAxis()) - _solution[i], 2);
        //}
        //
        Point3D center = new Point3D(0, 0, 0);
        Point3D test = new Point3D(0, 0, 0);
        Vector3D testVector = new Vector3D(0, 0, 0);
        Vector3D approx;
        
        double x = 0;
        double y = 0;
        double z = 0;

        double approxX = 0;
        double approxY = 0;
        double approxZ = 0;

        for (int i = 0; i < _grid.Elements.Length; i++)
        {
            center.X = _grid.Edges[_grid.Elements[i][3]].Point.X;
            center.Y = _grid.Edges[_grid.Elements[i][7]].Point.Y;
            center.Z = _grid.Edges[_grid.Elements[i][11]].Point.Z;
            
            test.X = _test.UValue(center, _timeGrid[itime], 0);
            test.Y = _test.UValue(center, _timeGrid[itime], 1);
            test.Z = _test.UValue(center, _timeGrid[itime], 2);
            testVector.UpdateVector(test.X, test.Y, test.Z); 
            
            approx = GetValue(center);

            x += Math.Pow(approx.X - test.X, 2);
            y += Math.Pow(approx.Y - test.Y, 2);
            z += Math.Pow(approx.Z - test.Z, 2);

            approxX += Math.Pow(approx.X, 2);
            approxY += Math.Pow(approx.Y, 2);
            approxZ += Math.Pow(approx.Z, 2);

            error += Math.Pow(testVector.GetLength() - approx.GetLength(), 2);
        }

        if (itime % 2 == 1)
        {
            Console.Write("\\rowcolor{gray!15} ");
        }

        _sum += Math.Sqrt(error / _grid.Elements.Length);
        Console.WriteLine($"${_timeGrid[itime]}$ & ${Math.Sqrt(error / _grid.Elements.Length)}$ \\\\ \\hline");
    }
    
    private int FindElement(Point2D point, Grid2D grid)
    {
        int rSteps, zSteps;

        if (point.X < grid.Nodes[0].X || point.Y < grid.Nodes[0].Y)
        {
            return -1;
        }

        for (rSteps = 1; rSteps < grid.SumSteps[0]; rSteps++)
        {
            if (point.X <= grid.Nodes[rSteps].X)
            {
                break;
            }
        }

        if (rSteps == grid.SumSteps[0] + 1)
        {
            return -1;
        }
        
        for (zSteps = 1; zSteps < grid.SumSteps[1] + 1; zSteps++)
        {
            if (point.Y <= grid.Nodes[zSteps * (grid.SumSteps[0] + 1)].Y)
            {
                var result = (zSteps - 1) * (grid.SumSteps[0]) + rSteps - 1;
                return result;
            }
        }
        return -1;
    }
    
    public double GetValue(Point2D point, Grid2D grid)
    {
        double result = 0;
        var kek = new Point2D(0, 0);

        var iElem = FindElement(point, grid);
        var elem = grid.Elements[iElem];

        kek.X = (point.X - grid.Nodes[elem[0]].X) / (grid.Nodes[elem[3]].X - grid.Nodes[elem[0]].X);
        kek.Y = (point.Y - grid.Nodes[elem[0]].Y) / (grid.Nodes[elem[3]].Y - grid.Nodes[elem[0]].Y);
        
        for (int i = 0; i < _basis2D.Size; i++)
        {
            result += _basis2D.GetPsi(i, kek) * _2DSolution[elem[i]];
        }
                

        return result;
    }

    public Vector3D GetValue(Point3D point)
    {
        var vector = new Vector3D(0, 0, 0);
        var kek = new Point3D(0, 0, 0);

        foreach (var elem in _grid.Elements)
        {
            if (point.X >= _grid.Edges[elem[0]].Point0.X && point.X <= _grid.Edges[elem[11]].Point1.X &&
                point.Y >= _grid.Edges[elem[0]].Point0.Y && point.Y <= _grid.Edges[elem[11]].Point1.Y &&
                point.Z >= _grid.Edges[elem[0]].Point0.Z && point.Z <= _grid.Edges[elem[11]].Point1.Z)
            {
                kek.X = (point.X - _grid.Edges[elem[0]].Point0.X) / _grid.Edges[elem[0]].Length;
                kek.Y = (point.Y - _grid.Edges[elem[4]].Point0.Y) / _grid.Edges[elem[4]].Length;
                kek.Z = (point.Z - _grid.Edges[elem[8]].Point0.Z) / _grid.Edges[elem[8]].Length;
                
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
                int stepsX = 8;
                int stepsY = 8;

                double hx = (_grid.Edges[_grid.Elements[ielem][11]].Point1.X -
                      _grid.Edges[_grid.Elements[ielem][0]].Point0.X) / stepsX;
                double hy = (_grid.Edges[_grid.Elements[ielem][11]].Point1.Y -
                      _grid.Edges[_grid.Elements[ielem][0]].Point0.Y) / stepsY;

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
        kek.Y = (point.Y - _grid.Edges[_grid.Elements[ielem][4]].Point0.Y) / _grid.Edges[_grid.Elements[ielem][4]].Length;
        kek.Z = (point.Z - _grid.Edges[_grid.Elements[ielem][8]].Point0.Z) / _grid.Edges[_grid.Elements[ielem][8]].Length;
        
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

    private double GetValueForRotAz(Point3D point)
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
                kek.Y = (point.Y - _grid.Edges[elem[4]].Point0.Y) / _grid.Edges[elem[4]].Length;
                kek.Z = (point.Z - _grid.Edges[elem[8]].Point0.Z) / _grid.Edges[elem[8]].Length;
                
                for (int i = 0; i < _basis.Size; i++)
                {
                    if (i is >= 4 and <= 7)
                    {
                        continue;
                    }
            
                    vector += _basis.GetDPsi(i, kek) * _solution[elem[i]];
                }
                
                break;
            }
        }

        return vector.Y - vector.X;
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
            if (i is >= 8 and <= 11)
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
}
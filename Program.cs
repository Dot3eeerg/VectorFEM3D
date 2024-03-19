using VectorFEM3D;

Grid grid = new Grid("GridParameters");
grid.BuildGrid();
grid.AccountBoundaryConditions();

//TimeGrid timeGrid = new TimeGrid("TimeGridParameters");
//timeGrid.BuildTimeGrid();

GeneratedTimeGrid timeGrid = new GeneratedTimeGrid("time_s.txt");

FEM fem = new FEM(grid, timeGrid);

fem.SetTest(new Test1(grid));

//fem.SetSolver(new BCGSTABLUSolver(1e-16, 1000));
fem.SetSolver(new BCGSTABSolver(1e-14, 5000));
//fem.SetSolver(new LUSolver());

fem.SetScheme(Scheme.Natural);
//fem.SetScheme(Scheme.Two_layer_Implicit);
//fem.SetScheme(Scheme.Three_layer_Implicit);
//fem.SetScheme(Scheme.Four_layer_Implicit);

fem.Compute();

//fem.GetValue(new Point3D(3, 2, -2));
//fem.GetValue(new Point3D(1.18, 1.18, 1.18));

var kek = fem.GetValue(new Point3D(1.18, 1.18, 1.18));
Console.WriteLine($"{kek.X}, {kek.Y}, {kek.Z}");

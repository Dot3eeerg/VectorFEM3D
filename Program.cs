using VectorFEM3D;

Grid grid = new Grid("GridParameters");
grid.BuildGrid();
grid.AccountBoundaryConditions();

TimeGrid timeGrid = new TimeGrid("TimeGridParameters");
timeGrid.BuildTimeGrid();

Grid2D grid2D = new Grid2D("GridParameters2D");
grid2D.BuildGrid(true);

//GeneratedTimeGrid timeGrid = new GeneratedTimeGrid("time_s.txt");

FEM fem = new FEM(grid, timeGrid);

fem.SetTest(new Test1(grid));

fem.SetSolver(new LOSLTSolver(1e-39, 5000));

fem.SetScheme(Scheme.Natural);
//fem.SetScheme(Scheme.Two_layer_Implicit);
//fem.SetScheme(Scheme.Three_layer_Implicit);
//fem.SetScheme(Scheme.Four_layer_Implicit);

fem.Convert2DSolution(grid2D);

fem.Compute();

namespace VectorFEM3D;

public enum ElementSide
{
    Left,
    Right,
    Bottom,
    Upper,
    Rear,
    Front
}

public enum Scheme
{
    Natural,
    Two_layer_Implicit,
    Three_layer_Implicit,
    Four_layer_Implicit
}
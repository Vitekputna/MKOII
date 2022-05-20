#pragma once
#include <vector>

typedef std::vector<std::vector<double>> field2;

class face
{
    public:
    std::vector<double> n = std::vector<double>(2,0);
    std::vector<double> s = std::vector<double>(2,0);
    std::vector<int> nodes = std::vector<int>(2,0);
    double S;

    face(std::vector<double> a, std::vector<double> b);
};

class cell
{
    public:
    std::vector<double> X,Y;
    std::vector<face> faces;
    std::vector<int> dirichlet_faces;
    std::vector<int> neumann_faces;
    std::vector<int> cell_faces;

    double x = 0,y = 0,V;

    cell(std::vector<std::vector<double>> nodes);
};

typedef std::vector<std::vector<cell>> cell_field2;

class mesh
{  
    public:  
    cell_field2 cells;

    mesh();
    void export_mesh(std::string name);

    //private:
    void load_mesh(std::string path,field2& nodes,field2& edges,field2& quads);
    void sort_mesh(field2 const& nodes,field2 const& edges,field2 const& quads);    
    std::vector<int> cell_boundary_check(std::vector<double> edge, std::vector<double> quad);
    void set_boundary(field2 const& edges, field2 const& quads, std::vector<cell> &unsorted_cells);
};
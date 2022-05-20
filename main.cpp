#include <iostream>
#include <vector>
#include "mesh.h"

int main()
{
    mesh K;
    //K.export_mesh("mesh.txt");

    for(int j = 0; j < K.cells.size(); j++)
    {
        for(int i = 0; i < K.cells[0].size(); i++)
        {
            std::cout << "cell " << i << " " << j << " -> dirichlet face ";
            for(int k = 0; k < K.cells[j][i].dirichlet_faces.size(); k++)
            {
                std::cout << K.cells[j][i].dirichlet_faces[k] << " ";
            }
            std::cout << std::endl;
        }
    }
}


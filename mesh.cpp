#include <fstream>
#include <vector>
#include <iostream>
#include "mesh.h"
#include "math.h"

face::face(std::vector<double> a, std::vector<double> b)
{
    double sx,sy;

    sx = b[0]-a[0];
    sy = b[1]-a[1];

    

    S = sqrt(sx*sx + sy*sy);

    n = {sy/S,-sx/S};
    s = {sx,sy};

    //std::cout << n[0] << " " << n[1] << " " << S << "\n";
}

cell::cell(std::vector<std::vector<double>> nodes)
{

    for(auto const& node : nodes)
    {
        X.push_back(node[0]);
        Y.push_back(node[1]);

        x += 0.25*node[0];
        y += 0.25*node[1];
    }

    this->faces.push_back(face(nodes[0],nodes[1]));
    this->faces.push_back(face(nodes[1],nodes[2]));
    this->faces.push_back(face(nodes[2],nodes[3]));
    this->faces.push_back(face(nodes[3],nodes[0]));

    V = 0.5*abs(faces[0].s[0] * faces[1].s[1] - faces[0].s[1] * faces[1].s[0]) + 0.5*abs(faces[2].s[0] * faces[3].s[1] - faces[2].s[1] * faces[3].s[0]);

    //std:: cout << V << "\n";
    
}

mesh::mesh()
{
    std::vector<std::vector<double>> nodes, edges, quads;

    load_mesh("block.mesh",nodes,edges,quads);
    sort_mesh(nodes,edges,quads);
    
}

void mesh::load_mesh(std::string path,field2& nodes,field2& edges,field2& quads)
{
    std::ifstream file(path);
	
	if(file.fail())
	{
		std::cout << "//////////////////////////////////// \n";
		std::cout << "File not found \n";
		std::cout << "//////////////////////////////////// \n";
		throw std::exception();
	}

	std::vector<std::string> text_vec;
	std::string text;

	// Reading as text file
	while (std::getline(file, text))
	{
		text_vec.push_back(text);
	}

	int num_nodes = stoi(text_vec[4]);

	// Reading nodes
	std::string word = "";
	double number;
	char k;
	std::vector<double> row;

	//std::cout << "NODES" << std::endl;
	for (unsigned int i = 5; i < num_nodes+5; i++)
	{
		for (unsigned int j = 0; j < text_vec[i].length(); j++)
		{
			if (text_vec[i][j] != 32)
			{
				while (text_vec[i][j] != 32 && j < text_vec[i].length())
				{
					k = text_vec[i][j];
					word.push_back(k);
					j++;
				}
				number = (std::stod(word));
				row.push_back(number);
				//std::cout << row.back() << " ";
				word = "";
			}
		}
		nodes.push_back(row);
		row.clear();
		//std::cout << "\n";
	}

	// reading edge 
	//std::cout << "EDGES" << std::endl;
	int edge_start = num_nodes + 6;
	int num_edges = stoi(text_vec[edge_start]);

	//std::cout << num_edges << "\n";

	for (unsigned int i = edge_start+1; i < num_edges + edge_start+1; i++)
	{
		for (unsigned int j = 0; j < text_vec[i].length(); j++)
		{
			if (text_vec[i][j] != 32)
			{
				while (text_vec[i][j] != 32 && j < text_vec[i].length())
				{
					k = text_vec[i][j];
					word.push_back(k);
					j++;
				}
				number = std::stod(word);
				row.push_back(number);
				//std::cout << row.back() << " ";
				word = "";
			}
		}
		edges.push_back(row);
		row.clear();
		//std::cout << "\n";
	}

	// reading quads
	//std::cout << "QUADS" << std::endl;
	int quads_start = num_edges + edge_start + 2;
	int num_quads = stoi(text_vec[quads_start]);

	//std::cout << num_quads << "\n";

	for (unsigned int i = quads_start+1; i < quads_start + num_quads +1; i++)
	{
		for (unsigned int j = 0; j < text_vec[i].length(); j++)
		{
			if (text_vec[i][j] != 32)
			{
				while (text_vec[i][j] != 32 && j < text_vec[i].length())
				{
					k = text_vec[i][j];
					word.push_back(k);
					j++;
				}
				number = std::stod(word);
				row.push_back(number);
				//std::cout << row.back() << " ";
				word = "";
			}
		}
		quads.push_back(row);
		row.clear();
		//std::cout << "\n";
	}
}

std::vector<int> mesh::cell_boundary_check(std::vector<double> edge, std::vector<double> quad)
{
    std::vector<bool> ans(2,false);
    std::vector<int> idx(2,0);

	for(unsigned int e = 0; e < edge.size()-1; e++)
	{
		for(unsigned int i = 0; i < quad.size()-1; i++)
		{
			if(quad[i] == edge[e])
			{
				ans[e] = true;
                idx[e] = i;
			}
		}
	}

	if(ans[0] && ans[1])
	{
		return idx;
	}
	else
	{
		return std::vector<int>{};
	}
}

void mesh::set_boundary(field2 const& edges, field2 const& quads, std::vector<cell> &unsorted_cells)
{
    std::vector<int> k;

	for(auto const& edge : edges)
	{
        for(unsigned int i = 0; i < unsorted_cells.size();i++)
        {

            k = cell_boundary_check(edge,quads[i]);

            if(k.size()) // find cell containing edge nodes
            {
                
                unsorted_cells[i].dirichlet_faces.push_back(k[0]);

            }
        }
	}
}

void mesh::sort_mesh(field2 const& nodes,field2 const& edges,field2 const& quads)
{

    // creating unsorted cells
	std::vector<cell> unsorted_cells;

	for (auto const& q : quads)
    {
        unsorted_cells.push_back(cell({nodes[q[0]-1],nodes[q[1]-1],nodes[q[2]-1],nodes[q[3]-1]}));     
    }

    set_boundary(edges,quads,unsorted_cells);

    int ny,nx;

	int n = quads.size();
	int N = nodes.size();

	nx = std::max((-1-n+N+sqrt((1+n-N)*(1+n-N) -4*n))/2,(-1-n+N-sqrt((1+n-N)*(1+n-N) -4*n))/2);
	ny = n/nx;

	std::cout << nx << " " << ny << "\n";

	std::vector<cell> row;

    //cells.push_back({});

	int k;
	for(unsigned int j = 0; j < ny; j++)
	{
        cells.push_back({});
		for(unsigned int i = 0; i < nx; i++)
		{
			k = (nx-1)*ny-i*ny+j;

			//row.push_back(unsorted_cells[k]);

            cells[j].push_back(unsorted_cells[k]);
		}
		
		//row.clear();
	}


}

void mesh::export_mesh(std::string name)
{
    std::ofstream f(name);
	f << "x coordinates:" << "\n";
	for (unsigned int i = 0; i < cells.size(); i++)
	{
		for (unsigned int j = 0; j < cells[0].size(); j++)
		{
			f << cells[i][j].x << " ";
		}
		f << "\n";
	}

	f << "\n";
	f << "y coordinates:" << "\n";
	for (unsigned int i = 0; i < cells.size(); i++)
	{
		for (unsigned int j = 0; j < cells[0].size(); j++)
		{
			f << cells[i][j].y << " ";
		}
		f << "\n";
	}

	f << "\n";
	f << "x normal component:" << "\n";
	for(unsigned int k = 0; k < 4;k++)
	{
		for (unsigned int i = 0; i < cells.size(); i++)
		{
			for (unsigned int j = 0; j < cells[0].size(); j++)
			{
				f << cells[i][j].faces[k].n[0] << " ";
			}
			f << "\n";
		}
	}

	f << "\n";
	f << "y normal component:" << "\n";
	for(unsigned int k = 0; k < 4;k++)
	{
		for (unsigned int i = 0; i < cells.size(); i++)
		{
			for (unsigned int j = 0; j < cells[0].size(); j++)
			{
				f << cells[i][j].faces[k].n[1] << " ";
			}
			f << "\n";
		}
	}

	f << "\n";
	f << "wall length:" << "\n";
	for(unsigned int k = 0; k < 4;k++)
	{
		for (unsigned int i = 0; i < cells.size(); i++)
		{
			for (unsigned int j = 0; j < cells[0].size(); j++)
			{
				f << cells[i][j].faces[k].S << " ";
			}
			f << "\n";
		}
	}

	f << "\n";
	f << "cell volume:" << "\n";
	for (unsigned int i = 0; i < cells.size(); i++)
	{
		for (unsigned int j = 0; j < cells[0].size(); j++)
		{
			f << cells[i][j].V << " ";
		}
		f << "\n";
	}
}




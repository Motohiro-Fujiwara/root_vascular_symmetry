/*
 *  cellform_vascul_cycle-initial2.cpp
 *  Created by Motohiro Fujiwara on 2022/12/06.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <iomanip>
using namespace std;

double uniform(void);
double rand_gamma( double theta, double kappa );
double Uniform(void){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}
double rand_normal(double mu, double sigma){
    double z=sqrt(-2.0*log(Uniform())) * sin(2.0*M_PI*Uniform());
    return mu + sigma*z;
}

const int pmax=10;//steady-state
const int qmax=10000;//eachhan13000//step wt10000 han11000
const int tmax=1;
const double tf=0.0416667;//1/24
const int ni=20;//gusunomi;
const int ni_2=ni*0.5;//;cell number
const double k_1=0.0003;//Area
const double k_2=0;//perimeter
const double k_3=0.00003;//apical
const double k_4=0.0001;//basal
const double k_5=0.0001;//lateral

class Vertex {
public:
	double coordinate[2];//x- y- coordinate 0 is origin
	double coordinate_h[2];//x- y- coordinate geometric center is origin
	vector<int> NeighborCell;
	vector<int> NeighborVertex;
	int Edge_Vertex[2];
	int Edge_Area[2];
	int T1_flag;
	double Fx, Fy;
};

class Cell  {
public:
	int i,j;
	vector<int> GroupVertex;
	vector<int> GroupVertex_divi;
	vector<Vertex> GroupVertex_h;
    vector<Vertex> GroupVertex_M;
	vector<int> GroupEdge;
	vector<int> GroupEdge_divi;
	vector<double> R_edge;
	vector<int> GroupVertex_T1;
	vector<int> GroupEdge_T1;
	vector<double> R_edge_T1;
	vector<double> af;
	
	double k_area;
	double k_peri;
	double k_api;
	double k_bas;
	double k_lat;
	int cell_type;
	int cell_lineage;
	int cell_defect;
	double Perimeter;
	double Area;
	double Area0;
    double Area0_0;
	double theta_mo;
	double Mxx;
	double Myy;
	double Mxy;
    double Sxx;
    double Syy;
    double Sxy;
    double S1;
    double S2;
    double Theta1;
    double Theta2;
	double time_division;
	int division_count;
	int division_xylem;
	int division_phloem;
	double centroid[2];//geometric center coordinate
	double cent_dr[2];//geometric center coordinate difference
	void firstset_cell();
};

class Tissue {
public:
	int i,imax,j,jmax,k,kmax,l,lmax;
	double dt;
	double E; //energy
	int cell_p1;	
	int cell_p2;
	double Area_divi;
	double Area_divimx;
	double Area_divise;
	double Area_divipr;
	double energy_bf;
	double energy_af;
	
	vector<Cell> Cellnum;
	vector<Vertex> Vertices;
	vector<Vertex> Edge;
	
	void firststage();
	void cellgroup();
	void VertexDynamics();
	void firstset_cell();
	void Geometory_cell();
	void T1_cell();
	int T1_count;
	void Peri_cell();
	void Division(int p);
	void mutation(int p);
	ofstream ofs[10000];
	
	void output(int p);
	void output_result();
	void output_barometer(int c);
	ofstream ofs_result[11];
	ofstream ofs_barometer[13];
};

void Tissue::firststage() {
	//Vertex
	i=0;
	ifstream ifs0("/datapass/initial2/vertex_initial.dat");
	string reading_line0;
	// read by line
	while (getline(ifs0, reading_line0)) {
		Vertices.push_back(Vertex());//vector add
		sscanf(reading_line0.data(), "%lf	%lf", &Vertices[i].coordinate[0], &Vertices[i].coordinate[1]);
		i++;
	}
	
	i=0;
	ifstream ifs1("/datapass/initial2/neigborcell_initial.dat");
	string reading_line1;
	// read by line
	while (getline(ifs1, reading_line1)) {
		Vertices[i].NeighborCell.push_back(int());//add
		Vertices[i].NeighborCell.push_back(int());//add
		Vertices[i].NeighborCell.push_back(int());//add
		sscanf(reading_line1.data(), "%d	%d	%d", &Vertices[i].NeighborCell[0], &Vertices[i].NeighborCell[1],&Vertices[i].NeighborCell[2]);
		i++;
	}
	
	i=0;
	ifstream ifs2("/datapass/initial2/neigborvertex_initial.dat");
	string reading_line2;
	// read by line
	while (getline(ifs2, reading_line2)) {
		Vertices[i].NeighborVertex.push_back(int());//add
		Vertices[i].NeighborVertex.push_back(int());//add
		Vertices[i].NeighborVertex.push_back(int());//add
		sscanf(reading_line2.data(), "%d	%d	%d", &Vertices[i].NeighborVertex[0], &Vertices[i].NeighborVertex[1],&Vertices[i].NeighborVertex[2]);
		i++;
	}
	
	//Edge
	i=0;
	ifstream ifs3("/datapass/initial2/edge-vertex_initial.dat");
	string reading_line3;
	// read by line
	while (getline(ifs3, reading_line3)) {
		Edge.push_back(Vertex());//add
		sscanf(reading_line3.data(), "%d	%d", &Edge[i].Edge_Vertex[0], &Edge[i].Edge_Vertex[1]);
		i++;
	}
	
	i=0;
	ifstream ifs4("/datapass/initial2/edge-cell_initial.dat");
	string reading_line4;
	// read by line
	while (getline(ifs4, reading_line4)) {
		sscanf(reading_line4.data(), "%d	%d", &Edge[i].Edge_Area[0], &Edge[i].Edge_Area[1]);
		i++;
	}
}

void Tissue::cellgroup() {
	//Vertex
	i=0;
	const char delimiter ='	';
	ifstream ifs5("/datapass/initial2/cell-vertex_initial.dat");
	string reading_line5;
	// read by line
	while (getline(ifs5, reading_line5)) {
		j=0;
		Cellnum.push_back(Cell());
		string separated_string;
		istringstream line_separater(reading_line5);
		while (getline(line_separater, separated_string, delimiter)) {
			Cellnum[i].GroupVertex.push_back(0);
			sscanf(separated_string.data(), "%d", &Cellnum[i].GroupVertex[j]);
			j++;
		}
		i++;
	}
	
	//Edge
	i=0;
	ifstream ifs6("/datapass/vascularsimu/initial2/cell-edge_initial.dat");
	string reading_line6;
	// read by line
	while (getline(ifs6, reading_line6)) {
		j=0;
		string separated_string;
		istringstream line_separater(reading_line6);
		while (getline(line_separater, separated_string, delimiter)) {
			Cellnum[i].GroupEdge.push_back(0);
			sscanf(separated_string.data(), "%d", &Cellnum[i].GroupEdge[j]);
			j++;
		}
		i++;
	}
	
}

void Tissue::firstset_cell() {
	//cell_type
	for (i=0; i<3; i++) {
		Cellnum[i].cell_type=1;//xylem
	}
	for (i=3; i<5; i++) {
		Cellnum[i].cell_type=2;//PSE
		//Cellnum[i].cell_type=3;//PSE_mutant
	}
	//6 meta-SE
	for (i=5; i<9; i++) {
		//Cellnum[i].cell_type=3;//procambium
        Cellnum[i].cell_type=7;//companion cell
	}
	for (i=9; i<18; i++) {
		Cellnum[i].cell_type=4;//pericycle
	}
	for (i=18; i<26; i++) {
		Cellnum[i].cell_type=5;//endodermis
	}
	Cellnum[26].cell_type=0;//boundary
	
	//cell lineage _defect
	imax=Cellnum.size();
	for (i=0; i<imax; i++) {
			Cellnum[i].cell_lineage=i;
			Cellnum[i].cell_defect=0;
	}
	
	//parameter_euler
	dt=3.0;
	
	imax=Cellnum.size();
	for (i=0; i<imax; i++) {
		Cellnum[i].k_area=k_1;
		Cellnum[i].k_peri=k_2;
		Cellnum[i].k_api=k_3;
		Cellnum[i].k_bas=0;
		Cellnum[i].k_lat=0;
	}
	
	for (i=0; i<3; i++) {
		Cellnum[i].k_area=k_1*1.0;
		Cellnum[i].k_peri=k_2*1;
		Cellnum[i].k_api=k_3*1;//xylem
	}
	for (i=3; i<5; i++) {
		Cellnum[i].k_area=k_1*1;
		Cellnum[i].k_peri=k_2*1;
		Cellnum[i].k_api=k_3*1;//PSE
	}
	for (i=5; i<9; i++) {
		Cellnum[i].k_area=k_1*1;
		Cellnum[i].k_peri=k_2*1;
		Cellnum[i].k_api=k_3*1;//procambium
	}
	for (i=9; i<18; i++) {
		Cellnum[i].k_area=k_1*1.0;
		Cellnum[i].k_peri=k_2*1;
		Cellnum[i].k_api=k_3*1;//pericycle
	}
	for (i=18; i<26; i++) {
		Cellnum[i].k_area=k_1*1.0;
		Cellnum[i].k_peri=k_2*1;
		Cellnum[i].k_api=k_3*1;//endodermis
	}	
	
	//back cell
	Cellnum[26].k_area=0.0;
	Cellnum[26].k_api=0.0;
	Cellnum[26].k_bas=0.0;
	Cellnum[26].k_lat=0.0;
	Cellnum[26].k_peri=0.0;//0.0001;
	
	imax=Edge.size();
	for (i=0; i<imax; i++) {
		Edge[i].k_au=k_au;
	}
	jmax=Cellnum[26].GroupVertex.size();	
	for (j=0; j<jmax; j++) {
			Edge[Cellnum[26].GroupEdge[j]].k_au=0;
	}
	for (i=0; i<imax; i++) {
		Edge[i].k_cy=k_cy;
	}
	jmax=Cellnum[26].GroupVertex.size();
	for (j=0; j<jmax; j++) {
			Edge[Cellnum[26].GroupEdge[j]].k_cy=0;
	}
	for (i=0; i<imax; i++) {
		Edge[i].k_aucy=k_aucy;
	}
	jmax=Cellnum[26].GroupVertex.size();
	for (j=0; j<jmax; j++) {
			Edge[Cellnum[26].GroupEdge[j]].k_aucy=0;
	}
	
	cout << imax << endl;//cell number
	
	//Geometory
	//Area
	imax=Cellnum.size();
	for (i=0; i<imax; i++) {
		Cellnum[i].Area0=0;							  
	}
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for (j=0; j<jmax; j++) {
			Cellnum[i].Area0+=0.5*(Vertices[Cellnum[i].GroupVertex[(j)%jmax]].coordinate[0]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]								  
							-Vertices[Cellnum[i].GroupVertex[(j)%jmax]].coordinate[1]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]);							  
		}
	}
	
	Cellnum[1].Area0=Cellnum[1].Area0*1.8;//MX
    Cellnum[0].Area0=Cellnum[0].Area0*1.6;//MX_divi
    Cellnum[2].Area0=Cellnum[2].Area0*2.0;//MX_divi
	Cellnum[3].Area0=Cellnum[3].Area0*2.5;//SE_divi
	Cellnum[5].Area0=Cellnum[5].Area0*1.4;//
    
	Area_divi=24*1.8;//ave50 procambium
	Area_divimx=24*3.2;//mx
	Area_divise=24*1.8;//pse
	Area_divipr=24*3.2;//peri
	Cellnum[26].Area0=0;	
	
	imax=Cellnum.size();
    for (i=0; i<imax; i++) {
        Cellnum[i].Area0_0=Cellnum[i].Area0;
    }
	
	for (i=9; i<18; i++) {
		   Cellnum[i].Area0_0=Cellnum[i].Area0*1.8;//pericycle
	   }
	   
	   for (i=18; i<26; i++) {
		   Cellnum[i].Area0_0=Cellnum[i].Area0*2.1;//endodermis
	   }
	
	//Edge
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for (j=0; j<jmax; j++) {
			Cellnum[i].R_edge.push_back(sqrt((Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])
											 *(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])
											 +(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])
											 *(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])));							  
		}
	}

	//centroid
	for (i=0; i<imax; i++) {
		for(k=0;k<2;k++){
			Cellnum[i].centroid[k]=0;
		}
	}
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for(j=0;j<jmax;j++){
			double gaiseki;
			gaiseki=(Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]								  
					 -Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]);
			for(k=0;k<2;k++){
				Cellnum[i].centroid[k]+=((Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[k])+(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[k]))*gaiseki;
			}
		}
		for(k=0;k<2;k++){
			Cellnum[i].centroid[k]=Cellnum[i].centroid[k]/(6.0*Cellnum[i].Area0);
		}
	}	

    //centroid-hosei
	imax=Cellnum.size();
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for (j=0; j<jmax; j++) {		
			Cellnum[i].GroupVertex_h.push_back(Vertex());
			Cellnum[i].GroupVertex_h.push_back(Vertex());
		}
	}	
	
	imax=Cellnum.size();
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex_h.size();	
		for (j=0; j<jmax; j++) {		
			Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[0]=0.0;
			Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[1]=0.0;
		}
	}
    
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[0]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Cellnum[i].centroid[0];
            Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[1]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Cellnum[i].centroid[1];
        }
    }
    
	//Mormento
    imax=Cellnum.size();
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].GroupVertex_M.push_back(Vertex());
            Cellnum[i].GroupVertex_M.push_back(Vertex());
        }
    }
    
    imax=Cellnum.size();
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex_h.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]=0.0;
            Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]=0.0;
        }
    }
    
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for (j=0; j<jmax; j++) {
			Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Cellnum[i].centroid[0];
			Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Cellnum[i].centroid[1];
		}
	}
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].af.push_back(Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
                                    -(Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]));
        }
    }
    for (i=0; i<imax; i++) {
        Cellnum[i].Mxx=0;
        Cellnum[i].Myy=0;
        Cellnum[i].Mxy=0;
    }
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].Mxx+=2*Cellnum[i].af[j%jmax]*(Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]
                                                     +Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
                                                     +Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1])*tf;
        }
    }
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].Myy+=2*Cellnum[i].af[j%jmax]*(Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]
                                                     +Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]
                                                     +Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0])*tf;
        }
    }
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].Mxy+=Cellnum[i].af[j%jmax]*(Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
                                                   +2*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]
                                                   +2*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
                                                   +Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1])*tf;
        }
    }
	
	//division initial time
	for (i=0; i<imax; i++) {
		Cellnum[i].time_division=-1;//int(100000*drand48());							  
	}
	//xylem pse
	/*for (i=0; i<9/*imax*/; i++) {
		Cellnum[i].time_division=1000+int(500*drand48());
	}*/
	//pericycle
	for (i=9; i<18; i++) {
		Cellnum[i].time_division=5000+int(3000*drand48());//8000+int(5000*drand48());//peri
	}
    //xylem
	Cellnum[1].time_division=-1;//int(12000+300*drand48());
	Cellnum[0].time_division=200+int(100*drand48());
	Cellnum[2].time_division=200+int(100*drand48());
    
	//PSE
    Cellnum[4].time_division=1000+int(500*drand48());//2000
    Cellnum[3].time_division=1000+int(500*drand48());//2100
    
	Cellnum[8].time_division=1000+int(500*drand48());//2100
	Cellnum[7].time_division=1600+int(500*drand48());//change
	Cellnum[6].time_division=1600+int(500*drand48());//change
   
	//pericycle
	Cellnum[9].time_division=7000+int(2000*drand48());//8000+int(5000*drand48());
    Cellnum[10].time_division=-1;//8000+int(5000*drand48());
    Cellnum[13].time_division=-1;//8000+int(5000*drand48());
     Cellnum[14].time_division=-1;//10000+int(5000*drand48());
	Cellnum[15].time_division=-1;

	Cellnum[26].time_division=-1;
	cell_p1=-1;	
	cell_p2=-1;
    
	cout <<  " time "<< endl;
	for (i=0; i<imax; i++) {
    cout <<     Cellnum[i].time_division << endl;
    }
	//division count
	for (i=0; i<imax; i++) {
		Cellnum[i].division_count=0;							  
	}
	//division count xylem
	for (i=0; i<imax; i++) {
		Cellnum[i].division_xylem=0;							  
	}
	
	//division count phloem
	for (i=0; i<imax; i++) {
		Cellnum[i].division_phloem=0;							  
	}
	
	//T1
	/*imax=Edge.size();
	for (i=0; i<imax; i++) {
		Edge[i].T1_flag=0;							  
	}
		T1_count=0;*/
}


void Tissue::VertexDynamics() {
	imax=Vertices.size();
	for (i=0; i<imax; i++) {
		Vertices[i].Fx=0;
		Vertices[i].Fy=0;
	}
	imax=Cellnum.size();
	double arrow;
	double expantion;
	int vertex1,vertex2;
	int area1,area2;
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();
		for (j=0; j<jmax; j++) {
			//area_elastic
			Vertices[Cellnum[i].GroupVertex[j]].Fx+=-Cellnum[i].k_area*(Cellnum[i].Area-Cellnum[i].Area0)*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[1]);//area_fx
			Vertices[Cellnum[i].GroupVertex[j]].Fy+=-Cellnum[i].k_area*(Cellnum[i].Area-Cellnum[i].Area0)*(-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]+Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[0]);//area_fy
			
			//perimeter
			Vertices[Cellnum[i].GroupVertex[j]].Fx+=-0.5*(Cellnum[i].k_peri*Cellnum[i].Perimeter*((Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0])/Cellnum[i].R_edge[j%jmax]
														  -(Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])/Cellnum[i].R_edge[(j-1+jmax)%jmax]));//perimeter_fx
			Vertices[Cellnum[i].GroupVertex[j]].Fy+=-0.5*(Cellnum[i].k_peri*Cellnum[i].Perimeter*((Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1])/Cellnum[i].R_edge[j%jmax]
														  -(Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])/Cellnum[i].R_edge[(j-1+jmax)%jmax]));//perimeter_fy
						
			//Edge
			expantion=1.0;
			//xylem-xylem
			/*if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type
				 && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==1){
				expantion=0.8;//0.5
			}*/
			//xylem-procambium
			 for (k=0; k<2;k++){
			if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==1){
			if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==3 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==6){
			expantion=1.2;
			}
			}
			}
            //xylem-pericycle
              for (k=0; k<2;k++){
             if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==1){
             if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
             expantion=1.5;
             }
             }
             }
			
            //sieve element-companion cell
          /*  for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==2){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==7){
                        expantion=0.9;
                    }
                }
            }*/
            //sieve element-procambium cell
           /* for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==2){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==5){
                        expantion=0.9;
                    }
                }
            }*/
            
            //sieve element-pericycle
             for (k=0; k<2;k++){
             if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==2){
             if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
             expantion=1.5;//1.5
             }
             }
             }
            
            //opc cell defect
			for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==6 && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_defect==1){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
                        expantion=1.5;//1.5
                    }
                }
            }
			
            /*// OPC-pericycle
            for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==6){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
                        expantion=1.5;//1.2
                    }
                }
            }*/
            
			//pericycle-pericycle
			/*if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type
				&& Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4){
				expantion=0.8;//0.8
			}
			
			//endodermis-endodermis
			if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type
			 && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==5){
				expantion=0.8;//0.8
			 }*/
			
            //pericycle-endodermis
			for (k=0; k<2;k++){
				if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==4){
					if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==5){
						expantion=1.5;//1.5
					}
				}
			}
			
			vertex1=Edge[Cellnum[i].GroupEdge[j]].Edge_Vertex[0];
			vertex2=Edge[Cellnum[i].GroupEdge[j]].Edge_Vertex[1];
			Vertices[vertex1].Fx+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[vertex1].coordinate[0]-Vertices[vertex2].coordinate[0])/Cellnum[i].R_edge[j%jmax]);
			Vertices[vertex1].Fy+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[vertex1].coordinate[1]-Vertices[vertex2].coordinate[1])/Cellnum[i].R_edge[j%jmax]);
			Vertices[vertex2].Fx+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[vertex2].coordinate[0]-Vertices[vertex1].coordinate[0])/Cellnum[i].R_edge[j%jmax]);
			Vertices[vertex2].Fy+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[vertex2].coordinate[1]-Vertices[vertex1].coordinate[1])/Cellnum[i].R_edge[j%jmax]);
		
		}
	}
	
	//perimeter_boundary
	i=26;
	jmax=Cellnum[i].GroupVertex.size();
	for (j=0; j<jmax; j++) {
		Vertices[Cellnum[i].GroupVertex[j]].Fx+=-0.5*(Cellnum[i].k_peri*Cellnum[i].Perimeter*((Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0])/Cellnum[i].R_edge[j%jmax]
													-(Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])/Cellnum[i].R_edge[(j-1+jmax)%jmax]));//perimeter_fx
		Vertices[Cellnum[i].GroupVertex[j]].Fy+=-0.5*(Cellnum[i].k_peri*Cellnum[i].Perimeter*((Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1])/Cellnum[i].R_edge[j%jmax]
													-(Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])/Cellnum[i].R_edge[(j-1+jmax)%jmax]));//perimeter_fy
	}
	//Vertices[15].Fy+=0.001;
	/*jmax=Cellnum[27].GroupVertex.size();	
	for (j=0; j<jmax; j++) {
	Vertices[Cellnum[27].GroupVertex[j]].Fx=0.0;
	Vertices[Cellnum[27].GroupVertex[j]].Fy=0.0;
	}*/
    
    
    //stress tensor
     imax=Cellnum.size();
     for (i=0; i<imax; i++) {
         Cellnum[i].Sxx=0;
         Cellnum[i].Syy=0;
         Cellnum[i].Sxy=0;
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].GroupVertex_h[j].Fx=0;
            Cellnum[i].GroupVertex_h[j].Fy=0;
     }
     }
    
    //1cell-force
    imax=Cellnum.size();
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            //area_elastic
            Cellnum[i].GroupVertex_h[j].Fx+=-Cellnum[i].k_area*(Cellnum[i].Area-Cellnum[i].Area0)*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[1]);//area_fx
           Cellnum[i].GroupVertex_h[j].Fy+=-Cellnum[i].k_area*(Cellnum[i].Area-Cellnum[i].Area0)*(-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]+Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[0]);//area_fy
            
            //perimeter
           Cellnum[i].GroupVertex_h[j].Fx+=-0.5*(Cellnum[i].k_peri*Cellnum[i].Perimeter*((Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0])/Cellnum[i].R_edge[j%jmax]
                                                                                                  -(Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])/Cellnum[i].R_edge[(j-1+jmax)%jmax]));//perimeter_fx
            Cellnum[i].GroupVertex_h[j].Fy+=-0.5*(Cellnum[i].k_peri*Cellnum[i].Perimeter*((Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1])/Cellnum[i].R_edge[j%jmax]
                                                                                                  -(Vertices[Cellnum[i].GroupVertex[(j-1+jmax)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])/Cellnum[i].R_edge[(j-1+jmax)%jmax]));//perimeter_fy
            
            //Edge-1cell-force
            expantion=1.0;
            //xylem-xylem
           /* if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type
                && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==1){
				expantion=0.8;
            }*/
            //xylem-procambium
			for (k=0; k<2;k++){
			if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==1){
			if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==3 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==6){
			expantion=1.2;
			}
			}
			}
            //xylem-pericycle
            for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==1){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
                        expantion=1.5;
                    }
                }
            }
			
            //sieve element-companion cell
          /*  for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==2){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==7){
                        expantion=0.9;
                    }
                }
            }*/
            
            //sieve element-procambium cell
           /* for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==2){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==5){
                        expantion=0.9;
                    }
                }
            }*/
            
            //sieve element-pericycle
            for (k=0; k<2;k++){
                if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==2){
                    if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
                        expantion=1.5;//1.5
                    }
                }
            }
			
			 //opc cell defect
				for (k=0; k<2;k++){
					if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==6 && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_defect==1){
						if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
							expantion=1.5;//1.2
						}
					}
				}
		
            //OPC-pericycle
				/*for (k=0; k<2;k++){
					if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==6){
						if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
							expantion=1.5;//1.2
							}
						}
					}*/
            
            //pericycle-pericycle
            if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type
                && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4){
                expantion=0.8;
            }
            //pericycle-procambium
           /*  for (k=0; k<2;k++){
             if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==4){
             if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==5){
             expantion=1.0;
             }
             }
             }*/
            
            //endodermis-endodermis
            if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type
                && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==5){
                expantion=0.8;
            }
			
            //pericycle-endodermis
			for (k=0; k<2;k++){
				if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==4){
					if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==5){
						expantion=1.5;//1.5
					}
				}
			}
			
            vertex1=Edge[Cellnum[i].GroupEdge[j]].Edge_Vertex[0];
            vertex2=Edge[Cellnum[i].GroupEdge[j]].Edge_Vertex[1];
            Cellnum[i].GroupVertex_h[j].Fx+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0])/Cellnum[i].R_edge[j%jmax]);
           Cellnum[i].GroupVertex_h[j].Fy+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1])/Cellnum[i].R_edge[j%jmax]);
            Cellnum[i].GroupVertex_h[(j+1)%jmax].Fx+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])/Cellnum[i].R_edge[j%jmax]);
            Cellnum[i].GroupVertex_h[(j+1)%jmax].Fy+=-0.5*(Cellnum[i].k_api*expantion*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])/Cellnum[i].R_edge[j%jmax]);
        }
    }
        
        
    imax=Cellnum.size();
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].Sxx+=1/(3*Cellnum[i].Area)*(Cellnum[i].GroupVertex_h[j].Fx*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[0]+Cellnum[i].GroupVertex_h[(j+1)%jmax].Fx*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[0])
                            +1/(6*Cellnum[i].Area)*(Cellnum[i].GroupVertex_h[j%jmax].Fx*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[0]+Cellnum[i].GroupVertex_h[(j+1)%jmax].Fx*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[0]);
            Cellnum[i].Syy+=1/(3*Cellnum[i].Area)*(Cellnum[i].GroupVertex_h[j%jmax].Fy*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[1]+Cellnum[i].GroupVertex_h[(j+1)%jmax].Fy*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[1])
                            +1/(6*Cellnum[i].Area)*(Cellnum[i].GroupVertex_h[j%jmax].Fy*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[1]+Cellnum[i].GroupVertex_h[(j+1)%jmax].Fy*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[1]);
            Cellnum[i].Sxy+=1/(6*Cellnum[i].Area)*(Cellnum[i].GroupVertex_h[j%jmax].Fx*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[1]+Cellnum[i].GroupVertex_h[j%jmax].Fy*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[0]
                                                   +Cellnum[i].GroupVertex_h[(j+1)%jmax].Fx*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[1]+Cellnum[i].GroupVertex_h[(j+1)%jmax].Fy*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[0])
                            +1/(12*Cellnum[i].Area)*(Cellnum[i].GroupVertex_h[j%jmax].Fx*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[1]+Cellnum[i].GroupVertex_h[j%jmax].Fy*Cellnum[i].GroupVertex_h[(j+1)%jmax].coordinate_h[0]+
                                                     Cellnum[i].GroupVertex_h[(j+1)%jmax].Fx*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[1]+Cellnum[i].GroupVertex_h[(j+1)%jmax].Fy*Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[0]);
        }
        
        Cellnum[i].S1=0.5*((Cellnum[i].Sxx+Cellnum[i].Syy)+sqrt((Cellnum[i].Sxx-Cellnum[i].Syy)*(Cellnum[i].Sxx-Cellnum[i].Syy)+4*Cellnum[i].Sxy*Cellnum[i].Sxy));
        Cellnum[i].S2=0.5*((Cellnum[i].Sxx+Cellnum[i].Syy)-sqrt((Cellnum[i].Sxx-Cellnum[i].Syy)*(Cellnum[i].Sxx-Cellnum[i].Syy)+4*Cellnum[i].Sxy*Cellnum[i].Sxy));
        
        Cellnum[i].Theta1=atan((Cellnum[i].S1-Cellnum[i].Sxx)/Cellnum[i].Sxy);
        Cellnum[i].Theta2=atan((Cellnum[i].S2-Cellnum[i].Sxx)/Cellnum[i].Sxy);
        
        //Cellnum[i].Theta1=0.5*atan(2*Cellnum[i].Sxy/(Cellnum[i].Sxx-Cellnum[i].Syy));
    }
    
    //Euler method
    imax=Vertices.size();
    for (i=0; i<imax; i++) {
        Vertices[i].coordinate[0]=Vertices[i].coordinate[0]+dt*Vertices[i].Fx;
        Vertices[i].coordinate[1]=Vertices[i].coordinate[1]+dt*Vertices[i].Fy;
    }
    
	//energy
	energy_af=0;
	imax=Cellnum.size();
	for (i=0; i<imax; i++) {
		energy_af+=Cellnum[i].k_area*(Cellnum[i].Area-Cellnum[i].Area0)*(Cellnum[i].Area-Cellnum[i].Area0);
		jmax=Cellnum[i].GroupVertex.size();
	   for (j=0; j<jmax; j++) {
		   expantion=1.0;
			  for (k=0; k<2;k++){
			 if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==1){
			 if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==3 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==6){
			 expantion=1.2;
			 }
			 }
			 }
			 //xylem-pericycle
			   for (k=0; k<2;k++){
			  if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==1){
			  if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
			  expantion=1.5;
			  }
			  }
			  }
			 
			 //sieve element-pericycle
			  for (k=0; k<2;k++){
			  if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==2){
			  if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
			  expantion=1.5;//1.5
			  }
			  }
			  }
			 
			 //opc cell defect
			 for (k=0; k<2;k++){
				 if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==6 && Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_defect==1){
					 if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==4){
						 expantion=1.5;//1.5
					 }
				 }
			 }
			 
			 //pericycle-endodermis
			 for (k=0; k<2;k++){
				 if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==4){
					 if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[(k+1)%2]].cell_type==5){
						 expantion=1.5;//1.5
					 }
				 }
			 }
		   
		energy_af+=0.5*(Cellnum[i].k_api*expantion*Cellnum[i].R_edge[j%jmax]);
		}
	}
	//cout <<(energy_af-energy_bf)*3333<<" "<<energy_af*3333 <<"	"<< energy_bf*3333<< " "<<endl;
	
	energy_bf=energy_af;
	
     /*cout << Cellnum[2].GroupVertex_h[0].Fx << endl;
    cout << Cellnum[2].S1 << endl;
     cout << Cellnum[2].S2 << endl;
      cout <<  Cellnum[2].Theta1<< endl;
    cout << endl;*/
	
}//

void Tissue::Geometory_cell() {
	imax=Cellnum.size();
	//Area
	for (i=0; i<imax; i++) {
			Cellnum[i].Area=0;							  
	}
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for (j=0; j<jmax; j++) {
			Cellnum[i].Area+=0.5*(Vertices[Cellnum[i].GroupVertex[(j)%jmax]].coordinate[0]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]								  
								  -Vertices[Cellnum[i].GroupVertex[(j)%jmax]].coordinate[1]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]);							  
		}
	}
	//Edge
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].R_edge.size();	
		for (j=0; j<jmax; j++) {
			Cellnum[i].R_edge[j%jmax]=(sqrt((Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])
											*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0])
											+(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])
											*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])));							  
		}
	}	
	//Perimeter
	for (i=0; i<imax; i++) {
		Cellnum[i].Perimeter=0;							  
	}
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for (j=0; j<jmax; j++) {
			Cellnum[i].Perimeter+=Cellnum[i].R_edge[j%jmax];							  
		}
	}
	
	//centoroid
	for (i=0; i<imax; i++) {
		for(k=0;k<2;k++){
			Cellnum[i].centroid[k]=0;
		}
	}
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
	for(j=0;j<jmax;j++){
		double gaiseki;
		gaiseki=(Vertices[Cellnum[i].GroupVertex[(j)%jmax]].coordinate[0]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]								  
				 -Vertices[Cellnum[i].GroupVertex[(j)%jmax]].coordinate[1]*Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]);
		for(k=0;k<2;k++){
			Cellnum[i].centroid[k]+=((Vertices[Cellnum[i].GroupVertex[(j)%jmax]].coordinate[k])+(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[k]))*gaiseki;
		}
	}
	for(k=0;k<2;k++){
		Cellnum[i].centroid[k]=Cellnum[i].centroid[k]/(6.0*Cellnum[i].Area);
	}
	}
    
    
    //centroid-hosei
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[0]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Cellnum[i].centroid[0];
            Cellnum[i].GroupVertex_h[j%jmax].coordinate_h[1]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Cellnum[i].centroid[1];
        }
    }
    
	//Mormento
	for (i=0; i<imax; i++) {
		jmax=Cellnum[i].GroupVertex.size();	
		for (j=0; j<jmax; j++) {		
			Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]-Cellnum[i].centroid[0];
			Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]=Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1]-Cellnum[i].centroid[1];
		}
	}
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].af[j%jmax]=Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
            -(Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]);
        }
    }
    for (i=0; i<imax; i++) {
        Cellnum[i].Mxx=0;
        Cellnum[i].Myy=0;
        Cellnum[i].Mxy=0;
    }
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].Mxx+=2*Cellnum[i].af[j%jmax]*(Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]
                                                     +Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
                                                     +Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1])*tf;
        }
    }
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].Myy+=2*Cellnum[i].af[j%jmax]*(Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]
                                                     +Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]
                                                     +Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0])*tf;
        }
    }
    for (i=0; i<imax; i++) {
        jmax=Cellnum[i].GroupVertex.size();
        for (j=0; j<jmax; j++) {
            Cellnum[i].Mxy+=Cellnum[i].af[j%jmax]*(Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
                                                   +2*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1]
                                                   +2*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[1]
                                                   +Cellnum[i].GroupVertex_M[(j+1)%jmax].coordinate_h[0]*Cellnum[i].GroupVertex_M[j%jmax].coordinate_h[1])*tf;
        }
    }

}

void Tissue::Division(int p) {
if (p<9000) {//echhan12000 wt9000 han10000
	imax=Cellnum.size();
	for (i=0; i<imax; i++) {
		if (Cellnum[i].time_division==p) {
					
			/*lmax=Cellnum[i].GroupEdge.size();
			 for (j=0; j<lmax; j++) {	
					for (k=0; k<2; k++) {
						if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[k]].cell_type==4) {		
											Cellnum[i].division_count+=1;	
						}
					}
			 }*/
		
		//perocycle neighbor division
		int pericycle_edge;
		pericycle_edge=0;
			
		if (Cellnum[i].cell_type!=3 && Cellnum[i].cell_type!=6) {
		jmax=Cellnum[i].GroupEdge.size();
		for (j=0; j<jmax; j++) {	
		if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){
			pericycle_edge+=1;
			if (Cellnum[i].time_division>=9000) {
				//Cellnum[i].cell_defect=1;
				pericycle_edge=0;
			}
		}
		} //pericycle neighbor
		} //if celltype
			
		//OPC division
		/*if (Cellnum[i].cell_type==6 && (Cellnum[i].cell_lineage==8)) {//8
			pericycle_edge+=1;
			Cellnum[i].cell_defect=1;
		}*/
			
		//IPC division
		//if (Cellnum[i].cell_type==3 && Cellnum[i].cell_defect==0) {
		/*if (Cellnum[i].cell_type==3 && Cellnum[i].cell_lineage==5) {
			jmax=Cellnum[i].GroupEdge.size();
			for (j=0; j<jmax; j++) {
			if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==1 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==1){
				pericycle_edge+=1;
				Cellnum[i].cell_defect=1;
				//if (Cellnum[i].division_count>2) {//1 time division
				//pericycle_edge=0;
				//}
			}
			}
		}*/
			
			if (Cellnum[i].cell_type==1) {
				pericycle_edge+=1;
				if (Cellnum[i].division_count>1) {//1 time division //x6-1
				pericycle_edge=0;
				}
			}
				
		if (pericycle_edge>=1) {//pericycle limit 1
             cout << p<<"	"<<i << endl;
			double theta,r,aa,bb,ramuda1,gra1,gra2,heiho;
			int vertex_max=0;
			int edge_max=0;
			int j1,hmax;
			int i_2[2];
			int j_2=0;
			double ab[2];
			int j_cross[2];
			int flag=0;
			
			Cellnum[i].GroupVertex_divi.assign(Cellnum[i].GroupVertex.begin(), Cellnum[i].GroupVertex.end());
			Cellnum[i].GroupEdge_divi.assign(Cellnum[i].GroupEdge.begin(), Cellnum[i].GroupEdge.end());
			jmax=Cellnum[i].GroupVertex_divi.size();
						
			//division_area
			//heiho=(Cellnum[i].Mxx+Cellnum[i].Myy)*(Cellnum[i].Mxx+Cellnum[i].Myy)-4*(Cellnum[i].Mxx*Cellnum[i].Myy-Cellnum[i].Mxy*Cellnum[i].Mxy*Cellnum[i].Mxy);
			//if (heiho >= 0){
			ramuda1=0.5*((Cellnum[i].Mxx+Cellnum[i].Myy)+sqrt((Cellnum[i].Mxx+Cellnum[i].Myy)*(Cellnum[i].Mxx+Cellnum[i].Myy)-4*(Cellnum[i].Mxx*Cellnum[i].Myy-Cellnum[i].Mxy*Cellnum[i].Mxy)));
			//}
			//if (heiho < 0){
			//	ramuda1=0.5*(Cellnum[i].Mxx+Cellnum[i].Myy);
			//}
			gra1=(Cellnum[i].Mxx-ramuda1)/Cellnum[i].Mxy;
			gra2=Cellnum[i].Mxy/(Cellnum[i].Myy-ramuda1);
			theta=atan(gra1);
			
			theta=theta*rand_normal(1,0.1);
			
			if (Cellnum[i].cell_type==1) {
				//if (Cellnum[i].division_xylem<1) {	
					theta=2*M_PI*0.22;//0.21
				//}
			}
			if (Cellnum[i].cell_type==1 && Cellnum[i].division_count>0) {
				//if (Cellnum[i].division_xylem<1) {
					theta=2*M_PI*0.21;//0.21
				//}
			}
			
			if (Cellnum[i].cell_type==2) {
				theta=2*M_PI*-0.03;
              /*  if (i==5) {
                theta=2*M_PI*0.1;
            }*/
			}
			
			/*if (Cellnum[i].cell_lineage==6 && Cellnum[i].division_count==0) {
				//if (Cellnum[i].division_xylem<1) {
					theta=2*M_PI*-0.21;//0.21 -0.18
				//}
			}
			if (Cellnum[i].cell_lineage==7 && Cellnum[i].division_count==0) {
				//if (Cellnum[i].division_xylem<1) {
					theta=2*M_PI*0.15;//0.21
				//}
			}*/
			
				/*if (Cellnum[i].cell_type==3 && Cellnum[i].cell_lineage==5) {
					theta=2*M_PI*-0.03*rand_normal(1,0.1);
				}*/
			
			/*if (Cellnum[i].cell_type==6) {
				jmax=Cellnum[i].GroupEdge.size();	
				for (j=0; j<jmax; j++) {	
					if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==2 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==2){
						int vertex1,vertex2;
						vertex1=Edge[Cellnum[i].GroupEdge[j]].Edge_Vertex[0];
						vertex2=Edge[Cellnum[i].GroupEdge[j]].Edge_Vertex[1];
						gra1=(Vertices[vertex1].coordinate[1]-Vertices[vertex2].coordinate[1])/(Vertices[vertex1].coordinate[0]-Vertices[vertex2].coordinate[0]);
						theta=atan(gra1);
						//theta=2*M_PI*drand48();
					}
				}
            
				//if (Cellnum[i].division_phloem<1) {
					//theta=2*M_PI*-0.03;
				//	Cellnum[i].division_phloem+=1;
				//}
				
			}*/
            
			/*if (Cellnum[i].cell_type==3) {
				theta=2*M_PI*drand48();
			}*/
			
			Cellnum[i].cent_dr[0]=cos(theta);
			Cellnum[i].cent_dr[1]=sin(theta);
			
			/*cout << (Cellnum[i].Mxx+Cellnum[i].Myy)*(Cellnum[i].Mxx+Cellnum[i].Myy)
			-4*(Cellnum[i].Mxx*Cellnum[i].Myy-Cellnum[i].Mxy*Cellnum[i].Mxy*Cellnum[i].Mxy) << endl;
			cout << "b1" << endl;
			cout << ramuda1 << endl;
			cout << Cellnum[i].Mxy << endl;
			cout << gra1 << endl;
			cout << gra2 << endl;
			cout << theta << endl;*/
			
			jmax=Cellnum[i].GroupVertex.size();	
			for (j=0; j<jmax; j++) {
				aa=Cellnum[i].cent_dr[0]*(Cellnum[i].centroid[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])
				-Cellnum[i].cent_dr[1]*(Cellnum[i].centroid[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]);
				bb=Cellnum[i].cent_dr[0]*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1])
				-Cellnum[i].cent_dr[1]*(Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0]);
				r=aa/bb;
				
				if (r>0&&r<1) {
					j1=j;
					ab[0]=Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[0];
			        ab[1]=Vertices[Cellnum[i].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[i].GroupVertex[j%jmax]].coordinate[1];
					
					vertex_max=Vertices.size();
					Vertices.push_back(Vertex());
                    Vertices[vertex_max].coordinate[0]=Vertices[Cellnum[i].GroupVertex[j]].coordinate[0]+r*ab[0];
                    Vertices[vertex_max].coordinate[1]=Vertices[Cellnum[i].GroupVertex[j]].coordinate[1]+r*ab[1];
					
					//cout << Vertices[vertex_max].coordinate[0]<< " " <<Vertices[vertex_max].coordinate[1] << endl;
					
					hmax=Cellnum[i].GroupVertex_divi.size();
					Cellnum[i].GroupVertex_divi.push_back(int());
					for (int h=hmax-1; h>j1+flag; h--) {
						Cellnum[i].GroupVertex_divi[(h+1)]=Cellnum[i].GroupVertex_divi[h];
					}
					
					Cellnum[i].GroupVertex_divi[j1+flag+1]=vertex_max;
					j_cross[flag]=j1;
					
					
					for (int k=0; k<2; k++) {
						if (Edge[Cellnum[i].GroupEdge[(j1)%jmax]].Edge_Area[k]!=i) {
							i_2[flag]=Edge[Cellnum[i].GroupEdge[(j1)%jmax]].Edge_Area[k];
							//cout << i << endl;
							//cout << k << endl;
							//cout << i_2[flag] << endl;
						}
					}
					//cout << "b1-2-2" << endl;
					/*hmax=Cellnum[i_2[flag]].GroupEdge.size();
					 for (int h=0; h<hmax; h++) {
					 cout << Cellnum[i_2[flag]].GroupEdge[h] << endl;
					 }
					 cout << endl;*/
					//cout << "b1-2-3" << endl;
					hmax=Cellnum[i_2[flag]].GroupVertex.size();
					for (int h=0; h<hmax; h++) {
						if (Cellnum[i].GroupVertex[j1%jmax]==Cellnum[i_2[flag]].GroupVertex[h]) {
							j_2=h;
						}
					}
					
					//cout << "b1-3" << endl;
					
					Cellnum[i_2[flag]].GroupVertex.push_back(int());
					for (int h=hmax-1; h>j_2-1; h--) {
						Cellnum[i_2[flag]].GroupVertex[(h+1)]=Cellnum[i_2[flag]].GroupVertex[h];
					}
					Cellnum[i_2[flag]].GroupVertex[j_2]=vertex_max;
					
					Cellnum[i_2[flag]].R_edge.push_back(double());//edge length add		
					Cellnum[i_2[flag]].af.push_back(double());
					Cellnum[i_2[flag]].GroupVertex_h.push_back(Vertex());
                    Cellnum[i_2[flag]].GroupVertex_M.push_back(Vertex());
					edge_max=Edge.size();
					Edge.push_back(Vertex());
			
					hmax=Cellnum[i_2[flag]].GroupEdge.size();
					Cellnum[i_2[flag]].GroupEdge.push_back(int());
					for (int h=hmax-1; h>j_2-1; h--) {
						Cellnum[i_2[flag]].GroupEdge[(h+1)]=Cellnum[i_2[flag]].GroupEdge[h];
					}
					
					//cout << "b1-4" << endl;
					
					Cellnum[i_2[flag]].GroupEdge[j_2]=edge_max;
					//cout <<edge_max<<endl;
					//cout << endl;
					
					/*hmax=Cellnum[i_2[flag]].GroupEdge.size();
					for (int h=0; h<hmax; h++) {
					 cout << Cellnum[i_2[flag]].GroupEdge[h] << endl;
					 }
					 cout << endl;*/
					
					flag+=1;
				}
			}
			edge_max=Edge.size();
			Edge[edge_max-2].Edge_Area[0]=i_2[0];
			Edge[edge_max-2].Edge_Area[1]=i;
			Edge[edge_max-2].Edge_Vertex[0]=vertex_max-1;
			Edge[edge_max-2].Edge_Vertex[1]=Cellnum[i].GroupVertex_divi[j_cross[0]];
			Edge[edge_max-2].T1_flag=0;
			//cout << i_2[0] << endl;
			Edge[edge_max-1].Edge_Area[0]=imax;
			Edge[edge_max-1].Edge_Area[1]=i_2[1];
			Edge[edge_max-1].Edge_Vertex[0]=vertex_max;
			Edge[edge_max-1].Edge_Vertex[1]=Cellnum[i].GroupVertex_divi[j_cross[1]+1];
			Edge[edge_max-1].T1_flag=0;
			//cout << i_2[1] << endl;
			Edge.push_back(Vertex());
			Edge[edge_max].Edge_Area[0]=i;
			Edge[edge_max].Edge_Area[1]=imax;
			Edge[edge_max].Edge_Vertex[0]=vertex_max-1;
			Edge[edge_max].Edge_Vertex[1]=vertex_max;			
			Edge[edge_max].T1_flag=0;
			
			/*jmax=Cellnum[i].GroupVertex_divi.size();
			 for (j=0; j<jmax; j++) {
			 cout << Cellnum[i].GroupVertex_divi[j] <<endl;
			 }*/
			
			//cout << "b2" << endl;
			
			//cross point 2ko
			if (flag==2) {
				//new cell
				Cellnum.push_back(Cell());//youso add
				imax=Cellnum.size();
				for (j=j_cross[0]+1; j<(j_cross[1]+1)+2; j++) {
					Cellnum[imax-1].GroupVertex.push_back(Cellnum[i].GroupVertex_divi[j]);
				}
				Cellnum[imax-1].cell_type=Cellnum[i].cell_type;
				Cellnum[imax-1].cell_lineage=Cellnum[i].cell_lineage;
				Cellnum[imax-1].cell_defect=Cellnum[i].cell_defect;
				Cellnum[imax-1].division_count=Cellnum[i].division_count;
				Cellnum[imax-1].k_area=Cellnum[i].k_area;
				Cellnum[imax-1].k_peri=Cellnum[i].k_peri;
				Cellnum[imax-1].k_api=Cellnum[i].k_api;
				Cellnum[imax-1].k_bas=0.0;
				Cellnum[imax-1].k_lat=0.0;
				Edge[edge_max-2].k_au=k_au;
				Edge[edge_max-1].k_cy=k_cy;
				Edge[edge_max].k_aucy=k_aucy;
                
				//Cellnum[i].Area0=Area_divi;
				//Cellnum[imax-1].Area0=Area_divi;
                
                Cellnum[imax-1].Area0_0=Cellnum[i].Area0_0*rand_normal(1,0.1);
                
				if (Cellnum[i].cell_type==1) {
					Cellnum[i].Area0_0=Area_divimx*rand_normal(1,0.1);
					Cellnum[imax-1].Area0_0=Area_divimx*rand_normal(1,0.1);
				}
				
				if (Cellnum[i].cell_type==2 || Cellnum[i].cell_type==7) {
					Cellnum[i].Area0_0=Area_divise*rand_normal(1,0.1);
					Cellnum[imax-1].Area0_0=Area_divise*rand_normal(1,0.1);
				}
				
				if (Cellnum[i].cell_type==4) {
					Cellnum[i].Area0_0=Area_divipr*rand_normal(1,0.1);
					Cellnum[imax-1].Area0_0=Area_divipr*rand_normal(1,0.1);
				}
                
                //Cellnum[imax-1].Area0=Cellnum[imax-1].Area0_0*1.0;
                
                Cellnum[imax-1].Area0=Cellnum[i].Area*0.5;
                Cellnum[i].Area0=Cellnum[i].Area*0.5;
				
				//cout << "b2-2" << endl;
				
				jmax=Cellnum[imax-1].GroupVertex.size();	
				for (j=0; j<jmax; j++) {
					Cellnum[imax-1].R_edge.push_back(sqrt((Vertices[Cellnum[imax-1].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[imax-1].GroupVertex[j%jmax]].coordinate[0])
														  *(Vertices[Cellnum[imax-1].GroupVertex[(j+1)%jmax]].coordinate[0]-Vertices[Cellnum[imax-1].GroupVertex[j%jmax]].coordinate[0])
														  +(Vertices[Cellnum[imax-1].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[imax-1].GroupVertex[j%jmax]].coordinate[1])
														  *(Vertices[Cellnum[imax-1].GroupVertex[(j+1)%jmax]].coordinate[1]-Vertices[Cellnum[imax-1].GroupVertex[j%jmax]].coordinate[1])));							  
				}
				jmax=Cellnum[imax-1].GroupVertex.size();	
				for (j=0; j<jmax; j++) {
					Cellnum[imax-1].GroupVertex_h.push_back(Vertex());
                    Cellnum[imax-1].GroupVertex_M.push_back(Vertex());
					Cellnum[imax-1].af.push_back(double());
				}								
				edge_max=Edge.size();
				for (j=j_cross[0]; j<(j_cross[1]); j++) {
					Cellnum[imax-1].GroupEdge.push_back(Cellnum[i].GroupEdge_divi[j]);
				}
				Cellnum[imax-1].GroupEdge.push_back(edge_max-2);
				Cellnum[imax-1].GroupEdge.push_back(edge_max-1);
				
				/*jmax=Cellnum[i].GroupVertex.size();
				 for (j=0; j<jmax; j++) {
				 cout << Vertices[Cellnum[i].GroupVertex[j]].coordinate[0]<< " " <<Vertices[Cellnum[i].GroupVertex[j]].coordinate[1] << endl;
				 }
				 cout << endl;*/
				
				int i_divi=0;
				int jmax_new;
				int han,han1,han2;
				jmax=Cellnum[i].GroupVertex_divi.size();
				jmax_new=Cellnum[imax-1].GroupVertex.size();
				han=jmax+2-jmax_new;
				han1=Cellnum[i].GroupVertex.size()-han;
				han2=han-Cellnum[i].GroupVertex.size();
				if (han1>0) {
					for (int k=0; k<han1; k++) {
						Cellnum[i].GroupVertex.pop_back();
						Cellnum[i].R_edge.pop_back();
						Cellnum[i].GroupEdge.pop_back();
						Cellnum[i].GroupVertex_h.pop_back();
                        Cellnum[i].GroupVertex_M.pop_back();
						Cellnum[i].af.pop_back();
					}
				}
				if (han2>0) {
					for (int k=0; k<han2; k++) {
						Cellnum[i].GroupVertex.push_back(int());
						Cellnum[i].R_edge.push_back(double());	
						Cellnum[i].GroupEdge.push_back(int());
						Cellnum[i].GroupVertex_h.push_back(Vertex());
                        Cellnum[i].GroupVertex_M.push_back(Vertex());
						Cellnum[i].af.push_back(double());
					}
				}
				
				for (j=0; j<j_cross[0]+2; j++) {
					Cellnum[i].GroupVertex[i_divi]=Cellnum[i].GroupVertex_divi[j];
					i_divi+=1;
				}
				for (j=(j_cross[1]+1)+1; j<jmax; j++) {
					Cellnum[i].GroupVertex[i_divi]=Cellnum[i].GroupVertex_divi[j];
					i_divi+=1;
				}
				
				jmax=Cellnum[i].GroupEdge_divi.size();
				i_divi=0;
				for (j=0; j<j_cross[0]; j++) {
					Cellnum[i].GroupEdge[i_divi]=Cellnum[i].GroupEdge_divi[j];
					i_divi+=1;
				}
				Cellnum[i].GroupEdge[i_divi]=edge_max-3;
				i_divi+=1;
				Cellnum[i].GroupEdge[i_divi]=edge_max-1;
				i_divi+=1;
				for (j=j_cross[1]; j<jmax; j++) {
					Cellnum[i].GroupEdge[i_divi]=Cellnum[i].GroupEdge_divi[j];
					i_divi+=1;
				}
				/*jmax=Cellnum[imax-1].GroupVertex.size();
				 for (j=0; j<jmax; j++) {
				 cout << Vertices[Cellnum[imax-1].GroupVertex[j]].coordinate[0]<< " " <<Vertices[Cellnum[imax-1].GroupVertex[j]].coordinate[1] << endl;
				 }
				 cout << endl;*/
				
			}
			
			//cout << "b3" << endl;
			
			//vertex
			imax=Cellnum.size();
			vertex_max=Vertices.size();
			jmax=Cellnum[i].GroupVertex_divi.size();
			Vertices[vertex_max-2].NeighborCell.push_back(imax-1);
			Vertices[vertex_max-2].NeighborCell.push_back(i);
			Vertices[vertex_max-2].NeighborCell.push_back(i_2[0]);
			Vertices[vertex_max-2].NeighborVertex.push_back(Cellnum[i].GroupVertex_divi[j_cross[0]%jmax]);
			Vertices[vertex_max-2].NeighborVertex.push_back(Cellnum[i].GroupVertex_divi[(j_cross[0]+2)%jmax]);
			Vertices[vertex_max-2].NeighborVertex.push_back(vertex_max-1);		
			Vertices[vertex_max-1].NeighborCell.push_back(i);
			Vertices[vertex_max-1].NeighborCell.push_back(imax-1);
			Vertices[vertex_max-1].NeighborCell.push_back(i_2[1]);
			Vertices[vertex_max-1].NeighborVertex.push_back(Cellnum[i].GroupVertex_divi[(j_cross[1]+1)%jmax]);
			Vertices[vertex_max-1].NeighborVertex.push_back(Cellnum[i].GroupVertex_divi[(j_cross[1]+1+2)%jmax]);
			Vertices[vertex_max-1].NeighborVertex.push_back(vertex_max-2);
			jmax=Cellnum[i].GroupVertex_divi.size();
			for (int k=0; k<3; k++) {
				if (Vertices[Cellnum[i].GroupVertex_divi[j_cross[0]%jmax]].NeighborVertex[k]==Cellnum[i].GroupVertex_divi[(j_cross[0]+2)%jmax]) {					
					Vertices[Cellnum[i].GroupVertex_divi[j_cross[0]%jmax]].NeighborVertex[k]=vertex_max-2;
				}
			}
			for (j=j_cross[0]+2; j<(j_cross[1]+1)+1; j++) {
				for (int k=0; k<3; k++) {
					if (Vertices[Cellnum[i].GroupVertex_divi[j%jmax]].NeighborCell[k]==i) {					
						Vertices[Cellnum[i].GroupVertex_divi[j%jmax]].NeighborCell[k]=imax-1;
					}
				}
			}
			for (int k=0; k<3; k++) {
				if (Vertices[Cellnum[i].GroupVertex_divi[(j_cross[0]+2)%jmax]].NeighborVertex[k]==Cellnum[i].GroupVertex_divi[j_cross[0]%jmax]) {					
					Vertices[Cellnum[i].GroupVertex_divi[(j_cross[0]+2)%jmax]].NeighborVertex[k]=vertex_max-2;
				}
			}
			for (int k=0; k<3; k++) {
				if (Vertices[Cellnum[i].GroupVertex_divi[(j_cross[1]+1)%jmax]].NeighborVertex[k]==Cellnum[i].GroupVertex_divi[(j_cross[1]+1+2)%jmax]) {					
					Vertices[Cellnum[i].GroupVertex_divi[(j_cross[1]+1)%jmax]].NeighborVertex[k]=vertex_max-1;
				}
			}
			for (int k=0; k<3; k++) {
				if (Vertices[Cellnum[i].GroupVertex_divi[(j_cross[1]+1+2)%jmax]].NeighborVertex[k]==Cellnum[i].GroupVertex_divi[(j_cross[1]+1)%jmax]) {					
					Vertices[Cellnum[i].GroupVertex_divi[(j_cross[1]+1+2)%jmax]].NeighborVertex[k]=vertex_max-1;
				}
			}
			
			//Edge
			jmax=Cellnum[i].GroupVertex_divi.size();
			hmax=Cellnum[i].GroupEdge_divi.size();
			for (int k=0; k<2; k++) {
				if (Edge[Cellnum[i].GroupEdge_divi[j_cross[0]%hmax]].Edge_Vertex[k]==Cellnum[i].GroupVertex_divi[(j_cross[0])%jmax]) {					
					Edge[Cellnum[i].GroupEdge_divi[j_cross[0]%hmax]].Edge_Vertex[k]=vertex_max-2;
				}
			}
			for (int k=0; k<2; k++) {
				if (Edge[Cellnum[i].GroupEdge_divi[(j_cross[1])%hmax]].Edge_Vertex[k]==Cellnum[i].GroupVertex_divi[(j_cross[1]+1)%jmax]) {					
					Edge[Cellnum[i].GroupEdge_divi[(j_cross[1])%hmax]].Edge_Vertex[k]=vertex_max-1;
				}
			}
			for (j=j_cross[0]; j<j_cross[1]; j++) {
				for (int k=0; k<2; k++) {
					if (Edge[Cellnum[i].GroupEdge_divi[j%hmax]].Edge_Area[k]==i) {					
						Edge[Cellnum[i].GroupEdge_divi[j%hmax]].Edge_Area[k]=imax-1;
					}
				}
			}
		//cout << "b4" << endl;
			
            //division timinng
            Cellnum[i].time_division=-1;//p+int(9000+2000*drand48());    //3500
            Cellnum[imax-1].time_division=-1;//p+int(9000+2000*drand48());//3500
            
            //xylem
			/*if (Cellnum[i].cell_type==1 && Cellnum[i].cell_lineage==2 && Cellnum[i].division_count==0) {
				jmax=Cellnum[i].GroupEdge.size();
				for (j=0; j<jmax; j++) {
				if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){
				Cellnum[i].time_division=p+int(5000*rand_normal(1,0.2));//-1;//8000
				}
				}
				jmax=Cellnum[imax-1].GroupEdge.size();
				for (j=0; j<jmax; j++) {
				if (Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[1]].cell_type==4){
				Cellnum[imax-1].time_division=p+int(5000*rand_normal(1,0.2));
				}
				}
				
			}*/
			
           /* if (Cellnum[i].cell_type==1) {
             //Cellnum[i].time_division=p+int(12000+2000*drand48());//-1;
             //Cellnum[imax-1].time_division=p+int(12000+2000*drand48());//-1;
             
             //Cellnum[i].time_division=p+int(15000+300*drand48());
             //Cellnum[imax-1].time_division=p+int(15000+300*drand48());
        if (p<qmax*0.1) {
                double judge1,judge2;
                judge1=drand48();
            if(judge1<0.0) {
             jmax=Cellnum[i].GroupEdge.size();
             for (j=0; j<jmax; j++) {
             if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){
             Cellnum[i].time_division=p+int(9000+300*drand48());//-1;
             }
             }
            }
                  
            judge2=drand48();
            if(judge2<0.0) {
             jmax=Cellnum[imax-1].GroupEdge.size();
             for (j=0; j<jmax; j++) {
             if (Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[1]].cell_type==4){
             Cellnum[imax-1].time_division=p+int(9000+300*drand48());
             }
             }
            }
              }
            }*/
            
            //procambium
            if (/*Cellnum[i].cell_type==3 || Cellnum[i].cell_type==6 ||*/ Cellnum[i].cell_type==7) {
                //Cellnum[i].time_division=p+int(4500+4500*drand48());//p+int(3000+1000*drand48());
               // Cellnum[imax-1].time_division=p+int(4500+4500*drand48());//p+int(3000+1000*drand48());
                Cellnum[i].cell_type=3;//pro
                 Cellnum[imax-1].cell_type=3;//pro
                
                jmax=Cellnum[i].GroupEdge.size();
                 kmax=Cellnum[i].GroupEdge.size();
                for (j=0; j<jmax; j++) {
                    if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==2 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==2){
                        for (k=0; k<kmax; k++) {
                        if(Cellnum[Edge[Cellnum[i].GroupEdge[k]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[k]].Edge_Area[1]].cell_type==4){//4
							Cellnum[i].time_division=p+int(2150*rand_normal(1,0.2));//2150
                            Cellnum[i].cell_type=7;//companion cell
                            
                            /*if (Cellnum[i].centroid[1]>Cellnum[1].centroid[1]*1){
                                Cellnum[i].time_division=p+int(2000+500*drand48());
                            }*/
                        }
                        }
                    }
                }
                jmax=Cellnum[imax-1].GroupEdge.size();
                 kmax=Cellnum[imax-1].GroupEdge.size();
                for (j=0; j<jmax; j++) {
                    if (Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[0]].cell_type==2 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[1]].cell_type==2){
                        for (k=0; k<kmax; k++) {
                        if(Cellnum[Edge[Cellnum[imax-1].GroupEdge[k]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[k]].Edge_Area[1]].cell_type==4){//4
                             Cellnum[imax-1].time_division=p+int(2150*rand_normal(1,0.2));
                            Cellnum[imax-1].cell_type=7;//companion cell
                            
                            /*if (Cellnum[i].centroid[1]>Cellnum[1].centroid[1]*1){
                                Cellnum[imax-1].time_division=p+int(2000+500*drand48());
                            }*/
                        }
                        }
                    }
                }
        
                
           /* if (p<qmax*0.4) {
                double judge1,judge2;
                judge1=drand48();
                if(judge1<0.0) {
                    jmax=Cellnum[i].GroupEdge.size();
                    for (j=0; j<jmax; j++) {
                       if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){
                     Cellnum[i].time_division=p+int(3000+1000*drand48());//3000
                    Cellnum[i].cell_type=7;//companion cell
                       }
                }
                }
                judge2=drand48();
                if(judge2<0.0) {
                    jmax=Cellnum[imax-1].GroupEdge.size();
                    for (j=0; j<jmax; j++) {
                        if (Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[1]].cell_type==4){
                Cellnum[imax-1].time_division=p+int(3000+1000*drand48());
                Cellnum[imax-1].cell_type=7;//companion cell
                        }
                    }
                }

                if (Cellnum[i].centroid[0]>Cellnum[1].centroid[0]*1 && Cellnum[i].centroid[1]<Cellnum[1].centroid[1]*1){
                    jmax=Cellnum[i].GroupEdge.size();
                    for (j=0; j<jmax; j++) {
                    if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){
                        Cellnum[i].time_division=p+int(3000+1000*drand48());//3000
                        Cellnum[i].cell_type=7;//companion cell
                   }
                    }
                    jmax=Cellnum[imax-1].GroupEdge.size();
                    for (j=0; j<jmax; j++) {
                    if (Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[1]].cell_type==4){
                        Cellnum[imax-1].time_division=p+int(3000+1000*drand48());
                        Cellnum[imax-1].cell_type=7;//companion cell
                    }
                    }
                }
                     
                 }*/
            }
            
            //phloem　procambiumの後にする
            /*if (i==4 || i==cell_p1) {//right-5
             Cellnum[i].time_division=-1;
             Cellnum[imax-1].time_division=-1;//p+int(3000+100*drand48());
             cell_p1=imax-1;
             Cellnum[i].cell_type=3;
             if (Cellnum[imax-1].division_count<2) {
             //Cellnum[imax-1].Area0=Area_divise*2;
             Cellnum[imax-1].Area0=Area_divise*1.0;
             }
             }
             if (i==3) {//right-4
             Cellnum[i].time_division=-1;//p+int(3000+100*drand48());
             Cellnum[imax-1].time_division=p+int(3000+100*drand48());//-1;
             Cellnum[imax-1].cell_type=3;
             if (Cellnum[i].division_count<2) {
             //Cellnum[i].Area0=Area_divise*2;
             Cellnum[i].Area0=Area_divise*1.0;
             }
             }*/
            
            //PSE
            if (Cellnum[i].cell_type==2) {
                //Cellnum[i].time_division=p+int(4500+4500*drand48());//-1;
                //Cellnum[imax-1].time_division=p+int(4500+4500*drand48());;
                Cellnum[i].cell_type=3;
                Cellnum[imax-1].cell_type=3;
                
                jmax=Cellnum[i].GroupEdge.size();
                for (j=0; j<jmax; j++) {
                    if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){//xylem-1 pericycle-4
                        Cellnum[i].cell_type=2;
                        Cellnum[i].time_division=p+int(4150*rand_normal(1,0.2));//p+int(7500+100*drand48());//-1;4150
                    }
                }
                jmax=Cellnum[imax-1].GroupEdge.size();
                for (j=0; j<jmax; j++) {
                    if (Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[1]].cell_type==4){//4
                        Cellnum[imax-1].cell_type=2;
                        Cellnum[imax-1].time_division=p+int(4150*rand_normal(1,0.2));//p+int(7500+100*drand48());
                    }
                }
            }
           
			//OPC
			if (Cellnum[i].cell_type==3){
			jmax=Cellnum[i].GroupEdge.size();
				for (j=0; j<jmax; j++) {
				if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){//4
					Cellnum[i].cell_type=6;
				}
				}
			}
			if (Cellnum[imax-1].cell_type==3){
			jmax=Cellnum[imax-1].GroupEdge.size();
				for (j=0; j<jmax; j++) {
				if(Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[imax-1].GroupEdge[j]].Edge_Area[1]].cell_type==4){//4
					Cellnum[imax-1].cell_type=6;
				}
				}
			}
			//OPC defect
            if (Cellnum[i].cell_type==6 && (Cellnum[i].cell_lineage==8)) {
                Cellnum[i].time_division=p+int(4150*rand_normal(1,0.2));//3150
            }
			if (Cellnum[imax-1].cell_type==6 && (Cellnum[i].cell_lineage==8)) {
                Cellnum[imax-1].time_division=p+int(4150*rand_normal(1,0.2));
            }
			
			//IPC defect
            //if (Cellnum[i].cell_type==3) {
			if (Cellnum[i].cell_type==3 && Cellnum[i].cell_lineage==5) {
                Cellnum[i].time_division=p+int(4150*rand_normal(1,0.2));//4150
            }
			//if (Cellnum[imax-1].cell_type==3) {
			if (Cellnum[imax-1].cell_type==3 && Cellnum[imax-1].cell_lineage==5) {
                Cellnum[imax-1].time_division=p+int(4150*rand_normal(1,0.2));//4150
            }
            
            //pericyle
            if (Cellnum[i].cell_type==4) {
                Cellnum[i].time_division=p+int(15000+5000*drand48());
                Cellnum[imax-1].time_division=p+int(15000+5000*drand48());
            }
            
            /*if (Cellnum[i].cell_type==1) {
             Cellnum[i].time_division=p+int(7000+300*drand48());
             Cellnum[imax-1].time_division=p+int(7000+300*drand48());
             }*/
            
           /* if (p>qmax*0.8) {
                Cellnum[i].time_division=-1;
                Cellnum[imax-1].time_division=-1;
            }*/
            
			Cellnum[i].Auxin=Cellnum[i].Auxin;	
			Cellnum[imax-1].Auxin=Cellnum[i].Auxin;
			Cellnum[i].Cytokinin=Cellnum[i].Cytokinin;	
			Cellnum[imax-1].Cytokinin=Cellnum[i].Cytokinin;
			
				 //cout << endl;
				 Cellnum[i].division_xylem+=1;
				 Cellnum[i].division_count+=1;
				 Cellnum[imax-1].division_count+=1;
				}//division_count
			 }//time settei
		}//i
	}//stop
}//kansu

void Tissue::T1_cell() {
	int vertex1,vertex2,vertex3,vertex4;
	double length_edge;
	int cell01,cell02;
	int cell1,cell2;
	int cell3,cell4;
	int j,j1,j2,j3,j4,jmax,imax,hmax;
	
	imax=Edge.size();
	for (i=0; i<imax; i++) {
		vertex1=Edge[i].Edge_Vertex[0];
		vertex2=Edge[i].Edge_Vertex[1];
		
			length_edge=sqrt((Vertices[vertex1].coordinate[0]-Vertices[vertex2].coordinate[0])
						 *(Vertices[vertex1].coordinate[0]-Vertices[vertex2].coordinate[0])
						 +(Vertices[vertex1].coordinate[1]-Vertices[vertex2].coordinate[1])
						 *(Vertices[vertex1].coordinate[1]-Vertices[vertex2].coordinate[1]));
		
		if (length_edge<0.5) {
			Edge[i].T1_flag+=1;
			//cout << "a0" << endl;
		}
		
		if (Edge[i].T1_flag>1000) {
			cell01=Edge[i].Edge_Area[0];
			cell02=Edge[i].Edge_Area[1];
			
			for (int k=0; k<3; k++) {
				if (Vertices[vertex1].NeighborVertex[k]==vertex2) {
					vertex3=Vertices[vertex1].NeighborVertex[(k+1)%3];
						 }
			}	
			for (int k=0; k<3; k++) {
				if (Vertices[vertex2].NeighborVertex[k]==vertex1) {
					vertex4=Vertices[vertex2].NeighborVertex[(k+1)%3];
				}
			}
			
			//cell1 cell2
			jmax=Cellnum[cell01].GroupVertex.size();
			for (j=0; j<jmax; j++) {
				if (Cellnum[cell01].GroupVertex[j]==vertex3) {
					cell1=cell01;
				}
				if (Cellnum[cell01].GroupVertex[j]==vertex4) {
					cell2=cell01;
				}
			}
			jmax=Cellnum[cell02].GroupVertex.size();
			for (j=0; j<jmax; j++) {
				if (Cellnum[cell02].GroupVertex[j]==vertex3) {
					cell1=cell02;
				}
				if (Cellnum[cell02].GroupVertex[j]==vertex4) {
					cell2=cell02;
				}
			}			
			
			Cellnum[cell1].GroupVertex_T1.assign(Cellnum[cell1].GroupVertex.begin(), Cellnum[cell1].GroupVertex.end());
			Cellnum[cell1].GroupEdge_T1.assign(Cellnum[cell1].GroupEdge.begin(), Cellnum[cell1].GroupEdge.end());
			Cellnum[cell1].R_edge_T1.assign(Cellnum[cell1].R_edge.begin(), Cellnum[cell1].R_edge.end());
			
			Cellnum[cell2].GroupVertex_T1.assign(Cellnum[cell2].GroupVertex.begin(), Cellnum[cell2].GroupVertex.end());
			Cellnum[cell2].GroupEdge_T1.assign(Cellnum[cell2].GroupEdge.begin(), Cellnum[cell2].GroupEdge.end());
			Cellnum[cell2].R_edge_T1.assign(Cellnum[cell2].R_edge.begin(), Cellnum[cell2].R_edge.end());

			//cell1
			jmax=Cellnum[cell1].GroupVertex.size();
			for (j=0; j<jmax; j++) {
				if (Cellnum[cell1].GroupVertex[j]==vertex1) {
					j1=j;
				}
			}
			
			for (j=jmax-1; j>j1; j--) {
				Cellnum[cell1].GroupVertex[j-1]=Cellnum[cell1].GroupVertex_T1[j];
			}
			/*for (j=jmax-1; j>j1; j--) {
				Cellnum[cell1].R_edge[j-1]=Cellnum[cell1].R_edge_T1[j];
			}*/
			for (j=jmax-1; j>j1; j--) {
				Cellnum[cell1].GroupEdge[j-1]=Cellnum[cell1].GroupEdge_T1[j];
			}
			Cellnum[cell1].GroupVertex.pop_back();
			Cellnum[cell1].R_edge.pop_back();
			Cellnum[cell1].GroupEdge.pop_back();
			Cellnum[cell1].GroupVertex_h.pop_back();
			Cellnum[cell1].GroupVertex_M.pop_back();
			Cellnum[cell1].af.pop_back();
			
			//cell2
			jmax=Cellnum[cell2].GroupVertex.size();
			for (j=0; j<jmax; j++) {
				if (Cellnum[cell2].GroupVertex[j]==vertex2) {
					j2=j;
				}
			}
			for (j=jmax-1; j>j2; j--) {
				Cellnum[cell2].GroupVertex[j-1]=Cellnum[cell2].GroupVertex_T1[j];
			}
			/*for (j=jmax-1; j>j2; j--) {
				Cellnum[cell2].R_edge[j-1]=Cellnum[cell2].R_edge_T1[j];
			}*/
			for (j=jmax-1; j>j2; j--) {
				Cellnum[cell2].GroupEdge[j-1]=Cellnum[cell2].GroupEdge_T1[j];
			}			
			Cellnum[cell2].GroupVertex.pop_back();
			Cellnum[cell2].R_edge.pop_back();
			Cellnum[cell2].GroupEdge.pop_back();
			Cellnum[cell2].GroupVertex_h.pop_back();
			Cellnum[cell2].GroupVertex_M.pop_back();
			Cellnum[cell2].af.pop_back();
			
			//cell3 cell4
			for (int k=0; k<3; k++) {
				if (Vertices[vertex1].NeighborCell[k]!=cell1 && Vertices[vertex1].NeighborCell[k]!=cell2) {
					cell3=Vertices[vertex1].NeighborCell[k];
				}
			}
			for (int k=0; k<3; k++) {
				if (Vertices[vertex2].NeighborCell[k]!=cell1 && Vertices[vertex2].NeighborCell[k]!=cell2) {
					cell4=Vertices[vertex2].NeighborCell[k];
				}
			}
			Cellnum[cell3].GroupVertex_T1.assign(Cellnum[cell3].GroupVertex.begin(), Cellnum[cell3].GroupVertex.end());
			Cellnum[cell3].GroupEdge_T1.assign(Cellnum[cell3].GroupEdge.begin(), Cellnum[cell3].GroupEdge.end());
			Cellnum[cell3].R_edge_T1.assign(Cellnum[cell3].R_edge.begin(), Cellnum[cell3].R_edge.end());
			
			Cellnum[cell4].GroupVertex_T1.assign(Cellnum[cell4].GroupVertex.begin(), Cellnum[cell4].GroupVertex.end());
			Cellnum[cell4].GroupEdge_T1.assign(Cellnum[cell4].GroupEdge.begin(), Cellnum[cell4].GroupEdge.end());
			Cellnum[cell4].R_edge_T1.assign(Cellnum[cell4].R_edge.begin(), Cellnum[cell4].R_edge.end());
			
			//cell3
			jmax=Cellnum[cell3].GroupVertex.size();
			Cellnum[cell3].GroupVertex.push_back(int());
			Cellnum[cell3].R_edge.push_back(double());	
			Cellnum[cell3].GroupEdge.push_back(int());
			Cellnum[cell3].GroupVertex_h.push_back(Vertex());
			Cellnum[cell3].GroupVertex_M.push_back(Vertex());
			Cellnum[cell3].af.push_back(double());
			
			for (int j=0; j<jmax; j++) {
				if (Cellnum[cell3].GroupVertex[j]==vertex1) {
					j3=j;
				}
			}
			for (int j=jmax-1; j>j3; j--) {
				Cellnum[cell3].GroupVertex[j+1]=Cellnum[cell3].GroupVertex_T1[j];
			}
			/*for (int j=jmax-1; j>j3-1; j--) {
				Cellnum[cell3].R_edge[j+1]=Cellnum[cell3].R_edge_T1[j];
			}*/
			for (int j=jmax-1; j>j3-1; j--) {
				Cellnum[cell3].GroupEdge[j+1]=Cellnum[cell3].GroupEdge_T1[j];
			}
			Cellnum[cell3].GroupVertex[j3+1]=vertex2;
			Cellnum[cell3].GroupEdge[j3]=i;
			
			jmax=Cellnum[cell3].GroupVertex.size();
			for (int k=0; k<2; k++) {
				if (Edge[Cellnum[cell3].GroupEdge[(j3+1)%jmax]].Edge_Vertex[k]==vertex1) {
					Edge[Cellnum[cell3].GroupEdge[(j3+1)%jmax]].Edge_Vertex[k]=vertex2;
				}
			}
	
			//cell4
			jmax=Cellnum[cell4].GroupVertex.size();
			Cellnum[cell4].GroupVertex.push_back(int());
			Cellnum[cell4].R_edge.push_back(double());	
			Cellnum[cell4].GroupEdge.push_back(int());
			Cellnum[cell4].GroupVertex_h.push_back(Vertex());
			Cellnum[cell4].GroupVertex_M.push_back(Vertex());
			Cellnum[cell4].af.push_back(double());
			
			for (int j=0; j<jmax; j++) {
				if (Cellnum[cell4].GroupVertex[j]==vertex2) {
					j4=j;
				}
			}
			for (int j=jmax-1; j>j4; j--) {
				Cellnum[cell4].GroupVertex[j+1]=Cellnum[cell4].GroupVertex_T1[j];
			}
			/*for (int j=jmax-1; j>j4-1; j--) {
				Cellnum[cell4].R_edge[j+1]=Cellnum[cell4].R_edge_T1[j];
			}*/
			for (int j=jmax-1; j>j4-1; j--) {
				Cellnum[cell4].GroupEdge[j+1]=Cellnum[cell4].GroupEdge_T1[j];
			}	
			Cellnum[cell4].GroupVertex[j4+1]=vertex1;	
			Cellnum[cell4].GroupEdge[j4]=i;
			
			jmax=Cellnum[cell4].GroupVertex.size();
			for (int k=0; k<2; k++) {
				if (Edge[Cellnum[cell4].GroupEdge[(j4+1)%jmax]].Edge_Vertex[k]==vertex2) {
					Edge[Cellnum[cell4].GroupEdge[(j4+1)%jmax]].Edge_Vertex[k]=vertex1;
				}
			}
			
			//配置換えによる番地変更　残り
			for (int k=0; k<3; k++) {
				if (Vertices[vertex1].NeighborCell[k]==cell1) {
					Vertices[vertex1].NeighborCell[k]=cell4;
				}
			}
			
			for (int k=0; k<3; k++) {
				if (Vertices[vertex2].NeighborCell[k]==cell2) {
					Vertices[vertex2].NeighborCell[k]=cell3;
				}
			}
			
			jmax=Cellnum[cell2].GroupVertex.size();
			for (int k=0; k<3; k++) {
				if (Vertices[vertex1].NeighborVertex[k]==vertex2) {
					Vertices[vertex1].NeighborVertex[(k+1)%3]=Vertices[vertex1].NeighborVertex[(k+2)%3];
					Vertices[vertex1].NeighborVertex[(k+2)%3]=Cellnum[cell2].GroupVertex[(j2+jmax-1)%jmax];
				}
			}
			
			for (int k=0; k<3; k++) {
				if (Vertices[Cellnum[cell2].GroupVertex[(j2+jmax-1)%jmax]].NeighborVertex[k]==vertex2) {
					Vertices[Cellnum[cell2].GroupVertex[(j2+jmax-1)%jmax]].NeighborVertex[k]=vertex1;
				}
			}
			
			jmax=Cellnum[cell1].GroupVertex.size();
			for (int k=0; k<3; k++) {
				if (Vertices[vertex2].NeighborVertex[k]==vertex1) {
					Vertices[vertex2].NeighborVertex[(k+1)%3]=Vertices[vertex2].NeighborVertex[(k+2)%3];
					Vertices[vertex2].NeighborVertex[(k+2)%3]=Cellnum[cell1].GroupVertex[(j1+jmax-1)%jmax];
				}
			}
					
			for (int k=0; k<3; k++) {
				if (Vertices[Cellnum[cell1].GroupVertex[(j1+jmax-1)%jmax]].NeighborVertex[k]==vertex1) {
					Vertices[Cellnum[cell1].GroupVertex[(j1+jmax-1)%jmax]].NeighborVertex[k]=vertex2;
				}
			}

			Edge[i].Edge_Area[0]=cell3;
			Edge[i].Edge_Area[1]=cell4;
			
			Edge[i].T1_flag=-1000;
			T1_count+=1;
			cout << "T1" <<"	"<<i<<endl;
		} // T1_flag
	} //edge
	
} //T1 close


void Tissue::mutation(int p) {
    imax=Cellnum.size();
    for (i=0; i<imax; i++){
        //  Cellnum[i].Area0+=Cellnum[i].Area0_0*0.00004;
        if (Cellnum[i].cell_type==1){
            //Cellnum[i].Area0+=Cellnum[i].Area0_0*0.00004;
            if ( Cellnum[i].Area0<=Cellnum[i].Area0_0){
                //Cellnum[i].Area0+=0.00015*Cellnum[i].Area0_0;
				Cellnum[i].Area0+=0.00025*24*1.8;
            }
        }
        
        if (Cellnum[i].cell_type==2||Cellnum[i].cell_type==7||Cellnum[i].cell_type==3||Cellnum[i].cell_type==6){
            //Cellnum[i].Area0+=Cellnum[i].Area0_0*0.00004;
            if ( Cellnum[i].Area0<=Cellnum[i].Area0_0){
				//Cellnum[i].Area0+=0.00025*Cellnum[i].Area0_0;
				Cellnum[i].Area0+=0.00025*24*1.8;
            }
        }
        if (Cellnum[i].cell_type==4){
            //Cellnum[i].Area0+=Cellnum[i].Area0_0*0.00004;
            if ( Cellnum[i].Area0<=Cellnum[i].Area0_0){
               // Cellnum[i].Area0+=0.0001*Cellnum[i].Area0_0;
				Cellnum[i].Area0+=0.0001*24*3.2;
            }
        }
        if (Cellnum[i].cell_type==5){
            //Cellnum[i].Area0+=Cellnum[i].Area0_0*0.00004;
            if ( Cellnum[i].Area0<=Cellnum[i].Area0_0){
                Cellnum[i].Area0+=0.0001*Cellnum[i].Area0_0;
				
            }
        }
        
        
        /*  if (Cellnum[i].cell_type==0 || Cellnum[i].cell_type==2){
         Cellnum[i].Area0=Cellnum[i].Area0_0;
         }*/
    }
    
	//cell-type setting
	imax=Cellnum.size();
	for (i=0; i<imax; i++){
		jmax=Cellnum[i].GroupEdge.size();
		   if (Cellnum[i].cell_type==7){
			   Cellnum[i].cell_type=3;
			   for (j=0; j<jmax; j++) {
			   if((Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==2 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==2)){//4
				   Cellnum[i].cell_type=7;
			   }
			   }
		   }
		   }
	
		imax=Cellnum.size();
		for (i=0; i<imax; i++){
		jmax=Cellnum[i].GroupEdge.size();
		if (Cellnum[i].cell_type==3 || Cellnum[i].cell_type==6){
			Cellnum[i].cell_type=3;
			for (j=0; j<jmax; j++) {
			if(Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==4){//4
				Cellnum[i].cell_type=6;
			}
			}
		}
		}
	
	
        /*imax=Cellnum.size();
        for (i=0; i<imax; i++){
            if (Cellnum[i].cell_type==5){
                Cellnum[i].Area0+=Cellnum[i].Area0*0.000001;
            }
        }*/

   /* imax=Cellnum.size();
    for (i=0; i<imax; i++){
        if (Cellnum[i].cell_type==2){
            Cellnum[i].Area0+=Cellnum[i].Area0*0.000001;
        }
    }*/
    // Cellnum[i].cell_type=2;
    
    //cell-collappse
	if (p>9000&&p<9800) {
		//Cellnum[1].k_area=k_1*0.005;
		/*if (Cellnum[1].Area>Cellnum[1].Area0_0*0.001) {
            //Cellnum[1].k_area=Cellnum[1].k_area*0;
            Cellnum[1].Area0-=0.002*Cellnum[1].Area0_0;
            Cellnum[1].k_api=k_3*0;
		}*/
		
		/*if (Cellnum[11].Area>Cellnum[11].Area0_0*0.001) {
            //Cellnum[1].k_area=Cellnum[1].k_area*0;
            //Cellnum[27].Area0-=0.001*Cellnum[27].Area0_0;
			Cellnum[11].Area0-=0.005*Cellnum[11].Area0;
			 //Cellnum[27].Area0-=0.001*Cellnum[27].Area0_0;
			//Cellnum[27].k_area=k_1*0;
            //Cellnum[27].k_api=k_3*0;
		}*/
		
		//Cellnum[0].k_api=k_3*0;
		//Cellnum[1].Area0=Cellnum[1].Area0*0.995;
		//Cellnum[1].k_bas=k_4*-5;
		//Cellnum[1].k_lat=k_5*-5;
		//Cellnum[2].k_lat=k_5*-5;
	}
		
	/*if (p>110000) {
		Cellnum[2].k_bas=k_4;		
	}
	if (p>120000) {
		Cellnum[2].k_lat=k_5;		
	}	*/
}
void Tissue::output(int a) {
	stringstream filename;
	filename << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_cellform" << a << ".dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs[a].open(filename.str());
	jmax=Cellnum.size();
	for (j=0; j<26; j++) {
		kmax=Cellnum[j].GroupVertex.size();//
		for (k=0; k<kmax; k++) {
			ofs[a] <<	Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] << 
			"	" << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1] << "	";
		}
		ofs[a] << endl;
	}
	for (j=27; j<jmax; j++) {
		kmax=Cellnum[j].GroupVertex.size();//
		for (k=0; k<kmax; k++) {
			ofs[a] <<	Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] << 
			"	" << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1] << "	";
		}
		ofs[a] << endl;
	}
	ofs[a].close();
	
	//Auxin-CYtokinin
	/*stringstream filename01;
	filename01 << "/Users/fujiwara/Desktop/vascul_model/data_form-a-c/zb_au-cy" << a << ".dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs[a].open(filename01.str());
	jmax=Cellnum.size();//
	ofs[a] <<	Auxin_ave	<< "	" <<	Cytokinin_ave;
	ofs[a] << endl;
	for (j=0; j<27; j++) {
		ofs[a] <<	Cellnum[j].Auxin	<< "	" <<	Cellnum[j].Cytokinin;
		ofs[a] << endl;
	}
	for (j=28; j<jmax; j++) {
		ofs[a] <<	Cellnum[j].Auxin	<< "	" <<	Cellnum[j].Cytokinin;
		ofs[a] << endl;
	}
	ofs[a].close();*/
	
	stringstream filename02;
	filename02 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_celltype" << a << ".dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs[a].open(filename02.str());
	jmax=Cellnum.size();
	for (j=0; j<26; j++) {
		ofs[a] <<	Cellnum[j].cell_type<<"	"<<	Cellnum[j].cell_lineage<<"	"<<	Cellnum[j].cell_defect;
		ofs[a] << endl;
	}
	for (j=27; j<jmax; j++) {
		ofs[a] <<	Cellnum[j].cell_type<<"	"<<	Cellnum[j].cell_lineage<<"	"<<	Cellnum[j].cell_defect;
		ofs[a] << endl;
	}
	ofs[a].close();
	
    stringstream filename03;
    filename03 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_stresstensor" << a << ".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs[a].open(filename03.str());
    jmax=Cellnum.size();
    for (j=0; j<26; j++) {
        ofs[a] <<    Cellnum[j].S1 <<" " << Cellnum[j].S2 <<" " <<   Cellnum[j].Theta1 <<" " << Cellnum[j].Theta2;
        ofs[a] << endl;
    }
    for (j=27; j<jmax; j++) {
       ofs[a] <<    Cellnum[j].S1 <<" " << Cellnum[j].S2 <<" " <<   Cellnum[j].Theta1 <<" " << Cellnum[j].Theta2;
        ofs[a] << endl;
    }
    ofs[a].close();
    
    stringstream filename04;
    filename04 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_centroid" << a << ".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs[a].open(filename04.str());
    jmax=Cellnum.size();
    for (j=0; j<26; j++) {
        ofs[a] <<    Cellnum[j].centroid[0] <<" " << Cellnum[j].centroid[1];
        ofs[a] << endl;
    }
    for (j=27; j<jmax; j++) {
        ofs[a] <<    Cellnum[j].centroid[0] <<" " << Cellnum[j].centroid[1];
        ofs[a] << endl;
    }
    ofs[a].close();
    
	//celltype-area
	 stringstream filename05;
    //filename05 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_cellarea" << a << ".dat";// 0tukerutoki setw(4) << setfill('0') <<
	//filename05 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_cellarea.dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs[a].open("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_cellarea.dat",ios::app);
    jmax=Cellnum.size();//
    for (j=0; j<26; j++) {
		   ofs[a] <<    Cellnum[j].cell_type <<","<<j <<" " << Cellnum[j].Area<<" ";
	   }
	   for (j=27; j<jmax; j++) {
		   ofs[a] <<    Cellnum[j].cell_type <<","<< j-1 <<" " << Cellnum[j].Area<<" ";
	   }
		ofs[a] << endl;
	   ofs[a].close();
    
}

void Tissue::output_result() {
	stringstream filename0;
	filename0 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/vertex_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[0].open(filename0.str());
	jmax=Vertices.size();//
	for (j=0; j<jmax; j++) {
		ofs_result[0] <<	Vertices[j].coordinate[0] << 
		"	" << Vertices[j].coordinate[1] << "	";
		ofs_result[0] << endl;
	}
	ofs_result[0].close();
	
	stringstream filename1;
	filename1 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/cell-vertex_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[1].open(filename1.str());
	jmax=Cellnum.size();//
	for (j=0; j<jmax; j++) {
		kmax=Cellnum[j].GroupVertex.size();
		for (k=0; k<kmax; k++) {
			ofs_result[1] << Cellnum[j].GroupVertex[k] << "	" ;
		}
		ofs_result[1] << endl;
	}
	ofs_result[1].close();
	
	stringstream filename2;
	filename2 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/neigborcell_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[2].open(filename2.str());
	jmax=Vertices.size();//
	for (j=0; j<jmax; j++) {
		kmax=Vertices[j].NeighborCell.size();//
		for (k=0; k<kmax; k++) {
			ofs_result[2] << Vertices[j].NeighborCell[k] << "	"  ;
		}
		ofs_result[2] << endl;
	}
	ofs_result[2].close();
	
	stringstream filename3;
	filename3 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/neigborvertex_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[3].open(filename3.str());
	jmax=Vertices.size();//
	for (j=0; j<jmax; j++) {
		kmax=Vertices[j].NeighborVertex.size();//
		for (k=0; k<kmax; k++) {
			ofs_result[3] << Vertices[j].NeighborVertex[k] << "	"  ;
		}
		ofs_result[3] << endl;
	}
	ofs_result[3].close();
	
	stringstream filename4;
	filename4 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/cell-cordinate_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[4].open(filename4.str());
	jmax=Cellnum.size();//
	for (j=0; j<jmax; j++) {
		kmax=Cellnum[j].GroupVertex.size();
		for (k=0; k<kmax; k++) {
			ofs_result[4] <<Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] <<
			"	" << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1] << "	";
		}
		ofs_result[4] << endl;
	}
	ofs_result[4].close();
	
	stringstream filename5;
	filename5 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/cellArea_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[5].open(filename5.str());
	jmax=Cellnum.size();//
	for (j=0; j<jmax; j++) {
		ofs_result[5] << Cellnum[j].Area ;
		ofs_result[5] << endl;
	}
	ofs_result[5].close();
	
	stringstream filename6;
	filename6 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/cell-edge_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[6].open(filename6.str());
	jmax=Cellnum.size();//
	for (j=0; j<jmax; j++) {
		kmax=Cellnum[j].GroupEdge.size();
		for (k=0; k<kmax; k++) {
			ofs_result[6] << Cellnum[j].GroupEdge[k] << "	" ;
		}
		ofs_result[6] << endl;
	}
	ofs_result[6].close();
	
	stringstream filename7;
	filename7 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/edge-vertex_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[7].open(filename7.str());
	jmax=Edge.size();//
	for (j=0; j<jmax; j++) {
		kmax=2;//
		for (k=0; k<kmax; k++) {
			ofs_result[7] << Edge[j].Edge_Vertex[k] << "	"  ;
		}
		ofs_result[7] << endl;
	}
	ofs_result[7].close();
	
	stringstream filename8;
	filename8 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/edge-cell_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[8].open(filename8.str());
	jmax=Edge.size();//
	for (j=0; j<jmax; j++) {
		kmax=2;//
		for (k=0; k<kmax; k++) {
			ofs_result[8] << Edge[j].Edge_Area[k] << "	"  ;
		}
		ofs_result[8] << endl;
	}
	ofs_result[8].close();

	stringstream filename9;
	filename9 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/celltype_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_result[9].open(filename9.str());
	jmax=Cellnum.size();//
	for (j=0; j<jmax; j++) {
		ofs_result[9] << Cellnum[j].cell_type ;
		ofs_result[9] << endl;
	}
	ofs_result[9].close();
    
    stringstream filename10;
    filename10 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/result2/cell-centroid_result.dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_result[10].open(filename10.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
        ofs_result[10] << Cellnum[j].centroid[0]<<
        "   " << Cellnum[j].centroid[1]<<  "  ";
        ofs_result[10] <<endl;
    }
    ofs_result[10].close();
	
}

void Tissue::output_barometer(int c) {
	int xylem_num;
	int procambium_num;
	int pericycel_num;
	int procambium_xylem_num;
	int xylem_edge;
	int pericycle_edge;
	
	xylem_num=0;
	procambium_num=0;
	pericycel_num=0;
	procambium_xylem_num=0;
	xylem_edge=0;	
	pericycle_edge=0;
	
	//cell-number
	stringstream filename10;
	filename10 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/cell_number"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_barometer[c].open(filename10.str());
	
	imax=Cellnum.size();//
	for (i=0; i<imax; i++) {
		if (Cellnum[i].cell_type==1) {
			xylem_num+=1;
		}
		if (Cellnum[i].cell_type==2 ||Cellnum[i].cell_type==3 || Cellnum[i].cell_type==7 ||Cellnum[i].cell_type==6) {
			procambium_num+=1;
			jmax=Cellnum[i].GroupEdge.size();
			for (j=0; j<jmax; j++) {	
				if (Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[0]].cell_type==1 || Cellnum[Edge[Cellnum[i].GroupEdge[j]].Edge_Area[1]].cell_type==1){
					if (xylem_edge==0) {
						procambium_xylem_num+=1;
						xylem_edge+=1;
					}
				}
			}
			xylem_edge=0;
		}
		if (Cellnum[i].cell_type==4) {
			pericycel_num+=1;
		}
	}
		ofs_barometer[c] <<	xylem_num << endl;
		ofs_barometer[c] <<	procambium_num << endl;
		ofs_barometer[c] <<	pericycel_num << endl;
		ofs_barometer[c] <<	procambium_xylem_num << endl;
	
	ofs_barometer[c].close();
	
	//xylem-coordinate
	stringstream filename11;
	filename11 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/xylem-vertex"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_barometer[c].open(filename11.str());
	jmax=Cellnum.size();//
	for (j=0; j<jmax; j++) {
		if (Cellnum[j].cell_type==1) {
			ofs_barometer[c] << j << endl;
			kmax=Cellnum[j].GroupVertex.size();
			for (k=0; k<kmax; k++) {
				ofs_barometer[c] <<Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] <<
				"	" << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1]<< endl;
			}
		ofs_barometer[c] << endl;
		}
	}
	ofs_barometer[c].close();
	
	//xylem-center-area
	stringstream filename12;
	filename12 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/xylem-center-area"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
	ofs_barometer[c].open(filename12.str());
	jmax=Cellnum.size();//
	for (j=0; j<jmax; j++) {
		if (Cellnum[j].cell_type==1) {
			ofs_barometer[c] << j << endl;
			ofs_barometer[c] << Cellnum[j].Area <<
            "	" << Cellnum[j].centroid[0]<<"	" << Cellnum[j].centroid[1]<<
			"	" <<Cellnum[j].S1 <<"	" << Cellnum[j].S2 <<"	" <<Cellnum[j].Theta1 <<"	"<<Cellnum[j].Theta2<<endl;
			ofs_barometer[c] << endl;
		}
	}
	ofs_barometer[c].close();
    
    //pericycle-coordinate
    stringstream filename13;
    filename13 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/pericycle-vertex"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename13.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
        if (Cellnum[j].cell_type==4) {
            ofs_barometer[c] << j << endl;
            kmax=Cellnum[j].GroupVertex.size();
            for (k=0; k<kmax; k++) {
                ofs_barometer[c] <<Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] <<
                "    " << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1]<< endl;
            }
            ofs_barometer[c] << endl;
        }
    }
    ofs_barometer[c].close();
    
    //pericycle-center-area
    stringstream filename14;
    filename14 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/pericycle-center-area"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename14.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
        if (Cellnum[j].cell_type==4) {
            ofs_barometer[c] << j << endl;
            ofs_barometer[c] << Cellnum[j].Area <<
            "    " << Cellnum[j].centroid[0]<<
            "    " << Cellnum[j].centroid[1]<< endl;
            ofs_barometer[c] << endl;
        }
    }
    ofs_barometer[c].close();
    
    //endodermis-coordinate
    stringstream filename15;
    filename15 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/endodermis-vertex"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename15.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
        if (Cellnum[j].cell_type==5) {
            ofs_barometer[c] << j << endl;
            kmax=Cellnum[j].GroupVertex.size();
            for (k=0; k<kmax; k++) {
                ofs_barometer[c] <<Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] <<
                "    " << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1]<< endl;
            }
            ofs_barometer[c] << endl;
        }
    }
    ofs_barometer[c].close();
    
    //endodermis-center-area
    stringstream filename16;
    filename16 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/endodermis-center-area"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename16.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
        if (Cellnum[j].cell_type==5) {
            ofs_barometer[c] << j << endl;
            ofs_barometer[c] << Cellnum[j].Area <<
            "    " << Cellnum[j].centroid[0]<<
            "    " << Cellnum[j].centroid[1]<< endl;
            ofs_barometer[c] << endl;
        }
    }
    ofs_barometer[c].close();
    
    //procambium-boundary-coordinate
    stringstream filename17;
    filename17 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/procambium-boundary-vertex"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename17.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
          if (Cellnum[j].cell_type==2 ||Cellnum[j].cell_type==3 || Cellnum[j].cell_type==7||Cellnum[i].cell_type==6) {
            kmax=Cellnum[j].GroupVertex.size();
                 for (k=0; k<kmax; k++) {
                    if (Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[0]].cell_type==1 || Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[1]].cell_type==1){
                        if (xylem_edge==0) {
                            ofs_barometer[7] << j << endl;
                              for (k=0; k<kmax; k++) {
                                ofs_barometer[c] <<Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] <<
                                "    " << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1]<< endl;
                            }
                            xylem_edge+=1;
                            ofs_barometer[c] << endl;
                        }
                    }
                }
                xylem_edge=0;
            }
        }
    ofs_barometer[c].close();
    
    //procambium-boundary-center-area
    stringstream filename18;
    filename18 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/procambium-boundary-center-area"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename18.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
        if (Cellnum[j].cell_type==2 ||Cellnum[j].cell_type==3 || Cellnum[j].cell_type==7||Cellnum[i].cell_type==6) {
            kmax=Cellnum[j].GroupVertex.size();
            for (k=0; k<kmax; k++) {
                if (Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[0]].cell_type==1 || Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[1]].cell_type==1){
                    if (xylem_edge==0) {
                        ofs_barometer[c] << j << endl;
                        ofs_barometer[c] << Cellnum[j].Area <<
                        "	" << Cellnum[j].centroid[0]<<"	" << Cellnum[j].centroid[1]<<
						"	" <<Cellnum[j].S1 <<"	" << Cellnum[j].S2 <<"	" <<Cellnum[j].Theta1 <<"	"<<Cellnum[j].Theta2<<endl;
                        ofs_barometer[c] << endl;
                        
                        xylem_edge+=1;
                    }
            }
            }
            xylem_edge=0;
        }
    }
    ofs_barometer[c].close();
    
    //procambium-coordinate
    stringstream filename19;
    filename19 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/procambium-vertex"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename19.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
        if (Cellnum[j].cell_type==2 ||Cellnum[j].cell_type==3 || Cellnum[j].cell_type==7||Cellnum[i].cell_type==6) {
            ofs_barometer[c] << j << endl;
            kmax=Cellnum[j].GroupVertex.size();
            for (k=0; k<kmax; k++) {
                ofs_barometer[c] <<Vertices[Cellnum[j].GroupVertex[k]].coordinate[0] <<
                "    " << Vertices[Cellnum[j].GroupVertex[k]].coordinate[1]<< endl;
            }
            ofs_barometer[c] << endl;
        }
    }
    ofs_barometer[c].close();
    
    //procambium-center-area
    stringstream filename20;
    filename20 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/procambium-center-area"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename20.str());
    jmax=Cellnum.size();//
    for (j=0; j<jmax; j++) {
       if (Cellnum[j].cell_type==2 ||Cellnum[j].cell_type==3 || Cellnum[j].cell_type==7||Cellnum[i].cell_type==6) {
            ofs_barometer[c] << j << endl;
            ofs_barometer[c] << Cellnum[j].Area <<
            "    " << Cellnum[j].centroid[0]<<
            "    " << Cellnum[j].centroid[1]<< endl;
            ofs_barometer[c] << endl;
        }
    }
    ofs_barometer[c].close();
    
    //procambium-boundary-stress
    stringstream filename21;
    filename21 << "/Users/fujiwara/Desktop/vertexmodel/vascularsimu/barometer-initial2/procambium-boundary-stress"<<c<<".dat";// 0tukerutoki setw(4) << setfill('0') <<
    ofs_barometer[c].open(filename21.str());
    jmax=Cellnum.size();//
	//PSE-PSE
	for (j=0; j<jmax; j++) {
       if (Cellnum[j].cell_type==2) {
            ofs_barometer[c] << Cellnum[j].centroid[0]<<"	" << Cellnum[j].centroid[1]<< endl;
        }
    }
	
    for (j=0; j<jmax; j++) {
        if (Cellnum[j].cell_type==2 ||Cellnum[j].cell_type==3 || Cellnum[j].cell_type==7 ||Cellnum[i].cell_type==6) {
            kmax=Cellnum[j].GroupEdge.size();
            for (k=0; k<kmax; k++) {
                if (Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[0]].cell_type==1 || Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[1]].cell_type==1){
                        xylem_edge+=1;
				}
				if (Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[0]].cell_type==4 || Cellnum[Edge[Cellnum[j].GroupEdge[k]].Edge_Area[1]].cell_type==4){
						pericycle_edge+=1;
				}
            }
			if (xylem_edge>0 && pericycle_edge==0) {
					ofs_barometer[c] << Cellnum[j].centroid[0]<< "	" << Cellnum[j].centroid[1]<<"	" <<Cellnum[j].S1 <<"	" << Cellnum[j].S2 <<"	" << Cellnum[j].Theta1 <<"	" << Cellnum[j].Theta2;
					ofs_barometer[c] << endl;
			}
            xylem_edge=0;
			pericycle_edge=0;
        }
    }
    ofs_barometer[c].close();
	
	
}

int main()
{
	srand48((unsigned int)time(NULL));
	 //srand((unsigned int)time(NULL));
	int i,imax;
	int j,jmax;
	int k,kmax;
	int p,q,t,a=0;
	double b=0;
     double c=0;
	Cell cell;
	Vertex vert;
	Tissue tiss;
	
	tiss.firststage();
	tiss.cellgroup();
	// cout << "bbpb" <<endl;
	tiss.output(a);
	a++;
	
	/*imax=tiss.Vertices.size();//
	for (i=0; i<imax; i++) {
		tiss.Vertices[i].show_vertex();
	}
	jmax=tiss.Cellnum.size();//
	cout << jmax << endl;
	for (j=0; j<jmax; j++) {
		kmax=tiss.Cellnum[j].GroupVertex.size();//
		cout << kmax << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].GroupVertex[k] << " " ;
		}
		cout << endl;	
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].GroupEdge[k] << " " ;
		}
		cout << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Vertex[0] << " " << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Vertex[1] << "  ";
		}
		cout << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Area[0] << " " << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Area[1] << "  ";
		}		
		cout << endl;
	}
	
	jmax=tiss.Vertices.size();//
	for (j=0; j<jmax; j++) {
		kmax=tiss.Vertices[j].NeighborCell.size();//
		for (k=0; k<kmax; k++) {
			cout << tiss.Vertices[j].NeighborCell[k] << " "  ;
		}
		cout << endl;
	}
	cout << endl;*/
	//cout << "bbbb" <<endl;
	/*jmax=tiss.Vertices.size();//
	for (j=0; j<jmax; j++) {
		kmax=tiss.Vertices[j].NeighborVertex.size();//
		for (k=0; k<kmax; k++) {
			cout << tiss.Vertices[j].NeighborVertex[k] << " "  ;
		}
		cout << endl;
	}*/
	
	tiss.firstset_cell();
	
	// cout << "bbbbb" <<endl;
	
	/*imax=tiss.Vertices.size();//
	for (i=0; i<imax; i++) {
		tiss.Vertices[i].show_vertex();
	}*/
	
	/*jmax=tiss.Cellnum.size();//
	 for (j=0; j<jmax; j++) {
	 kmax=2;//
	 for (k=0; k<kmax; k++) {
	 cout << tiss.Cellnum[j].centroid[k] << " " << tiss.Cellnum[j].centroid[k] << " " ;
	 }
	 cout << endl;
	 }*/
	
	/*jmax=tiss.Cellnum.size();//
	 for (j=0; j<jmax; j++) {
	 cout << tiss.Cellnum[j].Area0;
	 cout << endl;
	 }*/
	
	/*jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		cout << tiss.Cellnum[j].time_division;
		cout << endl;
	}*/
	
	/*tiss.Geometory_cell();
	jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		cout << tiss.Cellnum[j].Mxx << " "  ;
		cout << endl;
	}
	cout << endl;
	
	jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		cout << tiss.Cellnum[j].Myy << " "  ;
		cout << endl;
	}
	cout << endl;
	jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		cout << tiss.Cellnum[j].Mxy << " "  ;
		cout << endl;
	}
	cout << endl;*/
    
    for (p=0; p<pmax*100; p++) {
        tiss.VertexDynamics();
        tiss.Geometory_cell();
    }
    for (p=0; p<50; p++) {
	tiss.output(a);
	a++;
    }
	cout <<"start "<<endl;
	//step by step
	for (t=0; t<tmax; t++){
	for (q=0; q<qmax; q++) {
		for (p=0; p<pmax; p++) {
		tiss.VertexDynamics();
		tiss.Geometory_cell();
		}
		cout <<"time"<<endl;
		tiss.Division(q+qmax*t);
		tiss.Geometory_cell();
		
		//tiss.T1_cell();
		//tiss.Geometory_cell();
		
		tiss.mutation(q+qmax*t);
		tiss.Geometory_cell();
		
		if ((int)b==tmax) {
		tiss.output(a);
		a++;
		b=0;
		}
		b+=0.25;
        
        if (q==2500*c) {
            tiss.output_barometer(c);
            c++;
        }
	}
	}
	
	b=0;
    for (q=0; q<qmax*0.01; q++) {
        for (p=0; p<pmax; p++) {
            tiss.VertexDynamics();
            tiss.Geometory_cell();
        }
        if (b==1) {
            tiss.output(a);
            a++;
            b=0;
        }
        b+=0.25;
    }
    
	tiss.output_result();
	tiss.output_barometer(c);
	/*imax=tiss.Vertices.size();//
	for (i=0; i<imax; i++) {
		tiss.Vertices[i].show_vertex();
	}*/
	
	/*jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		kmax=tiss.Cellnum[j].GroupVertex.size();//
		for (k=0; k<kmax; k++) {
			tiss.Vertices[tiss.Cellnum[j].GroupVertex[k]].show_cell();
		}
		cout << endl;
	}*/	
	
	/*jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		kmax=2;//
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].centroid[k] << " " << tiss.Cellnum[j].centroid[k] << " " ;
		}
		cout << endl;
	}*/

	/*jmax=tiss.Cellnum.size();//
	 for (j=0; j<jmax; j++) {
	 kmax=tiss.Cellnum[j].GroupVertex.size();//
	 for (k=0; k<kmax; k++) {
	 cout << tiss.Vertices[tiss.Cellnum[j].GroupVertex[k]].Fx << " " << tiss.Vertices[tiss.Cellnum[j].GroupVertex[k]].Fy << " " ;
	 }
	 cout << endl;
	 }*/

	/*jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		kmax=tiss.Cellnum[j].GroupVertex.size();//
		for (k=0; k<kmax; k++) {
			cout << tiss.Vertices[tiss.Cellnum[j].GroupVertex[k]].Fx << " " << tiss.Vertices[tiss.Cellnum[j].GroupVertex[k]].Fy << " " ;
		}
		cout << endl;
	}*/
	
	/*jmax=tiss.Cellnum.size();//
	for (j=0; j<jmax; j++) {
		kmax=tiss.Cellnum[j].R_edge.size();//
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].R_edge[k] << " "  ;
		}
		cout << endl;
	}*/

	/*jmax=tiss.Vertices.size();//
	for (j=0; j<jmax; j++) {
		cout << j << endl;
		kmax=tiss.Vertices[j].NeighborCell.size();//
		for (k=0; k<kmax; k++) {
			cout << tiss.Vertices[j].NeighborCell[k] << " "  ;
		}
		cout << endl;
		kmax=tiss.Vertices[j].NeighborVertex.size();//
		for (k=0; k<kmax; k++) {
			cout << tiss.Vertices[j].NeighborVertex[k] << " "  ;
		}
		cout << endl;		
	}
	cout << endl;*/
	
	/*jmax=tiss.Cellnum.size();//
	cout << jmax << endl;
	for (j=0; j<jmax; j++) {
		kmax=tiss.Cellnum[j].GroupVertex.size();//
		cout << kmax << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].GroupVertex[k] << " " ;
		}
			cout << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].GroupEdge[k] << " " ;
		}
		cout << endl;
	}*/
	
	/*jmax=tiss.Cellnum.size();//
	cout << jmax << endl;
	for (j=0; j<jmax; j++) {
		kmax=tiss.Cellnum[j].GroupVertex.size();//
		cout << kmax << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].GroupVertex[k] << " " ;
		}
		cout << endl;	
		for (k=0; k<kmax; k++) {
			cout << tiss.Cellnum[j].GroupEdge[k] << " " ;
		}
		cout << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Vertex[0] << " " << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Vertex[1] << "  ";
		}
		cout << endl;
		for (k=0; k<kmax; k++) {
			cout << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Area[0] << " " << tiss.Edge[tiss.Cellnum[j].GroupEdge[k]].Edge_Area[1] << "  ";
		}		
		cout << endl;
	}
	
	jmax=tiss.Vertices.size();//
	 for (j=0; j<jmax; j++) {
	 kmax=tiss.Vertices[j].NeighborCell.size();//
	 for (k=0; k<kmax; k++) {
	 cout << tiss.Vertices[j].NeighborCell[k] << " "  ;
	 }
	 cout << endl;
	 }
	 cout << endl;
	 cout << "bb" <<endl;
	 jmax=tiss.Vertices.size();//
	 for (j=0; j<jmax; j++) {
	 kmax=tiss.Vertices[j].NeighborVertex.size();//
	 for (k=0; k<kmax; k++) {
	 cout << tiss.Vertices[j].NeighborVertex[k] << " "  ;
	 }
	 cout << endl;
	 }
	
	jmax=tiss.Cellnum.size();//
	 for (j=0; j<jmax; j++) {
		 if (tiss.Cellnum[j].cell_type==1) {//
		 for (k=0; k<2; k++) {
		 cout << tiss.Cellnum[j].centroid[k] << " "  ;
		 }
			 cout << endl;
		 }
	 }*/
	
	
	//cout << tiss.T1_count<< endl;
	cout << "stop"<< endl;
	return 0;
}

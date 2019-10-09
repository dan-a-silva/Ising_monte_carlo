#ifndef MonteCarlo3D_HPP
#define MonteCarlo3D_HPP

#include <ctime>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <random>
#include <iostream>
#include <iomanip>
#include <map>
#include <unordered_map>
//#include <vector>


using namespace std;

typedef vector<int> Row;
typedef vector<Row> Matrix;
typedef vector<double> Roww;
typedef vector<Roww> Matrixx;
typedef vector<Matrix> Matrix3d;

typedef struct {
    //norm = x*x + y*y + z*z
    int x,y,z;
    double pref;//pref = exp(-pi*pi*s0*s0*norm)/(pi*norm)
    int type;
    double rad;
} xyz;

typedef struct xyzs xyzs;


struct xyzs {
    //norm = x*x + y*y + z*z
    int x,y,z;
    double energ;
    int spin;
    xyzs **nbrs;
};


class MonteCarlo3d{

    Matrix old, newm;
    //  vector<xyz> nbr;
    Matrixx fccprimitve;
    int h, obs_mag ,inter, change_e, change_m, Length ;
    double T, del, mag,alph,obs_en ,total_m , total_e, total_sm;
    //vector<xyzs> nmbs;
    unordered_map <int,double> nmbs;
    vector<double> enrgies;
    default_random_engine e1;
  public:
//    void generate_coord(vector<xyz> & v2, vector<xyz> & nbr,int z, int choice, double aLatt);
    vector<int> return_neighbors ( int x, int y, int z, int L);
    double observable_e();
    double observable_m();
    double observable_sm();
    void initialize( int L, double h, int J,double d, double alpha, vector<xyzs> & spins, vector<xyzs*>& protons, int choice, unordered_map< int , double> evenns);
    int total_mag(vector<xyzs> & spins);
    int stag_mag(vector<xyzs> & spins);
    double total_energy(vector<xyzs> & spins, vector<xyzs*> protons,double alphz, int L);
    void status(vector<xyzs> & spins);
    void update(xyzs*  a , vector<xyzs*> & protons,int spiny, double dele );
    void propose(vector<xyzs> & spins , vector <xyzs*> & protons, double T, int L ,int sx, double alphz);
    double change(xyzs *a, vector<xyzs> & spins, vector<xyzs*> & protons , int spiny, int sx, int L,double alphz);
    double electro_en(xyzs* a,vector <xyzs*> & protons, int L);
    double dotproduct(vector<double> v1, vector<double> v2);
};

void MonteCarlo3d::initialize (int L, double h, int J,double d, double alpha, vector<xyzs> & spins, vector<xyzs*> & protons  ,int choice, unordered_map<int, double> evenns)
{
    //random_device r;               // Seed with a real random value, if available
    e1 = default_random_engine(time(0));
    int V = L*L*L;
    Length = L;
    mag = h;
    inter= J;
    del = d;
    alph = alpha;
    
    nmbs  = evenns;
    
    int j = 0;
    for(unsigned int ix=0; ix<L; ix++){
        for(unsigned int iy=0; iy<L; iy++){
            for(unsigned int iz=0; iz<L; iz++){
                spins[j].x =  ix;
                spins[j].y =  iy;
                spins[j].z =  iz;
//                (j < L ? spins[j].spin = 1 : spins[j].spin = -1);
                spins[j].spin = 1;
                j++;
        }
        }
    }
    
//    cout << "j = " << j << endl;

    
    
    for (size_t i = 0 ; i< spins.size(); i++)
    {
        spins[i].nbrs = new xyzs*[6];
        vector<int> x = return_neighbors( spins[i].x , spins[i].y, spins[i].z, L);
        for( size_t j =0; j< x.size(); j++){
            spins[i].nbrs[j] = &(spins[x[j]]);
        }
        if(spins[i].spin == 1) { protons.push_back(&spins[i]);}
    }
    total_m = total_mag(spins);
    total_sm = stag_mag(spins);
    total_e = total_energy(spins,protons,0,L);
    cout << "mag = "  << total_m << "\n";

}

vector<int>  MonteCarlo3d::return_neighbors ( int x, int y, int z, int L){
    vector<int> rtrn(6);
    int xn,yn=y,zn=z;
    xn = (x == 0 ? L-1 : x-1);
    rtrn[0] = xn * L * L + yn * L + zn;
    xn = (x == L-1 ? 0 : x+1);
    rtrn[1] = xn * L * L + yn * L + zn;
    xn =x;
    yn = (y == 0 ? L-1 : y-1);
    rtrn[2] = xn * L * L + yn * L + zn;
    yn = (y == L-1 ? 0 : y+1);
    rtrn[3] = xn * L * L + yn * L + zn;
    yn =y;
    zn = (z == 0 ? L-1 : z-1);
    rtrn[4] = xn * L * L + yn * L + zn;
    zn = (z == L-1 ? 0 : z+1);
    rtrn[5] = xn * L * L + yn * L + zn;
    
    return rtrn;
}


double MonteCarlo3d::observable_m(){
	return total_m;
}

double MonteCarlo3d::observable_sm(){
    return total_sm;
}

double MonteCarlo3d::observable_e(){
	return total_e;
}

int MonteCarlo3d::total_mag( vector<xyzs> & spins){

	int magn = 0;

    for (auto it = spins.begin(); it != spins.end(); ++it)
    {
        magn += it->spin;
    }
	obs_mag = magn;
	//std::cout<<magn<<" \n";
	return magn;
}

int MonteCarlo3d::stag_mag( vector<xyzs> & spins){

	int stag_magn = 0;
    
    for (auto it = spins.begin(); it != spins.end(); ++it)
    {
        stag_magn += (it->spin)*( (it->x + it->y + it->z )%2 ? -1:1 );
    }

    return stag_magn;
}

double MonteCarlo3d::total_energy(vector<xyzs> & spins, vector<xyzs*>  protons,double alphz, int L){

    double energ  = 0;
    double energy = 0;
    
    for (auto it = spins.begin(); it != spins.end(); ++it)
    {

        if(alphz != 0 && it->spin == 1){
            energy = electro_en( &(*it),protons,L);
        }
        
        energ += - inter* (it->spin)*(it->nbrs[0]->spin + it->nbrs[1]->spin + it->nbrs[2]->spin + it->nbrs[3]->spin + it->nbrs[4]->spin +it->nbrs[5]->spin) - mag*(it->spin) + del * (it->spin*it->spin)+ alphz*energy ;
    }

    obs_en = energ;
    return energ;
}

void MonteCarlo3d::status(vector<xyzs> &spins){
	 std::cout<<"The Lattice is : \n";
    for (auto it = spins.begin(); it != spins.end(); ++it)
    {
        cout << "x,y,z:spin" << it->x << ","<< it->y << ","<< it->z << ":" << it->spin << endl;
    }
	   
    std::cout<< "The energy is : " << obs_en << " the magnetization is : " <<  obs_mag << "\n";
}

double MonteCarlo3d::electro_en( xyzs * a ,vector<xyzs*> & protons, int L){
    vector<double> d(3);
    int d12;
    double en12 = 0;
    double less,more;
    int cnt;
    
    en12 = 0;
    
    if(a->spin == 1){
        for ( vector<xyzs*>::iterator it = protons.begin(); it != protons.end(); it++){
//        vector<xyzs*>::iterator it = protons.begin();
        if( *it != a ){
                //vector<double> tmp{(double)abs(it->x - a->x ), (double)abs(it->y - a->y) , (double)abs(it->z- a->z)};
            int xt,yt,zt;
            xt =  abs( (*it)->x - a->x );
            yt =abs((*it)->y - a->y) ;
            zt = abs((*it)->z- a->z);
            if (zt < xt) std::swap(zt,xt);
            if (zt < yt) std::swap(zt,yt);
            if (yt < xt) std::swap(xt,yt);
            int tmp  = xt* L*L + yt*L+ zt;
            
            en12 += nmbs[tmp];
//            }
//        }
    }
    }
    }
//    }
    return en12;
}

double MonteCarlo3d::dotproduct(vector<double> v1, vector<double> v2){
	double dott=0;
	for(unsigned int i=0; i<3; i++){
		dott += v1[i]*v2[i];
	}
	return dott;
}

double MonteCarlo3d::change( xyzs *a , vector<xyzs> & spins, vector<xyzs*> & protons ,int spiny, int sx, int L,double alphz){

    double change = 0;
    double energy = 0;
    double neighbors = a->nbrs[0]->spin + a->nbrs[1]->spin + a->nbrs[2]->spin + a->nbrs[3]->spin + a->nbrs[4]->spin + a->nbrs[5]->spin;
    if(alphz!=0){
        energy = electro_en(a,protons,L);
    }
    
    if(sx == 0){
        change = -inter*(spiny - a->spin)*neighbors - mag*(spiny - a->spin) - alphz*.5*energy*((spiny+1) - (a->spin + 1));
    }

    return change;
}

void MonteCarlo3d::update(xyzs * a , vector<xyzs*> & protons ,int spiny, double dele){
    total_m += (spiny- a->spin);
    total_sm += (spiny- a->spin) * ((a->x + a->y + a->z )%2 ? -1:1 );
    total_e += dele;
    a->spin = spiny;
    if( spiny == -1){
        if(protons.size() != 0){
            std::swap(protons.back(),a);
            protons.pop_back();
        }
    }
    
    //    (spiny == -1 ? protons.remove(a) : protons.push_back(a));
}

void MonteCarlo3d::propose(vector<xyzs> & spins, vector<xyzs*> & protons, double T, int L,int sx,double alphz){
    int spiny,accepted = 0;
    int ipick, ipick2,ipick3,ipick4;
    double ipick5;
    uniform_int_distribution <int> dist(0,L*L*L-1);
    uniform_real_distribution<double> acceptance(0.,1.);
    
//    for (auto it = spins.begin(); it != spins.end(); ++it)
    for( int i =0; i< (L*L*L) ;i++)
    {
        ipick4 = dist(e1);
        ipick5 = acceptance(e1);
        spiny = - spins[ipick4].spin;
        
        double dele = change(&spins[ipick4],spins,protons,spiny,sx,L,alphz);
        
        if (ipick5 < exp(-dele/T)){
            update(&spins[ipick4],protons,spiny, dele);
            accepted++;
        }
    }
}

//    ipick = ipick4 %L;
//    ipick4 = (ipick4-ipick)/L ;
//    ipick2 = ipick4 %L;
//    ipick3 = (ipick4-ipick2)/L ;

//                if (choice == 1){
//                    spins[l].x = i;
//                    spins[l].y = j;
//                    spins[l].z = k;
//                    spins[l].spin = pow(-1,i+j+k);
//                    l++;
//                }
//                if (choice == 2){
//                    spins.x = i;
//                    spins.y = j;
//                    spins.z = k;
//                    spins.spin = -1;
//                }
//                if (choice == 3){
//                    if ( (i+j+k)%2 == 1){
//                        spins.x = i;
//                        spins.y = j;
//                        spins.z = k;
//                        spins.spin = 0;
//                    }
//                    else{
//                        spins.x = i;
//                        spins.y = j;
//                        spins.z = k;
//                        spins.spin = -1;
//                    }
//                }
//                if (choice == 4){
//                    if ( (i+j+k)%2 == 1){
//                        spins.x = i;
//                        spins.y = j;
//                        spins.z = k;
//                        spins.spin = 0;
//                    }
//                    else{
//                        spins.x = i;
//                        spins.y = j;
//                        spins.z = k;
//                        spins.spin = 1;
//                    }
//                }

#endif

#include <string>
#include "monte3d.hpp"
//#include "binner.cxx"
#include <iterator>
#include <fftw3.h>



#define WARM 100000
#define MCS  400000

int main(int argc, char* argv[])
{
    if(argc != 2) return -1;
    
    double t0 = atof(argv[1]);
    time_t seconds;
    time(&seconds);
    srand((unsigned int) seconds);  // seed random numbers based on cpu time
    
    int L = 8;                  // The box length
    int V = L*L*L;              // Volume of box
    
//    double z = 8;
    
    int spin_mag = 0; // choice: 1 for spin 1 ; 0 for spin 1/2
    
    //  Here begins code that reads in data from a file which contains the tabulated
    //  euclidian distance between two charged particles vs the ewald coulomb energy between
    //  ************
    
    string su = to_string(L); //L
    ifstream is("energy_distance"+su+"sc.txt");
    istream_iterator<double> start(is), end;
    vector<double> numbers(start, end);
    cout << "Read " << numbers.size() << " numbers" << std::endl;
    
    ifstream is1("graph_two"+su+"sc.txt");
    istream_iterator<double> start1(is1), end1;
    vector<double> numbers1(start1, end1);
    cout<< "Read" << numbers1.size() << "numbers" << endl;
    
    // print the numbers to stdout
    cout << "numbers read in:\n";
    //    copy(numbers.begin(), numbers.end(),
    //         ostream_iterator<double>(cout, " " ));
    int i =0;
    for (auto it = numbers.begin(); it != numbers.end(); ++it,i++)
    {
        if(i%4==0){
            cout << *it <<"\t" << *(it+1) <<"\t" << *(it+2) <<"\t" << *(it+3) <<endl;
        }
    }
    
    int ik =0;
    for (auto it = numbers1.begin(); it != numbers1.end(); ++it,ik++)
    {
        if(ik%3==0){
            cout << *it <<"\t" << *(it+1) <<"\t" << *(it+2) << endl;
        }
    }
    
    cout << endl;
    vector<xyzs> even ;
    
    unordered_map<int, double> mp, mp1;
    
    
    //    xyzs temp;
    for(unsigned int i =0; i < numbers.size(); i++){
        if(i%4==0){
//            vector<double> temp{numbers[i],numbers[i+1],numbers[i+2]} ;
            //            temp.x = numbers[i];
            //            temp.y = numbers[i+1];
            //            temp.z = numbers[i+2];
            //            temp.energ =numbers[i+3];
            //            even.push_back(temp);
//            mp[temp] = numbers[i+3];
            
            int temp = numbers[i]*L*L + numbers[i+1]*L + numbers[i+2];
            mp[temp] = numbers[i+3];
        }
    }
    
    for(unsigned int i =0; i < numbers1.size(); i++){
        if(i%3==0){
            //            vector<double> temp{numbers[i],numbers[i+1],numbers[i+2]} ;
            //            temp.x = numbers[i];
            //            temp.y = numbers[i+1];
            //            temp.z = numbers[i+2];
            //            temp.energ =numbers[i+3];
            //            even.push_back(temp);
            //            mp[temp] = numbers[i+3];
            
            int temp = numbers1[i]*numbers1[i+1];
            mp1[temp] = numbers1[i+2];
        }
    }

    
    
    for (auto it = even.begin(); it != even.end(); ++it)
    {
        cout << it->x << "\t" << it->y << "\t" << it->z  << "\t" << it -> energ <<endl;
        
    }
    
    int alph = 0;
    
  
    //double t0 = .01;
    double h0 = 0;
    ofstream datum,pos,zx,zxx,static_struct3d_L,spinspxx;
    
    double J = 1;
    
    string ss = to_string(h0); //h
    
    
    string st = to_string(0); //d
    string sr = to_string(alph); //alph
    string sq = to_string(J); //j
    string tt = to_string(t0);
    
    string typee;

    
    
    
    //string tr = to_string(spin_mag);
    datum.open("observables_h0="+ss+"d="+st+"j="+sq+"L="+su+"alph="+sr+".txt");
    cout << "h0 : "+ ss<< endl;

    pos.open("lattice_L="+su+"_H="+ss+"_alph="+sr+".pdb");
    static_struct3d_L.open("ferro_ising3dstatic_struct_L="+su+"t="+tt+"_H="+ss+"_alph="+sr+".pdb");
    spinspxx.open("ising3dstatus_L="+su+"t="+tt+"_H="+ss+"_alph="+sr+".pdb");
    

    zx.open("static_struct3d_L="+su+"t="+tt+".txt");
    zxx.open("pair_corr3d_L="+su+"t="+tt+".txt");
    
    

    
    #pragma omp parallel for
    for( unsigned int q = 0 ; q< 40; q++)
    {
        vector<xyzs>  spins(V);
        vector<xyzs*> protons;
        double T = t0 + q*(.02);
        double H = h0;//- q*.004;
        MonteCarlo3d one;              // call MonteCarlo3d class
        one.initialize(L,H,J,0,alph,spins,protons,0,mp);
        //  cout<< abs(one.total_mag(spins))/V;
        
        double big_mag = 0;
        double big2_mag = 0;
        double big4_mag = 0;
        double big_energy = 0;
        double big_stag_mag = 0;
        double mag = 0;
        double energ = 0;
        double mag2 = 0;
        double mag4 = 0;
        double energ2 = 0;
        double mag_energy = 0;
        double mag2_energy = 0;
        double mag4_energy = 0;
        double alpha = alph;


#pragma omp critical
        {
            cout << "The magnetization is " << one.observable_m()/V <<"\n";
        }
        for (int j = 0; j< WARM; j++){
            
            one.propose(spins,protons, T,L,spin_mag,alpha);
            cout<< "\r c done: " << j*100.0/(MCS-WARM) << "%             "<<flush;
        }
        double nobins = 20;
        double binsize = (MCS-WARM)/nobins;
        //  Binner<double> bin_m(0,binsize,nobins);
        vector<double> bin_m(nobins,0),bin_m_average(nobins), bin_e_average(nobins), bin_m2(nobins), bin_susc(nobins,0), bin_spheat(nobins,0),bin_e(nobins,0),bin_binder(nobins,0),bin_mag_energy(nobins,0),bin_mag4(nobins,0),bin_mag2(nobins,0);
;   // observables from binning
        //for(int i =0; i<nobins;i++){ cout<< bin_m[i]<<endl; }
        
        vector<double>  energy_list(MCS-WARM), magnet_list(MCS-WARM),magnet2_list(MCS-WARM),magnet4_list(MCS-WARM),susceptibility_list(MCS-WARM), specific_heat_list(MCS-WARM),binder_list(MCS-WARM),deriv_m_list(MCS-WARM),deriv_m2_list(MCS-WARM),mag_energy_list(MCS-WARM);
        
        for ( int j = 0 ; j< MCS-WARM; j++){
//            cout <<"size of map : " <<mp.size() <<endl;
            
            one.propose(spins, protons, T,L,spin_mag,alpha);
            magnet_list[j] = abs(one.observable_m())/V;
            energy_list[j] = one.observable_e()/(2*V);
            bin_m[j%20] += magnet_list[j];
            //        cout<< "bin_m[j%20]" << bin_m[j%20] << endl;
            bin_e[j%20] += energy_list[j];
            magnet2_list[j] = magnet_list[j]*magnet_list[j];
            bin_mag2[j%20] += magnet2_list[j];
            mag2 += magnet2_list[j];
            magnet4_list[j] = magnet2_list[j]*magnet2_list[j];
            bin_mag4[j%20] += magnet4_list[j];
            mag4 += magnet4_list[j];
        

            energ2 += energy_list[j]*energy_list[j];
            
            
            mag_energy_list[j] = (magnet_list[j]*energy_list[j]);
            bin_mag_energy[j%20] += mag_energy_list[j];
            
            
            cout<< "\r done: " << j*100.0/(MCS-WARM) << "%             "<<flush;
            cout << endl;
            //
//#pragma omp critical
//            {
//            if(j == MCS-WARM-1){
//                double d;
//                int cnt=1;
//            for (auto it = spins.begin(); it != spins.end(); ++it)
//            {
//                if(it->spin  == 1){ typee="O";}
//                else{typee="N";}
//                pos <<"ATOM"<<setw(7)<<cnt<<setw(3)<<typee<< setw(4)<< "   MET X   1"<<setw(12) << it->x<<setw(8)<<it->y<<setw(8)<<it->z<<setw(2)<<""<<fixed<<setprecision(2)<<d<<endl;
//                cnt++;
//            }
//            }
//            }

        }
        
//
//        cout << "energ_mag = "<< deriv_m2_list[1]/(MCS-WARM) <<endl;

        
        fftw_complex *in, *out,*pair_corr, *stat_struct,spinz;
        
        //    spinz[0] = 1;
        //    spinz[1] = 1;
        fftw_plan p,x;
        
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L*L*L);
        out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L*L*L);
        stat_struct = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *L*L*L);
        pair_corr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *L*L*L);
        
        in[0][0] = 1;
        p = fftw_plan_dft_3d(L,L,L,in,out,-1,FFTW_ESTIMATE);
        
        int j = 0 ;
        for(unsigned int ix=0; ix<L; ix++){
            for(unsigned int iy=0; iy<L; iy++){
                for(unsigned int iz=0; iz<L; iz++){
                
//                    pick1 = acceptance(e1);
                   
                    in[iz+L*iy+L*L*ix][0] = spins[j].spin;
                    
                    in[iz+L*iy+L*L*ix][1] = 0;
                    j++;
                    }
                }
            }
        
        x = fftw_plan_dft_3d(L,L,L,stat_struct,pair_corr,1,FFTW_ESTIMATE);
        
        
        fftw_execute(p);
        for(unsigned int ix=0; ix<L; ix++){
            cout << "ix = " << ix << " out[ix] "<< out[ix][0] <<"  :  "<<out[ix][1] <<endl;
        }
        
        for(unsigned int ix=0; ix<L; ix++){
            for(unsigned int iy=0; iy<L; iy++){
                for(unsigned int iz=0; iz<L; iz++){
                    
                    
                    stat_struct[iz+L*iy+L*L*ix][0]= out[iz+L*iy+L*L*ix][0]*out[iz+L*iy+L*L*ix][0]+out[iz+L*iy+L*L*ix][1]*out[iz+L*iy+L*L*ix][1];
                    stat_struct[iz+L*iy+L*L*ix][1]= 0;
                    
                }
            }
        }
        
//        int cnt=1;

        
//        string typee;
//        double d;
//        int cnt=1;
//        for(unsigned int ix=0; ix<L; ix++){
//            for(unsigned int iy=0; iy<L; iy++){
//                for(unsigned int iz=0; iz<L; iz++){
//                    if( in[ix*L*L+L*iy+iz][0]== 1){ typee="O";}
//                    else if(in[ix*L*L+L*iy+iz][0] == -1){typee="N";}
//                    else{typee="C";}
//                    d = 1;
//                    spinspxx <<"ATOM"<<setw(7)<<cnt<<setw(3)<<typee<< setw(4)<< "   MET X   1"<<setw(12) << ix<<setw(8)<<iy<<setw(8)<<iz<<setw(2)<<""<<fixed<<setprecision(2)<<d<<endl;
//                    cnt++;
//                }
//            }
//        }
        
//        cnt =1;
        
//        for(unsigned int ix=0; ix<L; ix++){
//            for(unsigned int iy=0; iy<L; iy++){
//                for(unsigned int iz=0; iz<L; iz++){
////                    if( stat_struct[ix*L*L+L*iy+iz][0]/(L*L*L*L*L*L) == 1){ typee="O";}
////                    else if( stat_struct[ix*L*L+L*iy+iz][0]/(L*L*L*L*L*L) == 0){typee="N";}
////                    else{typee="C";}
//                    
//                    typee = "C";
//                    double d = stat_struct[ix*L*L+L*iy+iz][0]/(L*L*L);
//                    
//                    int tix = -ix ,tiy = -iy, tiz = -iz;
//                    //
//                    if(ix > L/2){
//                        tix = L - ix;
//                    }
//                    
//                    if(iy > L/2){
//                        tiy = L - iy;
//                    }
//                    
//                    if(iz > L/2){
//                        tiz = L - iz;
//                    }
//                    
//                    
//                    static_struct3d_L <<"ATOM"<<setw(7)<<cnt<<setw(3)<<typee<< setw(4)<< "   MET X   1"<<setw(12) << tix <<setw(8)<< tiy <<setw(8)<< tiz <<setw(2)<<""<<fixed<<setprecision(4)<<d<<endl;
//                    cnt++;
//                }
//            }
//        }
        
        fftw_execute(x);
        fftw_destroy_plan(p);
        fftw_destroy_plan(x);
        
        
        double z = 8;
        
//        dcomp static_struct;
        
        
        zx.close();
        zxx.close();
        
        
        fftw_free(in), fftw_free(out), fftw_free(pair_corr);
        
        
        cout<<endl;
        
        for(int i = 0 ; i<20 ; i++){
            bin_m[i] = bin_m[i]/binsize;
            bin_e[i] = bin_e[i]/binsize;
            bin_mag_energy[i] = bin_mag_energy[i]/binsize;
            
            bin_mag2[i] = bin_mag2[i]/binsize;
            bin_mag4[i] = bin_mag4[i]/binsize;
        }
        
        vector<double> bin_deriv_m(nobins,0),bin_deriv_m2(nobins,0);
        
        double bin_m_mean = 0;
        double bin_e_mean = 0;
        double bin_m2_mean= 0;
        double bin_susc_mean = 0;
        double bin_spheat_mean = 0;
        double bin_binder_mean = 0;
        double bin_mag_energy_mean = 0;
        double bin_deriv_m_mean = 0;
        double bin_deriv_m2_mean = 0;
        double bcumulant, bcumulant4 ,deriv_m,deriv_m2,deriv_binder2, susceptibility,susceptibility2,specific_heat ,m_error=0, e_error=0, susc_error=0,sphe_error=0,binder_error=0,deriv_m_error =0, deriv_m2_error=0;
        
        
        for(int i = 0 ; i<20; i++){
            bin_m_mean += bin_m[i]/nobins;
            bin_e_mean += bin_e[i]/nobins;
            bin_mag_energy_mean += bin_mag_energy[i]/nobins;

        }
        
        
        for( unsigned int j = 0; j < (MCS-WARM);  j++){
            susceptibility_list[j] = (magnet_list[j] - bin_m_mean)*(magnet_list[j] - bin_m_mean);
            bin_susc[j%20] += susceptibility_list[j];
            
            specific_heat_list[j] = (energy_list[j] - bin_e_mean)*(energy_list[j] - bin_e_mean);
            bin_spheat[j%20] += specific_heat_list[j];
            
            
//            binder_list[j] = 1 - magnet4_list[j]/(3*magnet2_list[j]*magnet2_list[j]);
//            bin_binder[j%20] += binder_list[j];
            
//            deriv_m_list[j] = mag_energy_list[j]/magnet_list[j] - energy_list[j];
//            deriv_m2_list[j] = mag_energy_list[j] - (energy_list[j]*magnet_list[j]);
//            
//            
//            bin_deriv_m[j%20] += deriv_m_list[j];
//            bin_deriv_m2[i] += deriv_m2_list[j];
            
            
        }
        
        for(int i = 0 ; i<20 ; i++){
            bin_susc[i] = bin_susc[i]/binsize;
            bin_spheat[i] = bin_spheat[i]/binsize;
            
            bin_binder[i] = 1- bin_mag4[i]/(3*bin_mag2[i]*bin_mag2[i]);
            
//            bin_binder[i] = bin_binder[i]/binsize;
            
            bin_deriv_m[i] = (bin_mag_energy[i]/bin_m[i]) - bin_e[i];
            bin_deriv_m2[i] = (bin_mag_energy[i]) - (bin_e[i]*bin_m[i]);
            
            
            
            
//            bin_deriv_m[i] += bin_deriv_m[i]/binsize;
//            bin_deriv_m2[i] += bin_deriv_m2[i]/binsize;
        }
        
        for(int i = 0 ; i<20; i++){
            bin_susc_mean += bin_susc[i]/nobins;
            bin_spheat_mean += bin_spheat[i]/nobins;
            bin_binder_mean += bin_binder[i]/nobins;
            
           
        }
        

        for(int i = 0 ; i<20; i++){
            bin_susc[i] = bin_susc[i]/T;
            bin_spheat[i] = bin_spheat[i]/(T*T);

        }
        
        
        for(int i = 0 ; i<20; i++){
        
            bin_deriv_m_mean += bin_deriv_m[i]/(nobins);
            bin_deriv_m2_mean += bin_deriv_m2[i]/(nobins);
        }
        
        
        bin_susc_mean = bin_susc_mean/(T);
        bin_spheat_mean = bin_spheat_mean/(T*T);
        
        energ2 = energ2/(MCS-WARM);
        mag2 = mag2/(MCS-WARM);
        mag4 = mag4/(MCS-WARM);
        mag_energy = mag_energy/(MCS-WARM);
        mag2_energy = mag2_energy/(MCS-WARM);
        
        mag4_energy = mag4_energy/(MCS-WARM);
        bcumulant = 1 - mag4/(3*mag2*mag2);
        bcumulant4= 1 - mag2/(3*bin_m_mean*bin_m_mean);
        
        susceptibility = 0;
        e_error = 0;
        susceptibility2= 0;
        
        deriv_m = (bin_mag_energy_mean)/(bin_m_mean) - (bin_e_mean);
        deriv_m2 = (bin_mag_energy_mean) - (bin_e_mean*bin_m_mean);
        deriv_binder2 = V*(1-bcumulant)*((bin_e_mean)- 2*(mag_energy)/bin_m_mean - mag2_energy/mag2);
        
        
        for( unsigned int j = 0; j < 20;  j++){
            m_error = pow((bin_m[j] - bin_m_mean),2);
        }
        
        for( unsigned int j = 0; j < 20;  j++){
            e_error += pow((bin_e[j] - bin_e_mean),2);
        }
        for( unsigned int j = 0; j < 20;  j++){
            susc_error += pow((bin_susc[j] - bin_susc_mean),2);
        }
        for( unsigned int j = 0; j < 20;  j++){
            sphe_error += pow((bin_spheat[j] - bin_spheat_mean),2);
        }
        for( unsigned int j = 0; j < 20;  j++){
            binder_error += pow((bin_binder[j] - bin_binder_mean),2);
        }
        for( unsigned int j = 0; j < 20;  j++){
            deriv_m_error += pow((bin_deriv_m[j] - bin_deriv_m_mean),2);
        }
        for( unsigned int j = 0; j < 20;  j++){
            deriv_m2_error += pow((bin_deriv_m2[j] - bin_deriv_m2_mean),2);
        }
        
        
        susceptibility2 =  mag2 - (bin_m_mean*bin_m_mean);
        specific_heat   =  energ2 - (bin_e_mean*bin_e_mean);
        
        
        m_error = sqrt(m_error)/20;
        e_error = sqrt(e_error)/20;
        susc_error = sqrt(susc_error)/20;
        sphe_error = sqrt(sphe_error)/20;
        susceptibility2 = susceptibility2/(T);
//        susceptibility = V*bin_susc_mean/(T)
        specific_heat = specific_heat/(T*T);
        
//        susc_error = sqrt()
        #pragma omp critical
        {
//            cout << "T = " <<T << "total_energy"<< one.total_energy(spins,,0,L)<<endl;
            datum << T << "\t"<< bin_m_mean << "\t" << m_error<< "\t"  << bin_susc_mean  << "\t"<< susc_error << "\t"<< bin_e_mean << "\t" << e_error<< "\t" << bin_spheat_mean  << "\t"<< sphe_error<< "\t" << bcumulant << "\t" << bin_binder_mean<< "\t" <<binder_error << "\t"<< bin_deriv_m_mean<< "\t"<<deriv_m_error << "\t" <<deriv_m2  << "\t" << deriv_m2_error<< "\t" <<deriv_binder2 <<endl;
        }
        
    }
    datum.close();
    pos.close();
    static_struct3d_L.close();
    spinspxx.close();

    
}

// Author: Yang Bo
// Date: 2014-2-5
//


#include <stdio.h>
#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include "expectvalue.h"
#include <time.h>
#include <sstream>
//#include "lanczos.h"
//#include <boost/lexical_cast.hpp>
//#include <boost/format.hpp>
#include <chrono>
#include "Eigen/Eigenvalues"



using Eigen::MatrixXd;
using Eigen::Matrix4f; 
using namespace std;



typedef unsigned long long INT;


int main() {

    cout<<"Code version Mar 9, 2020"<<endl;
    cout<<"This routine diagonalise within a selected hilbert space for selected two-body or three-body pseudopotentials on the sphere or disk"<<endl;
    cout<<"For sphere choose 1, for disk choose 2"<<endl;
    int g_opt;cin>>g_opt;
    cout<<"Is it (1) a two-body interaction or (2) a one-body background potential or (3) three-body interaction or (4) sum_i (L_i^dagger)^2 (L_i)^2 or (5) overlap or (6) L^dagger L just for the diagonals"<<endl;
    int hopt;cin>>hopt;
    time_t rawtime;
    struct tm * timeinfo;

    cout<<"give the number of orbitals"<<endl;
    cin>>n_orb;
    cout<<"give the number of electrons"<<endl;
    cin>>n_el;
    cout<<"give the dimension"<<endl; 
    int dim;cin>>dim;
    string nameroot;
    cout<<"give the wavefunction filename root (note they should all have the same ordered basis)\n (The file names are of the form <name root>_<index>)"<<endl; cin >> nameroot;

 /*   for (int i=0;i<dim;i++) {

		fflush(stdout); 
        string h_file;cin>>h_file;
        vector<INT> ev;vector<double> coe;
        readwf(ev,coe,h_file,n_el);
        if (evec.size()==0) evec.push_back(ev);
        e_coeff.push_back(coe);
        printf("\rProgress %d out of %d %d %d %d %d          ", i,dim,evec.size(),e_coeff.size(),ev.size(),coe.size());
        ev.clear();coe.clear();


    }*/

    for (int i=0;i<dim;i++) {
            fflush(stdout); 
            string h_file;//cin>>h_file;
            vector<INT> ev;vector<double> coe;
            h_file=nameroot+"_"+to_string(i); 
            cout<<"filename "<<h_file<<endl;
            readwf(ev,coe,h_file,n_orb);
            if (evec.size()==0) evec.push_back(ev);
            e_coeff.push_back(coe);
            printf("\rProgress %d out of %d %d %d %d %d          ",dim,evec.size(),e_coeff.size(),ev.size(),coe.size());
            ev.clear();coe.clear();
    }

    vector<vector<vector<int>>> evec_bin;
    vector<vector<vector<int>>> evec_occ;

    bool ccheck=true;

    for (int i=0;i<evec.size();i++) {

        vector<vector<int>> ev_bin,ev_occ;
        for (int j=0;j<evec[i].size();j++) {
            vector<int> vec,occ;
            if (hopt!=5) dec_bin(evec[i][j],vec);
            if (hopt==5&&j==0) dec_bin(evec[i][j],vec);

            if (ccheck) {
                ccheck=false;
                int e_num=0; double e_mom=0;
                double s=(vec.size()-1)/2.0;
                for (int k=0;k<vec.size();k++) {
                    if (vec[k]==1) {
                        e_num++;e_mom+=k-s;
                    }
                }
                cout<<endl;
                cout<<"Total number of electrons is "<<e_num<<"; L_z sector is "<<e_mom<<endl;
            }
            if (hopt!=5) ev_bin.push_back(vec);
           // bin_occ(vec,occ);
         //   ev_occ.push_back(occ);
            vec.clear();occ.clear();
    
        }

        
        if (hopt==6&&i==0) {
            evec_bin.push_back(ev_bin);
        } else if (hopt!=5) {
            evec_bin.push_back(ev_bin);
        }
       // evec_occ.push_back(ev_occ);
        ev_bin.clear();ev_occ.clear();
    }

//===================================================================================================
    vector<vector<double>> vmn;
    vector<double> vm1m2m3m4;


    if (hopt==1) {
	    for (int i=0;i<n_orb;i++) {
		    SPHERE_PSEUD.push_back(0);
	    }
        cout<<"Give the interaction in the form of pseudoptentials V_a (only isotropic ones are considered)"<<endl;
        vector<int> a;vector<double> coeff;
        read_pp(a,coeff);
	    for (int ii=0;ii<a.size();ii++) {
		// the loop is over all the pseudopotentials
		    mm=nn=a[ii];
		    SPHERE_PSEUD[a[ii]]=coeff[ii];
        }
	    if (g_opt==1) two_body_sphere(vm1m2m3m4);
        if (g_opt==2) two_body_disk(vm1m2m3m4,a,coeff);

    }

    if (hopt==2) {
        cout<<"The background potential is randomly generated, as an n_orb by n_orb matrix. The matrix elements lie between -1 and 1"<<endl;
        srand (time(NULL));
        vector<double> temp;
        for (int j=0;j<n_orb;j++) {
            for (int i=0;i<n_orb;i++) {
                //temp.push_back((rand() % 10000000)/10000000.0);
                temp.push_back(1);
            }
            vmn.push_back(temp);
            temp.clear();
        }

        /*for (int i=0;i<n_orb;i++) {
            for (int j=0;j<n_orb;j++) {
                cout<<vmn[i][j]<<endl;
            }
        }*/
    }

    if (hopt==3) {
        //threebodyH();
    }

    if (hopt==4) {
	   // if (g_opt==1) b4_sphere(vm1m2m3m4);
       // if (g_opt==2) b4_disk(vm1m2m3m4,a,coeff);
    }

    int minput;



//=========================================================

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  printf ( "Current local time and date: %s", asctime (timeinfo) );
  cout<<endl;

//=================================================================



    vector<int> row,col;
    vector<double> h;
    double ltt=(n_orb-1.0)/2.0;

    if (hopt==6) {
        for (int i=0;i<e_coeff.size();i++) {
            lminus(evec[0],evec_bin[0],e_coeff[i],i);
        }
        return 0;
    }

    int check=0;
    if (hopt==1) {
        hamiltonian_construction(evec[0],vm1m2m3m4);
        cout<<"Basis matrix constructed"<<endl;
    }

//=========================================================

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  printf ( "Current local time and date: %s", asctime (timeinfo) );
  cout<<endl;

//=================================================================

    for (int i=0;i<e_coeff.size();i++) {
        for (int j=i;j<e_coeff.size();j++) {

            printf("\rMatrix Element %d %d          ", i,j);
		    fflush(stdout);

            double sum=0;
            if (hopt==5) {
                for (int ii=0;ii<evec[0].size();ii++) {
                    sum+=e_coeff[i][ii]*e_coeff[j][ii];     
                }
                row.push_back(i);col.push_back(j);h.push_back(sum);
                continue;
            } 

            if (hopt==1) {
                row.push_back(i);col.push_back(j);
                if (check==0) {
                    h.push_back(matvec1(e_coeff[j],e_coeff[i]));check=1;
                } else {
                    h.push_back(matvec(e_coeff[j],e_coeff[i]));
                }
                continue;         
            }           



            for (int ii=0;ii<evec[0].size();ii++) {

                if (fabs(e_coeff[i][ii])<1e-15) continue;

                for (int jj=0;jj<evec[0].size();jj++) {
                    if (fabs(e_coeff[j][jj])<1e-15) continue;
                    double value=0;
                    if (hopt==1) {
                    } else {
                        if (hopt==2) value=onebody(evec_bin[0],evec_bin[0],ii,jj,vmn);
                    }
                    if (hopt==4&&evec[0][jj]==evec[0][ii]) {
                        for (int k=0;k<evec_bin[0].size();k++) {
                            if (evec_bin[0][jj][k]==1) {
                                double lm=ltt-k;
                                value+=(ltt+lm)*(ltt-lm+1)*(ltt+lm-1)*(ltt-lm+2);
                            }
                        }
                    }
                    if (value==0) continue;
                    sum+=value*e_coeff[i][ii]*e_coeff[j][jj];
                }
            }

            row.push_back(i);col.push_back(j);h.push_back(sum);            

        }
    }

//=========================================================

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  printf ( "Current local time and date: %s", asctime (timeinfo) );
  cout<<endl;

//=================================================================

    cout<<"display matrix?"<<endl;
    int dmatrix;cin>>dmatrix;
    if (dmatrix==1) {
        for (int i=0;i<row.size();i++) {
           if (fabs(h[i])>1e-12&&row[i]!=col[i]) cout<<row[i]<<" "<<col[i]<<" "<<h[i]<<endl;
           if (fabs(h[i]-1)>1e-12&&row[i]==col[i]) cout<<row[i]<<" "<<col[i]<<" "<<h[i]<<endl;

        }
    }




//=========================================================

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  printf ( "Current local time and date: %s", asctime (timeinfo) );
  cout<<endl;

//=================================================================

    if (h.size()==1) {

        cout<<"The energy is "<<h[0]<<endl;
        vector<double> eigen;
        eigen.push_back(h[0]);
        write(eigen,"eigen");


    } else {

	    int hsize=row.size();

	    for (int i=0;i<hsize;i++) {
		    if (row[i]!=col[i]) {
			    row.push_back(col[i]);
			    col.push_back(row[i]);
			    h.push_back(h[i]);
		    }
	    }
	    sort3(row,h,col);

        if (hopt==2) {

            double shift;

            cout<<"The effective matrix is as follows:"<<endl;
            for (int i=0;i<row.size();i++) {
                cout<<row[i]<<" "<<col[i]<<" "<<h[i]<<endl;
                if (row[i]==0&&col[i]==0) shift=h[i];
            }

            for (int i=0;i<row.size();i++) {
                if (row[i]==col[i]) h[i]-=shift;
            }
        }

        MatrixXd m(dim,dim);
        for (int i=0;i<dim;i++) {
		    for (int j=0;j<dim;j++) {
			    m(i,j)=0;
		    }
        }

        for (int i=0;i<row.size();i++) {
		    m(row[i],col[i])=h[i];
   	    }

        Eigen::SelfAdjointEigenSolver<MatrixXd> es;
        es.compute(m);

	    h.clear();row.clear();col.clear();

        cout<<"energy spectrum?"<<endl;
        int cc;cin>>cc;
        vector<double> eigen;
        for (int i=0;i<dim;i++) eigen.push_back(es.eigenvalues().col(0)[i]);
        if (cc==1) {
            write(eigen,"eigen");
        }

        /*for (int i=0;i<dim;i++) {
		    cout<<es.eigenvalues().col(0)[i]<<endl;
        }*/

        cout<<"give output name"<<endl;
        string output;cin>>output;

        for (int k=0;k<dim;k++) {
            vector<double> groundstate;
            for (int j=0;j<dim;j++) {
			    groundstate.push_back(es.eigenvectors().col(k)[j]);
    	    }

            if (hopt==5||hopt==1) {
                vector<INT> g_evec=evec[0];

                vector<double> g_coeff;
                for (int i=0;i<evec[0].size();i++) {
                    double sss=0;
                    for (int j=0;j<dim;j++) {
                        sss+=e_coeff[j][i]*groundstate[j];
                    }
                    g_coeff.push_back(sss);
                    sss=0;
                }
                double norm=0;
                for (int i=0;i<g_coeff.size();i++) norm+=g_coeff[i]*g_coeff[i];
                norm=sqrt(norm);
                for (int i=0;i<g_coeff.size();i++) g_coeff[i]=g_coeff[i]/norm;
                writewf(g_evec,g_coeff,output+to_string(k));
                g_evec.clear();g_coeff.clear();groundstate.clear();
                continue;
            }
 /*           vector<vector<double>> ee=e_coeff;

            for (int i=0;i<dim;i++) {

                for (int j=0;j<evec[0].size();j++) {
                    ee[i][j]=e_coeff[i][j]*groundstate[i];
                }            

            }

            vector<INT> g_evec;
            vector<double> g_coeff;

            for (int i=0;i<dim;i++) {
                for (int j=0;j<evec[0].size();j++) {
                    g_evec.push_back(evec[0][j]);
                    g_coeff.push_back(ee[i][j]);
                }
            }

            sort2b(g_evec,g_coeff);


            vector<INT> gg_evec;
            vector<double> gg_coeff;

            double sum=0;
            for (int i=0;i<g_evec.size()-1;i++) {

                sum+=g_coeff[i];
                if (g_evec[i+1]!=g_evec[i]) {
                    gg_evec.push_back(g_evec[i]);
                    gg_coeff.push_back(sum);
                    sum=0;
                }
            }

            sum+=g_coeff[g_coeff.size()-1];
            gg_evec.push_back(g_evec[g_evec.size()-1]);
            gg_coeff.push_back(sum);
            writewf(gg_evec,gg_coeff,output+to_string(k));
            gg_evec.clear();g_evec.clear();
            gg_coeff.clear();g_coeff.clear();
*/

        }

    }

}










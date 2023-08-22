/*
 *  cinc.h
 *  
 *
 *  Created by Bo Yang on 2/5/2014.
 *
 */

#include <vector>
#include "math.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <fstream>
#include <cstdio>
#include <random>

#define _USE_MATH_DEFINES

using namespace std;

typedef unsigned long long  INT;


int n_orb,n_el;
vector<vector<INT>> evec;
vector<vector<double>> e_coeff;
int mm,nn;
vector<vector<int>> binstate_vec,record;
vector<int> orbital;
vector<double> opotential;
vector<double> SPHERE_PSEUD;
vector<vector<int>> m_col,m_row;
vector<double> mdiag;
vector<vector<double>> moffdiag;

double* ccco3B;


int combi(int a,int b) {
    int c=1;
    for (int i=0;i<b;i++) {
        c=c*(a-i)/(i+1);
    }    
    return c;
}

void sortb(int n, int *ra,int *rb) {
    
    double l;
    int rra,rrb,ir,i,j;
    l=n/2+1;
    ir=n;
    
check1:
    
    if (l>1) {
        
        l--;
        rra=ra[int(l)-1];
        rrb=rb[int(l)-1];
        
    }
    else {
        
        rra=ra[int(ir)-1];
        rrb=rb[ir-1];
        ra[ir-1]=ra[0];
        rb[ir-1]=rb[0];
        ir--;
        
        if (ir==1) {
            
            ra[0]=rra;
            rb[0]=rrb;
            
            return;
            
        }
    }
    i=l;
    j=l+l;
    
check2:
    
    if (j<=ir) {
        
        if (j<ir) {
            
            if (ra[j-1]<ra[j]) j++;
            
        }
        
        if (rra<ra[j-1]) {
            
            ra[i-1]=ra[j-1];
            rb[i-1]=rb[j-1];
            
            i=j;
            j=j+j;
            
        }
        else {
            
            j=ir+1;
            
        }
        goto check2;
        
    }
    
    ra[i-1]=rra;
    rb[i-1]=rrb;
    
    goto check1;
    
    
}

void sortbb(int n, INT *ra,int *rb) {
    
    double l;
    int rrb,ir,i,j;
	INT rra;
    l=n/2+1;
    ir=n;
    
check1:
    
    if (l>1) {
        
        l--;
        rra=ra[int(l)-1];
        rrb=rb[int(l)-1];
        
    }
    else {
        
        rra=ra[int(ir)-1];
        rrb=rb[ir-1];
        ra[ir-1]=ra[0];
        rb[ir-1]=rb[0];
        ir--;
        
        if (ir==1) {
            
            ra[0]=rra;
            rb[0]=rrb;
            
            return;
            
        }
    }
    i=l;
    j=l+l;
    
check2:
    
    if (j<=ir) {
        
        if (j<ir) {
            
            if (ra[j-1]<ra[j]) j++;
            
        }
        
        if (rra<ra[j-1]) {
            
            ra[i-1]=ra[j-1];
            rb[i-1]=rb[j-1];
            
            i=j;
            j=j+j;
            
        }
        else {
            
            j=ir+1;
            
        }
        goto check2;
        
    }
    
    ra[i-1]=rra;
    rb[i-1]=rrb;
    
    goto check1;
    
    
}

void sort2(vector<int> &a,vector<int> &b) {
    
    int len=a.size();
    
    
    int *aa;
    
    aa= new int [len];
    
    for (int i=0;i<len;i++) {
        
        aa[i]=a[i];
        
    }
    
    int *iseq= new int [len];
    
    for (int i=0;i<len;i++) {
        
        iseq[i]=i;
        
    }
    
    sortb(len,aa,iseq);
    
    vector<int> bb;
    
    for (int i=0;i<len;i++) {
        
        a[i]=aa[i];
        bb.push_back(b[iseq[i]]);
        
    }	
    
    delete [] aa;
    
    for (int i=0;i<len;i++) {
        b[i]=bb[i];
    }
    
    bb.clear();
    
}

void sort2(int n, double *ra,int *rb) {
	
	double l,rra;
	int rrb,ir,i,j;
	l=n/2+1;
	ir=n;
	
check1:
	
	if (l>1) {
		
		l--;
		rra=ra[int(l)-1];
		rrb=rb[int(l)-1];
		
	}
	else {
		
		rra=ra[int(ir)-1];
		rrb=rb[ir-1];
		ra[ir-1]=ra[0];
		rb[ir-1]=rb[0];
		ir--;
		
		if (ir==1) {
			
			ra[0]=rra;
			rb[0]=rrb;
			
			return;
			
		}
	}
	i=l;
	j=l+l;
	
check2:
	
	if (j<=ir) {
		
		if (j<ir) {
			
			if (ra[j-1]<ra[j]) j++;
			
		}
		
		if (rra<ra[j-1]) {
			
			ra[i-1]=ra[j-1];
			rb[i-1]=rb[j-1];
			
			i=j;
			j=j+j;
			
		}
		else {
			
			j=ir+1;
			
		}
		goto check2;
		
	}
	
	ra[i-1]=rra;
	rb[i-1]=rrb;
	
	goto check1;
	
	
}

void sort2a(vector<int> &a,vector<double> &b) {
    
    int len=a.size();
    
    
    int *aa;
    
    aa= new int [len];
    
    for (int i=0;i<len;i++) {
        
        aa[i]=a[i];
        
    }
    
    int *iseq= new int [len];
    
    for (int i=0;i<len;i++) {
        
        iseq[i]=i;
        
    }
    
    sortb(len,aa,iseq);
    
    vector<double> bb;
    
    for (int i=0;i<len;i++) {
        
        a[i]=aa[i];
        bb.push_back(b[iseq[i]]);
        
    }	
    
    delete [] aa;
    
    for (int i=0;i<len;i++) {
        b[i]=bb[i];
    }
    
    bb.clear();
    
}

void sort2b(vector<INT> &a,vector<double> &b) {
    
    int len=a.size();
    
    
    INT *aa;
    
    aa= new INT [len];
    
    for (int i=0;i<len;i++) {
        
        aa[i]=a[i];
        
    }
    
    int *iseq= new int [len];
    
    for (int i=0;i<len;i++) {
        
        iseq[i]=i;
        
    }
    
    sortbb(len,aa,iseq);
    
    vector<double> bb;
    
    for (int i=0;i<len;i++) {
        
        a[i]=aa[i];
        bb.push_back(b[iseq[i]]);
        
    }	
    
    delete [] aa;
    
    for (int i=0;i<len;i++) {
        b[i]=bb[i];
    }
    
    bb.clear();
    
}

void sort3(vector<int> &a,vector<double> &b, vector<int> &c) {
    
    int len=a.size();
    
    
    int *aa;
    
    aa= new int [len];
    
    for (int i=0;i<len;i++) {
        
        aa[i]=a[i];
        
    }
    
    int *iseq= new int [len];
    
    for (int i=0;i<len;i++) {
        
        iseq[i]=i;
        
    }
    
    sortb(len,aa,iseq);
    
    vector<double> bb,cc;
    
    for (int i=0;i<len;i++) {
        
        a[i]=aa[i];
        bb.push_back(b[iseq[i]]);
        cc.push_back(c[iseq[i]]);
        
    }	
    
    delete [] aa;
    
    for (int i=0;i<len;i++) {
        b[i]=bb[i];c[i]=cc[i];
    }
    
    bb.clear();cc.clear();
    
}

void sort3(double *vec1,double *vec2,double *vec3,int len) {
	
	if (len==1) return;
	
	int *iseq= new int [len];
	double *temp1, *temp2;
	
	temp1= new double [len];
	temp2= new double [len];
	
	
	for (int i=1;i<=len;i++) {
		
		iseq[i-1]=i;
		
	}
	
	
	sort2(len,vec1,iseq);
	
	
	
	
	for (int i=1;i<=len;i++) {
		
		
		temp1[i-1]=vec2[iseq[i-1]-1];
		temp2[i-1]=vec3[iseq[i-1]-1];
		
	}
	
	
	delete [] iseq;
	
	for (int i=0;i<len;i++) {
		
		vec2[i]=temp1[i];
		vec3[i]=temp2[i];
		
	}
	
	delete [] temp1,temp2;
}

void sort3long(vector<INT> &a,vector<double> &b, vector<INT> &c) {
    
    int len=a.size();
    
    
    int *aa;
    
    aa= new int [len];
    
    for (int i=0;i<len;i++) {
        
        aa[i]=a[i];
        
    }
    
    int *iseq= new int [len];
    
    for (int i=0;i<len;i++) {
        
        iseq[i]=i;
        
    }
    
    sortb(len,aa,iseq);
    
    vector<double> bb,cc;
    
    for (int i=0;i<len;i++) {
        
        a[i]=aa[i];
        bb.push_back(b[iseq[i]]);
        cc.push_back(c[iseq[i]]);
        
    }	
    
    delete [] aa;
    
    for (int i=0;i<len;i++) {
        b[i]=bb[i];c[i]=cc[i];
    }
    
    bb.clear();cc.clear();
    
}





void dec_bin(INT dec,vector<int>& bin) {

	bin.clear();
    for (int i=0;i<n_orb;i++) {
        if (dec>0) {

		bin.push_back(dec%2);
        	dec=(dec-dec%2)/2;
	} else {
		bin.push_back(0);
	}
    }
}

int bin_dec(vector<int>& bin) {
    
    int dec=0;
    for (int i=0;i<bin.size();i++) {
        if (bin[i]==1) dec+=pow(2,i);
    }
    
    return dec;
}

void bin_occ(vector<int>& bin, vector<int>& occ) {
    
    int dec=0;
    for (int i=0;i<bin.size();i++) {
        if (bin[i]==1) occ.push_back(i);
    }
}

void readwf(vector<INT> &evec,vector<double> &coeff,string filename,int len) {
	
	
	int vec[len];
	double pnorm;
	
	FILE *file;
	
	char *fn;
	
	fn = new char [filename.size()+1];
	strcpy (fn,filename.c_str());

	file = fopen (fn,"r");

	delete [] fn;
	
	string suffix,svec;
	char *csvec;
	
	
	
	
	int s_len;
	fscanf(file,"%d",&s_len);

	while (feof(file) ==0) {
		
		INT index;
		fscanf(file, "%llu\n", &index);		

		
		evec.push_back(index);
		
		double coefficient;
		
		fscanf(file, "%lf\n",&coefficient);
		//cout.precision(20);
		//cout<<coefficient<<endl;
		
		coeff.push_back(coefficient);
	}
	
	assert(s_len==coeff.size());
	fclose(file);
	
}



void read_pp(vector<int>& a, vector<double>& coeff) {

    bool check=true;

    do {
	int aa;double cc;
	cin>>aa;
	if (aa!=0) {
		cin>>cc;
	}

	if (aa!=0) {
		a.push_back(aa);coeff.push_back(cc);
	} else {
		check=false;
	}
    } while (check);
}


double binom(int n, int m) // compute n choose m
{
  assert(n >= 0);
  assert(m >= 0);
  if(n >= m){
    double f = 1.0;
    for (int i = m; i >= 1; i--) {
      f *= (n + i - m) / static_cast<double>(i);
    }
    return f;
  }
  else
    return 0.0;
}

int abs ( const int x)
{
        return x < 0 ? -x : x ;
}


double cgcoef(int jj, int n1, int n2)
{
//-------------------------------------------------------------
// THIS ROUTINE CALCULATES THE CLEBSCH-GORDON
//  COEFFICIENT <L,M1,L,M2;J,M> WHERE 2L+1 = NORB 
//  VALID ARGUMENTS NORB > 0
//                  N1, N2 = 1,2,...,NORB
//                  JJ = 0,1,...NORB-1
	assert(n_orb>0);
	assert(n1 <= n_orb && n1 > 0);
	assert(n2 <= n_orb && n2 > 0);
	assert(jj < n_orb && jj >= 0);
	double cg = 0; 
	int mtot = n1 + n2 - n_orb - 1 ;
	if(jj < abs(mtot)) return cg;
	double fac, aa, b1, b2, b3, b4, factor;
	b1 = binom(n_orb+jj,2*jj+1);
	b2 = binom(2*jj,jj+n1+n2-n_orb-1);
	b3 = binom(n_orb - 1, n1 - 1);
	b4 = binom(n_orb - 1, n2 - 1);
	aa = sqrt(b1*b2*b3*b4);
	b1 = binom(n_orb-1, jj);
	fac = b1 / aa;
	int k1, k2, k3, l1, l2, l3, llmin, llmax, kmin, kmax, ll;
	k1 = 0;
      	k2 = n_orb-n1-jj;
      	k3 = n2-1-jj;
      	l1 = n_orb - 1 - jj;
      	l2 = n_orb - n1;
      	l3 = n2 - 1; 
      llmin = k1;
      if(k2>llmin) llmin = k2;
      if(k3>llmin) llmin = k3;

      llmax = l1;
      if(l2<llmax) llmax = l2;
      if(l3<llmax) llmax = l3;
      int isign = 1;
      if(llmin % 2 == 1) isign = -1; 
 
	factor = 0;
	kmin = llmin + 1;
	kmax = llmax + 1;
	for(int k = kmin; k<=kmax; k++)
	{
		ll = k - 1;
     	 	b1 =  binom(l1-k1,l1-ll);
      		b2 =  binom(l2-k2,l2-ll);
      		b3 =  binom(l3-k3,l3-ll);
      		factor += isign*b1*b2*b3;	
		isign = -isign; 
	}
	cg = fac * factor;
	return cg;
	 
}

  

double h_element_sphere(int m1,int m2,int m3,int m4) {
    
   if ((m1 + m2) != (m3 + m4)){
    return 0.0;
  }
     	int n1 = m1+1;
      	int n2 = m2+1;
      	int n3 = m3+1;
      	int n4 = m4+1;
     	int m =  m1+m2-n_orb+1;
	double temp = 0;
	for(int j = 0; j < n_orb; j++)   
	{

	if((abs(m) <= j) && ((n_orb-j)%2==0))
	{	
      	double cg1 =  cgcoef(j,n1,n2);
      	double cg2 =  cgcoef(j,n4,n3); 
	temp = temp + cg1 * cg2 * SPHERE_PSEUD[n_orb-1-j];
	}
	} 
	return -temp; 
}

double h_element_disk(int m1,int m2,int m3,int m4,vector<int>& a,vector<double>& coeff) {
    if ((m1 + m2) != (m3 + m4)){//Only consider isotropic cases
        return 0.0;
    }
    double sum=0;

    for (int i=0;i<a.size();i++) {

        mm=nn=a[i];
        double factor= pow(1/sqrt(2), m1+m2+m3+m4)*tgamma(m1+m2-mm+1)*sqrt(tgamma(mm+1)*tgamma(nn+1))
	    *sqrt(tgamma(m1+1)*tgamma(m2+1)*tgamma(m3+1)*tgamma(m4+1));

        for (int i=0;i<m2+1;i++) {
	        if (i>mm) break;
	        for (int j=0;j<m4+1;j++) {
		        if (j>nn) break;
		        for (int ii=0;ii<m1+1;ii++) {
			        if (ii+i>mm) break;
			        for (int jj=0;jj<m3+1;jj++) {
				        if (jj+j>nn) break;
				        if (ii+i==mm&&jj+j==nn) {
					        if ((i+j)%2==1) {
						        sum+=-factor/(tgamma(i+1)*tgamma(m2-i+1)*tgamma(ii+1)*tgamma(m1-ii+1)*tgamma(j+1)
						        *tgamma(m4-j+1)*tgamma(jj+1)*tgamma(m3-jj+1));
					        } else {
						        sum+=factor/(tgamma(i+1)*tgamma(m2-i+1)*tgamma(ii+1)*tgamma(m1-ii+1)*tgamma(j+1)
						        *tgamma(m4-j+1)*tgamma(jj+1)*tgamma(m3-jj+1));
					        }
				        }

			        }
		        }
	        }
        }
    }
    return sum;
}

int sign_check(int m1,int m2, int m3, int m4, vector<int>& state1,vector<int>& state2) {

	int sign=1;
	if (m1>m3) {
		for (int i=m3+1;i<m1;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}
	if (m1<m3) {
		for (int i=m1+1;i<m3;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}

	}

	if (m2<m4) {

		for (int i=m2+1;i<m4;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}

	if (m2>m4) {

		for (int i=m4+1;i<m2;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}
	return sign;
}

int sign_check_3b(int n3,int n2, int n1, int n4, int n5, int n6,vector<int>& state1,vector<int>& state2) {

    int sign=1;
	if ((n6 - n5) > 1) {
	    for (int jl = n5 + 1; jl < n6; jl++) {
			if (state2[jl] == 1) sign = -sign;
	    }
	}
   	if(n4 > 0){
		for (int jl = 0; jl < n4; jl++){
			if (state2[jl] == 1) sign = -sign;
		}
	}	  
	int isg = 1;
	if ((n1 - n2) > 1) {
	    for (int ik = n2 + 1; ik < n1; ik++) {
			if ((ik != n4) && (ik != n5) && (ik != n6) && (state2[ik] == 1)) isg = -isg;
	    }
	}

	if (n3 > 0) {
	    for (int ik = 0; ik < n3; ik++) {
			if ((ik != n4) && (ik != n5) && (ik != n6) && (state2[ik] == 1)) isg = -isg;
	    }
	}


/*	if (m1>m4) {
		for (int i=m4+1;i<m1;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}
	if (m1<m4) {
		for (int i=m1+1;i<m4;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}

	}

	if (m2<m5) {

		for (int i=m2+1;i<m5;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}

	if (m2>m5) {

		for (int i=m5+1;i<m2;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}

    if (m3<m6) {

		for (int i=m3+1;i<m6;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}

	if (m3>m6) {

		for (int i=m6+1;i<m3;i++) {
			if (state1[i]==1) sign=sign*(-1);
		}
	}*/
	return sign*isg;
}

void two_body_sphere(vector<double>& vm1m2m3m4) {

	for (int i=0;i<pow(n_orb,4);i++) {
		vm1m2m3m4.push_back(0);
	}

	for (int i=0;i<n_orb;i++) {
		for (int j=i+1;j<n_orb;j++) {
			for (int k=0;k<n_orb;k++) {
				for (int kk=k+1;kk<n_orb;kk++) {
					if ((i+j)!=(k+kk)) {
						continue;
					}

					vm1m2m3m4[i*pow(n_orb,3)+j*pow(n_orb,2)+k*n_orb+kk]+=h_element_sphere(i,j,k,kk);

				}
			}
		}
	}
}	

void b4_sphere(vector<double>& vm1m2m3m4) {

	for (int i=0;i<pow(n_orb,4);i++) {
		vm1m2m3m4.push_back(0);
	}

	for (int i=0;i<n_orb;i++) {
		for (int j=i+1;j<n_orb;j++) {
			for (int k=0;k<n_orb;k++) {
				for (int kk=k+1;kk<n_orb;kk++) {
					if ((i+j)!=(k+kk)) {
						continue;
					}

					vm1m2m3m4[i*pow(n_orb,3)+j*pow(n_orb,2)+k*n_orb+kk]+=h_element_sphere(i,j,k,kk);

				}
			}
		}
	}
}

void two_body_disk(vector<double>& vm1m2m3m4,vector<int>& a,vector<double>& coeff) {

	for (int i=0;i<pow(n_orb,4);i++) {
		vm1m2m3m4.push_back(0);
	}

	for (int i=0;i<n_orb;i++) {
		for (int j=i+1;j<n_orb;j++) {
			for (int k=0;k<n_orb;k++) {
				for (int kk=k+1;kk<n_orb;kk++) {
					if ((i+j)!=(k+kk)) {
						continue;
					}

					vm1m2m3m4[i*pow(n_orb,3)+j*pow(n_orb,2)+k*n_orb+kk]+=h_element_disk(i,j,k,kk,a,coeff);

				}
			}
		}
	}
}

void b4_disk(vector<double>& vm1m2m3m4,vector<int>& a,vector<double>& coeff) {

	for (int i=0;i<pow(n_orb,4);i++) {
		vm1m2m3m4.push_back(0);
	}

	for (int i=0;i<n_orb;i++) {
		for (int j=i+1;j<n_orb;j++) {
			for (int k=0;k<n_orb;k++) {
				for (int kk=k+1;kk<n_orb;kk++) {
					if ((i+j)!=(k+kk)) {
						continue;
					}

					vm1m2m3m4[i*pow(n_orb,3)+j*pow(n_orb,2)+k*n_orb+kk]+=h_element_disk(i,j,k,kk,a,coeff);

				}
			}
		}
	}
}

void setzero(double *vec,int x) {
	for (int i1=0;i1<x;i1++) {
		vec[i1]=0;
	}
}
void sumtospin(INT a,int *vec,int dim,int stat){
	
	
	
	INT b=a-1;
	INT b1;
	int dim2;
	
	for (int i=0;i<dim;i++) {
		
		vec[i]=0;
		
	}
	
	dim2=dim;
	
	if (stat==1) dim2=dim*dim;
	
	int tvec[dim2];
	
	for (int i1=1;i1<dim2+1;i1++) {
		
		if (log(b)*(1/log(2.0))<dim2-i1-1) {
			
			b1=b;
			tvec[dim2-i1]=0;
		}
		else {
			
			b1=b%(INT (pow(2,dim2-i1)));
			
			
			tvec[dim2-i1]=(b-b1)/(pow(2,(dim2-i1)));
			b=b1;
			
		}
		
	}
	
	if (stat==1) {
		
		int count=0,ptr=0;
		
		for (int i=0;i<dim2;i++) {
			
			if (tvec[i]==1) {
				
				vec[i+ptr]++;
				
				ptr--;
				
			}
			
			if (i+ptr==dim) break;
			
		}
		
	} else {
		
		for (int i=0;i<dim;i++) {
			
			vec[i]=tvec[i];
			
		}
		
	}
	
}
void lower(int *vec,INT index,double *list,int len,int n,int op,int stat) {
	
	
	//list[1:len] stores the coefficients, list[len+1:2len] stores the coordinates
	
	//for (int i=0;i<len;i++) {
	
	//cout<<vec[i];
	
	//}
	
	//cout<<endl;
	
	
	
	int count1=0;
	
	for (int i3=0;i3<len-1;i3++) {
		
		
		
		if (vec[i3]!=0) {
			
			if (vec[i3+1]==1&&stat==-1) continue;
			count1++;
			
			/*
			 !------------------------------
			 ! using the normalization of angular
			 ! momentum ladder operator
			 ! L^+ -> sqrt((j-m)(j+m+1))
			 ! L^- -> sqrt((j+m)(j-m+1))
			 !-------------------------------*/
			
			if (op==0) {
				
				int fact=1;
				
				for (int i=1;i<=vec[i3-1];i++) {
					
					fact=fact*i;
					
				}
				list[count1-1]=sqrt(vec[i3])*sqrt(vec[i3+1]+1)*sqrt((i3+1)*(n-i3));
				
				//cout<<"cof  "<<sqrt(vec[i3-1])*sqrt(vec[i3-2]+1)*sqrt((i3-1)*(n-i3+2))<<"  "<<n<<"  "<<i3<<endl;
			}
			else {
				list[count1-1]= i3-1;
			}
			
			if (stat==-1) list[len+count1-1]= index-pow(2,i3)+pow(2,i3+1);
			
			if (stat==1) {
				
				int tvec[len];
				
				for (int i=0;i<len;i++) {
					
					tvec[i]=vec[i];
					
				}
				tvec[i3-2]++;
				tvec[i3-1]--;
				
				int pointer=0;
				
				
				for (int i=0;i<len;i++) {
					
					if (tvec[i]!=0) {
						
						for (int i2=1;i2<=tvec[i];i2++) {
							
							list[len+count1-1]=list[len+count1-1]+pow(2,i+pointer);
							
							pointer++;
							
						}
						
					}
					
				}
				
				//cout<<"evec  "<<list[len+count1-1]<<endl;
				
			}
			
			
		}
		
	} 
	
	
}
void spmatmul(vector<double> &n1,vector<double> &n1r,vector<double> &n1c,\
			  vector<double> &n2,vector<double> &n2r,vector<double> &n2c,int len1,\
			  int len2,vector<double> &output,vector<double> &outputr,vector<double> &outputc) {
	
	
	double *m1 = new double[len1];
	double *m1c = new double[len1];
	double *m2 = new double[len2];
	double *m2r = new double[len2];
	double *m1r = new double[len1];
	double *m2c = new double[len2];
	
	
	for (int i=0;i<len1;i++) {
		
		m1[i]=n1[i];
		m1r[i]=n1r[i];
		m1c[i]=n1c[i];
		
	}
	
	for (int i=0;i<len2;i++) {
		
		m2[i]=n2[i];
		m2r[i]=n2r[i];
		m2c[i]=n2c[i];
		
	}
	
	sort3(m2c,m2,m2r,len2);
	sort3(m1r,m1,m1c,len1);
	
	n1.clear();n1r.clear();n1c.clear();n2.clear();n2r.clear();n2c.clear();
	
	
	double *m1rtemp,*m2ctemp;
	m1rtemp = new double[len1+1];
	m2ctemp = new double[len2+1];
	
	
	for (int i=0;i<len1;i++) {
		
		m1rtemp[i]=m1r[i];
		
	}
	
	
	m1rtemp[len1]=0;
	
	for (int i=0;i<len2;i++) {
		
		m2ctemp[i]=m2c[i];
		
	}
	m2ctemp[len2]=0;
	
	
	int check1=0,check2=0,flag=0;
	double workspace=0;
	
	float percentage,co=0;
	
	for (int i1=0;i1<len2;i1++) {
		
		percentage=100*i1/len2;
		
		
		
		if (m2ctemp[i1+1]!= m2ctemp[i1]) {
			
			for (int i2=0;i2<len1;i2++) {
				
				if (m1rtemp[i2+1]!= m1rtemp[i2]) {
					
					flag=0;
					
					
					
					for (int i3=check1;i3<i1+1;i3++) {
						
						for (int i4=check2;i4<i2+1;i4++) {
							
							
							if (m1c[i4]==m2r[i3]) {
								
								
								
								flag=1;
								workspace=workspace+m1[i4]*m2[i3];
								break;
								
							}
						}
					}
					
					if (flag==1&&fabs(workspace)>1e-17) {
						
						outputr.push_back(m1r[i2]);
						outputc.push_back(m2c[i1]);
						output.push_back(workspace);
						workspace=0;
						
					}
					check2=i2+1;
				}
				
			}
			
			check1=i1+1;
			check2=0;
			
		}
		
	}
	
	delete [] m1rtemp;
	delete [] m2ctemp;
	delete [] m1;
	delete [] m2;
	delete [] m1r;
	delete [] m1c;
	delete [] m2r;
	delete [] m2c;
	
	
}
void lminus(vector<INT> &vec, vector<vector<int>> &binvec,vector<double> &coeff,int ind) {

    vector<INT> index;vector<double> index_coeff;
	for (int i=0;i<binvec.size();i++) {
		
		for (int j=0;j<n_orb-1;j++) {
            if (binvec[i][j]==1&&binvec[i][j+1]!=1) {
                index.push_back(vec[i]-pow(2,j)+pow(2,j+1));				
				index_coeff.push_back(sqrt((j+1)*(n_orb-1-j))*coeff[i]);				
			}
        }
    }

    sort2b(index,index_coeff);

    double sum=0,ss=0;
    for (int i=0;i<index.size()-1;i++) {
        sum+=index_coeff[i];
        if (index[i+1]!=index[i]) {
            ss+=sum*sum;
            sum=0;
        }
    }
    sum+=index_coeff[index_coeff.size()-1];
    ss+=sum*sum;
    cout<<ind<<" "<<ss<<endl;
	
}


double h_evaluate(vector<double>& vm1m2m3m4,vector<int>& bin1,vector<int>& bin2,int b) {

	double sum=0;
				
    int count=0;
    for (int k=0;k<bin1.size();k++) {
            if (bin1[k]==1&&bin2[k]!=1) count++;
	        if (count>2) return sum;
    }
    vector<int> m1m2m3m4;

    if (count==2) {

        for (int k=0;k<bin1.size();k++) {
            if (bin1[k]==1&&bin2[k]==0) m1m2m3m4.push_back(k);
        }
        for (int k=0;k<bin2.size();k++) {
            if (bin2[k]==1&&bin1[k]==0) m1m2m3m4.push_back(k);
        }

		if (b!=0&&m1m2m3m4[2]+m1m2m3m4[3]<m1m2m3m4[0]+m1m2m3m4[1]) {
			int t1=m1m2m3m4[2],t2=m1m2m3m4[3];
			m1m2m3m4[2]=m1m2m3m4[0];m1m2m3m4[3]=m1m2m3m4[1];
			m1m2m3m4[0]=t1;m1m2m3m4[1]=t2;
		}
                    
        sum+=vm1m2m3m4[m1m2m3m4[0]*pow(n_orb,3)+m1m2m3m4[1]*pow(n_orb,2)+m1m2m3m4[2]*n_orb+m1m2m3m4[3]]*sign_check(m1m2m3m4[0],m1m2m3m4[1],m1m2m3m4[2],m1m2m3m4[3],bin1,bin2);
        m1m2m3m4.clear();
    }

	if (count==1) {

		    vector<int> index_s;

                    for (int k=0;k<bin1.size();k++) {
                        if (bin1[k]==1&&bin2[k]==1) index_s.push_back(k);
                    }

		    for (int kk=0;kk<index_s.size();kk++) {
			m1m2m3m4.push_back(index_s[kk]);
			for (int k=0;k<bin1.size();k++) {
                            	if (bin1[k]==1&&bin2[k]==0) m1m2m3m4.push_back(k);
                        }
		        m1m2m3m4.push_back(index_s[kk]);
                        for (int k=0;k<bin2.size();k++) {
                            	if (bin2[k]==1&&bin1[k]==0) m1m2m3m4.push_back(k);
                        }

		        if (b!=0&&m1m2m3m4[2]+m1m2m3m4[3]<m1m2m3m4[0]+m1m2m3m4[1]) {
			    	int t1=m1m2m3m4[2],t2=m1m2m3m4[3];
			    	m1m2m3m4[2]=m1m2m3m4[0];m1m2m3m4[3]=m1m2m3m4[1];
			    	m1m2m3m4[0]=t1;m1m2m3m4[1]=t2;
		        }
			sum+=vm1m2m3m4[m1m2m3m4[0]*pow(n_orb,3)+m1m2m3m4[1]*pow(n_orb,2)+m1m2m3m4[2]*n_orb+m1m2m3m4[3]]*sign_check(m1m2m3m4[0],m1m2m3m4[1],m1m2m3m4[2],m1m2m3m4[3],bin1,bin2);

			m1m2m3m4.clear();
		    }
		    index_s.clear();

	}

                
        if (count==0) {
                    
                  vector<int> index;
                  for (int k=0;k<bin1.size();k++) {
                        if (bin1[k]==1) index.push_back(k);
                  }
                  for (int k=0;k<index.size();k++) {
                        for (int kk=k+1;kk<index.size();kk++) {
                            	m1m2m3m4.push_back(index[k]);m1m2m3m4.push_back(index[kk]);
                            	m1m2m3m4.push_back(index[k]);m1m2m3m4.push_back(index[kk]);
                            	sum+=vm1m2m3m4[m1m2m3m4[0]*pow(n_orb,3)+m1m2m3m4[1]*pow(n_orb,2)+m1m2m3m4[2]*n_orb+m1m2m3m4[3]];m1m2m3m4.clear();
							
                        }
                  }
		  index.clear();

        }

	return sum;

}

void hamiltonian_construction(vector<INT> &basis,vector<double>& vm1m2m3m4) {

	vector<int> sstate;
	vector<INT> evec;
	vector<vector<int>> evec_bin;

	for (int i=0;i<basis.size();i++) {
		dec_bin(basis[i],sstate);
		evec_bin.push_back(sstate);
	}

	for (int i=0;i<basis.size();i++) {
		double value=h_evaluate(vm1m2m3m4,evec_bin[i],evec_bin[i],0);
		mdiag.push_back(value);
	}

	vector<int> holder;
	for (int i=0;i<basis.size();i++) {
		m_row.push_back(holder);
	}

    for (int i=0;i<basis.size();i++) {
		vector<double> holder1;
		vector<int> holder2;
		for (int j=i+1;j<basis.size();j++) {
			double value=h_evaluate(vm1m2m3m4,evec_bin[i],evec_bin[j],0);
			if (fabs(value)>1e-16) {
				holder1.push_back(value);holder2.push_back(j);
				m_row[j].push_back(i);
			}
		}
		moffdiag.push_back(holder1);m_col.push_back(holder2);
	}
	evec_bin.clear();
}

int binary_search(vector<int>& array,int n) {

	int s_size=array.size();
	int ss=0,s1=0,s2=s_size;
	while (array[ss]!=n) {
		if (n>array[ss]) {
			if (s2==ss+1) return ss+1;
			s1=ss;
			ss=(ss+s2)/2;
		}
		if (n<array[ss]) {
			if (s1==ss-1) return ss-1;
			s2=ss;
			ss=(s1+ss)/2;
		}

	}
	return ss;

}

double matvec1(vector<double>& v_in,vector<double> &v_out) {

	int h_dim=v_in.size();
	double ss,sum=0;
    for (int i=0;i<h_dim;i++) {  
        ss=v_in[i]*mdiag[i];
		for (int j=0;j<moffdiag[i].size();j++) {
			ss+=v_in[m_col[i][j]]*moffdiag[i][j];
		}
		vector<int> temp;     
		for (int j=0;j<m_row[i].size();j++) {
			int kk=binary_search(m_col[m_row[i][j]],i);
			temp.push_back(kk);
			ss+=v_in[m_row[i][j]]*moffdiag[m_row[i][j]][kk];
		}
		sum+=ss*v_out[i];
		record.push_back(temp);
    }

	return sum;

}

double matvec(vector<double>& v_in,vector<double> &v_out) {

	int h_dim=v_in.size();
	double ss,sum=0;
    for (int i=0;i<h_dim;i++) {  
        ss=v_in[i]*mdiag[i];
		for (int j=0;j<moffdiag[i].size();j++) {
			ss+=v_in[m_col[i][j]]*moffdiag[i][j];
		}
		for (int j=0;j<m_row[i].size();j++) {
			ss+=v_in[m_row[i][j]]*moffdiag[m_row[i][j]][record[i][j]];
		}
		sum+=ss*v_out[i];
    }

	return sum;

}

double onebody(vector<vector<int>>& evec_bin1,vector<vector<int>>& evec_bin2,INT s1,INT s2,vector<vector<double>>& vmn) {

	int count=0;
	double sum=0;
	int k1,k2;

    for (int k=0;k<evec_bin1[s1].size();k++) {
            if (evec_bin1[s1][k]==1&&evec_bin2[s2][k]!=1) count++;
			k1=k;
			if (evec_bin1[s1][k]==0&&evec_bin2[s2][k]==1) k2=k;
	        if (count>1) return sum;
    }

	if (count==0) {

		for (int i=0;i<evec_bin1[s1].size();i++) {
			sum+=vmn[i][0];
		}

		return sum;
	}

	int dk=k2-k1;

	int sign=1;

	if (dk>0) {

		for (int i=k1+1;i<k2;i++) {
			if (evec_bin1[s1][i]==1) sign=-sign;
		}

		return sign*vmn[k1][dk];
	} else {

		for (int i=k1-1;i>k2;i--) {
			if (evec_bin1[s1][i]==1) sign=-sign;
		}
		return sign*vmn[k1][-dk];
	}
	
}

/*extern "C"{ 

	void v123_(int *, int *, int *, int *, int *, int*, int *, double *);
} 

double three_bodyCG_sphere(int m1, int m2, int m3)
{

      int NORB=n_orb;
      double s = (NORB-1.0)/2.0;
      int lz1 = round((-s + m1)*2);
      int lz2 = round((-s + m2)*2);
      int lz3 = round((-s + m3)*2);
      int phi = NORB - 1;
      double result=0; 
      v123_(&m1, &lz1, &m2, &lz2, &m3, &lz3, &phi, &result);
      return result;
}

double three_body(int m1, int m2, int m3, int m4, int m5, int m6) {
	double res = 0.0;
	int stat = -1;
	if ((m1+m2+m3) != (m4+m5+m6))
	{ 
		res = 0.0; 
	}
	else
	{  
	res = 1.0 * three_bodyCG_sphere(m1, m2, m3) * three_bodyCG_sphere(m6, m5, m4) / 36.0 ;
	}
	return res; 
} 

void setccco3B(int norb, int m1, int m2, int m3, int m4, int m5, double* ccco, double cgd) {
  ccco[m1 * norb * norb * norb * norb + m2 * norb * norb * norb + m3 * norb * norb + m4 * norb + m5] = cgd;
} 

double getccco3B(int norb, int m1, int m2, int m3, int m4, int m5, double* ccco) {
  return ccco[m1 * norb * norb * norb * norb + m2 * norb * norb * norb + m3 * norb * norb + m4 * norb + m5];
}

void threebodyH() // 3BH on sphere
{   
	int stat = -1;
	int norb = n_orb; 
	double cg1;
    ccco3B = new double[norb * norb * norb * norb * norb];
    int q1 = 0;
 	for (int q4 = 0; q4 < norb; q4++) {
        for (int q5 = q4+1; q5 < norb; q5++) {
            for (int q6 = q5+1; q6 < norb; q6++) {
				for (int q3 = 0; q3 < norb; q3++) {
					for (int q2 = q3+1; q2 < norb; q2++) {
						q1 = q6 + q5 + q4 - q3 - q2;
						if ((q1 >= 0) && (q1 < norb) && (q1 > q2)) {
 							cg1 = three_body(q1,q2,q3,q4,q5,q6);
        			        cg1 = cg1 - three_body(q1,q3,q2,q4,q5,q6);
        			        cg1 = cg1 - three_body(q2,q1,q3,q4,q5,q6);
        			        cg1 = cg1 + three_body(q2,q3,q1,q4,q5,q6);
        			        cg1 = cg1 + three_body(q3,q1,q2,q4,q5,q6);
        			        cg1 = cg1 - three_body(q3,q2,q1,q4,q5,q6);
							setccco3B(norb, q1, q2, q3, q4, q5, ccco3B, cg1);  
						}
      		        }
    		    }
  		    }
 		}
	}	 	 
}

*/

/*double h3_evaluate(vector<vector<int>>& evec_bin1,vector<vector<int>>& evec_bin2,INT s1,INT s2,int b) {

    /*double sum=0;

    if (s1==s2) {

        for (int k = 0; k < n_el; k++) {
            n3 = evec_occ1[k];
            for (int i = k + 1; i < n_el; i++) {
	            n2 = evec_occ1[i];
	            for (int j = i + 1; j < NEL; j++) {
	                n1 = evec_occ1[j];
	                sum+= getccco3B(n_orb, n1, n2, n3, n3, n2, ccco3B);   
                }
            } 

        }
        
        return sum;

    }




	double sum=0;
				
    int count=0;
    for (int k=0;k<evec_bin1[s1].size();k++) {
            if (evec_bin1[s1][k]==1&&evec_bin2[s2][k]!=1) count++;
	        if (count>3) return sum;
    }
    int m1,m2,m3,m4,m5,m6;
    double value;


    if (count==3) {

        vector<int> m1m2m3,m4m5m6;

        for (int k=0;k<evec_bin1[s1].size();k++) {
            if (evec_bin1[s1][k]==1&&evec_bin2[s2][k]==0) m1m2m3.push_back(k);
        }
        for (int k=0;k<evec_bin2[s2].size();k++) {
            if (evec_bin2[s2][k]==1&&evec_bin1[s1][k]==0) m4m5m6.push_back(k);
        }

        m1=m1m2m3[0];m2=m1m2m3[1];m3=m1m2m3[2];m4=m4m5m6[0];m5=m4m5m6[1];m6=m4m5m6[2];  
        value=getccco3B(n_orb,m3,m2,m1,m4,m5,ccco3B)*sign_check_3b(m1,m2,m3,m4,m5,m6,evec_bin1[s1],evec_bin2[s2]);
        m1m2m3.clear();m4m5m6.clear();
        sum+=value;
        cout<<"check "<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<value<<endl;

    }

	if (count==2) {

        vector<int> m1m2,m3m4;
        for (int k=0;k<evec_bin1[s1].size();k++) {
            if (evec_bin1[s1][k]==1&&evec_bin2[s2][k]==0) m1m2.push_back(k);
        }
        for (int k=0;k<evec_bin2[s2].size();k++) {
            if (evec_bin2[s2][k]==1&&evec_bin1[s1][k]==0) m3m4.push_back(k);
        }
        m1=m1m2[0];m2=m1m2[1];m3=m3m4[0];m4=m3m4[1];
        for (int k=0;k<m1;k++) {
            if (evec_bin1[s1][k]==1&&evec_bin2[s2][k]==1) {
                if (k<m3) value=getccco3B(n_orb,m2,m1,k,k,m3,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                if (k>m3&&k<m4) value=getccco3B(n_orb,m2,m1,k,m3,k,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                if (k>m4) value=getccco3B(n_orb,m2,m1,k,m3,m4,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                cout<<"check1 "<<m1<<" "<<m2<<" "<<k<<" "<<m3<<" "<<m4<<" "<<k<<" "<<value<<endl;
                sum+=value;
            }
        }
        for (int k=m1+1;k<m2;k++) {
            if (evec_bin1[s1][k]==1&&evec_bin2[s2][k]==1) {
                if (k<m3) value=getccco3B(n_orb,m2,k,m1,k,m3,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                if (k>m3&&k<m4) value=getccco3B(n_orb,m2,k,m1,m3,k,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                if (k>m4) value=getccco3B(n_orb,m2,k,m1,m3,m4,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                cout<<"check2 "<<m1<<" "<<m2<<" "<<k<<" "<<m3<<" "<<m4<<" "<<k<<" "<<value<<endl;
                sum+=value;
            }


        }
        for (int k=m2+1;k<evec_bin2[s2].size();k++) {
            if (evec_bin1[s1][k]==1&&evec_bin2[s2][k]==1) {
                if (k<m3) value=getccco3B(n_orb,k,m2,m1,k,m3,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                if (k>m3&&k<m4) value=getccco3B(n_orb,k,m2,m1,m3,k,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                if (k>m4) value=getccco3B(n_orb,k,m2,m1,m3,m4,ccco3B)*sign_check(m1,m2,m3,m4,evec_bin1[s1],evec_bin2[s2]);
                cout<<"check3 "<<m1<<" "<<m2<<" "<<k<<" "<<m3<<" "<<m4<<" "<<k<<" "<<value<<endl;
                sum+=value;
            }
        }



    }
                
    if (count==0) {
                    
        vector<int> index;
        for (int k=0;k<evec_bin1[s1].size();k++) {
            if (evec_bin1[s1][k]==1) index.push_back(k);
        }
        for (int k=0;k<index.size();k++) {
            for (int kk=k+1;kk<index.size();kk++) {
                for (int kkk=kk+1;kkk<index.size();kkk++) {
                    value=getccco3B(n_orb,index[kkk],index[kk],index[k],index[k],index[kk],ccco3B);
                    cout<<"check "<<index[k]<<" "<<index[kk]<<" "<<index[kkk]<<" "<<index[k]<<" "<<index[kk]<<" "<<index[kkk]<<" "<<value<<endl;
                    sum+=value;
                }
							
            }
        }
		index.clear();

    }

	return sum;

}

*/
void write(vector<double>& y,string filename) {
	
	FILE *output;
	
	char *fn;
	
	fn = new char [filename.size()+1];
	strcpy (fn, filename.c_str());
	
	output = fopen(fn,"w");
	for (int i=0;i<y.size();i++) {
		
		fprintf(output,"%.15f\n",y[i]);
		
	}
	
	delete[] fn;
	fclose(output);
	
}

void writewf(vector<INT>& evec,vector<double>& coeff,string filename) {
	

	FILE *output;
	
	char *fn;
	
	fn = new char [filename.size()+1];
	strcpy (fn, filename.c_str());
	
	output = fopen(fn,"w");
	assert(evec.size()==coeff.size());
	fprintf(output,"%d\n",evec.size());
	vector<int> bin;
	for (int i=0;i<evec.size();i++) {

		if (fabs(coeff[i])<1e-12) coeff[i]=0;
		fprintf(output,"%llu\n",evec[i]);
		fprintf(output,"%.15f\n",coeff[i]);
		
	}
	
	delete[] fn;
	fclose(output);	

}
	



	
	
	
	
	
	
		
		
		


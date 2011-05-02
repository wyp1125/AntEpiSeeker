#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;
int iAntCount;
int iItCount;
int iLociCount;
int iLociModel;
int iEpiModel;
double alpha;
double rou;
double phe;
double pvalue;
int iTopModel;
int** loci_TopModel;
int iTopLoci;
int* loci_TopLoci;
double* eva_TopModel;
double* phe_TopLoci;
vector<vector<double> > sigepi;
int lcn; 
class SNP
{
public:
	int iSample;
	int iLoci;
	bool** data;
	char** SNPnames;
	int classvalues[2];
        void input_data(char* path);
        void destroy();
        void setpheromone(double level);
        double* pheromone;
        double* cdf;
        void display_pheromone();
};
void SNP::destroy()
{
int i;
for(i=0;i<iLoci;i++)
delete []SNPnames[i];
delete []SNPnames;
for(i=0;i<iSample;i++)
{
delete []data[i];
}
delete []data;
delete []pheromone;
delete []cdf;
}
void SNP::setpheromone(double level)
{
int i;
cdf[0]=0;
for(i=0;i<iLoci-1;i++)
{
pheromone[i]=level;
cdf[i+1]=double(i+1)/double(iLoci-1);
}
}
void SNP::display_pheromone()
{
int i;
for(i=10;i<111;i++)
{
std::cout<<setprecision(5)<<pheromone[i];
cout<<" ";
}
cout<<endl;
}
void SNP::input_data(char* path)
{
        cout<<"-----Reading data-----"<<endl;
	int i,j,temp;
        string line;
	ifstream in(path);
	getline(in,line);
	istringstream test(line);
	i=0;
	string word;
	while(!test.eof())
	{
		getline(test,word,',');
		i++;
	}
	iLoci=i;
	//cout<<SNPdata.iLoci<<endl;
	j=0;
	while(!in.eof())
	{       
		getline(in,line);
                //cout<<line<<endl;
		j++;
	}
	iSample=j-1;
	//cout<<iSample<<endl;
	in.close();
        pheromone=new double[iLoci-1];
        cdf=new double[iLoci];
        SNPnames=new char*[iLoci];
        for(i=0;i<iLoci;i++)
        SNPnames[i]=new char[20];
	data=new bool* [iSample];
	for(i=0;i<iSample;i++)
		data[i]=new bool [2*iLoci];
        classvalues[0]=0;
        classvalues[1]=0;
	ifstream in1(path);
        getline(in1,line);
	istringstream test1(line);
	i=0;
	while(!test1.eof())
	{
		if(i==iLoci)
                   break;
                getline(test1,word,',');
                strcpy(SNPnames[i],word.c_str());
		i++;
	}
        i=0;
	while(!in1.eof())
	{
                if(i==iSample)
                {
                    break;
                }
		getline(in1,line);
		istringstream values(line);
		j=0;
		while(!values.eof())
		{
                     getline(values,word,',');
                     istringstream int_iss(word);
                     int_iss>>temp;
                     if(j==2*iLoci-2)
                     {
                      classvalues[temp]++;
                     }
                     if(temp==0)
                     {
                      data[i][j++]=0;
                      data[i][j++]=0;
                      }
                     else if(temp==1)
                     {
                      data[i][j++]=0;
                      data[i][j++]=1;
                      }
                     else if(temp==2)
                     {
                      data[i][j++]=1;
                      data[i][j++]=0;      
                     }
                     else
                     {
                      data[i][j++]=1;
                      data[i][j++]=1;
                     }
		}
                i++;
	}
	in1.close();
        cout<<"-----Reading data completed!------"<<endl;
        cout<<"Number of loci: "<<iLoci-1<<endl;
        cout<<"Number of samples: "<<iSample<<endl;
}
SNP SNPdata;

#include<bits/stdc++.h>
#include<Windows.h>
#include<random>
#include<dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#define PI 3.14159265
#define tiaoshi  puts("what?")
using namespace std;
int tot=0;
map<string,int> contig2number,TNF2int;

static const std::string TN[] = { "GGTA", "AGCC", "AAAA", "ACAT", "AGTC", "ACGA", "CATA", "CGAA", "AAGT", "CAAA", "CCAG", "GGAC", "ATTA", "GATC", "CCTC",
		"CTAA", "ACTA", "AGGC", "GCAA", "CCGC", "CGCC", "AAAC", "ACTC", "ATCC", "GACC", "GAGA", "ATAG", "ATCA", "CAGA", "AGTA", "ATGA", "AAAT", "TTAA", "TATA",
		"AGTG", "AGCT", "CCAC", "GGCC", "ACCC", "GGGA", "GCGC", "ATAC", "CTGA", "TAGA", "ATAT", "GTCA", "CTCC", "ACAA", "ACCT", "TAAA", "AACG", "CGAG", "AGGG",
		"ATCG", "ACGC", "TCAA", "CTAC", "CTCA", "GACA", "GGAA", "CTTC", "GCCC", "CTGC", "TGCA", "GGCA", "CACG", "GAGC", "AACT", "CATG", "AATT", "ACAG", "AGAT",
		"ATAA", "CATC", "GCCA", "TCGA", "CACA", "CAAC", "AAGG", "AGCA", "ATGG", "ATTC", "GTGA", "ACCG", "GATA", "GCTA", "CGTC", "CCCG", "AAGC", "CGTA", "GTAC",
		"AGGA", "AATG", "CACC", "CAGC", "CGGC", "ACAC", "CCGG", "CCGA", "CCCC", "TGAA", "AACA", "AGAG", "CCCA", "CGGA", "TACA", "ACCA", "ACGT", "GAAC", "GTAA",
		"ATGC", "GTTA", "TCCA", "CAGG", "ACTG", "AAAG", "AAGA", "CAAG", "GCGA", "AACC", "ACGG", "CCAA", "CTTA", "AGAC", "AGCG", "GAAA", "AATC", "ATTG", "GCAC",
		"CCTA", "CGAC", "CTAG", "AGAA", "CGCA", "CGCG", "AATA" };
		
static const std::string TN_reverse_complement[] = {
    "TACC", "GGCT", "TTTT", "ATGT", "GACT", "TCGT", "TATG", "TTCG", "ACTT", "TTTG", "CTGG", "GTCC", "TAAT", "GATC", "GAGG",
    "TTAG", "TAGT", "GCCT", "TTGC", "GCGG", "GGCG", "GTTT", "GAGT", "GGAT", "GGTC", "TCTC", "CTAT", "TGAT", "TCTG", "TACT",
    "TCAT", "ATTT", "TTAA", "TATA", "CACT", "AGCT", "GTGG", "GGCC", "GGGT", "TCCC", "GCGC", "GTAT", "TCAG", "TCTA", "ATAT",
    "TGAC", "GGAG", "TTGT", "AGGT", "TTTA", "CGTT", "CTCG", "CCCT", "CGAT", "GCGT", "TTGA", "GTAG", "TGAG", "TGTC", "TTCC",
    "GAAG", "GGGC", "GCAG", "TGCA", "TGCC", "CGTG", "GCTC", "AGTT", "CATG", "AATT", "CTGT", "ATCT", "TTAT", "GATG", "TGGC",
    "TCGA", "TGTG", "GTTG", "CCTT", "TGCT", "CCAT", "GAAT", "TCAC", "CGGT", "TATC", "TAGC", "GACG", "CGGG", "GCTT", "TACG",
    "GTAC", "TCCT", "CATT", "GGTG", "GCTG", "GCCG", "GTGT", "CCGG", "TCGG", "GGGG", "TTCA", "TGTT", "CTCT", "TGGG", "TCCG",
    "TGTA", "TGGT", "ACGT", "GTTC", "TTAC", "GCAT", "TAAC", "TGGA", "CCTG", "CAGT", "CTTT", "TCTT", "CTTG", "TCGC", "GGTT",
    "CCGT", "TTGG", "TAAG", "GTCT", "CGCT", "TTTC", "GATT", "CAAT", "GTGC", "TAGG", "GTCG", "CTAG", "TTCT", "TGCG", "CGCG",
    "TATT"
};

void preprocess(){
	for(int i=0;i<136;i++)
		TNF2int[TN[i]]=TNF2int[TN_reverse_complement[i]]=i;
}

//Palindromic sequences
static const std::string TNP[] = { "ACGT", "AGCT", "TCGA", "TGCA", "CATG", "CTAG", "GATC", "GTAC", "ATAT", "TATA", "CGCG", "GCGC", "AATT", "TTAA", "CCGG", "GGCC" };

struct Sequence{
	string SEQ;
	vector<int> seq;
	int TNF[256];
	int number;
	string Realname;
	double sigma,average,length=0;
	int category;
	void initial(){memset(TNF,0,sizeof(TNF));}
	void transfer(){
		if(SEQ.size()>1000000){
			cout<<"oversize!";
		}
		for(int i=0;i<SEQ.size();i++){
			if(SEQ[i]=='a'||SEQ[i]=='A')
				seq.push_back(0);
			else if(SEQ[i]=='c'||SEQ[i]=='C')
				seq.push_back(1);
			else if(SEQ[i]=='g'||SEQ[i]=='G')
				seq.push_back(2);
			else seq.push_back(3);
		}
	}
	void statistic(){
		int LEN=SEQ.size();
		for(int i=0;i<=LEN-4;i++)
			TNF[TNF2int[SEQ.substr(i,4)]]++;
		for(int i=0;i<136;i++)
			length+=TNF[i]*TNF[i];
		length=sqrt(length);
	}
	Sequence(string h,string realname){
		this->SEQ=h;
		this->initial();
		this->transfer();
		this->statistic();
		this->number=++tot;
		this->Realname=realname;
		this->category=0;
		contig2number[realname]=tot;
//		cerr<<"struct realname: "<<realname<<" SEQ: "<<h<<'\n';
	}
};

double vector_compare(Sequence A,Sequence B){
//	cout<<"--------------------------------\n";
	double ans=0;
//	cout<<"SEQ1->";
//	for(int i=0;i<256;i++)cout<<A.TNF[i];
//	cout<<"\nSEQ2->";
//	for(int i=0;i<256;i++)cout<<B.TNF[i];
	for(int i=0;i<136;i++)
		ans+=A.TNF[i]*B.TNF[i];
//	cout<<"ans -> "<<ans<<'\n'<<"length: "<<length1<<" "<<length2<<'\n';
	double cosineangle=ans*1.0/(A.length*B.length);
//	cout<<"--------------------------------\n";
	return cosineangle;
}

vector<Sequence> Read_sequence_func(string filepath, int type=0){
	ifstream inputFile(filepath);
    vector<Sequence> DNA_DATA;
    if (!inputFile.is_open()) {
        cerr << "Could not open the file: " << filepath << std::endl;
        exit(1);
    }
    
    string line,name,fullline="";
    int number_of_sq=0;
    
    getline(inputFile,line);
    name=line.substr(1);
    istringstream iss(name);
    iss>>name;
    
    while (getline(inputFile, line)) {
        if(line[0]=='>'){
        	
        	Sequence TEMP(fullline,name);
    		DNA_DATA.push_back(TEMP);
    		number_of_sq++;
    		fullline="";
        		
        	name=line.substr(1);
        	istringstream iss(name);
        	iss>>name;
        	continue;
		}
		
		fullline+=line;
//		cerr<<"sq: "<<line<<'\n';
        
    }
    
    Sequence TEMP(fullline,name);
    DNA_DATA.push_back(TEMP);
    number_of_sq++;
    fullline="";
    
    cerr<<"total sequence number is: "<<number_of_sq<<'\n';
    inputFile.close();
    return DNA_DATA;
}

#define N 5000000

struct node{
	int u,v,w;
}e[N];
int nxt[N],first[N],cnt,MAXITERATION=3;
void buildgraph(int u,int v,int w){
//	cout<<u<<" link with "<<v<<'\n';
	e[++cnt].u=u;
	e[cnt].v=v;
	e[cnt].w=w;
	nxt[cnt]=first[u];
	first[u]=cnt;	
}
vector<Sequence> Meta,SCS;
void LPA(){
	std::vector<int> labels(tot+1);
    std::unordered_map<int, int> labelCounts;

    // Initialize labels (e.g., with node IDs or random values)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, tot - 1);
    for (int i = 0; i <=tot; ++i) {
        labels[i] = dis(gen);
    }
	labels[0]=1;
	for(int i=0;i<SCS.size();i++)
		labels[SCS[i].number]=1;
    for (int iter = 1; iter <= MAXITERATION; ++iter) {
    	cout<<"-------------------\n";
    	cout<<"For  iteration"<<iter<<" , the tags are: \n";
    	for (int i = 0; i <=tot; ++i) {
    	    cout<<labels[i]<<" ";
    	}
    	cout<<"\n\n";
        for (int i = 0; i <=tot; ++i) {
            std::vector<int> neighbors;
            for(int o=first[i];o;o=nxt[o])
            	neighbors.push_back(e[o].v);
            labelCounts.clear();

            // Count the occurrences of labels among neighbors
            for (int neighbor : neighbors) {
                int neighborLabel = labels[neighbor];
                labelCounts[neighborLabel]++;
            }

            // Choose the most frequent label as the new label for the current node
            int mostFrequentLabel = -1;
            int maxCount = 0;
            for (const auto& entry : labelCounts) {
                if (entry.second > maxCount) {
                    mostFrequentLabel = entry.first;
                    maxCount = entry.second;
                }
            }

            labels[i] = mostFrequentLabel;
        }
    }

    // Output the final labels
    std::cout << "Final labels after label propagation:" << std::endl;
    for (int i = 0; i<=tot; ++i) {
        std::cout << "Node " << i << ": Label " << labels[i] << std::endl;
    }
}

double sq(double x){return x*x;}
double PHI(double x){
	double X=x*sqrt(0.5);
	return erf(X)/2+0.5;
}
double abundancedistance(double mean1,double sigma1,double mean2,double sigma2){
	if(sigma2<sigma1)swap(sigma1,sigma2),swap(mean1,mean2);
	//normal_distribution<> norm1{mean1,sigma1},norm2{mean2,sigma2};
	double A=sq(sigma1)-sq(sigma2);
	double B=2*(sq(sigma2)*mean1-sq(sigma1)*mean2);
	double C=sq(sigma1)*sq(mean2)-sq(sigma2)*sq(mean1)-2*sq(sigma1)*sq(sigma2)*log(sigma1/sigma2);
	double k1_=(-B+sqrt(sq(B)-4*A*C))/(2*A);
	double k2_=(-B-sqrt(sq(B)-4*A*C))/(2*A);
	double k1=min(k1_,k2_),k2=max(k1_,k2_);
	//cout<<"k1->"<<k1<<" k2->"<<k2<<'\n';
	double k11,k12,k21,k22;
	k12=(k2-mean1)/sigma1;
	k11=(k1-mean1)/sigma1;
	k22=(k2-mean2)/sigma2;
	k21=(k1-mean2)/sigma2;
	double result=PHI(k12)-PHI(k11)+PHI(k21)-PHI(k22);
	return result;
} 

vector<vector<double>>MEAN,SIGMA;
int SAMPLE_SIZE=0;

double ADP(int i,int j, double c){
//	i=Meta[i].number,j=Meta[j].number;
	int flag=0;
	double P_ij=1;
	for(int n=0;n<SAMPLE_SIZE;n++){
		double mean_in=MEAN[n][i];
		double sigma_in=SIGMA[n][i];
		double mean_jn=MEAN[n][j];
		double sigma_jn=SIGMA[n][j];
	//	cout<<"mean_in->"<<mean_in<<"  sigma_in->"<<sigma_in<<"  mean_jn->"<<mean_jn<<"  sigma_jn->"<<sigma_jn<<'\n';
		int signal=(mean_in>c)|(mean_jn>c);
		flag+=signal;
		if(!signal)continue;
		P_ij*=abundancedistance(mean_in,sigma_in,mean_jn,sigma_jn);
		//cout<<"AD->"<<abundancedistance(mean_in,sigma_in,mean_jn,sigma_jn)<<'\n';
	}
	return pow(P_ij,1.0/flag);
}

double rad2deg(double rad){return rad*180*1.0/PI;}
void READ_PARTICULAR_STATISTIC_FILE(string folderpath,string filename,int nownumber,int size_meta){
	vector<double> newROW1(size_meta+1),newROW2(size_meta+1);
	string fullpath = folderpath + "\\" + filename;
	ifstream inputFile(fullpath);
	if(!inputFile.is_open()){
		cerr<<"unable to open the statistic file\n";
		exit(0);
	}
	string firstLine;
	getline(inputFile,firstLine);
	string current_string;
	int nowcount=0,refer;
	while(inputFile>>current_string){
		if(nowcount%5==0){
			string contig_name=current_string;
			refer=contig2number[contig_name];
//			cout<<"contig_name: "<<contig_name<<'\n';
//			cout<<"refer: "<<refer<<'\n';
			if(!refer){
				cerr<<"unmatched statistic data detected, now quit";
				exit(0);
			}
		}
		else if(nowcount%5==2){
			string contig_average=current_string;
//			cerr<<"average:"<<contig_average<<'\n';
			newROW1[refer]=stod(contig_average);
		}
		else if(nowcount%5==4){
			string contig_var=current_string;
//			cerr<<"variance:"<<contig_var<<'\n';
			newROW2[refer]=sqrt(stod(contig_var));
		}
		nowcount++;
		nowcount%=5;
	}
	MEAN.push_back(newROW1);
	SIGMA.push_back(newROW2);
}

void READSTAFILE(string folderpath,int size_meta){
	//cerr<<"folderpath: "<<folderpath<<'\n';
	DIR *dir;
	struct dirent* entry;
	dir=opendir(folderpath.c_str());
	if(dir==nullptr){
		std::cerr<<"unable to open the folder"<<endl;
		exit(721);
	}
	int nownumber=0;
	while((entry=readdir(dir))!=nullptr){
		struct stat filestat;
		string fullpath=folderpath+"\\"+entry->d_name;
		//cerr<<"fullpath:"<<fullpath<<'\n';
		if(stat(fullpath.c_str(),&filestat)==0){
			if(S_ISREG(filestat.st_mode)){
				cout<<"filename:"<<entry->d_name<<endl;
				READ_PARTICULAR_STATISTIC_FILE(folderpath,entry->d_name,nownumber,size_meta);
				SAMPLE_SIZE++;
				nownumber++;
			}
		}
	}
}

int main(int argc,char *argv[]){

    string fileF,fileS,Folder;
	double threshold,CUTOFF=0;
	
    for (int i=0;i<argc;i++) {
        std::string arg=argv[i];
        if (arg=="-meta"&&i+1<argc) {
            fileF=argv[i+1];
            i++;  // 跳过下一个参数，因为它是文件名
        } 
		else if(arg=="-single"&&i+1<argc) {
            fileS=argv[i+1];
            i++;  // 跳过下一个参数，因为它是文件名
        }
        else if(arg=="-h"||arg=="-help"){
        	cerr<<"Usage: "<<argv[0]<<" -meta <file1> -single <file2> -threshold <value> -folder <filepath> [-maxiter <value>]"<<endl;
        	cerr<<"-meta is contigs for metagenome file\n-single is for contig file derive from single cell sequencing\n";
        	cerr<<"-threshold should be input into degree within [0,90],not arc!\n";
        	cerr<<"-folder is for the mean and deviation file for contigs\n";
        	cerr<<"-maxiter is optional for deciding the number of iteration, the default would be 1\n";
		}
		else if(arg=="-threshold"&&i+1<argc) {
            threshold=stod(argv[i+1]);
            i++;  // 跳过下一个参数，因为它是文件名
        }
        else if(arg=="-maxiter"&&i+1<argc){
        	MAXITERATION=stoi(argv[i+1]);
        	i++;
		}
		else if(arg=="-folder"&&i+1<argc){
        	Folder=argv[i+1];
        	i++;
		}
    }
	if (argc<8) {
        cerr<<"Usage: "<<argv[0]<<" -meta <file1> -single <file2> -threshold <value> -folder <filepath> [-maxiter <value>]"<<endl;
        cerr<<"For detail, you can use -h parameter"<<endl;
        return 1;
    }
    // 执行基于参数的操作
    //std::cout << "Source file: " << sourceFile << std::endl;
    preprocess();
	cerr<<"preprocess completed"<<'\n';
    std::cout << "MetaFile specified after -f: " << fileF << std::endl;
    std::cout << "SingleFile specified after -s: " << fileS << std::endl;
    std::cout << "Folder specified after -folder: " << Folder << std::endl;
	Meta=Read_sequence_func(fileF);
	SCS=Read_sequence_func(fileS);
	cout<<"Metasize:"<<Meta.size()<<'\n';
	
//	READSTAFILE(Folder,Meta.size());
//	for(int i=0;i<256;i++)cerr<<Meta[1].TNF[i];
	for(int i=0;i<SCS.size();i++){
		buildgraph(0,SCS[i].number,0);//连接到超级点 
		buildgraph(SCS[i].number,0,0);
		SCS[i].category=1;
	}
	for(int i=0;i<Meta.size();i++){
		if(i%100==0)cout<<i<<'\n';
		for(int j=0;j<SCS.size();j++){
			double cosinee=vector_compare(Meta[i],SCS[j]);
			double angle=acos(cosinee);
//			cout<<"cos -> "<<cosinee<<'\n';
//			cout<<"angle -> "<<rad2deg(angle)<<'\n';
			if(rad2deg(angle)>threshold)continue;
			buildgraph(Meta[i].number,0,sin(angle));//连接单细胞和宏基因组像的点 
			buildgraph(0,Meta[i].number,sin(angle));
			break;
			Meta[i].category=1;
		}
	}
	for(int i=0;i<Meta.size();i++){
		if(i%100==0)cout<<i<<'\n';
		for(int j=i+1;j<Meta.size();j++){
			double cosinee=vector_compare(Meta[i],Meta[j]);
			double angle=acos(cosinee);
//			cout<<"cos -> "<<cosinee<<'\n';
//			cout<<"angle -> "<<rad2deg(angle)<<'\n';
			if(rad2deg(angle)>threshold)continue;
			buildgraph(Meta[i].number,Meta[j].number,sin(angle));//连接单细胞和宏基因组像的点 
			buildgraph(Meta[j].number,Meta[i].number,sin(angle));
		}
	}
//	cout<<"test->"<<abundancedistance(0,1,1,1.001)<<'\n';
	cerr<<"SAMPLESIZE -> "<<SAMPLE_SIZE<<'\n';
	for(int i=0;i<Meta.size();i++){
		for(int j=i+1;j<Meta.size();j++){
			double ADPij=ADP(Meta[i].number,Meta[j].number,CUTOFF);
//			cerr<<"ADP of "<<Meta[i].number<<" and "<<Meta[j].number<<" is: "<<ADPij<<'\n';
		}
	}
	LPA();
    return 0;
}

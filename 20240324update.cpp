#include<bits/stdc++.h>
#include<Windows.h>
#include<random>
#include<dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#define PI 3.14159265
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Dense>
#define tiaoshi  puts("what?")
using namespace std;
int tot=0;
map<string,int> contig2number,TNF2int;

//锟斤拷锟?
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

double cosconstant[181];
void preprocess(){
	for(int i=0;i<136;i++)
		TNF2int[TN[i]]=TNF2int[TN_reverse_complement[i]]=i;
	const double pi=acos(-1);
	for(int i=0;i<=180;i++)
		cosconstant[i]=cos(i*(pi/180.0));
}

template<class T, class S>
T random_unique(T begin, T end, S num_random) {
    S left = std::distance(begin, end);
    while (num_random--) {
        T r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

//Palindromic sequences
static const std::string TNP[] = { "ACGT", "AGCT", "TCGA", "TGCA", "CATG", "CTAG", "GATC", "GTAC", "ATAT", "TATA", "CGCG", "GCGC", "AATT", "TTAA", "CCGG", "GGCC" };

struct Sequence{
	string SEQ;
	vector<int> seq;
	double TNF[136];
	int number;
	string Realname;
	double sigma,average,length=0,logsize;
	int category;
	void initial(){memset(TNF,0,sizeof(TNF));}
	void transfer(){
		if(SEQ.size()>1000000){
			cerr<<"oversize!";
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
			length+=TNF[i]*TNF[i]*1.0;
		length=sqrt(length);
		for(int i=0;i<136;i++)
			TNF[i]/=length;
	}
	Sequence(string h,string realname){
		this->SEQ=h;
		this->initial();
		//this->transfer();
		this->statistic();
		this->number=++tot;
		this->Realname=realname;
		this->category=0;
		this->logsize=log10(SEQ.size());
		contig2number[realname]=tot;
//		cerr<<"struct realname: "<<realname<<" SEQ: "<<h<<'\n';
	}
};

int average_length(vector<Sequence> DNA){
	vector<int> LENGTH;
	for(auto it=DNA.begin(); it!=DNA.end(); ++it)
		LENGTH.push_back(it->SEQ.size());
	sort(LENGTH.begin(),LENGTH.end());
	return LENGTH[ceil(LENGTH.size()/4)];
}

void L2normal(array<double,136> &TNF){
	double length=0;
	for(int i=0;i<136;i++)
		length+=TNF[i]*TNF[i]*1.0;
	length=sqrt(length);
	for(int i=0;i<136;i++)
		TNF[i]/=length;
}

vector<vector<array<double,136>>> seq2seg(Sequence SINGLE,int averagelength, int slices=2){
	//骞惰鍖?浠ュ悗鍐?
	//姣忎釜鐗囨閮芥瀯閫犱竴涓悜閲忥紝涓€涓粍閲岄潰鏈夎嫢骞蹭釜鐗囨锛屼竴鏉″簭鍒楀悓鏃惰繕鏈塻lices涓粍
	vector<vector<array<double,136>>> SEG;//鍏堟瀯閫犱竴涓璞?
	for(int i=1;i<=slices;i++)SEG.push_back(vector<array<double,136>>());//鍏堝紑杈熺┖闂达紝鏂逛究鍚庨潰鐩存帴浣跨敤
	string h=SINGLE.SEQ;
	int len=h.size(),offset=averagelength/slices;
	for(int i=0;i<=len-4;i++){//寰幆姣忎釜鍥涙牳鑻烽吀搴忓垪
		int number=TNF2int[h.substr(i,4)];
		for(int j=0;j<slices;j++){
			int pos=i-j*offset;
			if(pos<0)break;
			pos=pos/averagelength;
			if(pos+1>SEG[j].size())SEG[j].push_back(array<double,136>{});
			SEG[j][pos][number]++;
		}
	}
	// for(int j=0;j<slices;j++){
	// 	for(auto it=SEG[j].begin();it!=SEG[j].end();++it){
	// 		L2normal(*it);
	// 	}
	// }
	return SEG;
}

int AVG_LENGTH,COVER_NUMBER;

typedef vector<vector<array<double,136>>> GROUP;

struct SEGMENT_TREE{
	#define lc (p<<1)
	#define rc (p<<1|1)
	vector<GROUP> TREE_SET;

	struct segment_tree{
		struct node{
			int l(0),r(0);
			array<double,136>TNFfrequency;
		};

		vector<node>t;
		vector<array<double,136>> slice;
		segment_tree(int n, const vector<array<double, 136>>& orig_slice): slice(orig_slice) { // 閫氳繃 move 鏉ユ帴鏀?orig_slice 鐨勫壇鏈?
        	t.resize(this->slice.size()*4); // 娉ㄦ剰锛氳繖閲屽簲璇ヤ娇鐢?this->slice.size()
        	build(1, 1, n);
			this->slice.clear();
    	}

		void pushup(int p){
    		for (int i = 0; i < 136; ++i) {
        		t[p].TNFfrequency[i] = t[lc].TNFfrequency[i] + t[rc].TNFfrequency[i];
			}
		}

		void build(int p,int l, int r){
			t[p].l=l,t[p].r=r;
			if(l==r){
				t[p].TNFfrequency=this->slice[l-1];
				L2normal(t[p].TNFfrequency);
				return;
			}
			int mid=l+r>>1;
			build(lc,l,mid);
			build(rc,mid+1,r);
			pushup(p);
			L2normal(t[p].TNFfrequency);
		}

		double vector_compare_segment(Sequence A, array<double,136> nowarr){
			double ans=0;
			for(int i=0;i<136;i++)
				ans+=A.TNF[i]*nowarr[i];
			return ans;
		}

		pair<int,double> query(Sequence A){
			int p=1;
			int estimate=A.SEQ.size()/AVG_LENGTH;
			double mincos=2;
			while(t[p].r-t[p].l+1>=estimate){
				if(t[p].l==0||t[p].r==0)break;
				if(t[p].l==t[p].r){
					
				}
				double left=vector_compare_segment(A,t[lc].TNFfrequency);
				double right=vector_compare_segment(A,t[rc].TNFfrequency);
				left>right?p=lc:p=rc;
				mincos=min(min(mincos,left),right);
			}
			return make_pair(p,mincos);
		}
	};

	SEGMENT_TREE(vector<Sequence> DNA){
		for(auto it=DNA.begin();it!=DNA.end();++it){
			TREE_SET.push_back(seq2seg(*it,AVG_LENGTH,COVER_NUMBER));
		}
		build_tree();
	}	
};

boost::numeric::ublas::matrix<double> matrixlization(vector<Sequence> T){
	boost::numeric::ublas::matrix<double> mat(T.size(),136);
	int i=0;
	for(auto it = T.begin(); it != T.end(); ++it, ++i){
		copy(it->TNF,it->TNF+136,&mat(i,0));
	}
	return mat;
}

boost::numeric::ublas::matrix<double> cosangle(boost::numeric::ublas::matrix<double> A, boost::numeric::ublas::matrix<double> B){
	auto T = boost::numeric::ublas::trans(B);
	return boost::numeric::ublas::prod(A,T);
}

double vector_compare(Sequence A,Sequence B,string type="cos"){
	double ans=0;
//	cerr<<"--------------------------------\n";
//	cerr<<"SEQ1->";
//	for(int i=0;i<256;i++)cerr<<A.TNF[i];
//	cerr<<"\nSEQ2->";
//	for(int i=0;i<256;i++)cerr<<B.TNF[i];
	for(int i=0;i<136;i++)
		ans+=A.TNF[i]*B.TNF[i];
//	cerr<<"ans -> "<<ans<<'\n'<<"length: "<<length1<<" "<<length2<<'\n';
	double cosineangle=ans;
	//*1.0/(A.length*B.length);
//	cerr<<"--------------------------------\n";
	if(type=="cos")return cosineangle;
	else return sqrt(2-2*cosineangle);
}

double TNF_DISTANCE(Sequence A,Sequence B){
	
	const double floor_prob=0.1;
	const double floor_preProb=log((1.0/floor_prob)-1.0);
	double LOG1=A.logsize,LOG2=B.logsize;
	double b, c; //parameters
	double d = vector_compare(A,B,"dis");
	double lw11 = min(LOG1,LOG2);
	double lw21 = max(LOG1,LOG2);
	double lw12 = lw11 * lw11;
	double lw13 = lw12 * lw11;
	double lw14 = lw13 * lw11;
	double lw15 = lw14 * lw11;
	double lw16 = lw15 * lw11;
	double lw17 = lw16 * lw11;
	double lw22 = lw21 * lw21;
	double lw23 = lw22 * lw21;
	double lw24 = lw23 * lw21;
	double lw25 = lw24 * lw21;
	double lw26 = lw25 * lw21;

	double prob;

	b = 46349.1624324381 + -76092.3748553155 * lw11 + -639.918334183 * lw21 + 53873.3933743949 * lw12 + -156.6547554844 * lw22 + -21263.6010657275 * lw13
			+ 64.7719132839 * lw23 + 5003.2646455284 * lw14 + -8.5014386744 * lw24 + -700.5825500292 * lw15 + 0.3968284526 * lw25 + 54.037542743 * lw16
			+ -1.7713972342 * lw17 + 474.0850141891 * lw11 * lw21 + -23.966597785 * lw12 * lw22 + 0.7800219061 * lw13 * lw23 + -0.0138723693 * lw14 * lw24
			+ 0.0001027543 * lw15 * lw25;
	c = -443565.465710869 + 718862.10804858 * lw11 + 5114.1630934534 * lw21 + -501588.206183097 * lw12 + 784.4442123743 * lw22 + 194712.394138513 * lw13
			+ -377.9645994741 * lw23 + -45088.7863182741 * lw14 + 50.5960513287 * lw24 + 6220.3310639927 * lw15 + -2.3670776453 * lw25 + -473.269785487 * lw16
			+ 15.3213264134 * lw17 + -3282.8510348085 * lw11 * lw21 + 164.0438603974 * lw12 * lw22 + -5.2778800755 * lw13 * lw23 + 0.0929379305 * lw14 * lw24
			+ -0.0006826817 * lw15 * lw25;

	//logistic model
	// prob = 1.0 / (1 + exp(-(b + c * d)));
	// if (prob >= .1)  //second logistic model
	double preProb = -(b + c * d); 
        // preProb <= LOG(9.0) yields prob > 0.1, so use second logistic model
	prob = preProb <= floor_preProb ? floor_prob: 1.0 / (1 + exp(preProb) );

	if (prob >= floor_prob) { //second logistic model
		b = 6770.9351457442 + -5933.7589419767 * lw11 + -2976.2879986855 * lw21 + 3279.7524685865 * lw12 + 1602.7544794819 * lw22 + -967.2906583423 * lw13
				+ -462.0149190219 * lw23 + 159.8317289682 * lw14 + 74.4884405822 * lw24 + -14.0267151808 * lw15 + -6.3644917671 * lw25 + 0.5108811613 * lw16
				+ 0.2252455343 * lw26 + 0.965040193 * lw12 * lw22 + -0.0546309127 * lw13 * lw23 + 0.0012917084 * lw14 * lw24 + -1.14383e-05 * lw15 * lw25;
		c = 39406.5712626297 + -77863.1741143294 * lw11 + 9586.8761567725 * lw21 + 55360.1701572325 * lw12 + -5825.2491611377 * lw22 + -21887.8400068324 * lw13
				+ 1751.6803621934 * lw23 + 5158.3764225203 * lw14 + -290.1765894829 * lw24 + -724.0348081819 * lw15 + 25.364646181 * lw25 + 56.0522105105 * lw16
				+ -0.9172073892 * lw26 + -1.8470088417 * lw17 + 449.4660736502 * lw11 * lw21 + -24.4141920625 * lw12 * lw22 + 0.8465834103 * lw13 * lw23
				+ -0.0158943762 * lw14 * lw24 + 0.0001235384 * lw15 * lw25;
		//prob = 1.0 / (1 + exp(-(b + c * d)));
		// prob = prob < .1 ? .1 : prob;
		preProb = -(b + c * d); // exp(preProb) <= 9 yields prob >= 0.1, so preProb <= LOG(9.0) to calculate, otherwise use the floor
		prob = preProb <= floor_preProb ? 1.0 / (1 + exp(preProb)) : floor_prob;
	}

	return prob;
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
//	cerr<<u<<" link with "<<v<<'\n';
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
    	cerr<<"-------------------\n";
    	cerr<<"For  iteration"<<iter<<" , the tags are: \n";
    	for (int i = 0; i <=tot; ++i) {
    	    cerr<<labels[i]<<" ";
    	}
    	cerr<<"\n\n";
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
    std::cerr << "Final labels after label propagation:" << std::endl;
    for (int i = 0; i<=tot; ++i) {
        std::cerr << "Node " << i << ": Label " << labels[i] << std::endl;
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
	//cerr<<"k1->"<<k1<<" k2->"<<k2<<'\n';
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
	//	cerr<<"mean_in->"<<mean_in<<"  sigma_in->"<<sigma_in<<"  mean_jn->"<<mean_jn<<"  sigma_jn->"<<sigma_jn<<'\n';
		int signal=(mean_in>c)|(mean_jn>c);
		flag+=signal;
		if(!signal)continue;
		P_ij*=abundancedistance(mean_in,sigma_in,mean_jn,sigma_jn);
		//cerr<<"AD->"<<abundancedistance(mean_in,sigma_in,mean_jn,sigma_jn)<<'\n';
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
//			cerr<<"contig_name: "<<contig_name<<'\n';
//			cerr<<"refer: "<<refer<<'\n';
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
				cerr<<"filename:"<<entry->d_name<<endl;
				READ_PARTICULAR_STATISTIC_FILE(folderpath,entry->d_name,nownumber,size_meta);
				SAMPLE_SIZE++;
				nownumber++;
			}
		}
	}
}


int main(int argc,char *argv[]){

    static string fileF,fileS,Folder;
	static double threshold,CUTOFF=0;
	
    for (int i=0;i<argc;i++) {
        std::string arg=argv[i];
        if (arg=="-meta"&&i+1<argc) {
            fileF=argv[i+1];
            i++;  // 锟斤拷锟斤拷锟斤拷一锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷为锟斤拷锟斤拷锟侥硷拷锟斤拷
        } 
		else if(arg=="-single"&&i+1<argc) {
            fileS=argv[i+1];
            i++;  // 锟斤拷锟斤拷锟斤拷一锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷为锟斤拷锟斤拷锟侥硷拷锟斤拷
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
            i++;  // 锟斤拷锟斤拷锟斤拷一锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷为锟斤拷锟斤拷锟侥硷拷锟斤拷
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
    
    // 执锟叫伙拷锟节诧拷锟斤拷锟侥诧拷锟斤拷
    //std::cerr << "Source file: " << sourceFile << std::endl;
    preprocess();
	cerr<<"preprocess completed"<<'\n';
    std::cerr << "MetaFile specified after -f: " << fileF << std::endl;
    std::cerr << "SingleFile specified after -s: " << fileS << std::endl;
    std::cerr << "Folder specified after -folder: " << Folder << std::endl;
	Meta=Read_sequence_func(fileF);
	SCS=Read_sequence_func(fileS);
	cerr<<"Metasize:"<<Meta.size()<<'\n';
	
	boost::numeric::ublas::matrix<double> MatrixMeta=matrixlization(Meta);
	boost::numeric::ublas::matrix<double> MatrixSCS=matrixlization(SCS);
	boost::numeric::ublas::matrix<double> MatrixCosMetaMeta=cosangle(MatrixMeta,MatrixMeta);
	boost::numeric::ublas::matrix<double> MatrixCosMetaSingle=cosangle(MatrixMeta,MatrixSCS);
	
	cerr<<MatrixCosMetaMeta<<'\n';
	
//	READSTAFILE(Folder,Meta.size());
//	for(int i=0;i<256;i++)cerr<<Meta[1].TNF[i];
	for(int i=0;i<SCS.size();i++){
		buildgraph(0,SCS[i].number,0);//锟斤拷锟接碉拷锟斤拷锟斤拷锟斤拷 
		buildgraph(SCS[i].number,0,0);
		SCS[i].category=1;
	}
	for(int i=0;i<Meta.size();i++){
		if(i%100==0)cerr<<i<<'\n';
		for(int j=0;j<SCS.size();j++){
			double cosinee=vector_compare(Meta[i],SCS[j]);
			double angle=acos(cosinee);
//			cerr<<"cos -> "<<cosinee<<'\n';
//			cerr<<"angle -> "<<rad2deg(angle)<<'\n';
			if(rad2deg(angle)>threshold)continue;
			buildgraph(Meta[i].number,0,sin(angle));//锟斤拷锟接碉拷细锟斤拷锟酵猴拷锟斤拷锟斤拷锟斤拷锟侥碉拷 
			buildgraph(0,Meta[i].number,sin(angle));
			break;
			Meta[i].category=1;
		}
	}
	for(int i=0;i<Meta.size();i++){
		if(i%100==0)cerr<<i<<'\n';
		for(int j=i+1;j<Meta.size();j++){
			double cosinee=vector_compare(Meta[i],Meta[j]);
			double angle=acos(cosinee);
//			cerr<<"cos -> "<<cosinee<<'\n';
//			cerr<<"angle -> "<<rad2deg(angle)<<'\n';
			if(rad2deg(angle)>threshold)continue;
			buildgraph(Meta[i].number,Meta[j].number,sin(angle));//锟斤拷锟接碉拷细锟斤拷锟酵猴拷锟斤拷锟斤拷锟斤拷锟侥碉拷 
			buildgraph(Meta[j].number,Meta[i].number,sin(angle));
		}
	}
//	cerr<<"test->"<<abundancedistance(0,1,1,1.001)<<'\n';
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

//-meta C:\Users\okura\Downloads\euk_test\euk_test\GCF_000227135.1_ASM22713v2_genomic.fna -single C:\Users\okura\Downloads\testdata\random_sequences.fasta -threshold 30 -folder C:\Users\okura\Downloads\testdata\statis

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
#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Dense>
#define tiaoshi  puts("what?")
using namespace std;
int tot=0;
int maxEdges=200;
double maxP=0.9;
double cutoff=0.90;
int debug=1;
int COVER_NUMBER=2;
double POLL_PARAMETER=0.75;
int MAXITERATION = 3;

map<string,int> contig2number,TNF2int;
vector<int>FILE_NUMBER;

#define N 5000000
//struct node{
//	int u,v;
//	double w;
//}e[N];
//int nxt[N],first[N],cnt,MAXITERATION=3,CNT;
//
//void buildgraph(int u,int v,int w){
////	cerr<<u<<" link with "<<v<<'\n';
//	e[++cnt].u=u;
//	e[cnt].v=v;
//	e[cnt].w=w;
////	nxt[cnt]=first[u];
////	first[u]=cnt;	
//}
//
//void BUILDGRAPH(int u,int v,int w){
//	E[++CNT].u=u;
//	E[CNT].v=v;
//	E[CNT].w=w;
//	nxt[CNT]=first[u];
//	first[u]=CNT;
//}



struct node {
    int u, v;
    double w;
};

vector<node> e,E; // 使用vector代替e数组
vector<int> nxt, First;
int cnt = 0, CNT = 0; // cnt和CNT初始化为0，因为vector的下标从0开始

void buildgraph(int u, int v, double w) {
    // 添加新的边信息
    e.push_back({u, v, w});
    
    if(debug){
    	cerr<<"\nu->"<<u<<" v->"<<v<<" w->"<<w<<'\n';
	}
    
}

void buildGraphWithAdjacencyList(int u, int v, double w) {
    // 对于邻接表的构建，同时需要更新E, nxt, 和first
    E.push_back({u, v, w});
    CNT = E.size() - 1; // 更新CNT为最新的边的索引
    nxt.push_back(First[u]); // nxt[CNT] = first[u];
    First[u] = CNT; // 更新first[u]为最新的边的索引
    
    if(debug){
    	cerr<<"\nU->"<<u<<" V->"<<v<<" W->"<<w<<'\n';
	}
    
}

vector<int> computeRanks(const vector<node>& arr) {
    vector<int> indices(arr.size());
    iota(indices.begin(), indices.end(), 0); // initial ranking array

    sort(indices.begin(), indices.end(), [&](int i, int j) { return arr[i].w > arr[j].w; });

    vector<int> ranks(arr.size());
    for (int i = 0; i < indices.size(); ++i) {
        ranks[indices[i]] = i + 1; // 直接赋值排名
    }

    return ranks;
}
//锟斤拷锟?
static const std::string TN[] = { 
		"GGTA", "AGCC", "AAAA", "ACAT", "AGTC", "ACGA", "CATA", "CGAA", "AAGT", "CAAA", "CCAG", "GGAC", "ATTA", "GATC", "CCTC",
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
	bool SINGLE_OR_NOT;
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
	Sequence(string h,string realname,int single,int category=0){
		this->SEQ=h;
		this->initial();
		//this->transfer();
		this->statistic();
		this->number=++tot;
		this->Realname=realname;
		this->category=category;
		this->logsize=log10(SEQ.size());
		this->SINGLE_OR_NOT=single;
		contig2number[realname]=tot;
//		cerr<<"struct realname: "<<realname<<" SEQ: "<<h<<'\n';
	}
};

int average_length(const vector<Sequence> &DNA){
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
		
//		if(debug&&averagelength==6&&slices==2){
//			cout<<number<<" ";
//		}
		
		for(int j=0;j<slices;j++){
			int pos=i-j*offset;
			if(pos<0)break;
			pos=pos/averagelength;
			if(pos+1>SEG[j].size())SEG[j].push_back(array<double,136>{});
			SEG[j][pos][number]++;
		}
	}
	
	if(debug&&averagelength==6&&slices==1){
		for(int i=0;i<SEG[0].size();i++){
			for(int j=0;j<136;++j){
//				cout<<SEG[1][i][j]<<" ";
				if(SEG[0][i][j]){
					cout<<TN[j]<<":"<<SEG[0][i][j]<<' ';
				}
			}
			cout<<'\n';
		}
	}
	
	// for(int j=0;j<slices;j++){
	// 	for(auto it=SEG[j].begin();it!=SEG[j].end();++it){
	// 		L2normal(*it);
	// 	}
	// }
	return SEG;
}

int AVG_LENGTH;

typedef vector<vector<array<double,136>>> GROUP;

struct SEGMENT_TREE{
	#define lc (p<<1)
	#define rc (p<<1|1)

	vector<GROUP> SEQ_SET;
	int start,category;

	struct segment_tree{
		struct Node{
			int l=0,r=0;
			array<double,136>TNFfrequency;
		};

		vector<Node>t;
		vector<array<double,136>> slice;
		int start,category;

		segment_tree(int n, const vector<array<double, 136>>& orig_slice): slice(orig_slice) { // 閫氳繃 move 鏉ユ帴鏀?orig_slice 鐨勫壇鏈?
        	t.resize(this->slice.size()*4); // 娉ㄦ剰锛氳繖閲屽簲璇ヤ娇鐢?this->slice.size()
        	build(1, 1, n);
        	
//			if(debug){
//        		puts("attention"); 
//			}
			 
			L2normal(t[1].TNFfrequency);
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
				
//				if(debug){
//					cout<<"p->"<<p<<" l->"<<l<<" r->"<<r<<'\n';
//					for(int i=0;i<136;i++)
//						if(t[p].TNFfrequency[i])cout<<TN[i]<<":"<<t[p].TNFfrequency[i]<<" ";
//					cout<<'\n';
//				}
				
				return;
			}
			int mid=l+r>>1;
			build(lc,l,mid);
			build(rc,mid+1,r);
			pushup(p);
			L2normal(t[lc].TNFfrequency);
			L2normal(t[rc].TNFfrequency);//miss
			
//			if(debug){
//				cout<<"\n\n\np->"<<p<<" l->"<<l<<" r->"<<r<<" lc->"<<lc<<" rc->"<<rc<<'\n';
//				cout<<"p:";
//				for(int i=0;i<136;i++)
//					if(t[p].TNFfrequency[i])cout<<TN[i]<<"->"<<t[p].TNFfrequency[i]<<" ";
//				cout<<'\n';
//				cout<<"lc:";
//				for(int i=0;i<136;i++)
//					if(t[lc].TNFfrequency[i])cout<<TN[i]<<"->"<<t[lc].TNFfrequency[i]<<" ";
//				cout<<'\n';
//				cout<<"rc:";
//				for(int i=0;i<136;i++)
//					if(t[rc].TNFfrequency[i])cout<<TN[i]<<"->"<<t[rc].TNFfrequency[i]<<" ";				
//				cout<<'\n';
//			}
			
		}

		double vector_compare_segment(const Sequence &A,const array<double,136> &nowarr){
			double ans=0;
			for(int i=0;i<136;i++)
				ans+=A.TNF[i]*nowarr[i];
			
//			if(debug)cout<<"ans->"<<ans<<'\n';
			
			return ans;
		}

		pair<int,double> query(Sequence A){
			int p=1;
			int estimate=A.SEQ.size()/AVG_LENGTH;
			double mincos=vector_compare_segment(A,t[1].TNFfrequency);
			while(t[p].r-t[p].l+1>=estimate){
				
//				if(debug){
//					cout<<p<<'\n';
//				}
				
				if(t[p].l==0||t[p].r==0)break;
				if(t[p].l==t[p].r){
					break;
				}
				double left=vector_compare_segment(A,t[lc].TNFfrequency);
				double right=vector_compare_segment(A,t[rc].TNFfrequency);
				left>right?p=lc:p=rc;
				mincos=max(max(mincos,left),right);
			}
			return make_pair(t[p].l,mincos);
		}
	};

	typedef vector<segment_tree> ONE_SEQ_TREE_SET;
	typedef vector<ONE_SEQ_TREE_SET> ALL_SEQ_TREE_SET;

	ALL_SEQ_TREE_SET SINGLE_SET;

	SEGMENT_TREE(vector<Sequence> DNA){
		start=DNA[0].number;
		category=DNA[0].category;
		for(auto it=DNA.begin();it!=DNA.end();++it){
			GROUP group=seq2seg(*it,AVG_LENGTH,COVER_NUMBER);
			ONE_SEQ_TREE_SET TEMP1;
			for(auto itt=group.begin();itt!=group.end();++itt){
				segment_tree T(itt->size(),*itt);
				TEMP1.push_back(T);
			}
			SEQ_SET.push_back(group);
			SINGLE_SET.push_back(TEMP1);
		}
	}

	void detect_relation(vector<Sequence> META, double cutoff){
		int number_of_seq=SINGLE_SET.size();
		for(auto it1=META.begin();it1!=META.end();it1++){
			int flag=0,oanji=0,high_credential=0,i=0;
			for(auto it2=SINGLE_SET.begin();it2!=SINGLE_SET.end();++it2,++i){//寰幆鍗曠粏鑳為噷姣忎竴鏉″簭鍒?
				double similarity=0;bool find_or_not=false;
				for(auto it3=it2->begin();it3!=it2->end();++it3){//姣忎竴鏉″簭鍒楃殑姣忎釜鍒嗙粍
					pair<int,double> where_and_similarity=it3->query(*it1);//鏌ヨ瀹忓熀鍥犵粍鐨勬煇鏉ontig鍦ㄥ垎缁勬爲閲岄潰鐨勬煡璇㈢粨鏋?
					similarity=max(similarity,where_and_similarity.second);
					if(where_and_similarity.second>cutoff){//濡傛灉宸茬粡澶т簬浜嗭紝灏变笉鐢ㄧ户缁煡涓嬮潰鐨勫垎缁勪簡锛岃瘉鏄庡凡缁忔湁鐩镐技鎬т簡
						find_or_not=true;
						flag++;//瀵规瘡鏉″簭鍒楄嚦澶氬彧浼氬姞涓€娆?
						break;
					}
				}
				if(flag>=number_of_seq*POLL_PARAMETER){//濡傛灉鎶曠エ鏁板凡缁忓崰澶氭暟浜嗭紝鍙互鑰冭檻杩為珮淇＄敤鍊艰竟
					high_credential=1;
				}
				if(flag>=number_of_seq*1.0*4/3){//濡傛灉鎶曠エ鏁板凡缁忕壒鍒珮浜嗭紝鐩存帴鍜屽崟缁嗚優褰掍负涓€绫伙紝涓嶅啀鍙桳PA鐨勫奖鍝?
					it1->category=category;
					it1->SINGLE_OR_NOT=1;
					oanji=1;
					break;
				}
				
//				if(debug)cerr<<i;
				
				if(find_or_not){
					if(high_credential){
						
//						if(debug){
//							cerr<<"\nstart->"<<start<<" i->"<<i<<'\n';
//						} 
						
						similarity=max(similarity,0.9);
						buildgraph(start+i,it1->number,similarity);
						buildgraph(it1->number,start+i,similarity);
					}
					else{
						
//						if(debug){
//							cerr<<"\nstart->"<<start<<" i->"<<i<<'\n';
//						} 
						
						buildgraph(start+i,it1->number,similarity);
						buildgraph(it1->number,start+i,similarity);
					}
				}
			}
		}
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
	
	
	/* to use this, the following code show the way, but whether it is time-consuming needs some practices*/
	
	
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

vector<Sequence> Read_sequence_func(string filepath, int type=0, int category=0){
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
        	
        	Sequence TEMP(fullline,name,type,category);
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
    
    Sequence TEMP(fullline,name,type,category);
    DNA_DATA.push_back(TEMP);
    number_of_sq++;
    fullline="";
    
    cerr<<"total sequence number is: "<<number_of_sq<<'\n';
    inputFile.close();
    return DNA_DATA;
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
            for(int o=First[i];o;o=nxt[o])
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

vector<string> GetFilesInDirectory(const string& directoryPath) {
    vector<string> files;
    DIR* dir = opendir(directoryPath.c_str()); // 使用C字符串打开目录

    if (dir != nullptr) {
        dirent* entity;
        while ((entity = readdir(dir)) != nullptr) { // 读取目录项
            // 构造完整的文件路径
            string fullPath = directoryPath + "\\" + string(entity->d_name);

            // 使用stat获取文件信息
            struct stat path_stat;
            stat(fullPath.c_str(), &path_stat);

            // 检查是否为普通文件
            if (S_ISREG(path_stat.st_mode)) {
                files.push_back(fullPath);
            }
        }
        closedir(dir); // 关闭目录句柄
    } else {
        cerr << "Cannot open directory: " << directoryPath << endl;
    }

    return files;
}

typedef  boost::numeric::ublas::matrix<double> MATRIX_DOUBLE;

int gen_tnf_graph_sample(double coverage, bool full, int nobs, const MATRIX_DOUBLE& matrix){
	
	std::vector<unsigned char> connected_nodes;
    connected_nodes.resize(nobs);
	int p = 999, pp = 1000;
	double cov = 0, pcov = 0;
    int round = 0;

	for (; p > 500; ) {
        round++;

		double cutoff = (double) p / 1000.;

		#pragma omp parallel for
		for (int i = 0; i < nobs; ++i) {
			int kk = nobs;
			for (int j = i + 1; j < kk; ++j) {
				double s = matrix(i,j);
				if (s >= cutoff) {
					connected_nodes[i] = 1;
					connected_nodes[j] = 1;
				}
			}
		}

		//cov = (double) connected_nodes.size() / _nobs;
		int counton = 0;
		#pragma omp parallel for reduction(+:counton)
                for (int i = 0; i < nobs; i++) {
                    if (connected_nodes[i] == 1) counton++;
                }
		cov = (double) counton / nobs;

		if (cov >= coverage) {
			//previous cov is closer to coverage then choose prev p instead current p
			if (cov - coverage > coverage - pcov) {
				p = pp;
				cov = pcov;
			}

			break;
		} else
		//	verbose_message("Preparing TNF Graph Building [pTNF = %2.1f; %d / %d (P = %2.2f%%) round %d]               \r", (double) p / 10., counton, nobs, cov * 100, round);

		pp = p; pcov = cov;

		if (p > 990) //99.9, 99.6, 99.3, 99.0
			p -= rand() % 3 + 1; //choose from 1,2,3
		else if (p > 900) //98.5, 98, 97.5, ... 90.0
			p -= rand() % 3 + 3; //choose from 3,4,5
		else //89, 88, 87, ..., 70
			p -= rand() % 3 + 9; //choose from 9,10,11

	}
	
	return p;
}

struct CompareEdge {
	constexpr bool operator() (std::pair<int, double> const & a, std::pair<int, double> const & b) const noexcept {
        return a.second > b.second;
    }
};

#include "tile.h"
void gen_tnf_graph(double cutoff, int nobs, const MATRIX_DOUBLE& matrix) {
	
	int TILE = 10;
	try {
		TILE = std::max((int) ( ( CacheSize() * 1024. ) / ( 2* sizeof(float) * 136 + maxEdges * ( 2 * sizeof(int) + 1 * sizeof(double) ) ) ), (int) 10);
	} catch (...) {}

//#pragma omp parallel for schedule(dynamic, 1) proc_bind(spread) reduction(merge_int: from) reduction(merge_int: to) reduction(merge_double: sTNF) 
//#pragma omp parallel for schedule(dynamic, 1) reduction(merge_int: from) reduction(merge_int: to) reduction(merge_double: sTNF) 
    for ( int ii= 0; ii < nobs; ii+=TILE ){
        std::vector<std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, CompareEdge> > edges(TILE);

        for ( int jj= 0; jj < nobs; jj+=TILE ) {
            for (int i = ii; i < ii+TILE && i<nobs; ++i){
                int que_index = i-ii;
                for (int j = jj; j < jj+TILE && j < nobs; ++j) {
                    if (i == j)
                        continue;
                    double sTNF = matrix(i, j);
                    if (sTNF > cutoff &&
                        (edges[que_index].size() < maxEdges || (edges[que_index].size() == maxEdges && sTNF > edges[que_index].top().second))) {
                        if (edges[que_index].size() == maxEdges)
                            edges[que_index].pop();
                        edges[que_index].push(std::make_pair(j, sTNF));
                    }
                }
            }
        }
        for ( int k = 0; k < TILE; ++k) {
            while (!edges[k].empty()) {
                std::pair<int, double> edge = edges[k].top();
                if ( (ii+k) < edge.first) {
//                    sTNF.push_back(edge.second);
//                    from.push_back((ii+k));
//                    to.push_back(edge.first);
                    buildgraph(edge.first,ii+k,edge.second);
                    buildgraph(ii+k,edge.first,edge.second);
                }
                edges[k].pop();
            }
        }
//        if (verbose) {
//            progress.track(TILE);
//            if (omp_get_thread_num() == 0 && progress.isStepMarker()) {
//                verbose_message("Building TNF Graph %s [%.1fGb / %.1fGb]                           \r",
//                                progress.getProgress(), getUsedPhysMem(), getTotalPhysMem() / 1024 / 1024);
//            }
//        }
    }

}

int label_propagation(vector<int>& membership, vector<int>& node_order) {
//	int no_of_nodes = g.getNodeCount();
//	int no_of_edges = g.getEdgeCount();
//
//        if (no_of_nodes == 0 || no_of_edges == 0) {
//                cerr << "There were " << no_of_nodes << " nodes and " << no_of_edges << " edges -- skipping label_propagation" << endl;
//                return 0;
//        }

//	if (g.sSCR.size() != no_of_edges) {
//		cerr << "sSCR != no_of_edges" << endl;
//		exit(1);
//	}

//	if (membership.size() != no_of_nodes) {
//		membership.resize(no_of_nodes);
//		std::iota(membership.begin(), membership.end(), 0);
//	}

	/* Do some initial checks */
//	if (*std::min_element(g.sSCR.begin(), g.sSCR.end()) < 0) {
//		cerr << "sSCR must be non-negative" << endl;
//		exit(1);
//	}

	std::unordered_map<int, std::unordered_set<int> > visited;
	std::unordered_set<int> blacklist;

	int nLeftMin = 2147483647;
	int attempt = 0;
	bool running = true;
	while (running) {
		running = false;

		int nLeft = 0;

		/* In the prescribed order, loop over the vertices and reassign labels */
		for (int i = 0; i < node_order.size(); i++) { //we reconsider all nodes regardless of its previous status, but is it better?
			int v1 = node_order[i];

			std::unordered_map<int, double> neighbor_scores; //sum of neighbors scores to cluster k
			std::unordered_map<int, int> neighbor_counts; //keep number of neighbors

//						std::vector<int>& ineis = g.incs[v1];
//						for (int j = 0; j < ineis.size(); j++) { //# of neighbors (edges connected to v1)
//							int edgeID = ineis[j];
//			
//							int_fast32_t k = membership[g.getOtherNode(edgeID, v1)]; //community membership of a neighbor (connected by j)
//			
//							if (neighbor_scores.find(k) == neighbor_scores.end()) {
//								neighbor_scores[k] = 0.;
//								neighbor_counts[k] = 0;
//							}
//							neighbor_scores[k] += LOG(1. - g.sSCR[edgeID]); //as p-value
//							neighbor_counts[k]++;
//						}
						
						
			for(int j=First[v1];j;j=nxt[j]){
				int k=membership[E[j].v];
				if (neighbor_scores.find(k) == neighbor_scores.end()) {
					neighbor_scores[k] = 0.;
					neighbor_counts[k] = 0;
				}
				if(fabs(E[j].w-1)<0.001)E[j].w=0.999;
				neighbor_scores[k] += log(1. - E[j].w); //as p-value
				neighbor_counts[k]++;
			}
			
			
			if (neighbor_scores.size() > 0) {
				for (auto &kv : neighbor_scores) {
					
//					if(debug){
//						cerr<<"kv.first:"<<kv.first<<" neighbor_counts[kv.first]:"<<neighbor_counts[kv.first]<<" kv.second:"<<kv.second<<'\n';
//					} 
					
					//Fisher's method to compare significance of different number of probs.
					boost::math::chi_squared chi_sqr_dist(2 * neighbor_counts[kv.first]);
					kv.second = boost::math::cdf(chi_sqr_dist, -2.0 * kv.second);
				}
				auto best_neighbor = std::max_element (neighbor_scores.begin(), neighbor_scores.end(),
					[] (const std::pair<int, double>& p1, const std::pair<int, double>& p2) { return p1.second < p2.second; });

				//however, if there was a clique (loop) out of >2 nodes
				int kPrev = membership[v1];
				if (kPrev != (int) best_neighbor->first && blacklist.find(v1) == blacklist.end()) {

					membership[v1] = best_neighbor->first;

					int kNext = membership[v1];
					if (visited.find(v1) == visited.end() || visited[v1].find(kNext) == visited[v1].end()) {
						//not have been assigned to the cls before
						nLeft++; //# of confirmation (that this choice is optimal) left
						running = true;
					} else {
						blacklist.insert(v1); //blacklist represents nodes that change cls in a circular form
					}
					visited[v1].insert(kNext);
				}
			}
		}

		if (nLeft < nLeftMin) {
			nLeftMin = nLeft;
			attempt = 0;
		} else {
			attempt++;
			if (attempt >= 10) {
				break;
			}
		}
		//cout << "nLeft: " << nLeft << " & attempt: " << attempt << endl;
	}

	return 0;
}

struct NodeIndex {
    node edge;
    int index;
};

bool compare(NodeIndex a, NodeIndex b) {
    return a.edge.w > b.edge.w;
}

vector<int> getRankPositions(const vector<node>& a) {
    int n = a.size();
    vector<NodeIndex> nodeIndices(n);

    // 填充结构体索引数组
    for (int i = 0; i < n; ++i) {
        nodeIndices[i] = {a[i], i};
    }

    // 根据结构体的w值进行排序
    sort(nodeIndices.begin(), nodeIndices.end(), compare);

    // 创建结果数组，并使用排序后的索引填充
    vector<int> b(n);
    for (int i = 0; i < n; ++i) {
        b[i] = nodeIndices[i].index;
    }

    return b;
}

void clustering( const vector<node>& e, const vector<Sequence> &Meta){
	vector<int> arr=getRankPositions(e);
	vector<int> node_order,mem;
	vector<double> p_schedule2;
	unordered_set<int> connected_nodes;
	for (int i = 1; i <= 10; ++i)
		p_schedule2.push_back(maxP / 10 * i);
		
	int nobs=Meta.size();
	mem.resize(tot+1);
	iota(mem.begin(),mem.end(),0);
	for(int i=0;i<nobs;i++){
		if(Meta[i].category)
			mem[i+1]=Meta[i].category;
	}
	
	for(int i=nobs+1;i<=tot;i++){
		mem[i]=lower_bound(FILE_NUMBER.begin(),FILE_NUMBER.end(),i-nobs)-FILE_NUMBER.begin()+nobs;
	}
	
//	if(debug){
//		for(auto & clu:mem){
//			cerr<<clu<<" ";
//		}
//		puts("");
//	}
	
	int which_p=0,nEdges=e.size();//undirect
	First.resize(tot+1);//this is a bug
	
	if(debug)cout<<"nEdges->"<<nEdges<<'\n';
	
	for(int i=0;i<nEdges;i+=2){
		int ii=e[arr[i]].u,jj=e[arr[i]].v;
		
		if(debug){
			cerr<<"ii->"<<ii<<" jj->"<<jj<<" weight->"<<e[arr[i]].w<<'\n'; 
			cerr<<"mem[ii]->"<<mem[ii]<<" mem[jj]->"<<mem[jj]<<'\n'; 
		}
		
		if (mem[ii] != mem[jj]) { // || which_p < 5 allow all edges from first 5 schedule
					
			if (connected_nodes.find(ii) == connected_nodes.end()) {
				if(ii<=nobs)node_order.push_back(ii);
				connected_nodes.insert(ii);
			}
			if (connected_nodes.find(jj) == connected_nodes.end()) {
				if(jj<=nobs)node_order.push_back(jj);
				connected_nodes.insert(jj);
			}

//			Similarity scr = g.sSCR[oSCR[i]];
//			if (scr > minS) {
//				g2.sSCR.push_back(scr);
//				g2.from.push_back(jj);
//				g2.to.push_back(ii);
//				g2.incs[ii].push_back(g2.from.size() - 1);
//				g2.incs[jj].push_back(g2.from.size() - 1);
//			} else {
//				i = nEdges - 1; //early stopping
//			}
			buildGraphWithAdjacencyList(ii,jj,e[arr[i]].w);
			buildGraphWithAdjacencyList(jj,ii,e[arr[i]].w);
		}
		
		if (E.size() > 0 && ((double) connected_nodes.size() / Meta.size() >= p_schedule2[which_p] || i >= nEdges - 1) ){
			//cout << "g2.sSCR.back(): " << g2.sSCR.back() << endl;
			label_propagation(mem, node_order);
			cerr<<"LPALPALPALPALPALPALPA\n";
//					if (debug) {
//						std::string osfileName("cluster.log." + boost::lexical_cast<std::string>(which_p));
//						std::ofstream os(osfileName);
//						if (!os) {
//							cerr << "[Error!] Failed to write to " << osfileName << endl;
//							return 1;
//						}
//						ClassMap _cls;
//						for (int i = 0; i < nobs; ++i) {
//							_cls[mems[i]].push_back(i);
//						}
//						for (int kk = 0; kk < nobs; ++kk) {
//							if (_cls[kk].size() > 1) {
//								os << kk << " : ";
//								ContigVector& vec = _cls[kk];
//								std::sort(vec.begin(), vec.end());
//								for (auto it2 = vec.begin(); it2 != vec.end(); ++it2) {
//									os << *it2 << ",";
//								}
//								os << endl;
//							}
//						}
//						os.close();
//						if (!os) {
//							cerr << "[Error!] Failed to write to " << osfileName << endl;
//							return 1;
//						}
//					}

			if (++which_p == p_schedule2.size())
				break;
				
//			if(debug){
//				for(auto & clu:mem){
//					cerr<<clu<<" ";
//				}
//				puts("");
//			}
		}
		
		if(debug){
			for(auto & clu:mem){
				cerr<<clu<<" ";
			}
			puts("");
		}
	
	}
}


int main(int argc,char *argv[]){

    preprocess();
	cerr<<"preprocess completed"<<'\n';
	
//	Sequence Test("AAAAAACCCCCCGGGGGGTTTTTT","test",0,0);
//	seq2seg(Test,6,2);
//	
//	while(1);
//  seq2seg has been checked

//	Sequence Test("AAAAAACCCCCCGGGGGGTTTTTT","test",0,0);
//	auto groupp=seq2seg(Test,6,1);
//	AVG_LENGTH=6;
//	SEGMENT_TREE::segment_tree T((groupp.begin())->size(),*(groupp.begin()));
////	cout<<(groupp.begin())->size()<<'\n';
//	Sequence comparetest("AAAACGT","test2",0,0);
//	pair<int,double> result=T.query(comparetest);
//	puts("here");
//	cout<<result.first<<" "<<result.second;
//	
//	while(1);
	

    static string fileF,fileS,Folder;
	static double threshold,CUTOFF=0;
	
    for (int i=0;i<argc;i++) {
        std::string arg=argv[i];
        if (arg=="-meta"&&i+1<argc) {
            fileF=argv[i+1];
            i++;  
        } 
		else if(arg=="-single"&&i+1<argc) {
            fileS=argv[i+1];
            i++;  
        }
        else if(arg=="-h"||arg=="-help"){
        	cerr<<"Usage: "<<argv[0]<<" -meta <file1> -single <file2> -threshold <value> -folder <filepath> [-maxiter <value>]"<<endl;
        	cerr<<"-meta is contigs for metagenome file \n-single is for contig file derive from single cell sequencing\n";
        	cerr<<"-threshold should be input into degree within [0,90],not arc!\n";
        	cerr<<"-folder is for the mean and deviation file for contigs\n";
        	cerr<<"-maxiter is optional for deciding the number of iteration, the default would be 1\n";
		}
		else if(arg=="-threshold"&&i+1<argc) {
            threshold=stod(argv[i+1]);
            i++; 
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
    std::cerr << "MetaFile specified after -f: " << fileF << std::endl;
    std::cerr << "SingleFile specified after -s: " << fileS << std::endl;
    std::cerr << "Folder specified after -folder: " << Folder << std::endl;
	vector<Sequence> Meta=Read_sequence_func(fileF);
	vector<string> files = GetFilesInDirectory(fileS);
	vector<vector<Sequence>> SCS;
	vector<SEGMENT_TREE> SCS_SEGMENT_TREE;

	int ii=Meta.size();
	AVG_LENGTH=average_length(Meta);
	
	if(debug){
		cerr<<"\nthe average length is: "<<AVG_LENGTH<<'\n';
	}

	FILE_NUMBER.push_back(0);

	for(const auto& file:files){
		
		if(debug){
			cerr<<"one of the names of the single is: "<<file<<'\n';
		}
		
		ii++;
		vector<Sequence> temp=Read_sequence_func(file,1,ii);
		
//		if(debug){
//			for(auto it=temp.begin();it!=temp.end();++it){
//				cerr<<"this is the sequence of single cell sequence: "<<it->SEQ<<", and type is:"<<it->SINGLE_OR_NOT<<'\n';
//			}
//		}

		FILE_NUMBER.push_back(FILE_NUMBER[FILE_NUMBER.size()-1]+temp.size());
		
		SCS.push_back(temp);
		SEGMENT_TREE T(temp);
		SCS_SEGMENT_TREE.push_back(T);
	}
	
//	if(debug){
//		for(auto &i:FILE_NUMBER){
//			cerr<<i<<" ";
//		}
//	}

	if(debug){
		cerr<<"\nstarting to clustering metagenome with single cell sequences\n\n";
	}
	
	for(auto it=SCS_SEGMENT_TREE.begin();it!=SCS_SEGMENT_TREE.end();++it){
		it->detect_relation(Meta,cutoff);
	}
	
	
//	if(debug)cerr<<"Metasize:"<<Meta.size()<<'\n';
	
	MATRIX_DOUBLE MatrixMeta=matrixlization(Meta);
	//boost::numeric::ublas::matrix<double> MatrixSCS=matrixlization(SCS[0]);
	MATRIX_DOUBLE MatrixCosMetaMeta=cosangle(MatrixMeta,MatrixMeta);
	//boost::numeric::ublas::matrix<double> MatrixCosMetaSingle=cosangle(MatrixMeta,MatrixSCS);
	
//	if(debug)cerr<<MatrixCosMetaMeta<<'\n';
	
	int pTNF=gen_tnf_graph_sample(maxP,false,Meta.size(),MatrixCosMetaMeta);
	
	if(debug){
		cout<<"the pTNF is:"<<pTNF<<'\n';
	}
	
	gen_tnf_graph(maxP,Meta.size(),MatrixCosMetaMeta);
	
	clustering(e,Meta);
	
	cerr<<"This is the end of the program";
	
//	READSTAFILE(Folder,Meta.size());
//	for(int i=0;i<256;i++)cerr<<Meta[1].TNF[i];
// 	for(int i=0;i<SCS.size();i++){
// 		buildgraph(0,SCS[i].number,0);//锟斤拷锟接碉拷锟斤拷锟斤拷锟斤拷 
// 		buildgraph(SCS[i].number,0,0);
// 		SCS[i].category=1;
// 	}
// 	for(int i=0;i<Meta.size();i++){
// 		if(i%100==0)cerr<<i<<'\n';
// 		for(int j=0;j<SCS.size();j++){
// 			double cosinee=vector_compare(Meta[i],SCS[j]);
// 			double angle=acos(cosinee);
// //			cerr<<"cos -> "<<cosinee<<'\n';
// //			cerr<<"angle -> "<<rad2deg(angle)<<'\n';
// 			if(rad2deg(angle)>threshold)continue;
// 			buildgraph(Meta[i].number,0,sin(angle));//锟斤拷锟接碉拷细锟斤拷锟酵猴拷锟斤拷锟斤拷锟斤拷锟侥碉拷 
// 			buildgraph(0,Meta[i].number,sin(angle));
// 			break;
// 			Meta[i].category=1;
// 		}
// 	}
// 	for(int i=0;i<Meta.size();i++){
// 		if(i%100==0)cerr<<i<<'\n';
// 		for(int j=i+1;j<Meta.size();j++){
// 			double cosinee=vector_compare(Meta[i],Meta[j]);
// 			double angle=acos(cosinee);
// //			cerr<<"cos -> "<<cosinee<<'\n';
// //			cerr<<"angle -> "<<rad2deg(angle)<<'\n';
// 			if(rad2deg(angle)>threshold)continue;
// 			buildgraph(Meta[i].number,Meta[j].number,sin(angle));//锟斤拷锟接碉拷细锟斤拷锟酵猴拷锟斤拷锟斤拷锟斤拷锟侥碉拷 
// 			buildgraph(Meta[j].number,Meta[i].number,sin(angle));
// 		}
// 	}
// //	cerr<<"test->"<<abundancedistance(0,1,1,1.001)<<'\n';
// 	cerr<<"SAMPLESIZE -> "<<SAMPLE_SIZE<<'\n';
// 	for(int i=0;i<Meta.size();i++){
// 		for(int j=i+1;j<Meta.size();j++){
// 			double ADPij=ADP(Meta[i].number,Meta[j].number,CUTOFF);
// //			cerr<<"ADP of "<<Meta[i].number<<" and "<<Meta[j].number<<" is: "<<ADPij<<'\n';
// 		}
// 	}
//	LPA();
    return 0;
}

// fasta.exe -meta C:\Users\okura\Downloads\euk_test\euk_test\GCF_000227135.1_ASM22713v2_genomic.fna -single C:\Users\okura\Downloads\testdata\single -threshold 30 -folder C:\Users\okura\Downloads\testdata\statis

// fasta.exe -meta C:\Users\okura\Downloads\testdata\meta.fasta -single C:\Users\okura\Downloads\testdata\single -threshold 30 -folder C:\Users\okura\Downloads\testdata\statis

//fasta.exe -meta C:\Users\okura\Downloads\testdata\Meta_GPT.fasta -single C:\Users\okura\Downloads\testdata\singlegpt -threshold 30 -folder C:\Users\okura\Downloads\testdata\statis


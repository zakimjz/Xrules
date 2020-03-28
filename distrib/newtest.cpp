#include<iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stack>
#include <string>
#include <functional>
#include <numeric>
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h> 

using namespace std;

#define microsec 1000000.0

class TimeTracker {
private:
   struct timeval start_time;
   struct timeval stop_time;
   bool  running;
public:
   TimeTracker() {
      running=false;
   }

   void Start() {
      gettimeofday(&start_time, (struct timezone *)0);
      running=true;
   }

  double Stop() {
     double st, en;
     if (!running) return(-1.0);
     else {
        gettimeofday(&stop_time, (struct timezone *)0);
        st = start_time.tv_sec + (start_time.tv_usec/microsec);
        en = stop_time.tv_sec + (stop_time.tv_usec/microsec);
        running=false;
        return (double)(en - st);
     }
  }
};
 
#define BranchIt -1 //-1 indicates a branch in the tree
#define Xlabel "classX"
#define Ylabel "classY"
#define Nolabel "NULL"
#define MAXRATIO 1000

#define wProportional 0
#define wEqual 1
#define wInverse 2
#define wCustom 3
#define avgconf 0
#define bestrule 1
#define bestk 2
 
//four possible values of a comparison
enum comp_vals {equals, subset, superset, notequal};
ofstream summary;

bool printmissed=false;
bool usemissedrules = false;
bool count_unique=true;
bool joinrules=false;
bool predict_defaultclass=false;
bool use_defaultclass=false;
double minconfX, minconfY; //need minconf to classify a test case X, Y
bool printijtest=false;
bool posijinput=true;
bool setrules=false; //for set mining
bool use_relclasssup = true; //use relative class sup
bool use_norm_ratio = false; //use raw ratio as default
bool omit_high_ratio = false; //omit rules with MAXRATIO
bool use_conf = false;// use conf or likelihood?

//vars for counting accuracy and coverage
int truex=0, truey=0;
int falsex=0, falsey=0;
int npredictx = 0, npredicty = 0;

//vars for default class prediction
int defaulttruex[4] = {0,0,0,0};
int defaulttruey[4] = {0,0,0,0};
int defaultfalsex[4] = {0,0,0,0};
int defaultfalsey[4] = {0,0,0,0};
double defaultweightx[4] = {0,0,0,0};
double defaultweighty[4] = {0,0,0,0};
string defaultlabel[4] = {Nolabel,Nolabel,Nolabel,Nolabel};

//vars for weight model
int weight_model=wInverse;
double weightx=0, weighty=0; //weights for each class N/N_i

string testf;
string rulef;
string deff;
string mrulef;

int scoringtype = avgconf;
int bestkval = 0;

class FileInfo{
public:
   FileInfo(int n=0, int x=0, int y=0):Ntrans(n), Xtrans(x), Ytrans(y){}
   
   int Ntrans;
   int Xtrans;
   int Ytrans;
};

FileInfo Train, Test, Default;


void parseargs(int argc, char**argv)
{
   if (use_norm_ratio){
      minconfX = minconfY = 0.5;
   }
   else{
      minconfX = minconfY = 1.0;
   }
   
   //read options
   if (argc < 2){
      cout << "newtest:\n";
      cout << "-a setrules?\n";
      cout << "-cx minconf\n";
      cout << "-cy minconf\n";
      cout << "-C use_conf\n";
      cout << "-d predictdefaultclass\n";
      cout << "-D usedefaultclass\n";
      cout << "-h omit_high_ratio?\n";
      cout << "-j join original and missed rules into one array\n";
      cout << "-m printmissed?\n";
      cout << "-M missedrulefile\n";
      cout << "-o printijtest\n";
      cout << "-r rulefilename\n";
      cout << "-s scoring type\n";
      cout << "-t testfileprefix\n";
      cout << "-w weightmodel\n";
      cout << "-x noposijinput\n";
      cout << "-n use_norm_ratio\n";
      cout << "-e supratio in DB?\n";
      cout << endl;
      exit(-1);
   }
   else{
      int i=1;
      char c;
      while(i < argc){
         c = argv[i][1];
         switch(c){
         case 'a':
            setrules = true;
            break;
         case 'c':
            c = argv[i][2];
            if (c == 'x') minconfX = atof(argv[++i]);
            else minconfY = atof(argv[++i]);
            break;
         case 'C':
            use_conf = true;
            break;
         case 'd':
            predict_defaultclass = true;
            ++i;
            deff = argv[i];
            break;
         case 'D':
            use_defaultclass = true;
            ++i;
            deff = argv[i];
            break;
         case 'e':
            use_relclasssup = false;
            break;
         case 'h':
            omit_high_ratio = true;
            break;
         case 'j':
            joinrules=true;
            break;
         case 'm':            
            printmissed = true;
            break;
         case 'M':
            usemissedrules=true;
            ++i;
            mrulef = argv[i];
            break;
         case 'n':
            use_norm_ratio = true;
            break;
         case 'o':
            printijtest=true;
            break;
         case 'r':
            ++i;
            rulef = argv[i];
            break;
         case 's':
            ++i;
            scoringtype = atoi(argv[i]);
            if (scoringtype == bestk){
               ++i;
               bestkval = atoi(argv[i]);
            }
            break;
         case 't':
            ++i;
            testf = argv[i];
            break;
         case 'w':
            ++i;
            weight_model = atoi(argv[i]);
            if (weight_model == wCustom){
               weightx = atof(argv[++i]);
               weighty = atof(argv[++i]);
            }
            break;
         case 'x':
            posijinput=false;
            break;
         default:
            cout << "wrong option\n";
            exit(-1);
            break;
         }
         ++i;
      }
   }
   if ((use_norm_ratio || use_conf) && (minconfX+minconfY < 1)){
      cerr << "ERROR: ensure that minconfX >= 1-minconfY and vice versa\n";
      exit(-1);
   }
}


 
template<class T>
class printnode: public unary_function<T, void>
{
public:
   void operator() (T x){ cout << x << ", "; }
   void operator() (T *x){ cout << *x << ", "; }
};

template <class T>
ostream & operator<<(ostream& fout, vector<T>& vec)
{
   //fout << "no ptr " << endl;
   for_each(vec.begin(), vec.end(), printnode<T>());
   return fout;
}

template <class T>
ostream & operator<<(ostream& fout, vector<T *>& vec)
{
   //fout << "no ptr " << endl;
   for_each(vec.begin(), vec.end(), printnode<T>());
   return fout;
}

//structure to store frequent subtrees with sup and ratio info for two classes
 
class FreqIt: public vector<int>{
public:
   FreqIt():totsup(0),supX(0), supY(0), ratioX(0.0), ratioY(0.0), 
      confX(0.0), confY(0.0), truelabel(Nolabel){}
   void print(bool printlist=false, bool printsc=true, 
              bool printij=false, ostream& fout=cout){
      int i;
      fout << (*this)[0];
      for (i=1 ; i < size(); ++i){
         fout << " " << (*this)[i];
      }
      
      if (printsc){
         fout << " - " << totsup << " " << supX << " " << supY
              << " " << ratioX << " " << ratioY
              << " " << confX << " " << confY << " " << truelabel;
      }
      else{
         fout << " - " << truelabel;
      }
      if (printij){
         fout << " " << confX << " " << confY << " " 
              << getlabel() << " " << posi << " " << posj 
              << " " << pdbfile;
      }

      if (printlist){
         fout << " -- " << rulelist;
      }
      fout << endl;
   };


   string getlabel(){
      if (confX > minconfX) return Xlabel;
      else if (confY > minconfY) return Ylabel;
      else return Nolabel;
      //if (confY > minconfY) return Ylabel;
      //else return Xlabel;
   }
   
   //use lots of globals! Caution!
   string predict(bool pred_default=true){
      //determine if test case belongs to class X or Y
      calc_conf();
      string plabel = getlabel();
      
      if(plabel == Nolabel){
         if (pred_default){
            //no class can be determined, use default class
            if (printmissed)
               print(false,false);
            
            if (predict_defaultclass){
               //simply increase the default class count 
               if (truelabel == Xlabel){
                  ++defaulttruex[0]; 
               }
               else ++defaulttruey[0];
            }
            else if (use_defaultclass){
               for (int w=0; w < 4; ++w){
                  //loop thru all cost models
                  if (truelabel == defaultlabel[w]){
                     //correct prediction
                     if (truelabel == Xlabel){
                        ++defaulttruex[w];
                        ++npredictx;
                     }
                     else{
                        ++defaulttruey[w];
                        ++npredicty;
                     }
                  }
                  else{
                     //incorrect prediction
                     if (truelabel == Xlabel){
                        ++defaultfalsex[w];
                        ++npredicty;
                     }
                     else{
                        ++defaultfalsey[w];
                        ++npredictx;
                     }
                  }         
               }
            }
         }
      }
      else if (plabel == truelabel){
         //correct prediction, increment true count and npredict
         if (truelabel == Xlabel){
            ++truex;
            ++npredictx;
         }
         else{
            ++truey;
            ++npredicty;               
         }
      }
      else{
         //incorrect prediction, increment false count & npedict for plabel
         if (truelabel == Xlabel){
            ++falsex;
            ++npredicty; //since plabel = Ylabel
         }
         else{
            ++falsey;
            ++npredictx; //since plabel = Xlabel
         }
      }      
      return plabel;
   }
   
   void calc_conf()
   {
      //called for the test data cases (not for rules). 
      confX = 0;
      confY = 0;
      int sz = rulelist.size();
      if (sz > 0){
         switch(scoringtype){
         case avgconf:
            //sz = rulelist.size();
            break;
         case bestrule:
            sz = 1;
            break;
         case bestk:
            sz = min(bestkval, sz);
            break;
         }
         double Xsum=0, Ysum=0;
         for (int i=0; i < sz; ++i){
            Xsum += rulelist[i]->ratioX;
            Ysum += rulelist[i]->ratioY;
         }
         confX = Xsum/sz;
         confY = Ysum/sz;
         //cout<< "calcconf SZ " << sz << " " << confX << " " << confY << endl;
      }
   }

   void addrule(FreqIt *rule)
   {
      rulelist.push_back(rule);
   }
   
   void clear_rules()
   {
      rulelist.clear();
   }
   
   void clear()
   {
      vector<int>::clear();
      clear_rules();
      totsup = supX = supY = 0;
      ratioX = ratioY = confX = confY = 0;
      truelabel = Nolabel;
   }

   static bool compare (const FreqIt *f1, const FreqIt *f2)
   {
      //comapre by conf, sup, size, items
      double r1 = max(f1->ratioX, f1->ratioY);
      double r2 = max(f2->ratioX, f2->ratioY);
      //compare confidence
      if (r1 == r2){
         r1 = max(f1->supX, f1->supY);
         r2 = max(f2->supX, f2->supY);
         //compare support
         if (r1 == r2){
            //compare size
            if (f1->size() == f2->size()){
               //compare items in lex order
               for (int i = 0; i < f1->size(); ++i){
                  if ((*f1)[i] < (*f2)[i]) return true;
                  else if ((*f1)[i] < (*f2)[i]) return false;
               }
               return false;
            }
            else return (f1->size() > f2->size());
         }
         else return (r1 > r2);
      }
      else return (r1 > r2);
   }
      
   int totsup;
   int supX; //frequency in class 1 and 2 DB
   int supY;
   double ratioX; //ratio in c1 and c2 DB
   double ratioY;
   double confX;
   double confY;
   string truelabel;
   vector<FreqIt *> rulelist; //matching freq rules for test case
   
   //vars for test file
   int cid;
   int tid;
   int nitems;
   int posi;
   int posj;
   string pdbfile;
};
 
typedef vector<FreqIt *> FreqAry;
 
FreqAry fary; //rules
FreqAry mary; //missed rules

 
void read_defaultinfo(const char *deff)
{
   int tx, dtx, ty, dty;
   ifstream fin(deff, ios::in);
   if (!fin){
      cerr << "cannot open rule file\n";
      exit(-1);
   }
 
   fin >> tx;
   fin >> dtx;
   fin >> Default.Xtrans;
   fin >> ty;
   fin >> dty;
   fin >> Default.Ytrans;
   
   double acx = ((double)dtx)/Default.Xtrans;
   double acy = ((double)dty)/Default.Ytrans;
   for (int w=0; w < 4; ++w){
      switch (w){
      case wProportional:
         defaultweightx[w] = ((double)Default.Xtrans)/
            (Default.Xtrans+Default.Ytrans);
         defaultweighty[w] = ((double)Default.Ytrans)/
            (Default.Xtrans+Default.Ytrans);
         break;
      case wEqual:
         defaultweightx[w] = defaultweighty[w] = 0.5;
         break;
      case wInverse:
         defaultweightx[w] = ((double)Default.Ytrans)/
            (Default.Xtrans+Default.Ytrans);
         defaultweighty[w] = ((double)Default.Xtrans)/
            (Default.Xtrans+Default.Ytrans);
         break;
      case wCustom:
         defaultweightx[w] = weightx;
         defaultweighty[w] = weighty;
         break;
      }
      defaultlabel[w] = ((defaultweightx[w]*acx > defaultweighty[w]*acy)? 
                         Xlabel:Ylabel);
   }  
   fin.close();
}
void write_defaultinfo(const char *ffile)
{
   ofstream fout(ffile, ios::out);
   if (!fout){
      cerr << "cannot open rule file\n";
      exit(-1);
   }
   
   //   for (int i=0; i < ary.size(); ++i)
   //   ary[i]->print(false, true, fout);
   
   fout << truex << " " << defaulttruex[0] << " "
        << Test.Xtrans <<  " " << truey << " " << defaulttruey[0] << " " 
        << Test.Ytrans << endl;
   fout.close();
}


void read_rules(const char *ffile, FreqAry &ary)
{
   const int lineSize=8192;
   const int wdSize=256;
 
   FreqIt *fit;
   int it; 
   char inBuf[lineSize];
   string inStr;
   int inSize;
 
   ifstream fin(ffile, ios::in);
   if (!fin){
      cerr << "cannot open rule file\n";
      exit(-1);
   }
 
   bool skipline;
   while(fin.getline(inBuf, lineSize)){
      inSize = fin.gcount();
      //istrstream ist(inBuf, inSize);
      istringstream ist(inBuf, istringstream::in);
      fit = NULL;
      ist >> inStr;
      if (inStr == "INFILE" || inStr == "ALG" 
          || inStr == "DBASE_MAXITEM" || inStr == "MINSUPPORT"
          || inStr == "TIME" || inStr == "NUMMAX" 
          || inStr[0] == 'F' || inStr[0] == '['){
         continue;
      }
      else if (inStr == "DBASE_NUM_TRANS"){
         ist >> inStr;
         ist >> Train.Ntrans;
      }
      else if (inStr == "DBASE_NUMCL_TRANS"){
         ist >> inStr;
         ist >> it;
         if (it == 0) ist >> Train.Xtrans;
         else ist >> Train.Ytrans;
      }
      else{
         fit = new FreqIt;
         while (inStr != "-"){
            it = atoi(inStr.c_str());
            fit->push_back(it);
            ist >> inStr;
         }

         while(fit->back() == BranchIt) fit->pop_back();
         ist >> fit->totsup;
         ist >> fit->supX;
         ist >> fit->supY;
         
         double xval, yval;
         if (use_relclasssup){
            xval = ((double)fit->supX)/Train.Xtrans;
            yval = ((double)fit->supY)/Train.Ytrans;
         }
         else{
            xval = ((double)fit->supX)/Train.Ntrans;
            yval = ((double)fit->supY)/Train.Ntrans;
         }

         if (use_conf){
            fit->ratioX = ((double)fit->supX)/fit->totsup;
            fit->ratioY = ((double)fit->supY)/fit->totsup;
         }
         else if (use_norm_ratio){
            fit->ratioX = xval/(xval+yval); //likelihood
            fit->ratioY = yval/(xval+yval);
         }
         else{
            if (yval > 0) fit->ratioX = xval/yval;
            else fit->ratioX = MAXRATIO;
            if (xval > 0) fit->ratioY = yval/xval;
            else fit->ratioY = MAXRATIO;
         }

         //cout << "RATIO " << fit->ratioX << " " << fit->ratioY << endl;
         if (omit_high_ratio && 
             (fit->ratioX == MAXRATIO || fit->ratioY == MAXRATIO)) continue;
         
         if (fit->ratioX > minconfX || fit->ratioY > minconfY)
            ary.push_back(fit);

         //cout << "RULE "; fit->print();
      }
   }
   cout << "RULE ARY SIZE " << ary.size() << endl;
   sort (ary.begin(), ary.end(), FreqIt::compare);
   //for (it = 0; it < ary.size(); ++it)
   //   ary[it]->print();
   fin.close();
}

bool read_next_test(istream &fin, FreqIt &fdata)
{
   static const int lineSize=8192;
   static const int wdSize=256;
   static char inBuf[lineSize];
   static char inStr[wdSize];

   string label;
   int inSize;
   int lbl;
   int i;
   
   fdata.clear();
   
   if(fin.getline(inBuf, lineSize)){
      inSize = fin.gcount();
      int it;
      //istrstream ist(inBuf, inSize);
      istringstream ist(inBuf, istringstream::in);
      
      ist >> lbl; //read class
      if (lbl==0){
         fdata.truelabel = Xlabel;
         ++Test.Xtrans;
      }
      else{
         fdata.truelabel=Ylabel;
         ++Test.Ytrans;
      }
      ++Test.Ntrans;
      ist >> fdata.cid; //skip cid
      //cout << inStr << endl << flush;
      ist >> fdata.tid; //skip tid
      ist >> fdata.nitems; //skip #items
      for(i= 0; i < fdata.nitems; ++i){
         ist >> it;
         fdata.push_back(it);
         //cout << it << " ";
      }

      while(fdata.back() == BranchIt) fdata.pop_back();

      if (posijinput){
         ist >> fdata.posi;
         ist >> fdata.posj;
         ist >> fdata.pdbfile;
      }
      
      //cout << "TEST "; fdata.print();
      return true;
   }
   return false; //means the data is over
}

 
comp_vals set_compare(FreqIt &fit1, FreqIt &fit2)
{
   int i,j;
   comp_vals res = notequal;
   
   if (fit1.size() <= fit2.size()){
      bool foundflg = false;
      for (i=0,j=0; i < fit1.size() && j < fit2.size() && !foundflg;){
         if (fit1[i] == fit2[j]){ 
            ++i; 
            ++j;
            if (i == fit1.size()) foundflg = true;
         }
         else if (fit1[i] < fit2[j]){ ++j; }
         else { break; } // not found
      }

      //cout << "came here " << fit1.size() << " " << fit2.size() << endl;
      if (foundflg){
         if (fit1.size() == fit2.size()) res = equals;
         else res = subset;
      }
      else res = notequal;
   }
   else{
      res = set_compare(fit2, fit1);
      if (res == subset) res = superset; //res cannot be equal or superset
      else res = notequal;
   }
   
   return res;
}

void check_subtree(FreqIt &fit, FreqIt &fit2,
                  int tpos, int ppos, int tscope,
                  stack<int> &stk, bool &foundflg)
{
   int i;
   int scope, ttscope;
   stack<int> tstk;

   scope = tscope;
   bool skip = false;
   if (fit[ppos] == BranchIt){
      skip = true;
      while(fit[ppos] == BranchIt){
         tstk.push(stk.top());
         stk.pop();
         ppos++;
      }
      tscope = tstk.top();
   }

   while (skip && scope >= tscope && tpos < fit2.size()){
      if (fit2[tpos] == BranchIt) scope--;
      else scope++;
      tpos++;
   }

   if (skip) tscope = stk.top();

   for (i=tpos; i < fit2.size() && !foundflg; i++){
      if (fit2[i] == BranchIt) scope--;
      else scope++;
      if (scope < tscope) break;
      if (fit2[i] == fit[ppos]){
         stk.push(scope);

         //for (int d = 0; d < stk.size(); d++) cout << "\t";

         if (ppos == fit.size()-1){
            //cout << ppos << " found at " << i << " " << scope << endl;
            //fit.dbsup++;
            foundflg = true;
         }
         else{
            //cout << ppos << " recurse at " << i << " " << scope << endl;
            check_subtree(fit, fit2, i+1, ppos+1, scope, stk, foundflg);
         }
         stk.pop();
      }
   }
   while(!tstk.empty()){
      stk.push(tstk.top());
      tstk.pop();
   }
}

comp_vals tree_compare(FreqIt &fit1, FreqIt &fit2)
{
   stack<int> stk;
   bool foundflg = false;

   comp_vals res = notequal;

   if (fit1.size() <= fit2.size()){
      check_subtree(fit1, fit2, 0, 0, 0, stk, foundflg);
      //cout << "came here " << fit1.size() << " " << fit2.size() << endl;
      if (foundflg){
         if (fit1.size() == fit2.size()) res = equals;
         else res = subset;
      }
      else res = notequal;
   }
   else{
      check_subtree(fit2, fit1, 0, 0, 0, stk, foundflg);
      if (foundflg) res = superset;
      else res = notequal;
   }

   return res;
}
 

void get_match_rules (FreqIt &testcase, FreqAry &rules)
{
   comp_vals cv;
   for (int j=0; j < rules.size(); j++){ 
      //compare pat1 vs pat2
      if (setrules) cv = set_compare (testcase, *rules[j]);
      else cv = tree_compare (testcase, *rules[j]);
      
      switch (cv){
      case equals: 
      case superset:
         //if rules matches testcase then add to test instance
         //cout << "EQSUP:  " << testcase << " --- " << *rules[j] << endl;
         testcase.addrule(rules[j]); 
         break;
      case subset: //in the following two case do nothing.
      case notequal:  
         //cout << "NOMATCH: " << testcase << " -- " << *rules[j] << endl;
         break;
      }
   }

   //testcase.print(true,true);
}



void print_covacc_use_default()
{
   double Xacc, Xcov, Yacc, Ycov, totacc, totcov;
   double wx, wy;
   for (int w=0; w < 4; ++w){
      if (w == wCustom && weight_model != wCustom) break;
      
      cout << "DEFAULT PREDICTIONS, LABEL= " 
           << defaultlabel[w] << " " << w << endl;
      wx = defaultweightx[w];
      wy = defaultweighty[w];
      
      cout << "class\t True\t False\n";
      cout << "X\t" << defaulttruex[w] << "\t" << defaultfalsex[w] << endl;
      cout << "Y\t" << defaulttruey[w] << "\t" << defaultfalsey[w] << endl;

      summary << "defaultw= " << w << " " << defaulttruex[w] << " " 
              << defaultfalsex[w] << " " << defaulttruey[w] << " " 
              << defaultfalsey[w] << " ";

      int numtx = truex+defaulttruex[w];
      int numfx = falsex+defaultfalsex[w];
      int numty = truey+defaulttruey[w];
      int numfy = falsey+defaultfalsey[w];

      Xacc = ((double) numtx)/(numtx+numfx);
      Yacc = ((double) numty)/(numty+numfy);
      Xcov = ((double) numtx)/Test.Xtrans;
      Ycov = ((double) numty)/Test.Ytrans;

      totacc = wx*Xacc+wy*Yacc;
      totcov = 1.0;
   
      cout << "TRUE/FALSE INFO FOR NPREDICT+DEFAULT\n";
      cout << "class\tTrue\tFalse\tAccuracy\tCoverage\n";
      cout << "X\t" << numtx << "\t" << numfx << "\t" 
           << Xacc << "\t" << Xcov << endl;      
      cout << "Y\t" << numty << "\t" << numfy << "\t" 
           << Yacc << "\t" << Ycov << endl;
      cout << "TOTAL\t" << Test.Xtrans << "\t" << Test.Ytrans << "\t"
           << totacc << "\t" << totcov << endl;
      
      summary << "full W=" << w << " " << Xacc << " " << Xcov << " " 
              << Yacc << " " << Ycov << " "
              << totacc << " " << totcov << " ";
   }
}


void processfreq(istream &fin)
{
   if (use_defaultclass) read_defaultinfo(deff.c_str());  

   FreqIt testdata;
   string plabel;
   while(read_next_test(fin, testdata)){      
      get_match_rules (testdata, fary);
      if (usemissedrules && !joinrules && !predict_defaultclass){
         plabel = testdata.predict(false);
         if (plabel == Nolabel){
            testdata.clear_rules();
            get_match_rules (testdata, mary);
            testdata.predict();
         }
      }
      else testdata.predict();
      if (printijtest){
         testdata.print(false, false, true);
      }     
   }
   
   double Xacc, Xcov, Yacc, Ycov, totacc, totcov;
   double wx, wy;
   Xacc = ((double)truex)/(truex+falsex);
   Xcov = ((double)truex)/Test.Xtrans;
   Yacc = ((double)truey)/(truey+falsey);
   Ycov = ((double)truey)/Test.Ytrans;
   totacc = ((double)(truex + truey))/(truex+truey+falsex+falsey);
   totcov = ((double)(truex + truey))/(Test.Xtrans+Test.Ytrans);

   cout << "TRUE/FALSE INFO FOR NPREDICT\n";
   cout << "class\tTrue\tFalse\tAccuracy\tCoverage\n";
   cout << "X\t" << truex << "\t" << falsex << "\t" 
        << Xacc << "\t" << Xcov << endl;      
   cout << "Y\t" << truey << "\t" << falsey << "\t" 
        << Yacc << "\t" << Ycov << endl;      
   cout << "NPRED\t" << npredictx << "\t" << npredicty << endl;
   cout << "TOTAL\t" << Test.Xtrans << "\t" << Test.Ytrans<< "\t"
        << totacc << "\t" << totcov << endl;      

   summary << "partial= " << Xacc << " " << Xcov << " " 
           << Yacc << " " << Ycov << " "
           << totacc << " " << totcov << " ";
   
   summary << "counts= " << truex << " " << falsex << " "
           << truey << " " << falsey << " ";

   if (use_defaultclass) print_covacc_use_default();
   else if (predict_defaultclass) write_defaultinfo(deff.c_str());
   
   summary << "ntrans= " << Test.Xtrans << " " << Test.Ytrans << " ";

   
}
 
int main(int argc, char **argv)
{
   TimeTracker tt;
   tt.Start();
   parseargs(argc,argv);
   summary.open("summary.out", ios::app);
   summary << "NEWTEST"
           << " -r " << rulef 
           << " -t " << testf 
           << " -m " << (usemissedrules?"true":"false")
           << " -M " << mrulef
           << " -s " << scoringtype  << " " << bestkval
           << " -cx " << minconfX
           << " -cy " << minconfY
           << " -w " << weight_model << " " << weightx << " " << weighty
           << " -j " << (joinrules?"true ":"false ")
           << " -C " << (use_conf?"true ":"false ")
           << " -d " << (predict_defaultclass?"true ":"false ")
           << " -D " << (use_defaultclass?"true ":"false ")
           << " -n " << (use_norm_ratio?"true ":"false ")
           << " -e " << (use_relclasssup?"true ":"false ");


   //read patterns from pat file
   read_rules (rulef.c_str(), fary);
   if (usemissedrules){
      if (joinrules) //join both rules in same array
         read_rules (mrulef.c_str(), fary); 
      else //keep original and missed separate
         read_rules (mrulef.c_str(), mary);
   }

   summary << "numrules = " << fary.size() << " " << mary.size() << " ";
   
   ifstream fin(testf.c_str(), ios::in);
   if (!fin){
      cerr << "cannot open test file\n";
      exit(-1);
   }

   processfreq(fin);

   fin.close();
   
   double tottime = tt.Stop();
   summary << tottime << endl;
   summary.close();
   
   exit(0);
}
 
 
 

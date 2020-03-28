#include<iostream>
#include <unistd.h>
#include <algorithm>
#include <stdio.h>
#include <list>
#include <limits.h>

using namespace std;

//headers
#include "xminer.h"
#include "timetrack.h"
#include "calcdb.h"
#include "eqclass.h"
#include "stats.h"
#include "hashtable.h"

#define VR_ITEM 999999999

//global vars
char infile[300];
Dbase_Ctrl_Blk *DCB;
Stats stats;

vector<double> MINSUP_PER(Max_Num_Classes,0.0);
vector<int> MINSUPPORT(Max_Num_Classes,-1);
int DBASE_MAXITEM;
int DBASE_NUM_TRANS;
vector<int> DBASE_NUMCL_TRANS(Max_Num_Classes,0);
int Num_Classes = 2; //default

//default flags
bool output = false; //don't print freq subtrees
bool output_idlist = false; //don't print idlist
bool count_unique = true; //count support only once per tree
bool use_fullpath = false; //use reduced scope to conserve mem
bool use_class_ratio=false; //relative sup of pattern wrt class
bool print_ratio=false; //print support ratio in DB or in class
int omit_item = VR_ITEM; //assume that VR_ITEM never occurs
int max_num_patterns = INT_MAX;
int num_patterns=0;

sort_vals sort_type = nosort; //default is to sort in increasing order
alg_vals alg_type = treeminer; //default is to find all freq patterns

prune_vals prune_type = noprune; //no prune
FreqHT FK; //to store freq subtrees for pruning
vector<int> clcnt[Max_Num_Classes];

vector<int> *ITCNT = NULL; //used for sorting F1
bool F1cmp(int x, int y){
   bool res = false;
   if ((*ITCNT)[x] < (*ITCNT)[y]) res = true;

   if (sort_type == incr) return res;
   else return !res;
}
   

void parse_args(int argc, char **argv)
{
   extern char * optarg;
   char c;
   int i,j;

   if (argc < 2){
      cout << "usage: treeminer\n";
      cout << "\t-i<infile>\n";
      cout << "\t-n<num_classes (default=2)> (specified before -s/-S)\n";
      cout << "\t-s<supports, one per class>\n";
      cout << "\t-S<absolute supports, one per class> (specify -s or -S)\n";
      cout << "\t-o (print patterns?)\n";
      cout << "\t-c<use relative support per class>\n";
      cout << "\t-C (print ratios?)\n";
      cout << "\t-r<max number of patterns to mine>\n";
      cout << "\t-b (binary_input?)\n";
      
   }
   else{
     for(i=1; i < argc; ++i){
       //     while ((c=getopt(argc,argv,"a:fi:lop:s:S:uz:"))!=-1){
       c=argv[i][1];
       //cout << "ARGVAL " << i << " -- " << c << endl;
       switch(c){
       case 'a':	 
          alg_type = (alg_vals) atoi(argv[++i]);
          break;
       case 'b':	 
          Dbase_Ctrl_Blk::binary_input = true;
          break;
       case 'c':
          ++i;
          use_class_ratio= (atoi(argv[i])?true:false);
          break;
       case 'C':
          print_ratio = true;
          break;
        case 'f':
          use_fullpath = true;
          break;
       case 'i': //input files
	 sprintf(infile,"%s",argv[++i]);
	 //cout << infile << endl;
	 break;
       case 'l': //print idlists along with freq subtrees
          output=true;
          output_idlist = true;
          break;
       case 'n':
          ++i;
          Num_Classes = atoi(argv[i]);
          if (Num_Classes > Max_Num_Classes){
             cout << "error: max classes can be " << Max_Num_Classes << endl;
             exit(-1);
          }
          break;
       case 'o': //print freq subtrees
	 output = true;
	 break;
       case 'p':
          prune_type = (prune_vals) atoi(argv[++i]);
          break;
       case 'r':
          max_num_patterns = atoi(argv[++i]);
          break;
       case 's': //support value for L2
          for (j=0; j < Num_Classes; ++j){
             MINSUP_PER[j] = atof(argv[++i]);
             //cout << "SUP = " << MINSUP_PER[j] << endl;
          }
          break;
       case 'S': //absolute support
          for (j=0; j < Num_Classes; ++j)
             MINSUPPORT[j] = atoi(argv[++i]);
	 break;
       case 'u': //count support multiple times per tree
	 count_unique = false;
	 break;
       case 'x':
          omit_item = atoi(argv[++i]);
          break;
       case 'z':
          sort_type = (sort_vals) atoi(argv[++i]);
          break;
       }               
     }
   }
}

void get_F1()
{
  TimeTracker tt;
  double te;

  int i, j, k, it;
  const int arysz = 10;
  
  vector<int> itcnt(arysz,0);
  
  for(j=0; j < Num_Classes; ++j) clcnt[j].resize(arysz,0);

  vector<int> flgs(arysz,-1);

  tt.Start();

  DBASE_MAXITEM=0;
  DBASE_NUM_TRANS = 0;
  
   while(DCB->get_next_trans())
   {
      //cout << "TRANS " << DCB->Cid << " " << DCB->Tid
      //     << " " << DCB->TransSz << " -- ";
      //for (i=0; i < DCB->TransSz; i++)
      //   cout << " " << DCB->TransAry[i];
      //cout << endl;
      for (i=0; i < DCB->TransSz; i++){
         it = DCB->TransAry[i];
         if (it == BranchIt) continue;
         
         if (it >= DBASE_MAXITEM){
            itcnt.resize(it+1,0);
            for (j=0; j < Num_Classes; ++j)
               clcnt[j].resize(it+1,0);
            flgs.resize(it+1,-1);
            DBASE_MAXITEM = it+1;
            //cout << "IT " << DBASE_MAXITEM << endl;
         }
         
         if (count_unique){
            if(flgs[it] == DCB->Cid) continue;
            else flgs[it] = DCB->Cid;
         }
         itcnt[it]++;
         clcnt[DCB->Class][it]++;
      }
      
      if (DCB->MaxTransSz < DCB->TransSz) DCB->MaxTransSz = DCB->TransSz;     
      DBASE_NUM_TRANS++;
      DBASE_NUMCL_TRANS[DCB->Class]++;
   }

   if (omit_item != VR_ITEM){ //ignore item "omit_item" totally even if freq
      itcnt[omit_item] = 0;
      for (j=0; j < Num_Classes; ++j) clcnt[j][omit_item] = 0; 
   }
   
   //for (i=0; i < DCB->TransSz; i++){
   //   it = DCB->TransAry[i];
   //   if (it != BranchIt){
   //      cout << it << " " << itcnt[it] << endl;
   //   }
   //}

   //set the value of MINSUPPORT
   for (j=0; j < Num_Classes; ++j){
      if (MINSUPPORT[j] == -1)
         MINSUPPORT[j] = (int) (MINSUP_PER[j]*DBASE_NUMCL_TRANS[j]+0.5);
      
      if (MINSUPPORT[j]<1) MINSUPPORT[j]=1;
   }
   
   cout<<"DBASE_NUM_TRANS : "<< DBASE_NUM_TRANS << endl;
   for (i=0; i < Num_Classes; ++i)
      cout<<"DBASE_NUMCL_TRANS : "<< i << " " << DBASE_NUMCL_TRANS[i] << endl;
   cout<<"DBASE_MAXITEM : "<< DBASE_MAXITEM << endl;
   for (j=0; j < Num_Classes; ++j)
      cout<<"MINSUPPORT : " << j << " " 
          << MINSUPPORT[j] << " (" << MINSUP_PER[j] << ")" << endl;

   //count number of frequent items
   DCB->NumF1 = 0;
   for (i=0; i < DBASE_MAXITEM; i++){
      bool flg=false;
      for (j=0; j < Num_Classes && !flg; ++j)
         if (clcnt[j][i] >= MINSUPPORT[j])
            flg=true;
      if (flg) DCB->NumF1++;
   }
   
      
   int *it_order = new int[DBASE_MAXITEM];
   for (i=0; i < DBASE_MAXITEM; i++)
      it_order[i] = i;
   
   if (sort_type != nosort){
      ITCNT = &itcnt;
      sort(&it_order[0], &it_order[DBASE_MAXITEM], F1cmp);
   }

   //construct forward and reverse mapping from items to freq items
   DCB->FreqIdx = new int [DCB->NumF1];
   DCB->FreqMap = new int [DBASE_MAXITEM];
   for (i=0,j=0; i < DBASE_MAXITEM && num_patterns < max_num_patterns; i++) {
      bool flg=false;
      for (k=0; k < Num_Classes && !flg; ++k)
         if (clcnt[k][it_order[i]] >= MINSUPPORT[k]) flg=true;
         
      if (flg){
         if (output){ //cout << i << " " << it_order[i] 
            cout << it_order[i]
                 << " - " << itcnt[it_order[i]];
            for (k=0; k < Num_Classes; ++k)
               cout << " " << clcnt[k][it_order[i]];
            
            if (print_ratio){
               if (use_class_ratio){
                  for (k=0; k < Num_Classes; ++k)
                     cout << " " << (1.0*clcnt[k][it_order[i]])/
                        DBASE_NUMCL_TRANS[k];
               }
               else{
                  for (k=0; k < Num_Classes; ++k)
                     cout << " " << (1.0*clcnt[k][it_order[i]])/
                        DBASE_NUM_TRANS;
               }
            }
            cout << endl;
         }
         DCB->FreqIdx[j] = it_order[i];
         DCB->FreqMap[it_order[i]] = j;
         ++j;
         
         ++num_patterns;
      }
      else DCB->FreqMap[it_order[i]] = -1;
   }
   
   //cout<< "F1 - " << DCB->NumF1 << " " << DBASE_MAXITEM << endl;  
   
   if (sort_type != nosort){
      ITCNT = NULL;
      delete [] it_order;
   }
   
   te = tt.Stop();
   stats.add(DBASE_MAXITEM, num_patterns, te);

}

list<Eqclass *> * get_F2()
{
  int i,j,k;
  int it1, it2;
  int scnt;
  
  TimeTracker tt;
  double te;

  tt.Start();

  list<Eqclass *> *F2list = new list<Eqclass *>;

  //itcnt2 an array of counts, flgs is array of flags
  int **itcnt2 = new int*[DCB->NumF1];
  int **clcnt2[Num_Classes];
  for (i=0; i < Num_Classes; ++i) clcnt2[i] = new int*[DCB->NumF1];
  int **flgs = new int*[DCB->NumF1];

  for (i=0; i < DCB->NumF1; i++){
    itcnt2[i] = new int [DCB->NumF1];
    for (j=0; j < Num_Classes; ++j) clcnt2[j][i] = new int [DCB->NumF1];
    flgs[i] = new int [DCB->NumF1];
    //cout << "alloc " << i << " " << itcnt2[i] << endl;
    for (j=0; j < DCB->NumF1; j++){
      itcnt2[i][j] = 0;
      for (k=0; k < Num_Classes; ++k) clcnt2[k][i][j] = 0;
      flgs[i][j] = -1;
    }
  }
    
   DCB->alloc_idlists();
   
   while(DCB->get_next_trans())
   {
      DCB->get_valid_trans();
      DCB->make_vertical();
      //DCB->print_trans();
      //count a pair only once per cid
      for (i=0; i < DCB->TransSz; i++){
         it1 = DCB->TransAry[i];
         if (it1 != BranchIt){
            scnt = 0;
            for (j=i+1; scnt >= 0 && j < DCB->TransSz; j++){
               it2 = DCB->TransAry[j];
               if (it2 != BranchIt){
		 scnt++;
		 if (count_unique){
		   if (flgs[it1][it2] == DCB->Cid) continue;
		   else flgs[it1][it2] = DCB->Cid;
		 }
		 //cout << "cnt " << it1 << " " << it2 << endl;
		 itcnt2[it1][it2]++;
		 clcnt2[DCB->Class][it1][it2]++;
               }
               else scnt--;
            }
         }
      }
   }                           
   
   int F2cnt=0;

   // count frequent patterns and generate eqclass
   Eqclass *eq;
   for (i=0; i < DCB->NumF1; i++) {
      eq = NULL;
      for (j=0; j < DCB->NumF1 && num_patterns < max_num_patterns; j++) {
         bool flg=false;
         //cout << "access " << i << " " << j << endl;
         for (k=0; k < Num_Classes && !flg; ++k)
            if (clcnt2[k][i][j] >= MINSUPPORT[k]) flg=true;
         
         if (flg){
            F2cnt++;
            if (eq == NULL){
               eq = new Eqclass;
               eq->prefix().push_back(i);
            }
	    eq->add_node(j,0);
            
	    if (output) {
	      cout << DCB->FreqIdx[i] << " " << DCB->FreqIdx[j] 
		   << " - " << itcnt2[i][j];
              for (k=0; k < Num_Classes; ++k)
                 cout << " " << clcnt2[k][i][j];

              if (print_ratio){
                 if (use_class_ratio){
                    for (k=0; k < Num_Classes; ++k)
                       cout << " " << 
                          (1.0*clcnt2[k][i][j])/DBASE_NUMCL_TRANS[k];
                 }
                 else{
                    for (k=0; k < Num_Classes; ++k)
                       cout << " " << (1.0*clcnt2[k][i][j])/DBASE_NUM_TRANS;
                 }
              }
              
              cout << endl;
            }

            ++num_patterns;
         }
      }   
      if (eq != NULL) F2list->push_back(eq);
   }
   
   for (i=0; i < DCB->NumF1; i++) {
     //cout << "dealloc " << i << " " << itcnt2[i] << endl;
     delete [] itcnt2[i];
     for (j=0; j < Num_Classes; ++j) delete [] clcnt2[j][i];
     //cout << "dealloc " << i << " " << flgs[i] << endl;
     delete [] flgs[i];
   }

   for (j=0; j < Num_Classes; ++j) delete [] clcnt2[j];
   delete [] itcnt2;
   delete [] flgs;
   
   
   //cout << "F2 - " << F2cnt << " " << DCB->NumF1 * DCB->NumF1 << endl;
   te = tt.Stop();
   stats.add(DCB->NumF1 * DCB->NumF1, F2cnt, te);

   return F2list;
}

static bool notfrequent (Eqnode &n){
  //cout << "IN FREQ " << n.sup << endl;
   bool flg=false;
   for (int i=0; i < Num_Classes && !flg; ++i)
      if (n.clsup[i] >= MINSUPPORT[i]) 
         flg=true; 
   if (flg) return false;
   else return true;
}


void check_ins(idlist *l1, idlist *l2, Eqnode *ins, 
                int st1, int st2, int en1, int en2, int clid){
   int i;
   static idnode *n1, *n2;
   static ival_type cmpval;

   bool found_flg = false;
   int pos1 = st1; //temporary position holder for n1 ival
   bool next2; //should we go to next n2 ival?

   //for each ival in n2, find the closest parent in n1
   int tpos = ins->tlist.size();
   while(st2 < en2){
      n1 = &(*l1)[pos1];
      n2 = &(*l2)[st2];

      next2 = true; //by default we go to next n2 ival
      cmpval = ival::compare(n1->itscope, n2->itscope);
      switch (cmpval){
      case sup: 
         //n2 was under some n1

         if (n1->path_equal(*n2)){
            if (en1-st1 > 1 || use_fullpath){
               ins->tlist.push_back(idnode(n2->cid, n1->parscope, 
                                           n1->itscope.lb, n2->itscope));
            }
            else{
               ins->tlist.push_back(idnode(n2->cid,n2->itscope));
            }
            if (!count_unique){
               ins->sup++;
               ins->clsup[clid]++;
            }
            
            found_flg = true;
         }
         
         next2 = false;
         break;
      case before: 
         //check next n1 ival for same n2 ival
         next2 = false; 
         break;
      }

      if (next2 || pos1+1 == en1){ //go to next n2 ival
         pos1 = st1;
         st2++;
      }
      else pos1++; //go to next n1 ival
   }
   
   if (found_flg && count_unique){
      ins->sup++;
      //cout << "CLASSID " << clid << endl;
      ins->clsup[clid]++;
   }
   
}

   
void check_outs(idlist *l1, idlist *l2, Eqnode *outs,
                int st1, int st2, int en1, int en2, int clid){
   static idnode *n1, *n2;
   static ival_type cmpval;
   int i;
   
   bool found_flg = false;
   bool next2;
   int pos1 = st1;
   
   //for each n2 ival find if there is a sibling to the left
   int tpos = outs->tlist.size();
   while(st2 < en2){
      n1 = &(*l1)[pos1];
      n2 = &(*l2)[st2];
      next2 = true;
      cmpval = ival::compare(n1->itscope, n2->itscope);
      switch (cmpval){
      case sup:
         next2 = false;
         break;
      case before: 
         //n1 is before n2. Check if n1.par is subset of or equal to n2.par
         if (n1->path_equal(*n2)){
            if (en1 - st1 > 1 || use_fullpath){
               outs->tlist.push_back(idnode(n2->cid, n1->parscope, 
                                            n1->itscope.lb, n2->itscope));
            }
            else{
               outs->tlist.push_back(idnode(n2->cid,n2->itscope));
            }
            if (!count_unique){
               outs->sup++;
               outs->clsup[clid]++;
            }
            
            found_flg = true;
         }
         next2 = false;
         break;
      }
      if (next2 || pos1+1 == en1){ //go to next n2 ival
         pos1 = st1;
         st2++;
      }
      else pos1++; //go to next n1 ival
   }
   if (found_flg && count_unique){
      outs->sup++;
      outs->clsup[clid]++;
   }
   
}


void get_intersect(idlist *l1, idlist *l2, Eqnode *ins, Eqnode *outs)
{
   static idnode *n1, *n2;
   int clid;
   int i1 = 0, i2 = 0;
   int e1, e2;

   while (i1 < l1->size() && i2 < l2->size()){
      n1 = &(*l1)[i1];
      n2 = &(*l2)[i2];

      //look for matching cids
      if (n1->cid < n2->cid) i1++;
      else if (n1->cid > n2->cid) i2++;
      else{
         //cids match
         e1 = i1;
         e2 = i2;

         //check the cid end positions in it1 and it2
         while (e1 < l1->size() && (*l1)[e1].cid == n1->cid) e1++;
         while (e2 < l2->size() && (*l2)[e2].cid == n2->cid) e2++;

         //increment support if candidate found
         clid = Dbase_Ctrl_Blk::classidx[n1->cid];
         if (ins) check_ins(l1, l2, ins, i1, i2, e1, e2, clid);
         if (outs) check_outs(l1, l2, outs, i1, i2, e1, e2, clid);

         //restore index to end of cids
         i1 = e1;
         i2 = e2;
      }
   }
}


bool lexsmaller(vector<int> &subtree, vector<int> &cand)
{
   int i,j;

   for (i=0, j=0; i < subtree.size() && j < cand.size();){

      if (subtree[i] > cand[j]){
         if (cand[j] != BranchIt) return false;
         else{
            while (cand[j] == BranchIt) j++;
            if (subtree[i] > cand[j]) return false;
            else if (subtree[i] < cand[j]) return true;
            else return false;
         }
         
      }
      else if (subtree[i] < cand[j]){
         if (subtree[i] != BranchIt) return true;
         else{
            while(subtree[i] == BranchIt) i++;
            if (subtree[i] > cand[j]) return false;
            else if (subtree[i] < cand[j]) return true;
            else return true;
         }
      }
      else{
         i++;
         j++;
      }
   }
   return false;
}


Eqnode *test_node(int iter, Eqclass *eq, int val, int pos)
{
   Eqnode *eqn = NULL;
   
   //if noprune, return with a new Eqnode
   if (prune_type == noprune){
      eqn = new Eqnode(val,pos);
      return eqn;
   }
   
   //perform pruning

   //prune based on frequent subtree
   static vector<int> cand;
   static vector<int> subtree;
 
   int hval;
   int scope, scnt;

   //form the candidate preifx
   cand = eq->prefix();
   scnt = eq->get_scope(pos, scope); //what is the scope of node.pos

   while(scnt > scope){
      cand.push_back(BranchIt);
      scnt--;
   }
   cand.push_back(val);

   //check subtrees
   int omita, omitb;
   bool res = true;
   //omit i-th item (omita) and associated BranchIt (omitb)
   int i,j,k;

   for (i=iter-3; i >= 0; i--){
      //find pos for i-th item
      for (j=0,k=0; j < cand.size(); j++){
         if (cand[j] != BranchIt){
            if (k == i){
               omita = j;
               break;
            }
            else k++;
         }
      }
      
      //find pos for associated BranchIt
      scnt = 0;
      for(j++; j < cand.size() && scnt >= 0; j++){
         if (cand[j] == BranchIt) scnt--;
         else scnt++;
      }
      if (scnt >= 0) omitb = cand.size();
      else omitb = j-1;

      //cout << "OMIT " << i << " " << omita << " " << omitb << endl;

      hval = 0;
      subtree.clear();
      bool rootless = false;
      scnt = 0;
      for (k=0; k < cand.size() && !rootless; k++){
         if (k != omita && k != omitb){
            subtree.push_back(cand[k]);
            if (cand[k] != BranchIt){
               hval += cand[k];
               scnt++;
            }
            else scnt--;
            if (scnt <= 0) rootless = true;
            
         }
      }

      //cout << "LEXTEST " << subtree << " vs " << cand;
   
      //skip a rootless subtree
      if (!rootless && lexsmaller(subtree, cand)){
         //cout << " -- SMALLER ";
         res = FK.find(iter-1, subtree, hval);  
         //cout << ((res)? " ** FOUND\n":" ** NOTFOUND\n");
         if (!res) return NULL; //subtree not found!
      }
      //else cout << " -- GREATER " << endl;
   }
   
  
   if (res) eqn = new Eqnode(val,pos);
   else eqn = NULL;

   return eqn;
}


void enumerate_freq(Eqclass *eq, int iter)
{
   TimeTracker tt;
   Eqclass *neq;
   list<Eqnode *>::iterator ni, nj;
   Eqnode *ins, *outs;

   if (prune_type == noprune) eq->sort_nodes(); //doesn't work with pruning

   //cout << "FX " << *eq << endl;

   for (ni = eq->nlist().begin(); ni != eq->nlist().end() &&
           num_patterns < max_num_patterns; ++ni){
      neq = new Eqclass;
      neq->set_prefix(eq->prefix(),*(*ni));
      tt.Start();
      for (nj = eq->nlist().begin(); nj != eq->nlist().end() &&
              num_patterns < max_num_patterns; ++nj){ 
         if ((*ni)->pos < (*nj)->pos) continue;

         ins = outs = NULL;
         if ((*ni)->pos > (*nj)->pos){
            outs = test_node(iter, neq, (*nj)->val, (*nj)->pos);
         }
         else{ 
            outs = test_node(iter, neq, (*nj)->val, (*nj)->pos);
            ins = test_node(iter, neq, (*nj)->val, neq->prefix().size()-1);
         }  

         //cout << "prefix " << neq->print_prefix() << " -- " 
         //     << *(*nj) << " " << outs_depth << endl;
         if (ins || outs){
            //if (ins) cout << "INS " << *ins << endl;
            //if (outs) cout << "OUTS " << *outs << endl;
            get_intersect(&(*ni)->tlist, &(*nj)->tlist, ins, outs);
         }
         
         if (outs){
            stats.incrcand(iter-1);
            //cout << "OUTS " << *outs;
            if (notfrequent(*outs)) delete outs;
            else if (num_patterns < max_num_patterns){
               ++num_patterns;
               neq->add_node(outs);
               stats.incrlarge(iter-1);
            }
         }
         if (ins){
            // cout << "INS " << *ins;
            stats.incrcand(iter-1);
            if (notfrequent(*ins)) delete ins;
            else if (num_patterns < max_num_patterns){
               neq->add_node(ins);
               ++num_patterns;
               stats.incrlarge(iter-1);
            }
         }
      }
      stats.incrtime(iter-1, tt.Stop());
      if (!neq->nlist().empty()){
         if (output) cout << *neq;
         if (prune_type == prune) FK.add(iter,neq);

         if (num_patterns < max_num_patterns)   
            enumerate_freq(neq, iter+1);
      }
      delete neq;
   }
}

void form_f2_lists(Eqclass *eq)
{
   list<Eqnode *>::iterator ni;
   idlist *l1, *l2;
   Eqnode *ins=NULL, *outs=NULL;
   int pit, nit;
   TimeTracker tt;
   
   tt.Start();
   pit = eq->prefix()[0];
   l1 = DCB->Idlists[pit];
   //cout << "IDLIST1 " << *l1;
   for (ni=eq->nlist().begin(); ni != eq->nlist().end(); ++ni){
      nit = (*ni)->val;
      l2 = DCB->Idlists[nit];
      ins = (*ni);
      //cout << "IDLIST2 " << *l2;
      //cout << "LISTS " << pit << " " << nit << " " << l1->size() 
      //     << " " << l2->size() << " " << ins->tlist.size() << endl;
      get_intersect(l1, l2, ins, outs);
      //cout << "f2prefix " << eq->prefix() << endl;
      //cout << "f2 " << *ins;
   }
   stats.incrtime(1,tt.Stop());
   //cout << *eq << endl;
}

void get_Fk(list<Eqclass *> &F2list){
   Eqclass *eq;

   while(!F2list.empty()){
      eq = F2list.front();
      form_f2_lists(eq);
      //cout << *eq;
      if (prune_type == prune) FK.add(2, eq);
      switch(alg_type){
      case treeminer:
         if (num_patterns < max_num_patterns)
            enumerate_freq(eq, 3); 
         break;
      case maxtreeminer:
         cout << "NOT IMPLEMENTED\n";
         break;
      }
      delete eq;
      F2list.pop_front();
   }
}

int main(int argc, char **argv)
{
   int i;
   TimeTracker tt;
   tt.Start(); 
   parse_args(argc, argv); 
  
   DCB = new Dbase_Ctrl_Blk(infile); 
   get_F1();
   list<Eqclass *> *F2list = get_F2();

   //DCB->print_vertical();
   get_Fk(*F2list);

   //for (i=2; i < stats.size(); i++){
   //   cout << "F" << i+1 << " - ";
   //   cout << stats[i].numlarge << " " << stats[i].numcand << endl;
   // }
   
   double tottime = tt.Stop();
   stats.tottime = tottime;

   cout << stats << endl;
   
   cout << "TIME = " << tottime << endl;

   //write results to summary file
   ofstream summary("summary.out", ios::app); 
   summary << "VTREEMINER ";
   switch(sort_type){
   case incr: summary << "INCR "; break;
   case decr: summary << "DECR "; break;
   default: break;
   }
   switch(prune_type){
   case prune: summary << "PRUNE "; break;
   deafult: break;
   }
   if (!count_unique) summary << "MULTIPLE ";
   if (use_fullpath) summary << "FULLPATH ";

   summary << infile << " " << DBASE_NUM_TRANS << " " << Num_Classes;
   for (i=0; i < Num_Classes; ++i)
      summary << " " << MINSUP_PER[i];
   for (i=0; i < Num_Classes; ++i)
      summary << " " << MINSUPPORT[i];
   summary << " " << stats << endl;
   summary.close();
   
   exit(0);
}




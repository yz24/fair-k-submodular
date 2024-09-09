#include <iostream>
#include <fstream>
#include <random>
#include <utility>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <time.h> // clock
#include <math.h>
#include <memory> // auto_ptr
#include "jlog.h"
#include "tools.h"
#include "mt19937ar.c"
#include <numeric>
#include <map>

using namespace std;
using namespace jlog_internal;
typedef long long LL;
typedef unsigned long long ULL;

vector<bool> X;
vector<bool> active;
int n;
int k;

// [z][u][i] = <v,p>
double simulate(vector<vector<vector<pair<int, double> > > > &es, vector<pair<int, int> > &S, int R, double delta) {
	Xorshift xs(0);
	LL sum = 0;
	for (int t = 0; t < R; t++) {
		for (int z = 0; z < k; z++) {
			vector<int> tmp;
			queue<int> Q;
			for (int i = 0; i < S.size(); i++) {
				if (S[i].second == z) {
					Q.push(S[i].first);
					X[S[i].first] = true;
				}
			}
			for (; !Q.empty();) {
				int u = Q.front();
				Q.pop();
				active[u] = true;
				tmp.push_back(u);
				for (int i = 0; i < es[z][u].size(); i++) {
					int v = es[z][u][i].first;
					double p = es[z][u][i].second;
					if (!X[v] && xs.nextDouble() < p) {
						X[v] = true;
						Q.push(v);
					}
				}
			}
			for (int i = 0; i < tmp.size(); i++) {
				X[tmp[i]] = false;
			}
		}
		int n1 = 0;
		for (int v = 0; v < n; v++) {
			if (active[v]) {
				n1++;
				active[v] = false;
			}
		}
		sum += n1;
	}
    double inf = 1.0 * sum / R;

    int random = rand()%2000;
    double div = inf * delta * 1000 / (random-1000);
    double inf1 = inf + div;

	return inf1;
}


void lazy_greedy( vector<vector<vector<pair<int, double> > > > &es, int budget, int beta, double delta) {
    cout<<"\n----------TS greedy with B="<<budget<<"----------"<<endl;
    bool TS = true;
    // < gain, < <node,item>, tick > >
    priority_queue<pair<double, pair<pair<int, int>, int> > > que;
    for (int v = 0; v < n; v++) {
        for (int z = 0; z < k; z++) {
            que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
        }
    }
    if(!TS){
        vector<int> budgets(k);
//    int B = 0;
        for (int j = 0; j < budgets.size(); j++) {
            budgets[j] = budget;
        }
    }


    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;
    vector<bool> used(n);
    vector<pair<int, int> > S;
    for (int j = 0; j < budget; j++) {
        printf(INFO "|S| = %d" DEF, j+1);
        double gain;
        pair<int, int> next;
        double start = get_current_time_sec();
        for (int num = 0;; num++) {
            double now = get_current_time_sec();

            pair<double, pair<pair<int, int>, int> > pp = que.top();
            que.pop();
            pair<int, int> s = pp.second.first;
            int last = pp.second.second;
//            if (used[s.first] || (!TS && budgets[s.second] == 0)) {
            if (used[s.first]) {
                continue;
            }

            if (last == j) {
                next = s;
                gain = pp.first;
                JLOG_ADD("num-of-eval", num);
                num_evaluations += num;
                break;
            }
            vector<pair<int, int> > SS(S);
            SS.push_back(s);
            candidate_f_value = simulate(es, SS, beta, delta);

            que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

        }
        current_f_value = candidate_f_value;
        S.push_back(next);
        used[next.first] = true;
//        if(!TS) budgets[next.second]--;
        if((j+1)%k ==0){
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;
            cout<<"seed set: ";
            vector<int> stats(k, 0);
            for(auto s: S){
                stats[s.second]++;
                cout<<"<"<<s.first<<", "<<s.second<<">, ";
            }
            cout<<endl;
            for(auto st: stats){
                cout<< st<<" ";
            }
            cout<<endl;
            double f_value = simulate(es, S, 500, 0);
            cout<<f_value<<endl;
        }

    }

    cout<<"seed set: ";
    vector<int> stats(k, 0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value = simulate(es, S, 500, 0);
    cout<<f_value<<endl;

}


void lazy_greedy_IS( vector<vector<vector<pair<int, double> > > > &es, vector<int> uppers, int budget, int beta, double delta) {
    cout<<"\n----------IS greedy with B="<<budget<<"----------"<<endl;
    bool TS = false;
    // < gain, < <node,item>, tick > >
    priority_queue<pair<double, pair<pair<int, int>, int> > > que;
    for (int v = 0; v < n; v++) {
        for (int z = 0; z < k; z++) {
            que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
        }
    }
    vector<int> budgets = std::move(uppers);
//    if(!TS){
//        vector<int> budgets(k);
////    int B = 0;
//        for (int j = 0; j < budgets.size(); j++) {
//            budgets[j] = budget;
//        }
//    }


    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;
    vector<bool> used(n);
    vector<pair<int, int> > S;
    for (int j = 0; j < budget; j++) {
//        printf(INFO "|S| = %d" DEF, j+1);
        double gain;
        pair<int, int> next;
        double start = get_current_time_sec();
        for (int num = 0;; num++) {
            double now = get_current_time_sec();

            pair<double, pair<pair<int, int>, int> > pp = que.top();
            que.pop();
            pair<int, int> s = pp.second.first;
            int last = pp.second.second;
//            if (used[s.first] || (!TS && budgets[s.second] == 0)) {
            if (used[s.first] || budgets[s.second]<=0) {
                continue;
            }

            if (last == j) {
                next = s;
                gain = pp.first;
//                JLOG_ADD("num-of-eval", num);
                num_evaluations += num;
                break;
            }
            vector<pair<int, int> > SS(S);
            SS.push_back(s);
            candidate_f_value = simulate(es, SS, beta, delta);

            que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

        }
        current_f_value = candidate_f_value;
        S.push_back(next);
        used[next.first] = true;
        budgets[next.second]--;
        if((j+1)%k ==0){
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;
            cout<<"seed set: ";
            vector<int> stats(k, 0);
            for(auto s: S){
                stats[s.second]++;
                cout<<"<"<<s.first<<", "<<s.second<<">, ";
            }
            cout<<endl;
            for(auto st: stats){
                cout<< st<<" ";
            }
            cout<<endl;
            double f_value = simulate(es, S, 500, 0);
            cout<<f_value<<endl;
        }

    }

    cout<<"seed set: ";
    vector<int> stats(k, 0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value = simulate(es, S, 500, 0);
    cout<<f_value<<endl;

}


bool extendable(vector<pair<int, int> > S, pair<int, int> s, vector<int> lowers, vector<int> uppers, int B){
    vector<int> supp(k, 0);
    supp[s.second]++;
    for(auto p : S){
        supp[p.second]++;
    }
    // check upper bounds, only check the s's type
    if(supp[s.second] > uppers[s.second]){
        return false;
    }
    // check lower bounds.
    int sum = 0;
    for(int i = 0; i<k; i++){
        sum += max(lowers[i], supp[i]);
    }
    if(sum > B){
        return false;
    }
    return true;
}


// [u][i] = <v, p[1..k]>
void fair_greedy(vector<vector<vector<pair<int, double> > > > &es, vector<int> lowers, vector<int> uppers, int B, int beta, double delta) {
    cout<<"\n---------- fair greedy with B="<<B<<"----------"<<endl;
	// < gain, < <node,item>, tick > >
	priority_queue<pair<double, pair<pair<int, int>, int> > > que;
	for (int v = 0; v < n; v++) {
		for (int z = 0; z < k; z++) {
			que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
		}
	}
//    vector<int> budgets(k);
//    int B = 0;
//    for (int j = 0; j < budgets.size(); j++) {
//        budgets[j] = budget;
//    }

    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;

	vector<bool> used(n);
	vector<pair<int, int> > S;
	for (int j = 0; j < B; j++) {
        pair<double, pair<pair<int, int>, int> > pp = que.top();
//        if(pp.first < 0){
//            break;
//        }
//		printf(INFO "|S| = %d" DEF, j+1);
		double gain;
		pair<int, int> next;
		for (int num = 0;; num++) {

			pair<double, pair<pair<int, int>, int> > pp = que.top();
			que.pop();
			pair<int, int> s = pp.second.first;
			int last = pp.second.second;

			if (used[s.first] || !extendable(S, s, lowers, uppers, B)) {
//                if(budgets[s.second] == 0)
//                    cout<<"type of "<<s.second<<" used up."<<endl;
				continue;
			}

			if (last == j) {
				next = s;
				gain = pp.first;
				JLOG_ADD("num-of-eval", num);
                num_evaluations += num;
				break;
			}
			vector<pair<int, int> > SS(S);
			SS.push_back(s);
//			double sigma = simulate(n, es, k, SS, beta);
//			double psigma = simulate(n, es, k, S, beta);
			candidate_f_value = simulate(es, SS, beta, delta);

			que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

		}
        current_f_value = candidate_f_value;
		S.push_back(next);
		used[next.first] = true;
//        if(!TS) budgets[next.second]--;

        if((j+1)%k==0){
            printf(INFO "|S| = %d" DEF, j+1);
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;
        }
	}
    cout<<"seed set: ";
    vector<int> stats(k,0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value = simulate(es, S, 500, 0);
    cout<<f_value<<endl;

}

// [u][i] = <v, p[1..k]>
void fair_threshold(vector<vector<vector<pair<int, double> > > > &es, vector<int> lowers, vector<int> uppers, int B, int beta, double eps, double delta) {
    cout<<"\n---------- fair threshold with eps="<<eps<<", B="<<B<<"----------"<<endl;
    // < gain, < <node,item>, tick > >
    map<pair<int, int>, int> evaluations;

    priority_queue<pair<double, pair<int, int>>> que;
    for (int v = 0; v < n; v++) {
        for (int z = 0; z < k; z++) {
            que.push(make_pair(1e12, make_pair(v, z)));
            evaluations[{v, z}] = 0;
        }
    }

    double current_f_value = 0.0;
    double candidate_f_value = 0.0;
    int num_evaluations = 0;

    vector<bool> used(n);
    vector<pair<int, int> > S;
    for (int j = 0; j < B; j++) {

        if(que.empty()){
            break;
        }

        pair<int, int> next;
        double gain;
        int num = 0;
        while (true) {
            pair<double, pair<int, int> > pp = que.top();
            que.pop();
            pair<int, int> s = pp.second;

            if (used[s.first] || !extendable(S, s, lowers, uppers, B)) {
                continue;
            }

            vector<pair<int, int> > SS(S);
            SS.push_back(s);
            candidate_f_value = simulate(es, SS, beta, delta);
            num ++;
            evaluations[s] ++;
            gain = candidate_f_value - current_f_value;

            if(gain >= que.top().first* (1-eps) ){
                S.push_back(s);
                next = s;
                used[s.first] = true;
                num_evaluations += num;

                break;
            }

            if(evaluations[s] <= int(log(B/(2*eps)))/eps){
                que.push(make_pair(gain, s));
            }
//          que.push(make_pair(candidate_f_value - current_f_value, make_pair(s, j)));

        }
        current_f_value = candidate_f_value;
        if((j+1)%k==0){
            printf(INFO "|S| = %d" DEF, j+1);
            JLOG_ADD("seed", next.first);
            JLOG_ADD("item", next.second);
            JLOG_ADD("gain", gain);
            JLOG_ADD("objective", current_f_value);
            JLOG_ADD("num-of-evaluation as of now", num_evaluations);
            cout<<endl;
        }

    }
    cout<<"seed set: ";
    vector<int> stats(k, 0);
    for(auto s: S){
        stats[s.second]++;
        cout<<"<"<<s.first<<", "<<s.second<<">, ";
    }
    cout<<endl;
    for(auto st: stats){
        cout<< st<<" ";
    }
    cout<<endl;
    double f_value = simulate(es, S, 500, 0);
    cout<<f_value<<endl;

}


vector<pair<int, int>> generate_random(int B, vector<int> lowers, vector<int> uppers){
//    vector<pair<int, int>> results(k);
    set<pair<int, int>> uniquePairs;
    vector<int> vec(n);
    vector<int> types(k);
    std::iota (std::begin(vec), std::end(vec), 0);
    std::iota (std::begin(types), std::end(types), 0);
    vector<int> type_count(k, 0);
    // Generate B unique pairs

    for(int i=0; i<k; i++){
        while(type_count[i]<lowers[i]){
            int randomVec = rand() % vec.size();
            if(uniquePairs.count({vec[randomVec], i}) > 0){
                continue;
            }
            uniquePairs.insert({vec[randomVec], i});
//            cout<<uniquePairs.size()<<": "<<vec[randomVec]<<" "<< i<<endl;

            type_count[i]++;
        }
    }

    int remaining = B - std::accumulate(lowers.begin(), lowers.end(), 0);
    while(remaining > 0){
        int randomVec = rand() % vec.size();
        int randomType = rand() % types.size();
        if(uniquePairs.count({vec[randomVec], types[randomType]}) > 0){
            continue;
        }
        uniquePairs.insert({vec[randomVec], types[randomType]});
//        cout<<uniquePairs.size()<<": "<<vec[randomVec]<<" "<< types[randomType]<<endl;
        remaining --;
        type_count[randomType] ++;
    }

    vector<pair<int, int>> results(uniquePairs.begin(), uniquePairs.end());

    return results;
}

vector<pair<int, int>> generate_random(int B){
//    vector<pair<int, int>> results(k);
    set<pair<int, int>> uniquePairs;
    vector<int> vec(n);
    vector<int> types(k);
    std::iota (std::begin(vec), std::end(vec), 0);
    std::iota (std::begin(types), std::end(types), 0);

    int c = 0;
    // Generate B unique pairs
    while (uniquePairs.size() < B) {
        int randomVec = rand() % vec.size();
        int randomType = rand() % types.size();
        uniquePairs.insert({vec[randomVec], types[randomType]});
//        cout<<c<<": "<<uniquePairs.size()<<", "<<vec[randomVec]<<" "<< types[randomType]<<endl;
        c++;
    }
    vector<pair<int, int>> results(uniquePairs.begin(), uniquePairs.end());

    return results;
}




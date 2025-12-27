#include <bits/stdc++.h>
#include <time.h>
#include <random>
#define N 100010
#define M 10010
using namespace std;
clock_t start, endt;
mt19937 myrand(time(0));
uniform_real_distribution<> dis(0.0, 1.0);
int n,m,c;
int h[M];
int hed[31], tal[31], tail[M];
int lst[N], nxt[N];
int Msg[N], User[N], Exe[N];
int Dedl[N];
int ExeTo[M], Using[31];
vector<vector<int>> info(M);
vector<int> realu;
int ins[N], Hd[31];
int lt[N], rt[N], EXE[N];
void Insert(int x, int y, int j){//insert task x behind y at core j
    int z;
    if(y) z = nxt[y], nxt[y] = x;
    else z = hed[j], hed[j] = x;
    if(z) lst[z] = x;
    else tal[j] = x;
    lst[x] = y, nxt[x] = z;
    
    if(Msg[x] == Msg[nxt[x]]){
        int tm = Exe[x];
        while(Msg[x] == Msg[nxt[x]]) x = nxt[x];
        EXE[x] += tm;
    }else{
        //multi-msg
        if(Msg[x] != Msg[y]){
            if(y && Msg[y] == Msg[z]){
                int u = y, v = z;
                int tm1 = 0, tm2 = 0;
                while(Msg[u] == Msg[lst[u]]) tm1 += Exe[u], u = lst[u];
                while(Msg[v] == Msg[nxt[v]]) tm2 += Exe[v], v = nxt[v];
                tm1 += Exe[u], tm2 += Exe[v];
                EXE[y] = tm1, EXE[v] = tm2;
                lt[x] = y, rt[x] = v, lt[y] = lt[v], lt[v] = x, rt[y] = x, EXE[x] = Exe[x];
                if(!lt[y]) Hd[j] = y;
                else rt[lt[y]] = y;
            }else{
                if(y) z = rt[y], rt[y] = x;
                else z = Hd[j], Hd[j] = x;
                if(z) lt[z] = x;
                lt[x] = y, rt[x] = z, EXE[x] = Exe[x];
            }
        }else{
            z = rt[y];
            if(z) lt[z] = x;
            if(!lt[y]) Hd[j] = x;
            else rt[lt[y]] = x;
            lt[x] = lt[y], rt[x] = z, lt[y] = rt[y] = 0, EXE[x] = Exe[x] + EXE[y], EXE[y] = 0;
        }
    }
}
void Erase(int x, int j){//erase task x from core j
    int y = lst[x], z = nxt[x];
    if(y) nxt[y] = z;
    else hed[j] = z;
    if(z) lst[z] = y;
    else tal[j] = y;
    if(Msg[x] != Msg[z]){
        z = rt[x];
        if(Msg[x] == Msg[y]){
            if(z) lt[z] = y;
            if(!lt[x]) Hd[j] = y;
            else rt[lt[x]] = y;
            rt[y] = z, lt[y] = lt[x], lt[x] = rt[x] = 0, EXE[y] = EXE[x] - Exe[x], EXE[x] = 0;
        }else{
            if(y && Msg[y] == Msg[z]){
                EXE[z] += EXE[y];
                EXE[y] = EXE[x] = 0;
                lt[z] = lt[y];
                if(lt[z]) rt[lt[z]] = z;
                else Hd[j] = z;
                lt[x] = rt[x] = 0, lt[y] = rt[y] = 0;
            }else{
                if(z) lt[z] = y;
                if(!y) Hd[j] = z;
                else rt[y] = z;
                lt[x] = rt[x] = 0, EXE[x] = 0;
            }
        }
    }else{
        int p = x;
        while(Msg[p] == Msg[nxt[p]]) p = nxt[p];
        EXE[p] -= Exe[x];
    }
    lst[x] = nxt[x] = 0;
}
int S1[31], S2[31];
int calc1(int j){//calculate increasing points
    if(j == 0) return 0;
    int Tim = 0, now = hed[j], sc = 0;
    while(now != 0){
        Tim += Exe[now];
        if(Tim <= Dedl[now]) sc++;
        now = nxt[now];
    }
    return sc;
}
int Calc1(int j){
    int now = hed[j], sc = 0;
    while(now != 0){
        if(Msg[now] == Msg[lst[now]]) sc++;
        now = nxt[now];
    }
    return sc;    
}

int rec[N], dat[N];
int calc2(int user, int j){//strategy of SA
    int sc1 = 0;
    for(auto i : info[user]){
        rec[i] = lst[i], dat[i] = ins[i] = 0;
        if(Msg[i] == Msg[lst[i]] || Msg[i] == Msg[nxt[i]]) sc1--;
        else if(lst[i] && Msg[lst[i]] == Msg[nxt[i]]) sc1++;
    }
    for(auto i : info[user]){
        Erase(i, j);
    }
    
    int Tim = Using[j], now = tal[j];
    int l = info[user].size() - 1;
    for(; l >= 0; l--){
        int i = info[user][l], y = now, curtm = Tim, antm = Tim;
        while(now != 0){
            //if(Tim <= Dedl[i] && !y) y = now, antm = Tim;
            if(Msg[now] == Msg[i]){
                if(!dat[i] || Tim <= Dedl[i]) dat[i] = now, curtm = Tim;
                if(Tim <= Dedl[i]){
                    break;
                }
            }
            Tim -= EXE[now], now = lt[now];
        }
        if(!dat[i] && l > 0) ins[i] = info[user][l-1], now = y, Tim = antm;
        if(dat[i]) ins[i] = dat[i], Tim = curtm - EXE[dat[i]], now = lt[dat[i]];
        Tim -= Exe[i];
    }
    
    for(auto i : info[user]){
        Insert(i, ins[i], j);
    }
    for(auto i : info[user]){
        if(Msg[i] == Msg[lst[i]] || Msg[i] == Msg[nxt[i]]) sc1++;
    }
    return sc1;
}
int Calc2(int user, int j, int k){
    int sc1 = 0;
    for(auto i : info[user]){
        rec[i] = lst[i], dat[i] = ins[i] = 0;
        if(Msg[i] == Msg[lst[i]] || Msg[i] == Msg[nxt[i]]) sc1--;
        else if(lst[i] && Msg[lst[i]] == Msg[nxt[i]]) sc1++;
    }

    for(auto i : info[user]){
        Erase(i, j);
    }
    if(j) Using[j] -= ExeTo[user];
    if(k) Using[k] += ExeTo[user];
    int Tim = Using[k], now = tal[k];
    
    int l = info[user].size() - 1;
    for(; l >= 0; l--){
        int i = info[user][l], y = now, curtm = Tim, antm = Tim;
        while(now != 0){
            //if(Tim <= Dedl[i] && !y) y = now, antm = Tim;
            if(Msg[now] == Msg[i]){
                if(!dat[i] || Tim <= Dedl[i]) dat[i] = now, curtm = Tim;
                if(Tim <= Dedl[i]){
                    break;
                }
            }
            Tim -= EXE[now], now = lt[now];
        }
        if(!dat[i] && l > 0) ins[i] = info[user][l-1], now = y, Tim = antm;
        if(dat[i]) ins[i] = dat[i], Tim = curtm - EXE[dat[i]], now = lt[dat[i]];
        Tim -= Exe[i];
    }
    
    for(auto i : info[user]){
        Insert(i, ins[i], k);
    }
    for(auto i : info[user]){
        if(Msg[i] == Msg[lst[i]] || Msg[i] == Msg[nxt[i]]) sc1++;
    }
    return sc1;
}
//i --- z - k --- j --> k --- j - i ---- z
void SA(){//SA for adjusting one user in one core
    while(true){
        int rept = min(1e4, 1e5 * m / n);
            
        for(int l = 1; l <= rept; l++){
            endt = clock();
            if((double)(endt - start) / CLOCKS_PER_SEC > 3.9) return;
            int x = myrand() % realu.size();
            x = realu[x];
            int j = h[x];
            int sc, SC, sc1, score;
            sc1 = calc2(x, j);
            sc = S1[j];
            SC = calc1(j);
            score = SC - sc + sc1;
            bool ok = (score >= 0);
            if(ok){
                S1[j] = SC;
            }else{
                for(auto i : info[x]){
                    Erase(i, j);
                }
                for(auto i : info[x]){
                    Insert(i, rec[i], j);
                }
            }
        }
    }
}

int rl;

double urg;
void SAA(){//SA for adjusting one user to another core
    double tm = 3.0;
    if(m <= 10) tm = 2.5;
    if(m <= 5) tm = 2.0;
    while(true){
        endt = clock();
        if((double)(endt - start) / CLOCKS_PER_SEC > tm) return;
        int rept = min(1e4, 1e5 * m / n);
        for(int l = 1; l <= rept; l++){
            int x = myrand() % realu.size();
            x = realu[x];
            int j = h[x], k = myrand() % (m+1);
            if(urg > 1.01 + 0.06 * log(1.0 * n * m)) k = myrand() % m + 1;
            int sc, sc1, sc2, SC, score, sc3;
            sc = S1[j];
            if(j != k) sc += S1[k];
            sc3 = Calc2(x, j, k);
            SC = sc1 = calc1(j);
            if(j != k) sc2 = calc1(k), SC += sc2;
            score = SC - sc + sc3;
            if(score >= 0){
                S1[j] = sc1;
                if(j != k) S1[k] = sc2, h[x] = k;
            }else{
                for(auto i : info[x]){
                    Erase(i, k);
                }
                for(auto i : info[x]){
                    Insert(i, rec[i], j);
                }
                if(k) Using[k] -= ExeTo[x];
                if(j) Using[j] += ExeTo[x];
            }
        }
        
    }
}
//cos((pi / 2) / (1 + 2x^2))
//cout << 1.0 / sqrt(1 + 1400 * pow(1.0 * n / m, 0.5) / pow(msg, 1.1) / pow(max(0.001, (1.0 * n / sz - 0.999)), 1.45));
#define pi 3.14159265
void init(){
    int mx = 0, te = 0, msg = 0;
    for(int i = 1; i <= n; i++) msg = max(msg, Msg[i]), mx = max(mx, Dedl[i]), te += Exe[i];
    double urgnt = 1.0 * mx * m / te;//urgency: from 0.1 to 4.0
    urg = urgnt;
    //calculate capability score and affinity score with different amount of cores, you can only use part of the cores, sacrifice some capability score to purchase higher affinity score
    double eff = cos(pi / 2 / (1 + 2.0 * urg * urg));//cosine function
    double aff = 1 - 1.0 / sqrt(1 + 5000 * pow(n, 0.5) / pow(m, 0.9) / pow(msg, 1.1) / pow(max(0.001, (1.0 * n / realu.size() - 0.999)), 1.45));//fitting function
    double Mx = eff + aff;
    while(m > 1){
        double Urg = 1.0 * mx * (m-1) / te;
        double Eff = cos(pi / 2 / (1 + 2.0 * Urg * Urg));
        double Aff = 1 - 1.0 / sqrt(1 + 5000 * pow(n, 0.5) / pow(m-1, 0.9) / pow(msg, 1.1) / pow(max(0.001, (1.0 * n / realu.size() - 0.999)), 1.45));
        if(Eff + Aff > Mx){
            Mx = Eff + Aff;
            eff = Eff, aff = Aff;
        }else break;
        m--;
    }
    sort(realu.begin(), realu.end(), [&](int x, int y){
        return 1.0 * ExeTo[x] / info[x].size() > 1.0 * ExeTo[y] / info[y].size();
    });
    //you can sacrifice part of users to purchase higher capability score
    for(auto x : realu){
        if(urgnt > 1.1 + 0.019 * log(1.0 * n * m)) break;
        te -= ExeTo[x], h[x] = -1;
        urgnt = 1.0 * mx * m / te;
    }
}

int Read(){
    char ch = getchar();
    while(ch < '0' || ch > '9') ch = getchar();
    int s = 0;
    while(ch >= '0' && ch <= '9'){
        s = (s * 10 + ch - '0');
        ch = getchar();
    }
    return s;
}

void input(){

    n = Read(), m = Read(), c = Read();
    rl = m;
    for(int i = 1; i <= n; i++){
        int msg, user, exe, dedl;
        msg = Read(), user = Read(), exe = Read(), dedl = Read();
        Msg[i] = msg, User[i] = user, Exe[i] = exe, Dedl[i] = dedl;
        ExeTo[user] += exe;
        info[user].push_back(i);
    }
    for(int i = 1; i <= 10000; i++){
        if(!info[i].empty()){
            realu.push_back(i);
        }
    }
    init();
}

void Print(int x){
    if(x > 9) Print(x / 10);
    putchar(char('0' + x % 10));
}

int solve(bool reborn){
    vector<vector<int>> task(m+1);
    auto calc = [&](){
        
        for(int i = 1; i <= n; i++){
            int msg, user, exe, dedl;
            msg = Msg[i], user = User[i], exe = Exe[i], dedl = Dedl[i];
            int j = h[user];
            int res = 0, now = tail[user] ? tail[user] : Hd[j];
            while(now && Msg[now] == Msg[nxt[now]]) now = nxt[now];
            while(now != 0){
                if(Msg[now] == msg){
                    res = now; break;
                }
                now = rt[now];
            }
            if(res){
                Insert(i, res, j);
            }else{
                Insert(i, tal[j], j);
            }
            tail[user] = i;
        }
        
        int r1 = 0;
        for(int j = 0; j <= m; j++){
            S1[j] = calc1(j), S2[j] = Calc1(j);
            r1 += S1[j] + S2[j];
        }
        return r1;
    };
    //greedy strategy
    if(m > 1){
        for(int i = 1; i <= n; i++){
            int msg, user, exe, dedl;
            msg = Msg[i], user = User[i], exe = Exe[i], dedl = Dedl[i];
            if(!h[user]){
                vector<tuple<double,int,int>> v;
                for(int j = 1; j <= m; j++){
                    int l = 0, k = 0;
                    while(k < task[j].size() && l < info[user].size()){
                        while(k < task[j].size()){
                            int now = task[j][k];
                            if(Msg[now] == Msg[info[user][l]]){
                                k++, l++;
                                break;
                            }
                            k++;
                        }
                    }
                    double s = l, timepoint = Using[j];
                    for(l = 0; l < info[user].size(); l++){
                        timepoint += Exe[info[user][l]];
                        if(Dedl[info[user][l]] >= timepoint) s += 0.95 + n / 1e5;
                    }
                    v.push_back({-s, Using[j], j});
                }
                sort(v.begin(), v.end());
                int j = get<2>(v.front());
                Using[j] += ExeTo[user], h[user] = j;
            }
            if(h[user] > 0) task[h[user]].push_back(i);
        }
    }else{
        for(auto x : realu) if(!h[x]) h[x] = 1;
    }
    for(int j = 0; j <= m; j++) Using[j] = 0;
    for(auto x : realu) if(h[x] > 0) Using[h[x]] += ExeTo[x];
    for(auto x : realu) if(h[x] == -1) h[x] = 0;
    
    calc();
    vector<vector<int>> ans(rl+1);
    for(int j = 0; j <= m; j++){
        int x = hed[j];
        while(x) ans[j].push_back(x), x = nxt[x];
    }
    SAA();
    SA();
    for(int j = 0; j <= rl; j++){
        int x = hed[j];
        if(j == 0){
            while(x){
                lt[x] = rt[x] = EXE[x] = 0;
                int y = nxt[x];
                Insert(x, tal[1], 1), x = y;
            }
        }
    }
    /*int r1 = 0;
    for(int j = 1; j <= rl; j++){
        S1[j] = calc1(j), S2[j] = Calc1(j);
        r1 += S1[j] + S2[j];
    }
    cout << r1 << endl;
    return 0;
    */for(int j = 1; j <= rl; j++){
        int now = hed[j];
        int count = 0;
        while(now){
            count ++, now = nxt[now];
        }
        Print(count), putchar(' ');
        now = hed[j];
        while(now){
            Print(Msg[now]), putchar(' '), Print(User[now]), putchar(' ');
            now = nxt[now];
        }
        putchar('\n');
    }
    return 0;
}
signed main()
{
    //freopen("data.txt", "r", stdin);
    
    start = clock();
    input();
    solve(false);
    
    return 0;
}


一般情况下，解决概率问题需要顺序循环，
而解决期望问题使用逆序循环，
如果定义的状态转移方程存在后效性问题，还需要用到 高斯消元 来优化。
概率 DP 也会结合其他知识进行考察，例如 状态压缩，树上进行 DP 转移等。

DP 求概率
这类题目采用顺推，也就是从初始状态推向结果。
同一般的 DP 类似的，难点依然是对状态转移方程的刻画，
只是这类题目经过了概率论知识的包装。
1. CodeForces 148 D Bag of mice
2. POJ3071 Football
3. CodeForces 768 D Jon and Orbs
DP 求期望题目在对具体是求一个值或是最优化问题上会对方程得到转移方式有一些影响，
但无论是 DP 求概率还是 DP 求期望，总是离不开概率知识和列出、化简计算公式的步骤，
在写状态转移方程时需要思考的细节也类似。
1. POJ2096 Collecting Bugs
2. HDU3853 LOOPS
3. HDU4035 Maze
4. 「NOIP2016」换教室
5. 「SCOI2008」奖励关


概率DP
https://ac.nowcoder.com/acm/contest/32708/G
http://poj.org/problem?id=3071


int n;
double p[1024][1024];
double dp[1024][1024]; //i轮后j队赢的概率

int main()
{
    IO;
    cout << fixed << setprecision(10);
    while(true)
    {
        cin >> n; if(n == -1) break;
        int len = 1 << n;
        for(int i = 1; i <= len;++i) {
            for(int j = 1; j <= len;++j) {
                cin >> p[i][j];
            }
        }

        memset(dp, 0, sizeof dp);
        //init
        for(int i = 1; i <= len;++i) dp[0][i] = 1.0;

        for(int i = 1; i <= n;++i) {
            for(int j = 1; j <= len;++j) {
                for(int k = 1; k <= len;++k) {
                    int x = ((j - 1) >> (i - 1)) ^ 1;
                    int y = (k - 1) >> (i - 1);
                    if(x == y) {
                        dp[i][j] += 1.0 * dp[i - 1][j] * dp[i - 1][k] * p[j][k];
                    }
                }
            }
        }    
        
        double res = 0.0;
        int res_i = 0;
        for(int i = 1; i <= len;++i) {
            //cout << dp[n][i] << '\n';
            if(dp[n][i] > res) {
                res = dp[n][i];
                res_i = i;
            }
        }
        cout << res_i << '\n';
    }

    return 0;
}


https://codeforces.com/problemset/problem/768/D
#include <bits/stdc++.h>

using namespace std;
#define IO ios::sync_with_stdio(false);cin.tie(0),cout.tie(0);

const int N = 1e4 + 10;
int k, q;
double p[N];
double dp[N][1024];
//第i天产生j钟龙晶的概率

int main()
{
    IO;
    cin >> k >> q;
    for(int i = 1; i <= q;++i) cin >> p[i];
    //init
    for(int i = 0;i <= k;++i) dp[0][i] = (i == 0)? 1.0 : 0;
    for(int i = 1; i < N;++i) {
        for(int j = 1; j <= k;++j) {
            dp[i][j] += 1. * dp[i - 1][j] * (1. * j / k) +
            dp[i - 1][j - 1] * (1. * (k - (j - 1)) / k);
        }
        
    }
    for(int j = 1; j <= q;++j) {
        for(int i = 1; i < N;++i) {
            if(dp[i][k] > 1. * p[j] / 2000) {
                cout << i << '\n';
                break; 
            }
        }    
    }

    return 0;
}

骰子
https://blog.csdn.net/qq_42754826/article/details/89043168

using i64 = long long;
i64 n, m;
i64 f[25][152];//投i次色子点数和为j的方案数

void solve()
{
    cin >> n >> m;

    for(int i = 1; i <= 6;++i) f[1][i] = 1;
    //f[0][0] = 1;
    for(int i = 2; i <= n;++i) {
        for(int j = 0;j <= n * 6;++j) {
            for(int k = 1; k <= 6;++k) {
                if(j >= k) f[i][j] += f[i - 1][j - k];
                //f[i + 1][j + k] += f[i][j];
            }
        }
    }
    //for(int j = 0; j <= n * 6 ;++j) cout << f[3][j] << ' ';

    i64 sum = pow(6, n);
    if(m == 0 || m == 1) cout << 1 << '\n';
    else if(m > 6 * n) cout << 0 << '\n';
    else {
        i64 res = 0;
        for(int i = m; i <= n * 6;++i) {
                res += f[n][i];
        }
        //cout << res << '\n';
        i64 g = __gcd(res, sum);
        cout << res / g << '/' << sum / g << '\n';
    } 
}


硬币/ Coins
https://atcoder.jp/contests/dp/tasks/dp_i

int n;
const int N = 3010;
double p[N];
double f[N][N];

void solve()
{
    cout << fixed << setprecision(10);

    cin >> n;
    for(int i = 1; i <= n;++i) cin >> p[i];

    f[1][0] = 1.0 - p[1];
    f[1][1] = p[1];
    for(int i = 2; i <= n;++i) {
        for(int j = 0; j <= i;++j) {
            f[i][j] += f[i - 1][j] * (1 - p[i]);
            if(j >= 1) f[i][j] += f[i - 1][j - 1] * p[i];
        }    
    }

    double res = 0.0;
    for(int i = 0; i <= n;++i) {
        //cout << f[n][i] << '\n';
        if(i > n / 2) res += f[n][i];
    }
    cout << res << '\n';
}


期望DP
Collecting Bugs
http://poj.org/problem?id=2096

#include <iostream>
#include <iomanip>
using namespace std;
#define IO ios::sync_with_stdio(false);cin.tie(0),cout.tie(0)

int n, s;
double dp[1010][1010];

int main()
{
    IO;
    cout << fixed << setprecision(4);

    cin >> n >> s;

    dp[n][s] = 0.0;
    for(int i = n; i >= 0;--i) {
        for(int j = s; j >= 0;--j) {
          if(i == n && j == s) continue;
          dp[i][j] += (1. * dp[i + 1][j + 1] * (n - i) * (s - j)
                    + dp[i + 1][j] * (n - i) * j 
                    + dp[i][j + 1] * i * (s - j) + n * s) / (n * s - i * j);  
        }
    }
    cout << dp[0][0] << '\n';
    return 0;
}

期望DP/记忆化搜索
https://atcoder.jp/contests/dp/tasks/dp_j
#include <iostream>
#include <iomanip>
#include <cstring>
using namespace std;
#define IO ios::sync_with_stdio(false);cin.tie(0),cout.tie(0)

int n;
const int N = 321;
int a[N];
double dp[N][N][N];//表示还剩下i,j,k个寿司这么多盘吃完所需要的期望值

double dfs(int c, int b, int a) {
    if(a == 0 && b == 0 && c == 0) {
        return 0.0;
    }
    if(dp[c][b][a] != 0.0) return dp[c][b][a];
    double cnt = 1.0 * n / (a + b + c);
    if(c) cnt += 1.0 * c / (a + b + c) * dfs(c - 1, b + 1, a);
    if(b) cnt += 1.0 * b / (a + b + c) * dfs(c, b - 1, a + 1);
    if(a) cnt += 1.0 * a / (a + b + c) * dfs(c, b, a - 1);
    return dp[c][b][a] = cnt; 
}

int main()
{
    IO;
    cout << fixed << setprecision(10);
    
    cin >> n;
    int one = 0, two = 0, three = 0;
    for(int i = 1; i <= n;++i) {
        cin >> a[i];
        if(a[i] == 1) ++one;
        else if(a[i] == 2) ++two;
        else if(a[i] == 3) ++three;
    }

    dfs(three, two, one);
    
    cout << dp[three][two][one] << '\n';
    return 0;
}
////////
网格DP
https://atcoder.jp/contests/dp/tasks/dp_y
#include <bits/stdc++.h>

using namespace std;
#define IO ios::sync_with_stdio(false);cin.tie(0),cout.tie(0)
using i64 = long long;
typedef long long ll;
i64 n, m, k;
const int M = 3010;
const int mod = 1e9 + 7;
pair<i64, i64> p[M];
i64 dp[M];

namespace CNM {//组合数板子
    const int N = 2e6 + 5;
    ll quick(ll x, ll n)
    {
        ll res = 1;
        while (n)
        {
            if (n & 1) res = (res*x) % mod;
            x = x * x%mod;
            n >>= 1;
        }
        return res;
    }
    ll inv(ll x) { return quick(x, mod - 2); }
    ll fac[N], invfac[N];
    void init()
    {
        fac[0] = 1;
        for (int i = 1; i < N; ++i) fac[i] = (fac[i - 1] * i) % mod;
        invfac[N - 1] = inv(fac[N - 1]);
        for (int i = N - 2; i >= 0; --i) invfac[i] = (invfac[i + 1] * (i + 1)) % mod;
    }
    ll C(int n, int m)
    {
        if (n < m || m < 0) return 0;
        return fac[n] * invfac[m] % mod*invfac[n - m] % mod;
    }
}


int main()
{
    IO;
    CNM::init();
    cin >> n >> m >> k;
    for(int i = 0; i < k;++i) {
        cin >> p[i].first >> p[i].second;
    }
    //排序
    sort(p, p + k);
    //计算从起点到达每个障碍的总方案数
    for(int i = 0; i < k;++i) {
        dp[i] = CNM::C(p[i].first + p[i].second - 2, p[i].first - 1);
        // cout << p[i].first + p[i].second - 2 << ' ' << p[i].first - 1 << '\n';
        // cout << dp[i] << '\n';
    }

    for(int i = 0; i < k;++i) {
        for(int j = 0; j < i;++j) {
            if(p[i].first >= p[j].first && p[i].second >= p[j].second) {
                i64 x = p[i].first - p[j].first + 1;
                i64 y = p[i].second - p[j].second + 1;
                dp[i] = dp[i] - dp[j] * CNM::C(x + y - 2, x - 1) % mod;
                dp[i] = (dp[i] + mod) % mod;
            }
        }
    }

    i64 ans = CNM::C(n + m - 2, n - 1);
    for(int i = 0; i < k;++i) {
        i64 x = n - p[i].first + 1;
        i64 y = m - p[i].second + 1;
        ans = ans - dp[i] * CNM::C(x + y - 2, x - 1) % mod;
        ans = (ans + mod) % mod;
    }

    cout << ans << '\n';

    return 0;
}

//博弈DP
https://atcoder.jp/contests/abc298/tasks/abc298_e

#include <bits/stdc++.h>

using namespace std;

int n, k;
bool dp[100010];
int a[110];

int main()
{
    ios::sync_with_stdio(false);cin.tie(0),cout.tie(0);

    cin >> n >> k;
    for(int i = 1; i <= n;++i) cin >> a[i];

    for(int i = 1; i <= k;++i) {
        bool ok = false;
        for(int j = 1; j <= n;++j) {
            if(a[j] <= i && dp[i - a[j]] == 0) {
                ok = true;
            }
        }
        dp[i] = ok;
        //转移下去全是先手必胜态则现在的状态必败
        //否则是先手必胜
    }

    if(dp[k]) cout << "First\n";
    else cout << "Second\n";

    return 0;
}

//记忆化搜索
#include <bits/stdc++.h>

using namespace std;

int n;
bool dp[100010][2];
int a[110];

int dfs(int k, int cnt) {
    if(dp[k][cnt] > 0) return dp[k][cnt];

    int res = 0;
    for(int i = 1; i <= n;++i) 
        if(a[i] <= k) res |= !dfs(k - a[i], !cnt);
    
    return dp[k][cnt] = res;
    
}

int main()
{
    ios::sync_with_stdio(false);cin.tie(0),cout.tie(0);
    int k;
    cin >> n >> k;
    for(int i = 1; i <= n;++i) cin >> a[i];

    if(dfs(k, 0)) cout << "First\n";
    else cout << "Second\n";

    return 0;
}

https://atcoder.jp/contests/abc298/tasks/abc298_e

const int mod = 998244353;
using i64 = long long;
i64 n, a, b, p, q;

void solve()
{
    cin >> n >> a >> b >> p >> q;
    auto inv_p = pow_mod(p, mod - 2, mod);
    auto inv_q = pow_mod(q, mod - 2, mod);

    vector<vector<i64> > first(n + 1, vector<i64>(n + 1));
    vector<vector<i64> > second(n + 1, vector<i64>(n + 1));

    for(int j = b; j < n; ++j) { //对于1和2 若1走到了n而2在b->n-1的时候，1赢的概率为1
        first[n][j] = second[n][j] = 1;
    }

    for(int i = n; i-- > a;) {
        for(int j = n; j-- > b; ) {
            for(int x = 1; x <= p;++x) {
                first[i][j] = (first[i][j] + second[min((i64)i + x, n)][j]) % mod;
            }
            first[i][j] = (first[i][j] * inv_p) % mod;

            for(int x = 1; x <= q;++x) {
                second[i][j] = (second[i][j] + first[i][min((i64)j + x, n)]) % mod;
            }
            second[i][j] = (second[i][j] * inv_q) % mod;
        }
    }

    cout << first[a][b] % mod << '\n';
}






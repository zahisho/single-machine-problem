#include <bits/stdc++.h>
#include <math.h>
#include <chrono>

#include <iostream>
#include <fstream>

using namespace std;

double hs[] = {0.2, 0.4, 0.6, 0.8};

struct job
{
  int p, a, b;
} jobs[1005];

const bool VALIDATE = 0;
const bool verbose = 0;

// Constants for GA

const int N_PAR = 10;        // Number of parents
const int N_CHD = 30;        // Number of children multiple of qty. of crossovers
const double P = 0.15;       // Probability for tournament (geometric distribution)
const double IMP_EPS = 0;    // EPS to use in improvement function
const int MAX_GEN = 1000;    // Maximum generations without improvement
const int TIME_PP = 30 * 60; // Limit time (s) by problem
const double MUTATION = 0.2; // Probability for a gen to mutate

const int CROSS_OVERS[] = {0, 1, 2, 3, 4, 5}; // Crossovers to apply
int n_crosses = sizeof(CROSS_OVERS) / sizeof(CROSS_OVERS[0]);

struct gen
{
  bool distribution[1005];
  int cost;

  bool operator<(const gen &g)
  {
    return cost < g.cost;
  }
} generation[N_PAR + N_CHD];

int ps,
    n;

int prev_job[1005], next_job[1005];
int group_a, group_b;
int expected_a[1005];
int expected_b[1005];
int ab_score[1005];
int process_order[1005];
int order_b[1005];
bool in_group_a[1005];
int best_cost;

ifstream in;
ofstream c_out;
ofstream i_out;
ofstream m_out;
ofstream log_out;

void random_code(char *s, int sz)
{
  for (int i = 0; i < sz - 1; i++)
  {
    s[i] = 'A' + (rand() % 26);
  }
  s[sz - 1] = 0;
  return;
}

void print_config()
{
  log_out << "CONFIG" << '\n';
  log_out << "\tN_PAR: " << N_PAR << '\n';
  log_out << "\tN_CHD: " << N_CHD << '\n';
  log_out << "\tP: " << P << '\n';
  log_out << "\tMUTATION: " << MUTATION << '\n';
  log_out << "\tCROSS_OVERS: ";
  for (int i = 0; i < n_crosses; i++)
  {
    if (i)
    {
      log_out << ", ";
    }
    log_out << CROSS_OVERS[i];
  }
  log_out << endl;
}

void print_group(int ini)
{
  for (int j = ini; j != -1;)
  {
    log_out << j << " ";
    j = next_job[j];
  }
  log_out << '\n';
}

int insert_a(int job)
{
  int cost_a = 0;
  int last_j = -1, lt = 0;
  bool inserted = false;
  for (int j = group_a; j != -1; j = next_job[j])
  {
    if (!inserted && jobs[job].p * jobs[j].b < jobs[j].p * jobs[job].b)
    {
      inserted = true;

      next_job[job] = j;
      prev_job[job] = prev_job[j];
      prev_job[j] = job;
      if (prev_job[job] == -1)
      {
        group_a = job;
      }
      else
      {
        next_job[prev_job[job]] = job;
      }

      lt += jobs[job].p;
      cost_a += lt * jobs[job].b;
    }
    lt += jobs[j].p;
    cost_a += lt * jobs[j].b;
    last_j = j;
  }
  if (!inserted)
  {
    if (group_a == -1)
    {
      group_a = job;
      prev_job[job] = -1;
      next_job[job] = -1;
    }
    else
    {
      next_job[last_j] = job;
      prev_job[job] = last_j;
      next_job[job] = -1;
    }

    lt += jobs[job].p;
    cost_a += lt * jobs[job].b;
  }
  in_group_a[job] = true;
  return cost_a;
}

int insert_b(int job)
{
  int cost_b = 0;
  int last_j = -1, lt = 0;
  bool inserted = false;
  for (int j = group_b; j != -1; j = next_job[j])
  {
    if (!inserted && jobs[job].p * jobs[j].a < jobs[j].p * jobs[job].a)
    {
      inserted = true;

      next_job[job] = j;
      prev_job[job] = prev_job[j];
      prev_job[j] = job;
      if (prev_job[job] == -1)
      {
        group_b = job;
      }
      else
      {
        next_job[prev_job[job]] = job;
      }

      cost_b += lt * jobs[job].a;
      lt += jobs[job].p;
    }
    cost_b += lt * jobs[j].a;
    lt += jobs[j].p;
    last_j = j;
  }
  if (!inserted)
  {
    if (group_b == -1)
    {
      group_b = job;
      prev_job[job] = -1;
      next_job[job] = -1;
    }
    else
    {
      next_job[last_j] = job;
      prev_job[job] = last_j;
      next_job[job] = -1;
    }

    cost_b += lt * jobs[job].a;
    lt += jobs[job].p;
  }
  in_group_a[job] = false;
  return cost_b;
}

void remove_a(int job)
{
  if (prev_job[job] == -1)
  {
    group_a = next_job[job];
  }
  else
  {
    next_job[prev_job[job]] = next_job[job];
  }
  if (next_job[job] != -1)
  {
    prev_job[next_job[job]] = prev_job[job];
  }
  in_group_a[job] = false;
}

void remove_b(int job)
{
  if (prev_job[job] == -1)
  {
    group_b = next_job[job];
  }
  else
  {
    next_job[prev_job[job]] = next_job[job];
  }
  if (next_job[job] != -1)
  {
    prev_job[next_job[job]] = prev_job[job];
  }
}

int cost_group_b()
{
  int cost = 0, lt = 0;
  for (int i = group_b; i != -1; i = next_job[i])
  {
    cost += lt * jobs[i].a;
    lt += jobs[i].p;
  }

  return cost;
}

int cost_group_a()
{
  int cost = 0, lt = 0;
  for (int i = group_a; i != -1; i = next_job[i])
  {
    lt += jobs[i].p;
    cost += lt * jobs[i].b;
  }

  return cost;
}

int calculate_cost()
{
  return cost_group_a() + cost_group_b();
}

void clean_groups()
{
  group_a = -1;
  group_b = -1;
  for (int i = 0; i < n; i++)
  {
    next_job[i] = -1;
    prev_job[i] = -1;
  }
}

void fill_groups(bool *arr)
{
  clean_groups();
  for (int i = 0; i < n; i++)
  {
    if (arr[i])
    {
      insert_a(i);
    }
    else
    {
      insert_b(i);
    }
  }
}

bool validate_group_b(int d)
{
  for (int i = group_b; i != -1; i = next_job[i])
  {
    d -= jobs[i].p;
  }
  return d >= 0;
}

bool validate_solution(int d)
{
  bool f_ans = true;
  bool ans[n];
  memset(ans, 0, n);
  for (int i = group_b; i != -1 && f_ans; i = next_job[i])
  {
    f_ans &= !ans[i] & !in_group_a[i];
    ans[i] = true;
  }
  if (!f_ans)
  {
    cout << "\tERROR 1" << endl;
    log_out << "\tERROR 1" << '\n';
  }
  for (int i = group_a; i != -1 && f_ans; i = next_job[i])
  {
    f_ans &= !ans[i] & in_group_a[i];
    ans[i] = true;
  }
  if (!f_ans)
  {
    cout << "\tERROR 2" << endl;
    log_out << "\tERROR 2" << '\n';
  }
  for (int i = 0; i < n && f_ans; i++)
  {
    f_ans &= ans[i];
  }
  if (!f_ans)
  {
    cout << "\tERROR 3" << endl;
    log_out << "\tERROR 3" << '\n';
  }
  return f_ans && validate_group_b(d);
}

bool long_validation(int d)
{
  bool ans = validate_solution(d);
  if (!ans)
  {
    cout << "\tERROR 4" << endl;
    log_out << "\tERROR 4" << '\n';
  }

  bool aux_group[n];
  for (int i = 0; i < n; i++)
  {
    aux_group[i] = in_group_a[i];
  }
  fill_groups(aux_group);
  for (int i = 0; i < n && ans; i++)
  {
    ans &= (aux_group[i] == in_group_a[i]);
  }
  if (!ans)
  {
    cout << "\tERROR 5" << endl;
    log_out << "\tERROR 5" << '\n';
  }

  return ans;
}

void shuffle_process_order()
{
  random_shuffle(&process_order[0], &process_order[n]);
}

/* Process order considering earliness vs tardiness score, take best between cost_a and cost_b */
int construction(int d)
{
  for (int i = 0; i < n; i++)
  {
    expected_a[i] = i;
    expected_b[i] = i;
    process_order[i] = i;
  }
  sort(expected_a, expected_a + n,
       [](int a, int b)
       { return jobs[a].p * jobs[b].b < jobs[b].p * jobs[a].b; });
  sort(expected_b, expected_b + n,
       [](int a, int b)
       { return jobs[a].p * jobs[b].a < jobs[b].p * jobs[a].a; });

  memset(ab_score, 0, n);

  for (int i = 0; i < n; i++)
  {
    ab_score[expected_b[i]] += i;
    ab_score[expected_a[i]] -= i;
  }

  for (int i = 0; i < n; i++)
  {
    order_b[expected_b[i]] = i;
  }

  sort(process_order, process_order + n,
       [](int a, int b)
       { return ab_score[a] > ab_score[b] || ab_score[a] == ab_score[b] &&
                                                 order_b[a] > order_b[b]; });
  clean_groups();
  int cost_a = 0;
  int cost_b = 0;
  for (int i = 0; i < n; i++)
  {
    int job = process_order[i];
    int ncost_a = insert_a(job);
    if (jobs[job].p <= d)
    {
      remove_a(job);
      int ncost_b = insert_b(job);
      if (ncost_b > ncost_a)
      {
        remove_b(job);
        insert_a(job);
      }
      else
      {
        cost_b = ncost_b;
        d -= jobs[job].p;
        ncost_a = cost_a;
      }
    }
    cost_a = ncost_a;
  }
  int cost = cost_a + cost_b;
  if (true && d && d < jobs[group_a].p)
  {
    int lt = d;
    cost_b = 0;
    for (int j = group_b; j != -1; j = next_job[j])
    {
      cost_b += lt * jobs[j].a;
      lt += jobs[j].p;
    }
    lt = -d;
    cost_a = 0;
    for (int j = group_a; j != -1; j = next_job[j])
    {
      lt += jobs[j].p;
      cost_a += lt * jobs[j].b;
    }
    cost = min(cost, cost_a + cost_b);
  }
  return cost;
}

int improvement(int cost, int d, double EPS)
{
  for (int j = group_b; j != -1; j = next_job[j])
  {
    d -= jobs[j].p;
  }
  while (true)
  {
    int n_cost = cost;
    int b = -1, a = -1;

    for (int i = 0; i < n; i++)
    {
      bool started_in_a = in_group_a[i];
      if (started_in_a)
      {
        remove_a(i);
        insert_b(i);
        if (d - jobs[i].p >= 0)
        {
          int aux = calculate_cost();
          if (aux < n_cost)
          {
            n_cost = aux;
            a = i;
            b = -1;
          }
        }
      }
      else
      {
        remove_b(i);
        insert_a(i);
        int aux = calculate_cost();
        if (aux < n_cost)
        {
          n_cost = aux;
          b = i;
          a = -1;
        }
      }
      for (int j = i + 1; j < n; j++)
      {
        if (started_in_a && !in_group_a[j] && d + jobs[j].p - jobs[i].p >= 0)
        {
          remove_b(j);
          insert_a(j);
          int aux = calculate_cost();
          if (aux < n_cost)
          {
            n_cost = aux;
            b = j;
            a = i;
          }
          remove_a(j);
          insert_b(j);
        }
        else if (!started_in_a && in_group_a[j] && d + jobs[i].p - jobs[j].p >= 0)
        {
          remove_a(j);
          insert_b(j);
          int aux = calculate_cost();
          if (aux < n_cost)
          {
            n_cost = aux;
            a = j;
            b = i;
          }
          remove_b(j);
          insert_a(j);
        }
      }
      if (started_in_a)
      {
        remove_b(i);
        insert_a(i);
      }
      else
      {
        remove_a(i);
        insert_b(i);
      }
    }

    if ((double)(cost - n_cost) / cost * 100 <= EPS || (a == -1 && b == -1))
    {
      break;
    }
    if (b != -1)
    {
      remove_b(b);
      insert_a(b);
      d += jobs[b].p;
    }
    if (a != -1)
    {
      remove_a(a);
      insert_b(a);
      d -= jobs[a].p;
    }
    cost = n_cost;
  }
  return cost;
}

void generate_random_solution(bool *arr, int d, double c)
{
  for (int i = 0; i < n; i++)
  {
    int j = process_order[i];
    arr[j] = (rand() % 100 < 100.0 / N_PAR * c);

    // force going to group A for time constraint
    arr[j] |= jobs[j].p > d;

    if (!arr[j])
    {
      d -= jobs[j].p;
    }
  }
  shuffle_process_order();
}

void cross1(bool *p1, bool *p2, bool *s, int d)
{
  for (int i = 0; i < n; i++)
  {
    int j = process_order[i];
    bool mutate = (rand() % 100) < (MUTATION * 100);
    s[j] = jobs[j].p > d || ((p1[j] || p2[j]) ^ mutate);
    if (!s[j])
    {
      d -= jobs[j].p;
    }
  }
  shuffle_process_order();
}

void cross2(bool *p1, bool *p2, bool *s, int d)
{
  for (int i = 0; i < n; i++)
  {
    int j = process_order[i];
    bool mutate = (rand() % 100) < (MUTATION * 100);
    s[j] = jobs[j].p > d || ((!p1[j] || !p2[j]) ^ mutate);
    if (!s[j])
    {
      d -= jobs[j].p;
    }
  }
  shuffle_process_order();
}

void cross3(bool *p1, bool *p2, bool *s, int d)
{
  for (int i = 0; i < n; i++)
  {
    int j = process_order[i];
    bool mutate = (rand() % 100) < (MUTATION * 100);
    s[j] = jobs[j].p > d || ((p1[j] && p2[j]) ^ mutate);
    if (!s[j])
    {
      d -= jobs[j].p;
    }
  }
  shuffle_process_order();
}

void cross4(bool *p1, bool *p2, bool *s, int d)
{
  for (int i = 0; i < n; i++)
  {
    int j = process_order[i];
    bool mutate = (rand() % 100) < (MUTATION * 100);
    s[j] = jobs[j].p > d || ((!p1[j] && !p2[j]) ^ mutate);
    if (!s[j])
    {
      d -= jobs[j].p;
    }
  }
  shuffle_process_order();
}

void cross5(bool *p1, bool *p2, bool *s, int d)
{
  for (int i = 0; i < n; i++)
  {
    int j = process_order[i];
    bool mutate = (rand() % 100) < (MUTATION * 100);
    s[j] = jobs[j].p > d || ((!p1[j] && p2[j]) ^ mutate);
    if (!s[j])
    {
      d -= jobs[j].p;
    }
  }
  shuffle_process_order();
}

void cross6(bool *p1, bool *p2, bool *s, int d)
{
  for (int i = 0; i < n; i++)
  {
    int j = process_order[i];
    bool mutate = (rand() % 100) < (MUTATION * 100);
    s[j] = jobs[j].p > d || ((p1[j] && !p2[j]) ^ mutate);
    if (!s[j])
    {
      d -= jobs[j].p;
    }
  }
  shuffle_process_order();
}

/* Create pointers to crosses functions to automatize calls */
void (*crosses[])(bool *, bool *, bool *, int) = {&cross1, &cross2, &cross3, &cross4, &cross5, &cross6};

default_random_engine generator;
geometric_distribution<int> distribution(P);

// geometric distribution
int run_tournament()
{
  int winner;
  do
  {
    winner = distribution(generator);
  } while (winner >= N_PAR);
  return winner;
}

void ga(int d)
{
  memcpy(generation[0].distribution, in_group_a, n);
  generation[0].cost = best_cost;

  if (verbose)
  {
    log_out << "\nPARENT " << 0 << " cost: " << generation[0].cost << "\nGROUP A\n";
    print_group(group_a);
    log_out << "GROUP B\n";
    print_group(group_b);
  }
  for (int i = 1; i < N_PAR; i++)
  {
    generate_random_solution(generation[i].distribution, d, i);
    fill_groups(generation[i].distribution);
    int aux_cost = calculate_cost();
    generation[i].cost = improvement(aux_cost, d, IMP_EPS);
    memcpy(generation[i].distribution, in_group_a, n);

    if (verbose)
    {
      log_out << "\nPARENT " << i << " cost: " << generation[i].cost << "\nGROUP A\n";
      print_group(group_a);
      log_out << "GROUP B\n";
      print_group(group_b);
    }
  }

  sort(generation, generation + N_PAR);
  if (generation[0].cost < best_cost)
  {
    cout << "\t\tIMPROVING " << best_cost << " > " << generation[0].cost << endl;
    log_out << "\t\tIMPROVING " << best_cost << " > " << generation[0].cost << '\n';
    best_cost = generation[0].cost;
  }

  int cycles = 0;
  int n_gen = 0;

  auto tini = std::chrono::high_resolution_clock::now();

  while (cycles++ < MAX_GEN)
  {
    n_gen++;
    if (verbose)
    {
      cout << "\tGENERATION: " << n_gen << endl;
      log_out << "\tGENERATION: " << n_gen << '\n';
    }

    int idx_crosses = 0;
    for (int i = N_PAR; i < N_PAR + N_CHD; idx_crosses++, i++)
    {
      int p1 = run_tournament();
      int p2 = run_tournament();

      if (verbose)
      {
        log_out << "\nSELECTED: p1 " << p1 << " cost: " << generation[p1].cost << " && p2 " << p2 << " cost: " << generation[p2].cost << "\n";
      }
      if (idx_crosses == n_crosses)
      {
        idx_crosses = 0;
      }
      int idx_aux = CROSS_OVERS[idx_crosses];
      crosses[idx_aux](generation[p1].distribution, generation[p2].distribution, generation[i].distribution, d);
      fill_groups(generation[i].distribution);

      int aux_cost = calculate_cost();

      if (verbose)
      {
        log_out << "\nCHILD " << i << " cost: " << aux_cost << "\nGROUP A\n";
        print_group(group_a);
        log_out << "GROUP B\n";
        print_group(group_b);
      }

      generation[i].cost = improvement(aux_cost, d, IMP_EPS);
      memcpy(generation[i].distribution, in_group_a, n);

      if (verbose)
      {
        log_out << "\nIMPROVED " << i << " cost: " << generation[i].cost << "\nGROUP A\n";
        print_group(group_a);
        log_out << "GROUP B\n";
        print_group(group_b);
      }
    }

    sort(generation, generation + N_PAR + N_CHD);
    if (generation[0].cost < best_cost)
    {
      cout << "\t\tIMPROVING " << best_cost << " > " << generation[0].cost << endl;
      log_out << "\t\tIMPROVING " << best_cost << " > " << generation[0].cost << '\n';
      cycles = 0;
      best_cost = generation[0].cost;
    }
    auto tend = std::chrono::high_resolution_clock::now();

    if (chrono::duration_cast<chrono::seconds>(tend - tini).count() >= TIME_PP)
    {
      cout << "Timeout!" << endl;
      log_out << "Timeout!" << '\n';
      cycles = MAX_GEN;
    }
  }
  fill_groups(generation[0].distribution);
}

int main(int argv, char **argvs)
{
  srand(time(NULL));
  char inFile[30], outFile[30], logFile[30], code[5];
  random_code(code, 5);
  cout << code << endl;

  sprintf(inFile, "dados/sch%s.txt", argvs[1]);
  sprintf(outFile, "out%s_%s.txt", argvs[1], code);
  sprintf(logFile, "log%s_%s.txt", argvs[1], code);

  in.open(inFile);
  c_out.open("c_" + string(outFile));
  i_out.open("i_" + string(outFile));
  m_out.open("m_" + string(outFile));
  log_out.open(logFile);

  cout << "INPUT " << inFile << endl;
  log_out << "INPUT " << inFile << '\n';
  cout << "OUTPUT " << outFile << endl;
  log_out << "OUTPUT " << outFile << '\n';

  print_config();

  in >> ps;
  auto ttini = std::chrono::high_resolution_clock::now();
  //
  while (ps--)
  {
    in >> n;

    int total_time = 0;

    for (int i = 0; i < n; i++)
    {
      in >> jobs[i].p >> jobs[i].a >> jobs[i].b;
      total_time += jobs[i].p;
    }

    for (int i = 0; i < 4; i++)
    {
      cout << "Problem: " << ps << " " << hs[i] << endl;
      log_out << "Problem: " << ps << " " << hs[i] << '\n';

      int d = total_time * hs[i];
      int cost = construction(d);
      // print constructive result
      c_out << cost << endl;
      best_cost = INT_MAX;
      best_cost = improvement(cost, d, IMP_EPS);
      // print improvement result
      i_out << best_cost << endl;
      ga(d);

      if (VALIDATE && !long_validation(d))
      {
        cout << "\tERROR validation\n";
        log_out << "\tERROR validation\n";
        m_out << -1 << endl;
        continue;
      }
      // print metaheuristic result
      m_out << best_cost << endl;
    }
  }
  auto ttend = chrono::high_resolution_clock::now();

  cout << "Time: " << chrono::duration_cast<chrono::microseconds>(ttend - ttini).count() << '\n';
  log_out << "Time: " << chrono::duration_cast<chrono::microseconds>(ttend - ttini).count() << '\n';

  c_out.close();
  i_out.close();
  m_out.close();
  log_out.close();

  return 0;
}

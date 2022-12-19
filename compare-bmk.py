
import numpy as np


if __name__ == "__main__":
  bmk_path = 'outputs/test6/m_out'  # '../benchmarks/bmk'  # 'outputs/test6/i_out'
  out_path = 'outputs/test-less-imp/i_out'

  print(out_path)

  ts = [10, 20, 50, 100, 200]  # , 500, 1000]
  hs = 4
  ds = []
  algo = 0
  for t in ts:
    bls = open(bmk_path+str(t)+'.txt', 'r').readlines()
    ols = open(out_path+str(t)+'.txt', 'r').readlines()

    diff = [0]*hs

    for i in range(40):
      fb = int(bls[i].strip())
      fo = int(ols[i].strip())
      diff[i % 4] += (fo-fb)/fb*10
    ds.append(diff)
    print('%s & %s & %s \\\\' %
          (t, ' & '.join(map(str, np.around(diff, 2))), "%.2f" % (sum(diff)/4)))
    #print("%.2f" % (sum(diff)/4))
    algo += sum(diff)/hs

  algo2 = 0
  for i in range(hs):
    algo2 += sum(map(lambda x: x[i], ds))/len(ts)
  print('\\textbf{Mean} & %s & %s \\\\' % (' & '.join(
      ["%.2f" % (sum(map(lambda x: x[i], ds))/len(ts)) for i in range(4)]), "%.2f" % (algo2 / 4)))
  print(algo/len(ts), algo2/4)

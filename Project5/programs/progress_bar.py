import sys

class progress_bar(object):
  """
  A command-line progress bar for use in large loops (as a sanity check).
  -----------------------------------------------------------------------
  Argument N = number of loop iterations. Usage:

  BAR = misclib.progress_bar(N)
  for i in range(N):
    # do something
    BAR.update()
  """
  def __init__(self,N,add_space=False):
    self.n = 0
    self.N = N
    self.add_space = add_space
    if add_space: print('')
    sys.stdout.write('\r['+' '*20+']    0 % ')
  def update(self):
    self.n += 1
    if self.n < self.N:
      sys.stdout.write('\r[{:20s}] {:4.0f} % '.format('='*int(20*self.n/self.N),100.*self.n/self.N))
    else:
      sys.stdout.write('\r['+'='*20+']  100 % Done!\n')
      if self.add_space: print('')


from memops.math.MathException import MathException

# private function
# a = n x n matrix
# piv = n vector
def _findPivot(a, piv):

  n = len(a)
  maxVal = 0
  maxRow = maxCol = -1
  for i in range(n):
    if piv[i] != 1:
      for j in range(n):
        if piv[j] == 0:
          m = abs(a[i][j])
          if m >= maxVal:
            maxVal = m
            maxRow = i
            maxCol = j
        elif piv[j] > 1:
          raise MathException('pivot > 1 for row = %d, col = %d' % (i, j))

  piv[maxCol] += 1

  return (maxRow, maxCol)

# private function
# a = n x n matrix
# b = n vector
def _interchangeVector(a, b, row, col):

  n = len(a)
  for j in range(n):
    a[row][j], a[col][j] = a[col][j], a[row][j]

  b[row], b[col] = b[col], b[row]

# private function
# a = n x n matrix
# b = n vector
def _pivotVector(a, b, col):

  pivInv = a[col][col]

  if pivInv == 0:
    raise MathException('pivInv = 0 for col = %d' % col)

  pivInv = 1.0 / pivInv
  n = len(a)

  a[col][col] = 1.0
  for j in range(n):
    a[col][j] *= pivInv

  b[col] *= pivInv

  for i in range(n):
    if i != col:
      x = float(a[i][col])
      a[i][col] = 0.0

      for j in range(n):
        a[i][j] -= x * a[col][j]

      b[i] -= x * b[col]

# private function
# a = n x n matrix
# rows = n vector
# cols = n vector
def _unscrambleVector(a, rows, cols):

  n = len(a)
  for j in range(n-1, -1, -1):
    row = rows[j]
    col = cols[j]
    if row != col:
      for i in range(n):
        a[i][row], a[i][col] = a[i][col], a[i][row]

# solves a x = b
# a, b input
# on output a --> a inverse and b --> x
# a = n x n matrix
# b = n vector
def gaussJordan(a, b):

  n = len(a)
  piv = n * [0]
  rows = n * [0]
  cols = n * [0]

  for i in range(n):
    (maxRow, maxCol) = _findPivot(a, piv)
    if maxRow != maxCol:
      _interchangeVector(a, b, maxRow, maxCol)

    rows[i] = maxRow
    cols[i] = maxCol

    _pivotVector(a, b, maxCol)

  _unscrambleVector(a, rows, cols)

if __name__ == '__main__':

  a = [[1, 2], [3, 4]]
  b = [1, 1]

  print 'before', a, b
  gaussJordan(a, b)
  print 'after', a, b


# Written by Tim J. Stevens
import numpy

def kMeans(data, k):

  indices = set()
  n = len(data)
  while len(indices) < k:
    indices.add(numpy.random.randint(0, n-1))

  centers = data[list(indices)]
  change = 1.0

  while change > 0.0001:


    clusters = [[] for x in range(k)]
    for vector in data:
      diffs = centers - vector
      dists = (diffs * diffs).sum(axis=1)
      closest = dists.argmin()
      clusters[closest].append(vector)

    change = 0
    
    for i, cluster in enumerate(clusters):
      cluster = numpy.array(cluster)
      center = cluster.sum(axis=0)/len(cluster)
      diff = center - centers[i]
      change += (diff * diff).sum()
      centers[i] = center
    
    #print change  
    
  return centers, clusters
  
if __name__ == '__main__':

  testDataA = numpy.random.random((1000,2)) # No clumps

  testDataB1 = numpy.random.normal(0.0, 2.0, (100,2))
  testDataB2 = numpy.random.normal(7.0, 2.0, (100,2))
  testDataB = vstack([testDataB1, testDataB2]) # Two clumps

  means, clusters = kMeans(testDataB, 2)

  from matplotlib import pyplot
  colors = ['#FF0000','#00FF00','#0000FF',
            '#FFFF00','#00FFFF','#FF00FF']

  for i, cluster in enumerate(clusters):
     x = [vec[0] for vec in cluster]
     y = [vec[1] for vec in cluster]

     color = colors[i % len(colors)]
     pyplot.scatter(x, y, c=color, marker='o')

  x = [vec[0] for vec in means]
  y = [vec[1] for vec in means]
  pyplot.scatter(x, y, s=40, c='black', marker='o')

  pyplot.show()


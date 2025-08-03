// Implementation of the Weisfeiler-Leman Algorithm for
// colored graphs.
//
// (c) 2017 Sven Reichard
//
// Input is taken from standard input.
// Format (white space is ignored):
// r n a11 a12 a13.... ann
// that is, the rank and the order followed by the rows of the adjacency matrix.
// All entries are assumed to be unsigned integers.
// The diagonal relations are separated at the beginning of the program.
//
// Example: Mobius ladder of order 8
// 2
// 8
// 0 1 0 0 1 0 0 1 
// 1 0 1 0 0 1 0 0 
// 0 1 0 1 0 0 1 0 
// 0 0 1 0 1 0 0 1 
// 1 0 0 1 0 1 0 0 
// 0 1 0 0 1 0 1 0 
// 0 0 1 0 0 1 0 1 
// 1 0 0 1 0 0 1 0 
// 
// Taken from: https://github.com/sven-reichard/stabilization/blob/master/weisfeiler.cpp


# include <iostream>
# include <iomanip>
# include <vector>
# include <set>
# include <map>
# include <list>
# include <iterator>
# include <algorithm>
# include <numeric>
# include <ctime>
# include <cstdlib>
# include <cassert>

using  std::vector;
using  std::map;
using  std::set;
using  std::list;
using  std::cout;
using  std::endl;
using  std::cin;

int order;
int rank;
typedef vector<vector<int> > Matrix;
Matrix matrix;
Matrix matrix2;
Matrix X;
Matrix Y;

vector<int> x_hash;
vector<int> y_hash;

class
counting_iterator
{
  int data;
public:
  counting_iterator(int value)
    : data(value)
  {};
  counting_iterator& operator++()
  {
    ++ data;
    return *this;
  };
  int operator *() const
  {
    return data;
  }
};

void
resizeMatrix(Matrix& m, int o=order)
{
  m.resize(o);
  //  for_each(m.begin(), m.end(), bind2nd(mem_fun(&vector<int>::resize),o));
  for (unsigned int i = 0; i < m.size(); i++)
    m[i].resize(o);
}

void expect(const std::string& s1, const std::string& s2)
{
  if (s1 != s2)
    {
      std::cerr << "expected " << s2 << ", got " << s1 << std::endl;
      exit(1);
    }
}

void read()
{
  //  cin >> rank;
  int dimension;
  std::string word;
  cin >> word >> order;
  expect(word, "order");
  cin >> word >> dimension;
  expect(word, "dimension");

  resizeMatrix(matrix);
  resizeMatrix(matrix2);
  resizeMatrix(X);
  resizeMatrix(Y);
  rank = 0;
  for (int i= 0; i < order; i++)
    for (int j = 0; j < order; j++)
      {
        cin >> matrix[i][j];
        if (matrix[i][j] > rank)
          rank = matrix[i][j];
      }
  rank ++;
};

void normalize()
{  
  set<int> diagonal;
  set<int> offDiagonal;
  for (int i = 0; i < order; i++)
    for (int j = 0; j < order; j++)
      if (i == j)
        diagonal.insert(matrix[i][j]);
      else
        offDiagonal.insert(matrix[i][j]);
  vector<int> diagonalColors(*diagonal.rbegin()+1);
  vector<int> offColors(*offDiagonal.rbegin()+1);
  rank = 0;
  for (set<int>::const_iterator iter = diagonal.begin();
       iter != diagonal.end(); iter++)
    diagonalColors[*iter] = rank++;
  for (set<int>::const_iterator iter = offDiagonal.begin();
       iter != offDiagonal.end(); iter++)
    offColors[*iter] = rank++;

  for (int i = 0; i < order; i++)
    for (int j = 0; j < order; j++)
      if (i == j)
        matrix[i][i] = diagonalColors[matrix[i][i]];
      else
        matrix[i][j] = offColors[matrix[i][j]];
};

void write()
{
  cout << order << endl;
    for (int i = 0; i < order; i++) cout << matrix[i][i] << " ";
    cout << endl;
  cout << "rank " << rank << endl;
};

void
symmetrize()
{
  vector<set<int> > values(rank);
  for (int x = 0; x < order; x++)
    for (int y = 0; y < order; y++)
      {
        int color = matrix[x][y];
        int otherColor = matrix[y][x];
        values[color].insert(otherColor);
      }
  vector<int> partialSum(values.size());
  partialSum[0] = 0;
  for (unsigned int i = 1; i < partialSum.size(); i++)
    partialSum[i] = partialSum[i-1] + values[i-1].size();
  for (int x = 0; x < order; x++)
    for (int y = 0; y < order; y++)
      {
        int oldColor = matrix[x][y];
        int pos = partialSum[oldColor];
        int value = matrix[y][x];
        for (set<int>::const_iterator iter = values[oldColor].begin();
             (iter != values[oldColor].end()) && (*iter != value);
             iter ++, pos ++);
        matrix2[x][y] = pos;

      }
  rank = partialSum.back() + values.back().size();
  swap(matrix, matrix2);
}

void
collectColors(vector<map<int, int> >& valueMap,
              const vector<set<int> >& values )
{
  rank = 0;
  for (unsigned int i = 0; i < valueMap.size(); i++)
    for (set<int>::const_iterator iter = values[i].begin();
         (iter != values[i].end());
         iter++, rank++)
      valueMap[i][*iter] = rank;
  rank = rank;
}

void
recolorMatrix2(vector<map<int, int> >& valueMap)
{

  for (int x = 0; x < order; x++)
    for (int y = 0; y < order; y++)
      {
        const int oldColor = matrix[x][y];
        const int value = matrix2[x][y];
        matrix2[x][y] = valueMap[oldColor][value];
      }
}

void
renormalizeColors(vector<set<int> >& values)
{
  for (int x = 0; x < order; x++)
    for (int y = 0; y < order; y++)
      {
        values[matrix[x][y]].insert(matrix2[x][y]);
      }
  vector<map<int, int> > valueMap(values.size());
  collectColors(valueMap, values);
  recolorMatrix2(valueMap);
}

void
selectRandomValues()
{
  x_hash.resize(rank);
  generate(x_hash.begin(), x_hash.end(), rand);
  y_hash.resize(rank);
  generate(y_hash.begin(), y_hash.end(), rand);
}

int
totalDegree(const vector<int>& color, int x)
{

  int result = 0;
  for (int y = 0; y < order; y++)
    result += x_hash[matrix[x][y]] * y_hash[color[y]];
  return result;
}

void
normalizeVector(std::vector<int>& vec)
{
  set<int> allValues(vec.begin(), vec.end());
  map<int, int> m;
  int c = 0;
  for (set<int>::const_iterator iter = allValues.begin(); iter != allValues.end(); iter++)
    m[*iter] = c++;
  for (int i = 0; i < order; i++)
    vec[i] = m[vec[i]];
};

vector<int>
degreePartition()
{
  vector<int> color(order);
  for (int i = 0; i < order; i++)
    color[i] = matrix[i][i];
  int numberOfCells = 0;
  while (true)
    {
      selectRandomValues();
      vector<int> newColor(order);
      normalizeVector(newColor);
      int newNumber = (*max_element(newColor.begin(), newColor.end())) + 1;
      if (newNumber == numberOfCells)
        break;
      numberOfCells = newNumber;
      swap(color, newColor);
    }
  return color;
}

void
substituteValues()
{
  for (int x = 0; x < order; x++)
    for (int y = 0; y < order; y++)
      {
        X[x][y] = x_hash[matrix[x][y]];
        Y[x][y] = y_hash[matrix[x][y]];
      }
}


bool WeisfeilerLemanStep()
{
  const int oldRank = rank;
  vector<set<int> > values(rank);
  selectRandomValues();
  substituteValues();
  for (int x = 0; x < order; x++)
    for (int y = 0; y < order; y++)
      {
        matrix2[x][y] = inner_product(X[x].begin(), X[x].end(), Y[y].begin(), 0);
        //values[matrix[x][y]].insert(matrix2[x][y]);
      }
  renormalizeColors( values );
  swap(matrix, matrix2);
  return oldRank != rank;
}

int main(void)
{
  srand(time(0));
  read();
  normalize();
  symmetrize();
  while (WeisfeilerLemanStep() );
  write();
  // cout << "rank: " << rank << endl; 
    for (int x = 0; x < order; x++)
    for (int y = 0; y < order; y++)
        cout << matrix[x][y] << " ";
    cout << endl;
    return 0;
};


// Local Variables:
// compile-command: "make"
// End:


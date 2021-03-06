// Programmers: Craig Olson, Stephen Chan, Danny Eiden
// Github: https://github.com/deiden26/AI_Assignment2


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <cstring>
#include <string>
#include <math.h>

using namespace std;

#define  consistencyMax 10000000

/*###### rowCol Class ######*/
class rowCol
{
public:
	int row;
	int col;
	rowCol();
	rowCol(int inputRow, int inputCol);
	void setRowCol(int r, int c);
};

/* Craig */
rowCol::rowCol() 
{
	row = 0;
	col = 0;
}
/* Danny */
rowCol::rowCol(int inputRow, int inputCol)
{
	row = inputRow;
	col = inputCol;
	return;
}
/* Craig */
void rowCol::setRowCol(int r, int c) {
	row = r;
	col = c;
}

/*###### Board Class ######*/

class Board
{
	int dim;
	int ** cells;
	long totalChecks;
	bool failed;
	int *** remainingValues; // For forward checking
	int ** index;
public:
	Board (int);
	~Board();
	string toString();
	void set_square_value(int,int,int);
	int get_square_value(int,int);
	static Board * fromFile(string);
	bool checkForVictory();
	int get_dim() {return dim;}
	rowCol findEmptySquare();
	bool FCheck(int r, int c, int x);
	bool isValidMove(rowCol position, int value);
	bool hasFailed() {return failed;}
	void setFailed(bool state);
	void printBoard();
	int numberOfConstraints(rowCol position);
	rowCol mostConstrainedFreeSquare();
  	int numberOfConstraining(rowCol position);
  	rowCol mrvMcvFreeSquare();
  	void leastConstrainingValues(rowCol position, int* values);
  	int getIndex(int r, int c) {return index[r-1][c-1];};
  	void setIndex(int r, int c, int i) {index[r-1][c-1] = i;};
  	int *** getRemainingValues() {return remainingValues;};
  	void setRemainingValues(int*** r) {remainingValues = r;};
	int* getRemainingValues(int r, int c) {return remainingValues[r-1][c-1];}
	void pushRemainingValues(int r, int c, int x);
	int popRemainingValues(int r, int c, int x);
	void resetRemainingValues(int r, int c, int x);
	void printRemainingValues(int r, int c);
	void findConstrainingValues(int r, int c);
	void initializeRemainingValues();
};



Board::Board(int d) {
	if(d > 62)
		throw ("Dimensions must be at most 62");
	dim = d;
	cells = new int*[dim];
	remainingValues = new int**[dim];
	index = new int*[dim];
	for(int i=0; i<dim;i++) {
		cells[i] = new int[dim];
		index[i] = new int[dim];
		remainingValues[i] = new int*[dim];
		for(int j=0; j<dim;j++) {
			cells[i][j] = 0;
			index[i][j] = dim;
			remainingValues[i][j] = new int[dim];
			for (int k=0; k<dim; k++) {
				remainingValues[i][j][k] = k+1;
			}
		}
	}
	failed = false;
	totalChecks = 0;
}

Board::~Board() {
	for(int i=0; i<dim;i++) {
		delete [] cells[i];
	}
	delete [] cells;
}

string Board::toString() {
	stringstream s;
	for(int i=0; i<dim;i++) {
		for(int j=0; j<dim; j++) {
			if(cells[i][j]==0)
				s << '_';
			else
				s << cells[i][j];
		}
		s << endl;
	}
	return s.str();
}

void Board::set_square_value(int row, int col, int val) {
	cells[row-1][col-1] = val;
}

int Board::get_square_value(int row, int col) {
	return cells[row-1][col-1];
}

Board * Board::fromFile(string inFileName) {
  string line;
  ifstream inFile (inFileName.c_str());
  Board * out;
  if (inFile.is_open()) {
	  getline (inFile,line);
	  int d = atoi(line.c_str());
	  out = new Board(d);
	  getline (inFile, line);
	  int numVals = atoi(line.c_str());
	  for(int i=0; i<numVals;i++) {
		string field;
		getline (inFile,field, '\t');
		int row = atoi(field.c_str());
		getline (inFile,field, '\t');
		int col = atoi(field.c_str());
		getline (inFile,field);
		int val = atoi(field.c_str());
		out->set_square_value(row, col, val);
	  }
  }
  inFile.close();
  return out;
}

bool Board::checkForVictory() {
	unsigned long victory = 0;
	//optimization: check if it's filled:
	for(int i=1; i<dim+1;i++)
		for(int j=1;j<dim+1;j++)
			if(this->get_square_value(i,j)==0)
				return false;
	for(int i=1; i<dim+1; i++)
		victory += 1 << i;
	//check rows and columns:
	for(int i=1; i<dim+1;i++) {
		unsigned long rowTotal = 0;
		unsigned long columnTotal = 0;
		for(int j=1; j<dim+1; j++) {
			rowTotal += 1 << this->get_square_value(i, j);
			columnTotal += 1 << this->get_square_value(j, i);
		}
		if(rowTotal!=victory||columnTotal!=victory)
			return false;
	}
	int dimsqrt = (int)(sqrt((double)dim));
	//check little squares:
	for(int i=0;i<dimsqrt;i++) {
		for(int j=0;j<dimsqrt;j++) {
			unsigned long squareTotal = 0;
			for(int k=1; k<dimsqrt+1;k++) {
				for(int m=1; m<dimsqrt+1;m++) {
					squareTotal += 1 << this->get_square_value(i*dimsqrt+k, j*dimsqrt+m);
				}
			}
			if(squareTotal != victory)
				return false;
		}
	}
	return true;
}
/* Danny */
rowCol Board::findEmptySquare()
{
	/* Search every element of the board  until	*/
	/* an element is found with value == 0		*/
	for(int i=1; i<dim+1;i++)
	{
		for(int j=1;j <dim+1;j++)
		{
			if(this->get_square_value(i,j)==0)
			{
				return rowCol(i, j);
			}
		}
	}
	/* returning [0,0] signifies that no	*/
	/* free square was found				*/
	return rowCol(0, 0);
}

/* Danny */
bool Board::isValidMove(rowCol position, int newValue)
{
	for (int i=1; i<dim+1 ;i++)
	{
		/* Check column for any square with value == newValue */
		if (this->get_square_value(i, position.col) == newValue)
		{
			return false;
		}
		/* Check row for any square with value == newValue */
		if (this->get_square_value(position.row, i) == newValue)
		{
			return false;
		}
	}

	/* Get the dimension of the subBoard */
	int dimsqrt = (int)(sqrt((double)dim));
	/* Get the row and column of the subBoard (zero indexed) */
	int subBoardRow = ceil((float)position.row/(float)dimsqrt)-1;
	int subBoardCol = ceil((float)position.col/(float)dimsqrt)-1;

	/* Check subBoard for any square with value == newValue */
	for(int k=1; k<dimsqrt+1;k++)
	{
		for(int m=1; m<dimsqrt+1;m++)
		{
			if (this->get_square_value(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m) == newValue)
			{
				return false;
			}
		}
	}
	/* If no important squares have been found	*/
	/* with value == newValue, return true	 	*/
	return true;
}

/* Danny */
void Board::setFailed(bool state)
{
	failed = state;
	return;
}

/* Danny */
void Board::printBoard()
{
	/* Print top, horizontal line of the board */
	cout << endl << endl;
	for (int k=0; k<dim*5; k++)
	{
		cout << "-";
	}
	cout << endl;

	int dimsqrt = (int)(sqrt((double)dim)); 

	/* Print the numbers of the board by row */
	for(int i=1; i<dim+1;i++)
	{
		for(int j=1; j<dim+1; j++)
		{
			if (this->get_square_value(i, j) == 0)
			{
				cout << "_" << setw(4);
			}
			else
			{
				cout << this->get_square_value(i, j) << setw(4);
			}
			/* print vertical lines of the board*/
			if (j % dimsqrt == 0)
			{
				cout << "|";
				if (j != dim)
				{
					cout << setw(4);
				}
			}
		}
		cout << endl;
		/* Print the horizontal lines of the board */
		if (i % dimsqrt == 0)
		{
			for (int k=0; k<dim*5; k++)
			{
				cout << "-";
			}
			cout << endl;
		}
	}
	return;
}

/* Danny */
int Board::numberOfConstraints(rowCol position)
{
  // need to make this a variable sized array depending on size of board
	int constraints[dim];
 	memset(constraints, 0, dim * sizeof(int));
	int constraintCount = 0;

	for (int i=1; i<dim+1 ;i++)
	{
		/* Check column for any square with value != 0 */
		if (this->get_square_value(i, position.col) != 0)
		{
			constraints[this->get_square_value(i, position.col) -1] = 1;
		}
		/* Check row for any square with value != 0 */
		if (this->get_square_value(position.row, i) != 0)
		{
			constraints[this->get_square_value(position.row, i) -1] = 1;
		}
	}

	/* Get the dimension of the subBoard */
	int dimsqrt = (int)(sqrt((double)dim));
	/* Get the row and column of the subBoard (zero indexed) */
	int subBoardRow = ceil((float)position.row/(float)dimsqrt)-1;
	int subBoardCol = ceil((float)position.col/(float)dimsqrt)-1;

	/* Check subBoard for any square with value != 0 */
	for(int k=1; k<dimsqrt+1;k++)
	{
		for(int m=1; m<dimsqrt+1;m++)
		{
			if (this->get_square_value(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m) != 0)
			{
				constraints[this->get_square_value(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m) -1] = 1;
			}
		}
	}

	/* Sum Constraints */
	for (int i=0; i < dim; i++)
	{
		constraintCount += constraints[i];
	}
  return constraintCount;

}

/* Danny */
rowCol Board::mostConstrainedFreeSquare()
{
  rowCol mostConstrained = rowCol(0, 0);
  int maxConstraints = 0;
	for(int i=1; i < dim+1;i++)
	{
		for(int j=1; j < dim+1;j++)
		{
      		if (this->get_square_value(i, j) != 0)
      		{
          		// skip squares that are already filled
          		continue;
      		}
			int currentConstraints = this->numberOfConstraints(rowCol(i,j));
			/* If the number of contraints == dim for any square, there is no valid move for it */
			/* Not allowing this quare to be returned will prevent backtracking */
			if (currentConstraints == dim)
			{
				mostConstrained = rowCol(i,j);
				return mostConstrained;
			}
			if (currentConstraints > maxConstraints)
			{
				maxConstraints = currentConstraints;
				mostConstrained = rowCol(i,j);
			}
		}
	}
	/* returning [0,0] signifies that no	*/
	/* free square was found				*/
	return mostConstrained;
}


/* Stephen */
int Board::numberOfConstraining(rowCol position)
{
  /* For each square, there are up to 
   * 2 * (dim - 1) + (dim - 1) * (dim - 1) free squares
   * that it could impose a constraint on upon assignment.
   * How many free squares are in the constraint scope of 
   * this square? */

  //int constraining[possibleFreeInScope];
  //memset(constraining, 0, possibleFreeInScope * sizeof(int));

	int constrainingCount = 0;

	for (int i=1; i<dim+1 ;i++)
	{
		/* Check column for any square with value == 0 */
		if (this->get_square_value(i, position.col) == 0)
		{
			constrainingCount += 1;
		}
		/* Check row for any square with value == 0 */
		if (this->get_square_value(position.row, i) == 0)
		{
			constrainingCount += 1;
		}
	}

	/* Get the dimension of the subBoard */
	int dimsqrt = (int)(sqrt((double)dim));
	/* Get the row and column of the subBoard (zero indexed) */
	int subBoardRow = ceil((float)position.row/(float)dimsqrt)-1;
	int subBoardCol = ceil((float)position.col/(float)dimsqrt)-1;

	/* Check subBoard for any square with value == 0 */
	for(int k=1; k<dimsqrt+1;k++)
	{
		for(int m=1; m<dimsqrt+1;m++)
		{
			// We don't want to double count any of the free squares that were already counted
			// in the row and column
			if ( subBoardRow * dimsqrt + k == position.row || subBoardCol * dimsqrt + m == position.col )
			{
				continue;
			}
			if (this->get_square_value(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m) == 0)
			{
				constrainingCount += 1;
			}
		}
	}

  return constrainingCount;
}

/* Craig */
void Board::pushRemainingValues(int r, int c, int x)
{
	int ind = index[r-1][c-1];
	if (ind != this->get_dim()) {
		remainingValues[r-1][c-1][ind] = x;
		ind++;
	} 
	index[r-1][c-1] = ind;
	return;
}

/* Craig */
int Board::popRemainingValues(int r, int c, int x)
{
	int ind = index[r-1][c-1];
	for (int i=0; i<ind; i++) {
		if (remainingValues[r-1][c-1][i] == x) {
			//remainingValues[r-1][c-1][i] = -1;
			for (; i < ind; i++) {
				remainingValues[r-1][c-1][i] = remainingValues[r-1][c-1][i+1];
			} 
			ind--;
			index[r-1][c-1] = ind;
			return x;
		}
	}
	return -1;
}

/* Craig */
void Board::resetRemainingValues(int r, int c, int x) {
	// Find the restraining values and then remove them from the list of remaining values
	rowCol row, col, sub;

	for (int i=1; i<dim+1 ;i++)
	{
		row.setRowCol(r, i);
		col.setRowCol(i, c);
		if (this->isValidMove(row, x))
			pushRemainingValues(r, i, x);
		if (this->isValidMove(col, x))
			pushRemainingValues(i, c, x);
		
	}

	/* Get the dimension of the subBoard */
	int dimsqrt = (int)(sqrt((double)dim));
	/* Get the row and column of the subBoard (zero indexed) */
	int subBoardRow = ceil((float)r/(float)dimsqrt)-1;
	int subBoardCol = ceil((float)c/(float)dimsqrt)-1;

	/* Check subBoard for any square with value == 0 */
	for(int k=1; k<dimsqrt+1;k++)
	{
		for(int m=1; m<dimsqrt+1;m++)
		{
	      // We don't want to double count any of the free squares that were already counted
	      // in the row and column
	      if ( subBoardRow * dimsqrt + k == r || subBoardCol * dimsqrt + m == c ) {
	        continue;
	      }
		  sub.setRowCol(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m);
		  if (this->isValidMove(sub, x))
		 	 pushRemainingValues(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m, x);
			
		}
	}
	return;
}

/* Craig */
void Board::printRemainingValues(int r, int c) {
	cout << "Remaining Values:\n\t";
	int ind = index[r-1][c-1];
	for (int i=0; i<ind; i++) {
		cout << remainingValues[r-1][c-1][i] << "\t";
	}
	cout << endl;
}
/* Craig */
void Board::findConstrainingValues(int r, int c) {
	// Find the restraining values and then remove them from the list of remaining values
	
	for (int i=1; i<dim+1 ;i++)
	{
		/* Check column for any square with value == 0 */
		if (this->get_square_value(i, c) != 0)
		{
			popRemainingValues(r, c, this->get_square_value(i, c));
		}
		/* Check row for any square with value == 0 */
		if (this->get_square_value(r, i) != 0)
		{
			popRemainingValues(r, c, this->get_square_value(r, i));
		}
	}

	/* Get the dimension of the subBoard */
	int dimsqrt = (int)(sqrt((double)dim));
	/* Get the row and column of the subBoard (zero indexed) */
	int subBoardRow = ceil((float)r/(float)dimsqrt)-1;
	int subBoardCol = ceil((float)c/(float)dimsqrt)-1;

	/* Check subBoard for any square with value == 0 */
	for(int k=1; k<dimsqrt+1;k++)
	{
		for(int m=1; m<dimsqrt+1;m++)
		{
	      // We don't want to double count any of the free squares that were already counted
	      // in the row and column
	      if ( subBoardRow * dimsqrt + k == r || subBoardCol * dimsqrt + m == c ) {
	        continue;
	      }
			if (this->get_square_value(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m) != 0)
			{
				popRemainingValues(r, c, this->get_square_value(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m));
			}
		}
	}
	
}

/* Craig */
void Board::initializeRemainingValues() {
	for (int i=1; i<=dim; i++) {
		for (int j=1; j<=dim; j++) {
			if (cells[i-1][j-1] == 0) {
				findConstrainingValues(i, j);
			}
		}
	}
}

/* Craig & Danny */
bool Board::FCheck(int r, int c, int x) {
	// Find the restraining values and then remove them from the list of remaining values
	
	for (int i=1; i<dim+1 ;i++)
	{
		if (i != c)
		{
			if (this->get_square_value(r, i) == 0)
			{
				popRemainingValues(r, i, x);
				if (this->index[r-1][i-1] == 0)
				{
					resetRemainingValues(r, c, x);
					return false;
				}
			}
		}
		if (i != r)
		{
			if (this->get_square_value(i, c) == 0)
			{
				popRemainingValues(i, c, x);
				if (this->index[i-1][c-1] == 0)
				{
					resetRemainingValues(r, c, x);
					return false;
				}
			}
		}
	}

	/* Get the dimension of the subBoard */
	int dimsqrt = (int)(sqrt((double)dim));
	/* Get the row and column of the subBoard (zero indexed) */
	int subBoardRow = ceil((float)r/(float)dimsqrt)-1;
	int subBoardCol = ceil((float)c/(float)dimsqrt)-1;

	/* Check subBoard for any square with value == 0 */
	for(int k=1; k<dimsqrt+1;k++)
	{
		for(int m=1; m<dimsqrt+1;m++)
		{
	      // We don't want to double count any of the free squares that were already counted
	      // in the row and column
	      if ( subBoardRow * dimsqrt + k == r || subBoardCol * dimsqrt + m == c )
	      {
	        continue;
	      }
		  if (this->get_square_value(subBoardRow * dimsqrt + k, subBoardCol * dimsqrt + m))
		  {
		  	continue;
		  }
		  popRemainingValues(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m, x);
		  if (this->index[subBoardRow*dimsqrt+k-1][subBoardCol*dimsqrt+m-1] == 0)
		  {
		  	resetRemainingValues(r, c, x);
		  	return false;
		  }
			
		}
	}
	return true;
}

/* Stephen Danny */
rowCol Board:: mrvMcvFreeSquare()
/* Chooses minimum remaining value square of the board.
 * Ties are broken using the constraining variable count */
{
  rowCol mostConstrained = rowCol(0, 0);
  int maxConstraints = 0;
	for(int i=1; i < dim+1;i++)
	{
		for(int j=1; j < dim+1;j++)
		{
			if (this->get_square_value(i, j) != 0)
			{
			  // skip squares that are already filled
			  continue;
			}
			rowCol current = rowCol(i, j);
			int currentConstraints = this->numberOfConstraints(rowCol(i,j));
			if (currentConstraints > maxConstraints)
			{
				maxConstraints = currentConstraints;
				mostConstrained = current;
			}
			if (currentConstraints == maxConstraints)
			{
				// break the tie!
				int mostConstrainedConstraining = 0;
				// Test for mostConstrained == (0,0) which isn't a real board position 
				if (mostConstrained.row != 0)
				{
					mostConstrainedConstraining= numberOfConstraining(mostConstrained);
				}
				int currentConstraining = numberOfConstraining(current);
				if (mostConstrainedConstraining < currentConstraining)
				{
					// then replace 
					mostConstrained = current;
				}
			}
		}
	}
	/* returning [0,0] signifies that no	*/
	/* free square was found				*/
	return mostConstrained;
}

/* Stephen Danny */
void Board::leastConstrainingValues(rowCol position, int* values)
{
	/* Returns a list of all values that can be used to fill the square at
	 * position in increasing order of how many constraints the value would cause */
	int totalConstraints[dim];
	int currentTotalConstraints;

	/* Get the dimension of the subBoard */
	int dimsqrt = (int)(sqrt((double)dim));
	/* Get the row and column of the subBoard (zero indexed) */
	int subBoardRow = ceil((float)position.row/(float)dimsqrt)-1;
	int subBoardCol = ceil((float)position.col/(float)dimsqrt)-1;

	for (int j=1; j < dim+1; j++)
	{
		if (this->isValidMove(position,j))
		{
			currentTotalConstraints = 0;
			/* Fill square with the test value */
			this->set_square_value(position.row, position.col, j);
			/* Count the number of constraints on every variable that the square constrains */
			for (int i=1; i<dim+1 ;i++)
			{
				/* Check column for any square with value == 0 */
				if (this->get_square_value(i, position.col) == 0)
				{
					currentTotalConstraints += this->numberOfConstraints(rowCol(i, position.col));
				}
				/* Check row for any square with value == 0 */
				if (this->get_square_value(position.row, i) == 0)
				{
					currentTotalConstraints += this->numberOfConstraints(rowCol(position.row, i));
				}
			}

			/* Check subBoard for any square with value == 0 */
			for(int k=1; k<dimsqrt+1;k++)
			{
				for(int m=1; m<dimsqrt+1;m++)
				{
					// We don't want to double count any of the free squares that were already counted
					// in the row and column
					if ( subBoardRow * dimsqrt + k == position.row || subBoardCol * dimsqrt + m == position.col )
					{
						continue;
					}
					if (this->get_square_value(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m) == 0)
					{
						currentTotalConstraints += this->numberOfConstraints(rowCol(subBoardRow*dimsqrt+k, subBoardCol*dimsqrt+m));
					}
				}
			}
			/* Store the constraint count in an array with indeces
			 * corresponding to possible square values */
			totalConstraints[j-1] = currentTotalConstraints;
		}
		else
		{
			/* If the move isn't allowed, give it the highest possible constraint count */
			totalConstraints[j-1] = 2147483646; //Highest 32bit int - 1 (Essential Infinity)
		}
	}


	/* This section is used to order potential values
	 * by the number of variables they would constrain */
	int lowestPosition;
	int lowestConstraint;
	/* At each itteration of the loop, find the lowest number in the array */
	for (int j=0; j<dim; j++)
	{
		lowestConstraint = 2147483647; //Highest 32bit int (Essential Infinity)
		/* Find the lowest number in the array */
		for (int i=0; i<dim; i++)
		{
	        if(totalConstraints[i]<lowestConstraint)
	        {
		        lowestConstraint=totalConstraints[i];
		    	lowestPosition = i;
		    }
		}
		/* Use the position of the lowest number to find what value it corresponds to */
		values[j] = lowestPosition+1;
		/* Set the number of this value to (basically) Infinity so it won't be selected again */
		totalConstraints[lowestPosition] = 2147483647; //Highest 32bit int (Essential Infinity)
	}
	/* Clear the value of the square that was being tested */
	this->set_square_value(position.row, position.col, 0);
	return;
}

/*###### BackTracking Search ######*/
/* Danny */
int recursiveBackTrackingSearch(Board *currentBoard, int consistencyCount)
{
	//currentBoard->printBoard(); //See board at every recursion
	/* Find a square to attempt to fill */
	rowCol emptySquare = currentBoard->findEmptySquare();
	/* If there isn't a free square, the search is complete */
	if (emptySquare.row == 0)
		return consistencyCount;
	/* Try all possible values for the given free square */
	for (int i=1; i<currentBoard->get_dim()+1; i++)
	{
		consistencyCount++;
		if (consistencyCount >= consistencyMax)
		{
			return consistencyCount;
		}
		//cout << i << " "; //See each attempted input
		if (currentBoard->isValidMove(emptySquare,i))
		{
			/* If you find a number that is allowed, fill it in and recurse */
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, i);
			consistencyCount = recursiveBackTrackingSearch(currentBoard, consistencyCount);

			if(!(currentBoard->hasFailed()))
			{
				return consistencyCount;
			}
			/* If the search failed, erase value and attempt with different value */
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, 0);
			currentBoard->setFailed(false);
		}
	}
	/* Fail if there isn't a valid number to put in the selected free square */
	currentBoard->setFailed(true);
	return consistencyCount;
}

int backTrackingSearch(Board *initialBoard)
{
	int numberOfConsistencyChecks = recursiveBackTrackingSearch(initialBoard, 0);
	return numberOfConsistencyChecks;
}

/*###### BackTracking Search W/ MRV+MCV ######*/
/* Stephen Danny */
int recursiveBackTrackingSearchMrvMcv(Board *currentBoard, int consistencyCount)
{
	//currentBoard->printBoard(); //See board at every recursion
	/* Find a square to attempt to fill */
	rowCol emptySquare = currentBoard->mrvMcvFreeSquare();
	/* If there isn't a free square, the search is complete */
	if (emptySquare.row == 0)
		return consistencyCount;
	/* Try all possible values for the given free square */
	for (int i=1; i<currentBoard->get_dim()+1; i++)
	{
		consistencyCount++;
		if (consistencyCount >= consistencyMax)
		{
			return consistencyCount;
		}
		//cout << i << " "; //See each attempted input
		if (currentBoard->isValidMove(emptySquare,i))
		{
			/* If you find a number that is allowed, fill it in and recurse */
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, i);
			consistencyCount = recursiveBackTrackingSearchMrvMcv(currentBoard, consistencyCount);

			if(!(currentBoard->hasFailed()))
			{
				return consistencyCount;
			}
			/* If the search failed, erase value and attempt with different value */
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, 0);
			currentBoard->setFailed(false);
		}
	}
	/* Fail if there isn't a valid number to put in the selected free square */
	currentBoard->setFailed(true);
	return consistencyCount;
}

int backTrackingSearchMrvMcv(Board *initialBoard)
{
	int numberOfConsistencyChecks = recursiveBackTrackingSearchMrvMcv(initialBoard, 0);
	return numberOfConsistencyChecks;
}
/*###### BackTracking Search W/ MRV+MCV+LCV ######*/
/* Stephen Danny */
int recursiveBackTrackingSearchMrvMcvLcv(Board *currentBoard, int consistencyCount)
{
	//currentBoard->printBoard(); //See board at every recursion
	/* Find a square to attempt to fill */
	rowCol emptySquare = currentBoard->mrvMcvFreeSquare();
	/* If there isn't a free square, the search is complete */
	if (emptySquare.row == 0)
		return consistencyCount;
	/* Get ordered list of values to try */
	int values[currentBoard->get_dim()];
	currentBoard->leastConstrainingValues(emptySquare, values);
	/* Try all possible values for the given free square */
	for (int i=0; i<currentBoard->get_dim(); i++)
	{
		int attemptValue = values[i];
		consistencyCount++;
		if (consistencyCount >= consistencyMax)
		{
			return consistencyCount;
		}
		//cout << attemptValue << " "; //See each attempted input
		if (currentBoard->isValidMove(emptySquare,attemptValue))
		{
			/* If you find a number that is allowed, fill it in and recurse */
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, attemptValue);
			consistencyCount = recursiveBackTrackingSearchMrvMcvLcv(currentBoard, consistencyCount);

			if(!(currentBoard->hasFailed()))
			{
				return consistencyCount;
			}
			/* If the search failed, erase value and attempt with different value */
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, 0);
			currentBoard->setFailed(false);
		}
	}
	/* Fail if there isn't a valid number to put in the selected free square */
	currentBoard->setFailed(true);
	return consistencyCount;
}

int backTrackingSearchMrvMcvLcv(Board *initialBoard)
{
	int numberOfConsistencyChecks = recursiveBackTrackingSearchMrvMcvLcv(initialBoard, 0);
	return numberOfConsistencyChecks;
}

/*###### BackTracking Search W/ Forward Checking ######*/
/* Craig Danny */
int recursiveBackTrackingSearchFCheck(Board *currentBoard, int consistencyCount)
{
	//currentBoard->printBoard(); //See board at every recursion
	/* Find a square to attempt to fill */
	rowCol emptySquare = currentBoard->findEmptySquare();
	rowCol* remainingValueLocations;
	/* If there isn't a free square, the search is complete */
	if (emptySquare.row == 0)
		return consistencyCount;
	/* Try all possible values for the given free square */
	for (int i=1; i<currentBoard->get_dim()+1; i++)
	{
		consistencyCount++;
		if (consistencyCount >= consistencyMax)
		{
			return consistencyCount;
		} 
		//cout << i << " "; //See each attempted input
		if (currentBoard->isValidMove(emptySquare,i) && currentBoard->FCheck(emptySquare.row, emptySquare.col, i))
		{
			/* If you find a number that is allowed, fill it in and recurse */
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, i);
			// Now forward check, but maintain all current remaining values
			// Need to save current remaining values and put them back, if necessary
			//currentBoard->FCheck(emptySquare.row, emptySquare.col, i);
			consistencyCount = recursiveBackTrackingSearchFCheck(currentBoard, consistencyCount);


			if(!(currentBoard->hasFailed()))
			{
				return consistencyCount;
			}
			/* If the search failed, erase value and attempt with different value */
			// Reset remainingValues
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, 0);
			currentBoard->resetRemainingValues(emptySquare.row, emptySquare.col, i);
			
			currentBoard->setFailed(false);
		}
	}
	/* Fail if there isn't a valid number to put in the selected free square */
	currentBoard->setFailed(true);
	return consistencyCount;
}

int backTrackingSearchFCheck(Board *initialBoard)
{
	// Write function that initializes remainingValues 
	initialBoard->initializeRemainingValues();
	int numberOfConsistencyChecks = recursiveBackTrackingSearchFCheck(initialBoard, 0);
	return numberOfConsistencyChecks;
}


/*###### main Function ######*/
/* Danny */
int main(int argc, char* argv[])
{
	/*~~~~ 4x4 board BackTracking ~~~~*/
	Board *inputBoard_4x4 = Board::fromFile("4x4.sudoku");

	cout << "4x4 Board BackTracking\n";
	int numberOfConsistencyChecks_4x4 = backTrackingSearch(inputBoard_4x4);

	if (inputBoard_4x4->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 4x4 board MRV+MCV BackTracking ~~~~*/
	inputBoard_4x4 = Board::fromFile("4x4.sudoku");

	cout << "4x4 Board MRV+MCV BackTracking\n";
	int numberOfConsistencyChecks_4x4_MRV = backTrackingSearchMrvMcv(inputBoard_4x4);
	
	if (inputBoard_4x4->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";
	/*~~~~ 4x4 board MRV+MCV+LVC BackTracking ~~~~*/
	inputBoard_4x4 = Board::fromFile("4x4.sudoku");

	cout << "4x4 Board MRV+MCV+LVC BackTracking\n";
	int numberOfConsistencyChecks_4x4_MRV_LVC = backTrackingSearchMrvMcvLcv(inputBoard_4x4);
	
	if (inputBoard_4x4->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 4x4 board Forward Checking BackTracking ~~~~*/
	inputBoard_4x4 = Board::fromFile("4x4.sudoku");

	cout << "4x4 Board Forward Checking BackTracking\n";
	int numberOfConsistencyChecks_4x4_FCheck = backTrackingSearchFCheck(inputBoard_4x4);
	
	if (inputBoard_4x4->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 9x9 board ~~~~*/
	Board *inputBoard_9x9 = Board::fromFile("9x9.sudoku");

	cout << "9x9 Board BackTracking\n";
	int numberOfConsistencyChecks_9x9 = backTrackingSearch(inputBoard_9x9);

	if (inputBoard_9x9->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 9x9 board MRV+MCV BackTracking ~~~~*/
	inputBoard_9x9 = Board::fromFile("9x9.sudoku");

	cout << "9x9 Board MRV+MCV BackTracking\n";
	int numberOfConsistencyChecks_9x9_MRV = backTrackingSearchMrvMcv(inputBoard_9x9);
	
	if (inputBoard_9x9->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 9x9 board MRV+MCV+LVC BackTracking ~~~~*/
	inputBoard_9x9 = Board::fromFile("9x9.sudoku");

	cout << "9x9 Board MRV+MCV+LVC BackTracking\n";
	int numberOfConsistencyChecks_9x9_MRV_LVC = backTrackingSearchMrvMcvLcv(inputBoard_9x9);
	
	if (inputBoard_9x9->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 9x9 board Forward Checking BackTracking ~~~~*/
	inputBoard_9x9 = Board::fromFile("9x9.sudoku");

	cout << "9x9 Board Forward Checking BackTracking\n";
	int numberOfConsistencyChecks_9x9_FCheck = backTrackingSearchFCheck(inputBoard_9x9);
	
	if (inputBoard_9x9->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 16x16 board ~~~~*/
	Board *inputBoard_16x16 = Board::fromFile("16x16.sudoku");

	cout << "16x16 Board BackTracking\n";
	int numberOfConsistencyChecks_16x16 = backTrackingSearch(inputBoard_16x16);

	if (inputBoard_16x16->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 16x16 board MRV+MCV BackTracking ~~~~*/
	inputBoard_16x16 = Board::fromFile("16x16.sudoku");

	cout << "16x16 Board MRV+MCV BackTracking\n";
	int numberOfConsistencyChecks_16x16_MRV = backTrackingSearchMrvMcv(inputBoard_16x16);
	
	if (inputBoard_16x16->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 16x16 board MRV+MCV+LVC BackTracking ~~~~*/
	inputBoard_16x16 = Board::fromFile("16x16.sudoku");

	cout << "16x16 Board MRV+MCV+LVC BackTracking\n";
	int numberOfConsistencyChecks_16x16_MRV_LVC = backTrackingSearchMrvMcvLcv(inputBoard_16x16);
	
	if (inputBoard_16x16->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 16x16 board Forward Checking BackTracking ~~~~*/
	inputBoard_16x16 = Board::fromFile("16x16.sudoku");

	cout << "16x16 Board Forward Checking BackTracking\n";
	int numberOfConsistencyChecks_16x16_FCheck = backTrackingSearchFCheck(inputBoard_16x16);
	
	if (inputBoard_16x16->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 25x25 board ~~~~*/
	Board *inputBoard_25x25 = Board::fromFile("25x25.sudoku");

	cout << "25x25 Board BackTracking\n";
	int numberOfConsistencyChecks_25x25 = backTrackingSearch(inputBoard_25x25);

	if (inputBoard_25x25->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 25x25 board MRV+MCV BackTracking ~~~~*/
	inputBoard_25x25 = Board::fromFile("25x25.sudoku");

	cout << "25x25 Board MRV+MCV BackTracking\n";
	int numberOfConsistencyChecks_25x25_MRV = backTrackingSearchMrvMcv(inputBoard_25x25);
	
	if (inputBoard_25x25->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 25x25 board MRV+MCV+LVC BackTracking ~~~~*/
	inputBoard_25x25 = Board::fromFile("25x25.sudoku");

	cout << "25x25 Board MRV+MCV+LVC BackTracking\n";
	int numberOfConsistencyChecks_25x25_MRV_LVC = backTrackingSearchMrvMcvLcv(inputBoard_25x25);
	
	if (inputBoard_25x25->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~~ 25x25 board Forward Checking BackTracking ~~~~*/
	inputBoard_25x25 = Board::fromFile("25x25.sudoku");

	cout << "25x25 Board Forward Checking BackTracking\n";
	int numberOfConsistencyChecks_25x25_FCheck = backTrackingSearchFCheck(inputBoard_25x25);
	
	if (inputBoard_25x25->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/*~~~ Performance Table ~~~*/
	cout << "\n------------------------------------------------------------------------------\n";
	cout << "| "  << "Problem" << setw(4) << "| " << "BackTracking" << setw(5) << "| " << "MRV + MCV" << setw (3) << "| ";
	cout << "MRV + MCV + LCV" << setw(3) << "| " << "Forward Checking" << setw(2) << "|" << endl;
	cout << "| -------------------------------------------------------------------------- |\n";
	/* 4x4 Row */
	cout << "| " << "4x4" << setw(8) << "| " << numberOfConsistencyChecks_4x4 << setw(14) << "| ";
	cout << numberOfConsistencyChecks_4x4_MRV << setw(10) << "| ";
	cout << numberOfConsistencyChecks_4x4_MRV_LVC << setw(16) << "| ";
	cout << numberOfConsistencyChecks_4x4_FCheck << setw(16) << "|" << endl;
	/* 9x9 Row */
	cout << "| " << "9x9" << setw(8) << "| " << numberOfConsistencyChecks_9x9 << setw(11) << "| ";
	cout << numberOfConsistencyChecks_9x9_MRV << setw(9) << "| ";
	cout << numberOfConsistencyChecks_9x9_MRV_LVC << setw(14) << "| ";
	cout << numberOfConsistencyChecks_9x9_FCheck << setw(12) << "|" << endl;
	/* 16x16 Row */
	cout << "| " << "16x16" << setw(6) << "| " << numberOfConsistencyChecks_16x16 << setw(9) << "| ";
	cout << numberOfConsistencyChecks_16x16_MRV << setw(6) << "| ";
	cout << numberOfConsistencyChecks_16x16_MRV_LVC << setw(12) << "| ";
	cout << numberOfConsistencyChecks_16x16_FCheck << setw(10) << "|" << endl;

	/* 25x25 Row */
	cout << "| " << "25x25" << setw(6) << "| " << numberOfConsistencyChecks_25x25 << setw(9) << "| ";
	cout << numberOfConsistencyChecks_25x25_MRV << setw(4) << "| ";
	cout << numberOfConsistencyChecks_25x25_MRV_LVC << setw(10) << "| ";
	cout << numberOfConsistencyChecks_25x25_FCheck << setw(10) << "|" << endl;
	cout << endl;

	return 1;
}

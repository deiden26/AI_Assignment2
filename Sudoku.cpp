// Sudoku.cpp : Basic class for holding a Sudoku board, reading a board from files, a writing a board to the screen
//

#include<iostream>
#include<fstream>
#include<sstream>

#include<string>
#include<math.h>

using namespace std;

/*###### rowCol Class ######*/
class rowCol
{
public:
	int row;
	int col;
	rowCol(int inputRow, int inputCol);
};

rowCol::rowCol(int inputRow, int inputCol)
{
	row = inputRow;
	col = inputCol;
	return;
}

/*###### Board Class ######*/

class Board
{
	int dim;
	int ** cells;
	long totalChecks;
	bool failed;
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
	bool isValidMove(rowCol position, int value);
	bool hasFailed() {return failed;}
	void setFailed(bool state);
	void printBoard();
};

Board::Board(int d) {
	if(d > 62)
		throw ("Dimensions must be at most 62");
	dim = d;
	cells = new int*[dim];
	for(int i=0; i<dim;i++) {
		cells[i] = new int[dim];
		for(int j=0; j<dim;j++)
			cells[i][j] = 0;
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
rowCol Board::findEmptySquare()
{
	for(int i=1; i<dim+1;i++)
	for(int j=1;j <dim+1;j++)
		if(this->get_square_value(i,j)==0)
			return rowCol(i, j);
	return rowCol(0, 0);
}

bool Board::isValidMove(rowCol position, int value)
{
	for (int i=1; i<dim+1 ;i++)
	{
		/* Check column */
		if (this->get_square_value(i, position.col) == value)
		{
			return false;
		}
		/* Check row */
		if (this->get_square_value(position.row, i) == value)
		{
			return false;
		}
	}

	int dimsqrt = (int)(sqrt((double)dim));
	int squareRow = ceil((float)position.row/(float)dimsqrt)-1;
	int squareCol = ceil((float)position.col/(float)dimsqrt)-1;

	for(int k=1; k<dimsqrt+1;k++)
	{
		for(int m=1; m<dimsqrt+1;m++)
		{
			if (this->get_square_value(squareRow*dimsqrt+k, squareCol*dimsqrt+m) == value)
			{
				return false;
			}
		}
	}
	return true;
}

void Board::setFailed(bool state)
{
	failed = state;
	return;
}

void Board::printBoard()
{
	cout << endl << endl;
	for (int k=0; k<dim*10; k++)
	{
		cout << "-";
	}
	cout << endl;
	int dimsqrt = (int)(sqrt((double)dim)); 
	for(int i=1; i<dim+1;i++)
	{
		for(int j=1; j<dim+1; j++)
		{
			cout << this->get_square_value(i, j) <<  "\t";
			if (j % dimsqrt == 0)
			{
				cout << "|\t";
			}
		}
		cout << endl;
		if (i % dimsqrt == 0)
		{
			for (int k=0; k<dim*10; k++)
			{
				cout << "-";
			}
			cout << endl;
		}
	}
	return;
}


/*###### BackTracking Search ######*/


void recursiveBackTrackingSearch(Board *currentBoard)
{
	//currentBoard->printBoard(); //See board at every recursion
	rowCol emptySquare = currentBoard->findEmptySquare();
	if (emptySquare.row == 0)
		return;
	for (int i=1; i<currentBoard->get_dim()+1; i++)
	{
		//cout << i << " "; //See each attempted input
		if (currentBoard->isValidMove(emptySquare,i))
		{
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, i);
			recursiveBackTrackingSearch(currentBoard);
			if(!(currentBoard->hasFailed()))
			{
				return;
			}
			currentBoard->set_square_value(emptySquare.row, emptySquare.col, 0);
			currentBoard->setFailed(false);
		}
	}
	currentBoard->setFailed(true);
	return;
}

void backTrackingSearch(Board *initialBoard)
{
	recursiveBackTrackingSearch(initialBoard);
	return;
}

/*###### main Function ######*/

int main(int argc, char* argv[])
{
	/* 4x4 board */
	Board *inputBoard_4x4 = Board::fromFile("4x4.sudoku");
	inputBoard_4x4->printBoard();
	backTrackingSearch(inputBoard_4x4);
	inputBoard_4x4->printBoard();

	if (inputBoard_4x4->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/* 9x9 board */
	Board *inputBoard_9x9 = Board::fromFile("9x9.sudoku");
	inputBoard_9x9->printBoard();
	backTrackingSearch(inputBoard_9x9);
	inputBoard_9x9->printBoard();

	if (inputBoard_9x9->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/* 16x16 board */
	Board *inputBoard_16x16 = Board::fromFile("16x16.sudoku");
	inputBoard_16x16->printBoard();
	backTrackingSearch(inputBoard_16x16);
	inputBoard_16x16->printBoard();

	if (inputBoard_16x16->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	/* 25x25 board */
	Board *inputBoard_25x25 = Board::fromFile("25x25.sudoku");
	inputBoard_25x25->printBoard();
	backTrackingSearch(inputBoard_25x25);
	inputBoard_25x25->printBoard();

	if (inputBoard_25x25->checkForVictory())
		cout << "Victory!\n";
	else
		cout << "Defeat\n";

	return 1;
}

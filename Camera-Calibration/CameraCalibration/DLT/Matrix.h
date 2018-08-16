#pragma 
#include <math.h>
#include <vector>
using namespace std;

#define PI 3.1415926535897931

typedef vector<vector<double>> mMatrix;


// �������溯��1,matrixΪ���뷽��NΪ����Ĵ�С,PΪ���ص������
mMatrix Inverse( mMatrix& matrix )
{
	double d = 0.0;
	double temp = 0.0;
	int m = matrix.size();//���������
	int n = matrix[0].size();//���������
	mMatrix resultMatrix(m);	
	for (int i = 0; i < m; i++)
	{	resultMatrix[i].resize(n);
	for (int j = 0; j < n; j++)
		resultMatrix[i][j] = matrix[i][j];
	}
	int *mainRow = new int[m];
	int *mainCol = new int[m];
	if (m == n)//��������������
	{
		for (int k = 0; k <= m - 1; k++)
		{
			d = 0.0;
			for (int i = k; i <= m - 1; i++)
			{
				for (int j = k; j <= m - 1; j++)
				{
					temp = abs(resultMatrix[i][j]);
					if (temp > d)
					{
						d = temp; mainRow[k] = i; mainCol[k] = j;
					}
				}
			}
			if (mainRow[k] != k)
			{
				for (int j = 0; j <= m - 1; j++)
				{
					temp = resultMatrix[mainRow[k]][j]; resultMatrix[mainRow[k]][j] = resultMatrix[k][j]; resultMatrix[k][j] = temp;
				}
			}
			if (mainCol[k] != k)
			{
				for (int i = 0; i <= m - 1; i++)
				{
					temp = resultMatrix[i][mainCol[k]]; resultMatrix[i][mainCol[k]] = resultMatrix[i][k]; resultMatrix[i][k] = temp;
				}
			}
			resultMatrix[k][k] = 1.0 / resultMatrix[k][k];
			for (int j = 0; j <= m - 1; j++)
			{
				if (j != k)
				{
					resultMatrix[k][j] *= resultMatrix[k][k];
				}
			}

			for (int i = 0; i <= m - 1; i++)
			{
				if (i != k)
				{
					for (int j = 0; j <= m - 1; j++)
					{
						if (j != k)
						{
							resultMatrix[i][j] -= resultMatrix[i][k] * resultMatrix[k][j];
						}
					}
				}
			}
			for (int i = 0; i <= m - 1; i++)
			{
				if (i != k)
				{
					resultMatrix[i][k] = -resultMatrix[i][k] * resultMatrix[k][k];
				}
			}
		}
		for (int k = m - 1; k >= 0; k--)
		{
			if (mainCol[k] != k)
			{
				for (int j = 0; j <= m - 1; j++)
				{
					temp = resultMatrix[k][j]; resultMatrix[k][j] = resultMatrix[mainCol[k]][j]; resultMatrix[mainCol[k]][j] = temp;
				}
			}
			if (mainRow[k] != k)
			{
				for (int i = 0; i <= m - 1; i++)
				{
					temp = resultMatrix[i][k]; resultMatrix[i][k] = resultMatrix[i][mainRow[k]]; resultMatrix[i][mainRow[k]] = temp;
				}
			}
		}
	}
	return resultMatrix;
}


/// ������˺���
bool Multiply( const mMatrix& matrix1, const mMatrix& matrix2, mMatrix& resultMatrix)
{
	resultMatrix.clear();
	int m1, n1, m2, n2;
	double temp;
	m1 = matrix1.size(); n1 = matrix1[0].size();
	m2 = matrix2.size(); n2 = matrix2[0].size();
	resultMatrix.resize(m1);
	for (int i = 0; i < m1; i++)
	{
		resultMatrix[i].resize(n2);
		for (int j = 0; j < n2; j++)
		{
			temp = 0;
			for (int k = 0; k < n1; k++)
			{
				temp = temp + matrix1[i][k] * matrix2[k][j];
			}
			resultMatrix[i][j] = temp;
		}
	}
	return false;
}

bool Add( const mMatrix& matrix1, const mMatrix& matrix2, mMatrix& resultMatrix,int flag )
{
	resultMatrix.clear();
	int m1, n1, m2, n2;

	m1 = matrix1.size(); n1 = matrix1[0].size();
	m2 = matrix2.size(); n2 = matrix2[0].size();
	resultMatrix.resize(m1);

	if (m1 == m2 && n1 == n2)//���������ά����ͬʱ���Ž�����ӻ�����ļ��㣬������������ֵȫΪ��
	{
		for (int i = 0; i < m1; i++)
		{
			resultMatrix[i].resize(n2);
			for (int j = 0; j < n2; j++)
			{
				resultMatrix[i][j] = matrix1[i][j] + flag * matrix2[i][j];//��flag=1ʱ��Ϊ������ӣ���flag=-1ʱ��Ϊ�������
			}
		}
		return true;

	}
	return false;
}

bool Transpose( const mMatrix& matrix, mMatrix& resultMatrix )
{
	resultMatrix.clear();
	int m, n;//��������д�С
	m = matrix.size(); n = matrix[0].size();
	resultMatrix.resize(n);
	for (int i = 0; i < n; i++)
	{
		resultMatrix[i].resize(m);
		for (int j = 0; j < m; j++)
		{
			resultMatrix[i][j] = matrix[j][i];//ת��
		}
	}
	return false;
}
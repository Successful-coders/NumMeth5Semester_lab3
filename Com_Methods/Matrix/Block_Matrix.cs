using System;
using System.IO;

namespace Com_Methods
{
    class Block_Matrix : IMatrix
    {
        //размер матрицы
        public int M { set; get; }
        public int N { set; get; }
        public int Size_Block { set; get; }
        //элементы матрицы
        public Matrix[][] Block { set; get; }

        //конструктор по умолчанию
        public Block_Matrix() { }

        //конструктор по размеру
        public Block_Matrix(int Row, int Column, int SIZE_BLOCK) 
        {
            M = Row;
            N = Column;
            Size_Block = SIZE_BLOCK;
            Block = new Matrix[Row][];
            for (int i = 0; i < Row; i++)
            {
                Block[i] = new Matrix[Column];
                for (int j = 0; j < Column; j++)
                    Block[i][j] = new Matrix(Size_Block, Size_Block);
            }
        }

        //конструктор матрицы по бинарному файлу
        public Block_Matrix(string PATH, int SIZE_BLOCK)
        {
            //чтение размера системы
            using (var Reader = new BinaryReader(File.Open(PATH + "Size.bin", FileMode.Open)))
            {
                M = Reader.ReadInt32();
                N = M;
            }
            //проверка размера
            if (M % SIZE_BLOCK != 0) throw new Exception("Block_Matrix: error in the block size...");
            //размер матрицы
            M /= SIZE_BLOCK;
            N /= SIZE_BLOCK;
            Size_Block = SIZE_BLOCK;
            Block = new Matrix[M][];

            //чтение матрицы
            using (var Reader = new BinaryReader(File.Open(PATH + "Matrix.bin", FileMode.Open)))
            {
                //считываем каждое значение из файла
                try
                {
                    for (int i = 0; i < M; i++)
                    {
                        //выделили место под блочную матрицу
                        Block[i] = new Matrix[N];

                        //выделяем место под каждый блок
                        for (int j = 0; j < N; j++) Block[i][j] = new Matrix(Size_Block, Size_Block);

                        //заполняем строки в блоках i-ой строки блочной матрицы 
                        for (int ii = 0; ii < Size_Block; ii++)
                        {
                            for (int j = 0; j < N; j++)
                                for (int k = 0; k < Size_Block; k++)
                                {
                                    Block[i][j].Elem[ii][k] = Reader.ReadDouble();
                                }
                        }

                        //диагональный блок необходимо преобразовать в LU-разложение
                        var LU_Decomp = new LU_Decomposition(Block[i][i]);
                        Block[i][i]   = LU_Decomp.LU;  
                    }
                }
                catch { throw new Exception("Block_Matrix: data file is not correct..."); }
            }
        }
    }
}

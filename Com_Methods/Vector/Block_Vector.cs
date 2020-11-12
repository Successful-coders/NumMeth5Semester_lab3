using System;
using System.IO;

namespace Com_Methods
{
    class Block_Vector : IVector
    {
        //размер вектора
        public int N { set; get; }
        public int Size_Block { set; get; }
        public Vector[] Block { set; get; }

        //конструктор по умолчанию
        public Block_Vector() { }

        //конструктор по размеру
        public Block_Vector(int Size, int SIZE_BLOCK)
        {
            N = Size;
            Size_Block = SIZE_BLOCK;
            Block = new Vector[Size];
            for (int i = 0; i < Size; i++)
                Block[i] = new Vector(Size_Block);
        }

        //конструктор вектора по бинарному файлу
        public Block_Vector(String PATH, int SIZE_BLOCK)
        {
            //чтение размера системы
            using (var Reader = new BinaryReader(File.Open(PATH + "Size.bin", FileMode.Open)))
            {
                N = Reader.ReadInt32();
            }

            if (N % SIZE_BLOCK != 0) throw new Exception("Block_Vector: error in the block size...");

            //размер вектора
            N /= SIZE_BLOCK;
            Size_Block = SIZE_BLOCK;
            Block = new Vector[N];

            //чтение вектора правой части
            using (var Reader = new BinaryReader(File.Open(PATH + "F.bin", FileMode.Open)))
            {
                // считываем данные из файла
                try
                {
                    for (int j = 0; j < N; j++)
                    {
                        Block[j] = new Vector(Size_Block);
                        for (int k = 0; k < Size_Block; k++)
                            Block[j].Elem[k] = Reader.ReadDouble();
                    }
                }
                catch { throw new Exception("Block_Vector: data file is not correct..."); }
            }
        }

        //вывод вектора на консоль
        public void Console_Write_Vector()
        {
            for (int i = 0; i < N; i++)
                for (int k = 0; k < Size_Block; k++)
                {
                    Console.WriteLine("X[{0}] = {1, -20}", i * Size_Block + k + 1, Block[i].Elem[k]);
                }
        }
    }
}

using System;


namespace Com_Methods
{
    //блочный метод релаксации
    class SOR_Method : IIteration_Solver
    {
        //максимальное число итераций
        public int Max_Iter { set; get; }
        //точность решения
        public double Eps { set; get; }
        //текущая итерация 
        public int Iter { set; get; }

        //конструктор
        public SOR_Method(int MAX_ITER, double EPS)
        {
            Max_Iter = MAX_ITER;
            Eps = EPS;
            Iter = 0;
        }

        //точечный метод релаксации: реализация
        public Vector Start_Solver(Matrix A, Vector F, double w)
        {
            //норма разности решений на двух итерациях
            double Norm_Xnew_Xold;

            //инициализация вектора результата
            var RES = new Vector(F.N);
            var RES_New = new Vector(F.N);

            //начальное приближение
            for (int k = 0; k < RES.N; k++) RES.Elem[k] = 0.0;

            //итерации метода релаксации
            do
            {
                Norm_Xnew_Xold = 0.0;
                //цикл по неизвестным
                for (int i = 0; i < RES.N; i++)
                {
                    //инициализация скобки (F - Ax)
                    double F_Ax = F.Elem[i];

                    //произведение нижнего треугольника матрицы на новые X
                    for (int j = 0; j < i; j++)
                    {
                        F_Ax -= A.Elem[i][j] * RES_New.Elem[j];
                    }

                    //произведение верхнего треугольника матрицы на предыдущие X
                    for (int j = i + 1; j < A.N; j++)
                    {
                        F_Ax -= A.Elem[i][j] * RES.Elem[j];
                    }

                    //формируем результат для i-ой компоненты и квадрат нормы
                    F_Ax /= A.Elem[i][i];
                    RES_New.Elem[i] = (1.0 - w) * RES.Elem[i] + w * F_Ax;
                    Norm_Xnew_Xold += Math.Pow(RES.Elem[i] - RES_New.Elem[i], 2);
                }

                //копирование полученного результата
                RES.Copy(RES_New);

                //норма разности решений
                Norm_Xnew_Xold = Math.Sqrt(Norm_Xnew_Xold);
                Iter++;

                Console.WriteLine("Iter {0,-10} {1}", Iter, Norm_Xnew_Xold/RES.Norm());
            }
            while (Norm_Xnew_Xold > Eps && Iter < Max_Iter);
            return RES;
        }

        //метод блочной релаксации: реализация
        public Block_Vector Start_Solver(Block_Matrix A, Block_Vector F, double w)
        {
            //норма разности решений на двух итерациях
            double Norm_Xnew_Xold;

            //вектор результата на текущей и следующей итерации
            var RES = new Block_Vector(A.N, A.Size_Block);
            var RES_New = new Block_Vector(A.N, A.Size_Block);

            //LU-решатель для обращения диагональных блоков
            var LU_Solver = new LU_Decomposition();

            //вспомогательный вектор для вычисления скобки (F - Ax)
            var F_Ax = new Vector(A.Size_Block);

            //итерации метода блочной релаксации
            do
            {
                Norm_Xnew_Xold = 0.0;
                //цикл по неизвестным
                for (int i = 0; i < RES.N; i++)
                {
                    //инициализация суммы
                    for (int k = 0; k < RES.Size_Block; k++) F_Ax.Elem[k] = F.Block[i].Elem[k];

                    //произведение блоков матрицы на новые X
                    for (int j = 0; j < i; j++)
                    {
                        var Current_Matrix_Block = A.Block[i][j];
                        var Current_Vector_Block = RES.Block[j];
                        for (int Row = 0; Row < Current_Matrix_Block.M; Row++)
                            for (int Column = 0; Column < Current_Matrix_Block.N; Column++)
                                F_Ax.Elem[Row] -= Current_Matrix_Block.Elem[Row][Column] * Current_Vector_Block.Elem[Column];
                    }

                    //произведение блоков матрицы на предыдущие X
                    for (int j = i + 1; j < A.N; j++)
                    {
                        var Current_Matrix_Block = A.Block[i][j];
                        var Current_Vector_Block = RES.Block[j];
                        for (int Row = 0; Row < Current_Matrix_Block.M; Row++)
                            for (int Column = 0; Column < Current_Matrix_Block.N; Column++)
                                F_Ax.Elem[Row] -= Current_Matrix_Block.Elem[Row][Column] * Current_Vector_Block.Elem[Column];
                    }

                    //решение СЛАУ с диагональной матрицей
                    LU_Solver.LU = A.Block[i][i];
                    F_Ax = LU_Solver.Start_Solver(F_Ax);

                    //формируем результат для i-ой компоненты
                    for (int k = 0; k < RES.Size_Block; k++)
                    {
                        double X_NEW = (1 - w) * RES.Block[i].Elem[k] + w * F_Ax.Elem[k];
                        Norm_Xnew_Xold += Math.Pow(X_NEW - RES.Block[i].Elem[k], 2);
                        RES.Block[i].Elem[k] = X_NEW;
                    }
                }

                //норма разности решений
                Norm_Xnew_Xold = Math.Sqrt(Norm_Xnew_Xold);
                Iter++;

                Console.WriteLine("Iter {0,-10} {1}", Iter, Norm_Xnew_Xold);
            }
            while (Norm_Xnew_Xold > Eps && Iter < Max_Iter);
            return RES;
        }
    }
}
